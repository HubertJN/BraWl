!----------------------------------------------------------------------!
! Wang-Landau module                                                   !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module wang_landau
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

contains

  subroutine wl_main(setup, wl_setup, my_rank)
    ! Rank of this processor
    integer, intent(in) :: my_rank
    integer :: ierror, num_proc

    ! Input file data types
    type(run_params) :: setup
    type(wl_params) :: wl_setup

    ! Loop integers and error handling variable
    integer :: i, j, ierr

    ! Temperature and temperature steps
    real(real64) :: temp, beta

    ! wl variables and arrays
    integer :: bins, resets
    real(real64) :: bin_width, energy_to_ry, wl_f, wl_f_prev, tolerance, flatness_tolerance
    real(real64) :: target_energy, min_val, flatness, acceptance, bins_buffer, bins_min, radial_min, radial_min_buffer
    real(real64), allocatable :: bin_edges(:), wl_hist(:), wl_logdos(:), bin_energy(:)
    logical :: first_reset, hist_reset, first_hist_reset, rho_saved

    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file

    ! window variables
    integer, allocatable :: window_indices(:, :)
    integer :: num_windows, num_walkers

    ! mpi variables
    real(real64) :: start, end, time_max, time_min
    integer :: mpi_bins, mpi_start_idx, mpi_end_idx, mpi_index, mpi_start_idx_buffer, mpi_end_idx_buffer, beta_index
    integer, allocatable :: intervals(:,:)
    real(real64) :: scale_factor, scale_count, wl_logdos_min, bin_overlap, beta_diff, beta_original, beta_merge
    real(real64), allocatable :: mpi_bin_edges(:), mpi_wl_hist(:), wl_logdos_buffer(:), wl_logdos_write(:)
    real(real64), allocatable :: rank_time(:), rank_time_buffer(:,:)

    ! radial density across energy
    real(real64), allocatable :: rho_of_E(:,:,:,:), rho_of_E_buffer(:,:,:,:)
    integer, allocatable :: radial_record(:, :)

    ! Minimum value in array to be considered
    min_val = wl_setup%tolerance*1e-1_real64

    ! Radial densities as a function of energy
    allocate(rho_of_E(setup%n_species, setup%n_species, setup%wc_range, wl_setup%bins))
    if (my_rank == 0) then
      allocate(rho_of_E_buffer(setup%n_species, setup%n_species, setup%wc_range, wl_setup%bins))
    end if
    rho_of_E = 0.0_real64

    ! Path to radial file and name of radial file
    radial_file = "radial_densities/rho_of_E.dat"

    ! Get number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierror)

    ! Check if number of MPI processes is divisible by number of windows
    bins = wl_setup%bins
    num_windows = wl_setup%num_windows
    num_walkers = num_proc/num_windows
    wl_setup%radial_samples = wl_setup%radial_samples/num_walkers
    if (MOD(num_proc, num_windows) /= 0) then
      if (my_rank == 0) then
        write (6, '(72("~"))')
        write (6, '(5("~"),x,"Error: Number of MPI processes not divisible by num_windows",x,6("~"))')
        write (6, '(72("~"))')
      end if
      call MPI_FINALIZE(ierror)
      call EXIT(0)
    end if

    ! Get start and end indices for energy windows
    allocate(window_indices(num_windows, 2))
    allocate(intervals(num_windows, 2))
    
    call divide_range(1, bins, num_windows, intervals)
    bin_overlap = wl_setup%bin_overlap

    window_indices(1, 1) = intervals(1,1)
    window_indices(1,2) = INT(intervals(1,2) + ABS(intervals(2,1)-intervals(2,2))*bin_overlap)
    do i = 2, num_windows-1
      window_indices(i, 1) = INT(intervals(i,1) - ABS(intervals(i-1,1)-intervals(i-1,2))*bin_overlap)
      window_indices(i, 2) = INT(intervals(i,2) + ABS(intervals(i+1,1)-intervals(i+1,2))*bin_overlap)
    end do
    window_indices(num_windows, 1) = INT(intervals(num_windows,1) - ABS(intervals(num_windows-1,1) &
                                    -intervals(num_windows-1,2))*bin_overlap)
    window_indices(num_windows,2) = intervals(num_windows,2)

    mpi_index = my_rank/num_walkers + 1
    mpi_start_idx = window_indices(mpi_index, 1)
    mpi_end_idx = window_indices(mpi_index, 2)
    mpi_bins = mpi_end_idx - mpi_start_idx + 1

    if (my_rank == 0) then
      print*, "Window Indices"
      print*, window_indices(:, 1)
      print*, window_indices(:, 2)
      print*, "Window Bins"
      print*, window_indices(:, 2) - window_indices(:, 1) + 1
    end if

    ! allocate arrays
    allocate(radial_record(intervals(mpi_index,2) - intervals(mpi_index,1) + 1, 2))
    allocate(bin_edges(bins + 1))
    allocate(wl_hist(bins))
    allocate(wl_logdos(bins))
    allocate(wl_logdos_buffer(bins))
    allocate(bin_energy(bins))

    ! MPI arrays
    allocate(mpi_bin_edges(mpi_bins + 1))
    allocate(mpi_wl_hist(mpi_bins))
    allocate(rank_time(num_windows))
    allocate(rank_time_buffer(num_windows,3))

    if (my_rank == 0) then
      allocate(wl_logdos_write(bins))
    end if

    !  Set temperature in appropriate units for simulation
    temp = setup%T*k_b_in_Ry
    beta = 1.0_real64/temp

    ! Conversion meV/atom to Rydberg
    energy_to_ry = setup%n_atoms/(eV_to_Ry*1000)

    ! Create energy bins and set mpi bins
    j = 1
    bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)*energy_to_ry
    do i = 1, bins + 1
      bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i - 1)*bin_width
    end do
    do i = 1, bins
      bin_energy(i) = wl_setup%energy_min*energy_to_ry + (i - 0.5)*bin_width
    end do
    do i = mpi_start_idx, mpi_end_idx + 1
      mpi_bin_edges(j) = bin_edges(i)
      j = j + 1
    end do

    ! Target energy for burn in
    target_energy = (mpi_bin_edges(1) + mpi_bin_edges(SIZE(mpi_bin_edges)))/2

    ! Initialize
    wl_hist = 0.0_real64; wl_logdos = 0.0_real64; wl_f = wl_setup%wl_f; wl_f_prev = wl_f
    flatness = 0.0_real64; first_reset = .False.; hist_reset = .True.; first_hist_reset = .False.
    mpi_wl_hist = 0.0_real64; wl_logdos_buffer = 0.0_real64; rank_time = 0.0_real64
    radial_record = wl_setup%radial_samples
    rho_saved = .False.

    ! Set up the lattice
    call initial_setup(setup, config)
    call lattice_shells(setup, shells, config)

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      write (6, '(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
      print *, "Number of atoms", setup%n_atoms
    end if
    call comms_wait()

    !---------!
    ! Burn in !
    !---------!
    call wl_burn_in(setup, config, target_energy, MINVAL(mpi_bin_edges), MAXVAL(mpi_bin_edges), &
                    num_proc, num_walkers, my_rank, mpi_index)
    print*, "Rank: ", my_rank, "Burn-in complete"
    call comms_wait()
    if (my_rank == 0) then
      write (*, *)
      write (6, '(27("-"),x,"Burn-in complete",x,27("-"),/)')
      write (*, *)
    end if

    !--------------------!
    ! Main Wang-Landau   !
    !--------------------!
    resets = 1
    tolerance = wl_setup%tolerance
    flatness_tolerance = wl_setup%flatness
    do while (wl_f > tolerance)
      if (hist_reset .eqv. .True.) then
        ! Start timer
        start = mpi_wtime()
        hist_reset = .False.
      end if
      acceptance = run_wl_sweeps(setup, wl_setup, config, temp, bins, bin_edges, mpi_start_idx, mpi_end_idx, &
                                 mpi_wl_hist, wl_logdos, wl_f, mpi_index, intervals, radial_record, rho_of_E)

      flatness = minval(mpi_wl_hist)/(sum(mpi_wl_hist)/mpi_bins)
      bins_min = count(mpi_wl_hist > min_val)/REAL(mpi_bins)
      radial_min = REAL(sum(radial_record(:,2)))/REAL(wl_setup%radial_samples*SIZE(radial_record(:,2)))

      if (first_hist_reset .neqv. .True.) then
        write (6, '(a,i0,a,f6.2,a,f6.2,a)') "Rank: ", my_rank, " | Bins visited: ", bins_min*100_real64, &
              "% | Flatness:", flatness*100_real64, "%"
      end if

      if (first_reset .neqv. .True.) then ! First reset after system had time to explore
        if (minval(mpi_wl_hist) > 10.0_real64) then
          first_reset = .True.
          mpi_wl_hist = 0.0_real64
          radial_record = 0
        end if
      end if
      !Check if we're flatness_tolerance% flat
      if (flatness > flatness_tolerance .and. minval(mpi_wl_hist) > 10) then
        !Reset the histogram
        mpi_wl_hist = 0.0_real64
        hist_reset = .True.
        first_hist_reset = .True.
        !Reduce f
        wl_f = wl_f*1/2
        ! End timer
        end = mpi_wtime()
        wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > 0.0_real64))!wl_logdos_min
        wl_logdos = ABS(wl_logdos * merge(0, 1, wl_logdos < 0.0_real64))
        !Average DoS across walkers
        do i = 1, num_windows
          if (mpi_index == i) then
            if (my_rank /= (i - 1)*num_walkers) then
              call MPI_Send(wl_logdos, bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, ierror)
              !print*, my_rank, "send", (i - 1)*num_walkers
            end if
            if (my_rank == (i - 1)*num_walkers) then
              do j = 1, num_walkers - 1
                call MPI_Recv(wl_logdos_buffer, bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                              MPI_STATUS_IGNORE, ierror)
                !print*, my_rank, "recv", (i - 1)*num_walkers + j
                wl_logdos = wl_logdos + wl_logdos_buffer
              end do
            end if
            if (my_rank == (i - 1)*num_walkers) then
              wl_logdos = wl_logdos/num_walkers
              do j = 1, num_walkers - 1
                call MPI_Send(wl_logdos, bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                              ierror)
                !print*, my_rank, "send", (i - 1)*num_walkers + j
              end do
            else
              call MPI_Recv(wl_logdos_buffer, bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE, ierror)
              !print*, my_rank, "recv", (i - 1)*num_walkers
            end if
          end if
        end do
        call comms_wait()
      end if

      if (hist_reset .eqv. .True.) then
        rank_time = 0.0_real64
        rank_time(mpi_index) = end - start
        call MPI_ALLREDUCE(end - start, time_max, 1, MPI_DOUBLE_PRECISION, &
        MPI_MAX, MPI_COMM_WORLD, ierror)

        call MPI_REDUCE(end - start, time_min, 1, MPI_DOUBLE_PRECISION, &
        MPI_MIN, 0, MPI_COMM_WORLD, ierror)

        call MPI_REDUCE(rank_time/num_walkers, rank_time_buffer(:,1), num_windows, MPI_DOUBLE_PRECISION, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        call MPI_REDUCE(rank_time, rank_time_buffer(:,3), num_windows, MPI_DOUBLE_PRECISION, &
        MPI_MAX, 0, MPI_COMM_WORLD, ierror)

        rank_time = time_max
        rank_time(mpi_index) = end - start

        call MPI_REDUCE(rank_time, rank_time_buffer(:,2), num_windows, MPI_DOUBLE_PRECISION, &
        MPI_MIN, 0, MPI_COMM_WORLD, ierror)

        call MPI_REDUCE(radial_min, radial_min_buffer, 1, MPI_DOUBLE_PRECISION, &
        MPI_SUM, 0, MPI_COMM_WORLD, ierror)
        radial_min_buffer = radial_min_buffer/num_proc

        call comms_wait()

        if (my_rank == 0) then
          wl_logdos_write = wl_logdos
          write (6, '(a,f20.18,a,f8.2,a)', advance='no') "Flatness reached for W-L F of: ", wl_f_prev, &
                  " | Radial samples: ", radial_min_buffer*100_real64, "%"
            write (*, *)
          do i=1, num_windows
            write (6, '(a,i0,a,f8.2,a,f8.2,a,f8.2,a)') "MPI Window: ", i, " | Avg. time: ", rank_time_buffer(i,1), &
            "s | Time min: ", rank_time_buffer(i,2), "s Time max: " , rank_time_buffer(i,3), "s"
          end do
          wl_f_prev = wl_f
          write (*, *)
        end if

        ! MPI send and recieve calls for combining window DoS
        do i = 2, num_windows
          if (my_rank == (i - 1)*num_walkers) then
            call MPI_Send(wl_logdos, bins, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
            call MPI_Send(mpi_start_idx, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
            call MPI_Send(mpi_end_idx, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
          end if
          if (my_rank == 0) then
            call MPI_Recv(wl_logdos_buffer, bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, 0, MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, ierror)
            call MPI_Recv(mpi_start_idx_buffer, 1, MPI_INT, (i - 1)*num_walkers, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_Recv(mpi_end_idx_buffer, 1, MPI_INT, (i - 1)*num_walkers, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            scale_factor = 0.0_real64
            scale_count = 0.0_real64
            beta_diff = wl_logdos_write(2) - wl_logdos_write(1)
            do j = 0,  window_indices(i-1, 2) - window_indices(i, 1) - 1
              beta_original = wl_logdos_write(mpi_start_idx_buffer + j + 1) - wl_logdos_write(mpi_start_idx_buffer + j)
              beta_merge = wl_logdos_buffer(mpi_start_idx_buffer + j + 1) - wl_logdos_buffer(mpi_start_idx_buffer + j)
              if (ABS(beta_original - beta_merge) < beta_diff) then
                beta_diff = ABS(beta_original - beta_merge)
                beta_index = mpi_start_idx_buffer + j
              end if
            end do
            scale_factor = wl_logdos_write(beta_index) - wl_logdos_buffer(beta_index)
            do j = beta_index, mpi_end_idx_buffer
              wl_logdos_write(j) = wl_logdos_buffer(j) + scale_factor
            end do
          end if
        end do

        !rho_saved
        if (rho_saved .eqv. .False. .and. INT(radial_min) == 100) then
          
          rho_saved = .True.
          call MPI_REDUCE(rho_of_E, rho_of_E_buffer, setup%n_species*setup%n_species*setup%wc_range*wl_setup%bins, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
          if (my_rank == 0) then
            print*, "Radial densities saved"
            rho_of_E_buffer = rho_of_E_buffer/wl_setup%radial_samples/num_walkers
            call ncdf_radial_density_writer_across_energy(radial_file, rho_of_E_buffer, shells, bin_energy, setup)
          end if
        end if

        if (my_rank == 0) then
          ! Write output files
          call ncdf_writer_1d("wl_dos_bins.dat", ierr, bin_edges)
          call ncdf_writer_1d("wl_dos.dat", ierr, wl_logdos_write)
          call ncdf_writer_1d("wl_hist.dat", ierr, wl_hist)
        end if
      end if
    end do

    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine wl_main

  integer function bin_index(energy, bin_edges, bins) result(index)
    integer, intent(in) :: bins
    real(real64), intent(in) :: energy
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64) :: bin_range

    bin_range = bin_edges(bins + 1) - bin_edges(1)
    index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  function run_wl_sweeps(setup, wl_setup, config, temp, bins, bin_edges, &
                         mpi_start_idx, mpi_end_idx, mpi_wl_hist, wl_logdos, wl_f, &
                         mpi_index, intervals, radial_record, rho_of_E) result(acceptance)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: bins, mpi_start_idx, mpi_end_idx, mpi_index
    integer, dimension(:,:), intent(inout) :: radial_record, intervals
    real(real64), dimension(:,:,:,:), intent(inout) :: rho_of_E
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64), dimension(:), intent(inout) :: mpi_wl_hist, wl_logdos
    real(real64), intent(in) :: wl_f, temp

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e, beta
    integer :: acceptance, i, ibin, jbin, iradial
    integer(int16) :: site1, site2

    ! Set inverse temp
    beta = 1.0_real64/temp

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped
    acceptance = 0.0_real64

    do i = 1, wl_setup%mc_sweeps*setup%n_atoms

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      call pair_swap(config, rdm1, rdm2)

      ! Calculate energy if different species
      if (site1 /= site2) then
        e_swapped = setup%full_energy(config)
      end if

      ibin = bin_index(e_unswapped, bin_edges, bins)
      jbin = bin_index(e_swapped, bin_edges, bins)

      ! Only compute energy change if within limits where V is defined and within MPI region
      if (jbin > mpi_start_idx - 1 .and. jbin < mpi_end_idx + 1) then

        ! Calculate radial density and add to appropriate location in array
        if (ibin > intervals(mpi_index,1) - 1 .and. ibin < intervals(mpi_index,2) + 1) then
          iradial = ibin - intervals(mpi_index,1) + 1
          if(i > MOD(radial_record(iradial, 1),wl_setup%mc_sweeps*setup%n_atoms) &
            .and. radial_record(iradial, 2) < wl_setup%radial_samples ) then
            radial_record(iradial, 1) = radial_record(iradial, 1) + INT(setup%n_atoms*0.05)
            radial_record(iradial, 2) = radial_record(iradial, 2) + 1
            rho_of_E(:,:,:,ibin) = rho_of_E(:,:,:,ibin) + radial_densities(setup, config, setup%wc_range, shells)
          end if
        end if

        ! Add change in V into diff_energy
        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (genrand() .lt. exp((wl_logdos(ibin) - wl_logdos(jbin)))) then
          acceptance = acceptance + 1
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
          jbin = ibin
        end if
        mpi_wl_hist(jbin - mpi_start_idx + 1) = mpi_wl_hist(jbin - mpi_start_idx + 1) + 1.0_real64
        wl_logdos(jbin) = wl_logdos(jbin) + wl_f
      else
        ! reject and reset
        call pair_swap(config, rdm1, rdm2)
      end if
    end do

  end function run_wl_sweeps

  subroutine wl_burn_in(setup, config, target_energy, min_e, max_e, &
                        num_proc, num_walkers, my_rank, mpi_index)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    real(real64), intent(in) :: target_energy, min_e, max_e
    integer, intent(in) :: num_proc, num_walkers, my_rank, mpi_index

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e
    integer(int16) :: site1, site2
    logical :: stop_burn = .False., flag
    integer :: rank, rank_index, request, ierr, status, con_sum

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)

    ! Non-blocking receive setup
    if (mpi_index /= num_proc) then
      call MPI_IRECV(stop_burn, 1, MPI_LOGICAL, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, request, ierr)
    end if
    do while (.True.)
      call MPI_TEST(request, flag, MPI_STATUS_IGNORE, ierr)
      if (stop_burn .eqv. .True. .and. mpi_index /= num_proc) then
        call MPI_RECV(config, SIZE(config), MPI_SHORT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        !print*, my_rank, setup%full_energy(config), min_e, max_e
        exit
      else if (e_unswapped > min_e .and. e_unswapped < max_e .and. mpi_index /= num_proc) then
        stop_burn = .True.
        call MPI_CANCEL(request, ierr)
        call MPI_REQUEST_FREE(request, ierr)
        do rank=0, num_proc-1
          rank_index = rank/num_walkers + 1 
          if (rank /= my_rank .and. rank_index == mpi_index) then
            !print*, mpi_index, my_rank, "Sending to", rank_index, rank
            call MPI_SEND(stop_burn, 1, MPI_LOGICAL, rank, 0, MPI_COMM_WORLD, ierr)
            call MPI_SEND(config, SIZE(config), MPI_SHORT, rank, 0, MPI_COMM_WORLD, ierr)
          end if
        end do
        exit
      else if (e_unswapped > min_e .and. e_unswapped < max_e .and. mpi_index == num_proc) then
        exit
      end if

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      if (site1 /= site2) then
        call pair_swap(config, rdm1, rdm2)
        e_swapped = setup%full_energy(config)

        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (e_swapped > target_energy .and. delta_e < 0) then
          e_unswapped = e_swapped
        else if (e_swapped < target_energy .and. delta_e > 0) then
          e_unswapped = e_swapped
        else if (genrand() .lt. 0.001_real64) then ! to prevent getting stuck in local minimum (should adjust this later to something more scientific instead of an arbitrary number)
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do
  end subroutine wl_burn_in

  subroutine divide_range(start, finish, num_intervals, intervals)
    integer, intent(in) :: num_intervals
    integer, intent(in) :: start, finish
    integer, intent(out) :: intervals(num_intervals, 2)
    real(real64) :: factor, index, power, b, n, g, idx, z
    integer :: i

    intervals(1,1) = start
    intervals(num_intervals,2) = finish

    power = 2
    b = finish
    n = num_intervals + 1
    g = 0.2
    idx = 3
    z = 2

    factor = (b-1.0_real64)/((n+1.0_real64)**power-1.0_real64)
    !factor = (1.0_real64 - b)/(EXP(g)-EXP(g*n))
    !factor = (idx - b)/(EXP(g*z)-EXP(g*n))

    do i = 2, num_intervals
      index = FLOOR(factor*(i**power-1)+1)
      !index = FLOOR(factor*EXP(g*i)+1-factor*EXP(g))
      !index = FLOOR(factor*EXP(g*i)+idx-factor*EXP(g*z))
      intervals(i-1,2) = INT(index)
      intervals(i,1) = INT(index + 1)
    end do
  end subroutine divide_range

end module wang_landau
