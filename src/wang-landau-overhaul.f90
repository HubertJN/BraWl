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

    ! Input file data types
    type(run_params) :: setup
    type(wl_params) :: wl_setup

    ! MPI variables
    integer :: ierror, mpi_process

    ! Convert input variables to simulation variables
    eV_per_atom_to_Ry = setup%n_atoms/(eV_to_Ry*1000) ! Conversion meV/atom to Rydberg
    setup%T = setup%T*k_b_in_Ry ! Convert temperature to simulation units
    wl_setup%energy_max = wl_setup%energy_max*eV_per_atom_to_Ry
    wl_setup%energy_min = wl_setup%energy_min*eV_per_atom_to_Ry

    ! Adjusting due to MPI processes saving data for the same bin
    wl_setup%radial_samples = wl_setup%radial_samples/num_walkers

    ! Get number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_process, ierror)

    ! Number of MPI processes per MPI window
    num_walkers = mpi_process/wl_setup%num_windows

    ! Obtain values for bin_edges array
    call energy_binning(wl_setup, bin_edges)

    ! Obtain non-uniform MPI window window_intervals
    call divide_range(wl_setup, window_intervals)

    ! Populate MPI arrays and indlude MPI window overlap
    call mpi_arrays(wl_setup, my_rank, window_intervals, &
                    window_indices, mpi_bins, mpi_index)

    ! Set up the lattice
    call initial_setup(setup, config)
    call lattice_shells(setup, shells, config)

    ! Burn-in MPI windows
    call burn_in(setup, wl_setup, config, my_rank, mpi_index, window_rank_index, bin_edges, &
                window_min_e, window_max_e, radial_record, rho_of_E, window_intervals)

    ! Perform initial pre-sampling
    do while(minval(mpi_wl_hist) < 10.0_real64)
      call sweep(setup, wl_setup, config, bin_edges, &
                mpi_start_idx, mpi_end_idx, mpi_wl_hist, wl_logdos, wl_f, &
                mpi_index, window_intervals, radial_record, rho_of_E)
    end do
    mpi_wl_hist = 0.0_real64

    ! Have one sweep between radial densities
    ! Program timing such that it excludes radial density calculation

    ! Main Wang-Landau sampling loop
    time_start = mpi_wtime()
    do while(condition)

      ! Perform sweeps
      call sweep(setup, wl_setup, config, bin_edges, wl_logdos, wl_f, &
                mpi_start_idx, mpi_end_idx, mpi_wl_hist, mpi_index, &
                window_intervals, radial_record, rho_of_E, &
                time_start)

      ! Calculate flatness
      flatness =  calc_flatness(mpi_wl_hist)

      if (flatness > wl_setup%flatness .and. minval(mpi_wl_hist) > 1) then
        time_end = mpi_wtime()
        rank_time = 0.0_real64
        rank_time(mpi_index) = time_end - time_start
        ! Average Density of States within MPI windows
        call dos_average(wl_setup, mpi_index, my_rank, num_walkers, &
                        window_rank_index, wl_logdos, wl_logdos_buffer)

        ! Combine Density of States across MPI windows
        call dos_combine(wl_setup, mpi_index, my_rank, num_walkers, &
                        mpi_start_idx, mpi_start_idx_buffer, mpi_end_idx, mpi_end_idx_buffer, &
                        window_indices, window_rank_index, &
                        wl_logdos, wl_logdos_buffer, wl_logdos_write)

        ! Save data with NetCDF 
        call save_data(my_rank, bin_edges, wl_logdos, wl_hist) ! Include radial density

        ! Adjust sampling parameters 
        call adjust_parameters(mpi_wl_hist, wl_f)

        ! MPI metadata and MPI window optimisation
        call mpi_metadata(...)
        call mpi_window_optimise(...)
        time_start = mpi_wtime()
      end if

    end do

  end subroutine wl_main

  ! Subroutine for finding the bin index
  integer function bin_index(energy, bin_edges, bins) result(index)
    integer, intent(in) :: bins
    real(real64), intent(in) :: energy
    real(real64), dimension(:), intent(in) :: bin_edges
    real(real64) :: bin_range

    bin_range = bin_edges(bins + 1) - bin_edges(1)
    index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  ! Subroutine that calculates energy range binning and populates bin_edges
  subroutine energy_binning(wl_setup, bin_edges, bin_energy)
    ! Subroutine input
    type(wl_params) :: wl_setup

    ! Subroutine output
    real(real64), dimension(:), intent(out) :: bin_edges, bin_energy

    ! Loop indices
    integer :: i

    ! Calculation variables
    real(real64) :: bin_width

    bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)
    do i = 1, wl_setup%bins + 1
      bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i - 1)*bin_width
    end do
    do i = 1, bins
      bin_energy(i) = wl_setup%energy_min*energy_to_ry + (i - 0.5)*bin_width
    end do

  end subroutine energy_binning

  ! Subroutine that generates non-uniform bin distribution for MPI windows
  ! and distributes available MPI processes
  subroutine divide_range(wl_setup, mpi_process, window_intervals, window_rank_index)
    ! Subroutine input
    type(wl_params) :: wl_setup
    integer, intent(in) :: mpi_process

    ! Subroutine output
    integer, dimension(:, :), intent(out) :: window_intervals, window_rank_index

    ! Loop indices
    integer :: i

    ! Calculation variables
    real(real64) :: a, b, power, factor
    ! Set first values outside of loop
    window_intervals(1,1) = 1
    window_intervals(wl_setup%num_windows,2) = wl_setup%bins

    ! values for equation to generate bin distribution
    ! factor derived from equation of form:
    ! y = Ax^B + C
    power = 2
    b = wl_setup%bins
    n = wl_setup%num_windows + 1 ! +1 since "edges" are being obtained for window_intervals array

    factor = (b-1.0_real64)/((n+1.0_real64)**power-1.0_real64)

    do i = 2, num_intervals
      window_intervals(i-1,2) = INT(FLOOR(factor*(i**power-1)+1))
      window_intervals(i,1) = INT(FLOOR(factor*(i**power-1)+1) + 1)
    end do

    ! Distribute MPI processes
    window_rank_index(1,1) = 0
    window_rank_index(wl_setup%num_windows, 2) = mpi_process
    do i = 2, wl_setup%num_windows
      window_rank_index(i,1) = INT(FLOOR(mpi_process/wl_setup%num_windows))*(i-1)
      window_rank_index(i-1,2) = window_rank_index(i,1) - 1
    end do
    do i = 1, MOD(mpi_process, wl_setup%num_windows)
      window_rank_index(i+1) = window_rank_index(i+1) + i
    end do
  end subroutine divide_range

  ! Subroutine to create overlapping indexing of MPI windows and associated variables
  subroutine mpi_arrays(wl_setup, my_rank, window_intervals, &
                        window_indices, mpi_bins, mpi_index)
    ! Subroutine input
    type(wl_params) :: wl_setup
    integer, intent(in) :: my_rank
    integer, dimension(:,:), intent(in) :: window_intervals

    ! Subroutine output
    integer, dimension(:,:), intent(out) :: window_indices
    integer, intent(out) :: mpi_bins, mpi_index

    ! Loop indices
    integer :: i

    window_indices(1, 1) = window_intervals(1,1)
    window_indices(1,2) = INT(window_intervals(1,2) + ABS(window_intervals(2,1)-window_intervals(2,2))*wl_setup%bin_overlap)
    do i = 2, wl_setup%num_windows-1
      window_indices(i, 1) = INT(window_intervals(i,1) - ABS(window_intervals(i-1,1)-window_intervals(i-1,2))*wl_setup%bin_overlap)
      window_indices(i, 2) = INT(window_intervals(i,2) + ABS(window_intervals(i+1,1)-window_intervals(i+1,2))*wl_setup%bin_overlap)
    end do
    window_indices(wl_setup%num_windows, 1) = INT(window_intervals(wl_setup%num_windows,1) - ABS(window_intervals(wl_setup%num_windows-1,1) &
                                    -window_intervals(wl_setup%num_windows-1,2))*wl_setup%bin_overlap)
    window_indices(wl_setup%num_windows,2) = window_intervals(wl_setup%num_windows,2)

    mpi_index = my_rank/wl_setup%num_windows + 1
    mpi_bins = window_indices(mpi_index, 2) - window_indices(mpi_index, 1) + 1

  end subroutine mpi_arrays

  ! Subroutine to burn-in MPI windows
  subroutine burn_in(setup, wl_setup, config, my_rank, mpi_index, window_rank_index, bin_edges, &
                     window_min_e, window_max_e, radial_record, rho_of_E, window_intervals)

    ! Subroutine input
    type(run_params) :: setup
    type(wl_params) :: wl_setup
    integer(int16), dimension(:, :, :, :) :: config
    integer, intent(in) :: my_rank, mpi_index, 
    integer, dimension(:,:), intent(in) :: window_intervals
    real(real64) :: window_min_e, window_max_e

    ! Subroutine input-output
    integer, dimension(:,:), intent(inout) :: radial_record
    integer, dimension(:,:,:,:), intent(inout) :: rho_of_E

    ! Subroutine internal
    integer :: bin, rank, request, ierror, status
    logical :: stop_burn_in, flag

    ! Init
    stop_burn_in = .False.

    ! Target energy
    target_energy = (window_min_e + window_max_e)/2.0_real64

    ! Initial configuration energy
    e_unswapped = congif%full_energy(config)

    ! Non-blocking MPI receive
    call MPI_IRECV(stop_burn_in, 1, MPI_LOGICAL, MPI_ANY_SOURCE, mpi_index, MPI_COMM_WORLD, request, ierror)

    do while(.True.)
      ! Check if MPI message received
      call MPI_TEST(request, flag, MPI_STATUS_IGNORE, ierror)

      ! Stop burn if other rank in window is burnt in
      ! or if burnt in send configuration to rest of window
      if (stop_burn_in .eqv. .True.) then
        call MPI_RECV(config, SIZE(config), MPI_SHORT, MPI_ANY_SOURCE, mpi_index, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror) ! Check if you can put an array of accepted values in the source variable
        exit
      else if (e_unswapped > window_min_e .and. e_unswapped < window_max_e) then
        stop_burn_in = .True.
        call MPI_CANCEL(request, ierror)
        call MPI_REQUEST_FREE(request, ierror)
        do rank=window_rank_index(mpi_index, 1), window_rank_index(mpi_index, 2)
          call MPI_SEND(stop_burn, 1, MPI_LOGICAL, rank, mpi_index, MPI_COMM_WORLD, ierror)
          call MPI_SEND(config, SIZE(config), MPI_SHORT, rank, mpi_index, MPI_COMM_WORLD, ierror)
        end do
        exit
      end if

      bin = bin_index(e_unswapped, bin_edges, wl_setup%bins)
      call radial_density_record(setup, config, radial_record, rho_of_E, &
                                intervals, mpi_index, bin)

      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      ! Only caclulate proceed for sites with different elements
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
        else    ! Loop indices    ! Loop indices
          call pair_swap(config, rdm1, rdm2)
        end if    ! Loop indices
      end if
    end do

  end subroutine burn_in

  subroutine radial_density_record(setup, config, radial_record, rho_of_E, &
                                   intervals, mpi_index, bin)
    ! Subroutine input
    type(run_params) :: setup
    integer(int16), dimension(:, :, :, :) :: config
    integer, intent(in) :: mpi_index, bin
    integer, dimension(:,:), intent(in) :: intervals

    ! Subroutine input-output
    integer, dimension(:,:), intent(inout) :: radial_record
    integer, dimension(:,:,:,:), intent(inout) :: rho_of_E

    radial_record(iradial, 1) = radial_record(iradial, 1) + 1

    if (bin > intervals(mpi_index,1) - 1 .and. bin < intervals(mpi_index,2) + 1) then
      iradial = bin - intervals(mpi_index,1) + 1
      if((radial_record(iradial, 1) >= setup%n_atoms) &
        .and. radial_record(iradial, 2) < wl_setup%radial_samples ) then
        radial_record(iradial, 1) = 0
        radial_record(iradial, 2) = radial_record(iradial, 2) + 1
        rho_of_E(:,:,:,bin) = rho_of_E(:,:,:,bin) + radial_densities(setup, config, setup%wc_range, shells)
      end if
    end if
  end subroutine radial_density_record

  subroutine sweep(setup, wl_setup, config, bin_edges, wl_logdos, wl_f, &
                  mpi_start_idx, mpi_end_idx, mpi_wl_hist, mpi_index, &
                  window_intervals, radial_record, rho_of_E, &
                  time_start)
    
    ! Subroutine input
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: mpi_start_idx, mpi_end_idx, mpi_index
    real(real64), intent(in) :: wl_f
    real(real64), dimension(:), intent(in) :: bin_edges

    ! Subroutine input-output
    integer, dimension(:,:), intent(inout) :: radial_record, intervals
    real(real64), dimension(:), intent(inout) :: mpi_wl_hist, wl_logdos
    real(real64), dimension(:,:,:,:), intent(inout) :: rho_of_E
    real(real64), intent(inout) :: time_start

    ! Subroutine internal
    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, delta_e, radial_start, radial_end
    integer :: i, ibin, jbin
    integer(int16) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped

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

      ibin = bin_index(e_unswapped, bin_edges, wl_setup%bins)
      jbin = bin_index(e_swapped, bin_edges, wl_setup%bins)

      radial_start = mpi_wtime()
      call radial_density_record(setup, config, radial_record, rho_of_E, &
                                intervals, mpi_index, jbin)
      radial_end = mpi_wtime()
      time_start = time_start - (radial_end - radial_start)

      ! Only compute energy change if within limits where V is defined and within MPI region
      if (jbin > mpi_start_idx - 1 .and. jbin < mpi_end_idx + 1) then

        ! Add change in V into diff_energy
        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (genrand() .lt. exp((wl_logdos(ibin) - wl_logdos(jbin)))) then
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
  end subroutine sweep

  function calc_flatness(mpi_wl_hist) result(flatness)
    ! Function input
    real(real64), dimension(:), intent(in) :: mpi_wl_hist
    ! Function internal
    real(real64) :: flatness

    flatness = minval(mpi_wl_hist)/(sum(mpi_wl_hist)/mpi_bins)
  end function calc_flatness

  subroutine average_dos(wl_setup, mpi_index, my_rank, num_walkers, &
                        window_rank_index, wl_logdos, wl_logdos_buffer)

    ! Input
    type(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: mpi_index, my_rank, num_walkers
    integer, dimension(:,:), intent(in) :: window_rank_index

    ! Input-Output
    real(real64), dimension(:), intent(inout) :: wl_logdos, wl_logdos_buffer

    ! Internal
    integer :: i, j, ierror

    do i = 1, wl_setup%num_windows
      if (mpi_index == i) then
        if (my_rank /= window_rank_index(i, 1)) then
          call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, ierror)
        else
          do j = 1, num_walkers - 1
            call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                          MPI_STATUS_IGNORE, ierror)
            wl_logdos = wl_logdos + wl_logdos_buffer
          end do
        end if
        if (my_rank == window_rank_index(i, 1)) then
          wl_logdos = wl_logdos - minval(num_walkers)
          do j = 1, num_walkers - 1
            call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers + j, i, MPI_COMM_WORLD, &
                          ierror)
          end do
        else
          call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, i, MPI_COMM_WORLD, &
                        MPI_STATUS_IGNORE, ierror)
        end if
      end if
    end do
  end subroutine average_dos

  subroutine combine_dos(wl_setup, mpi_index, my_rank, num_walkers, &
                        mpi_start_idx, mpi_start_idx_buffer, mpi_end_idx, mpi_end_idx_buffer, &
                        window_indices, window_rank_index, &
                        wl_logdos, wl_logdos_buffer, wl_logdos_write)

    ! Input
    type(wl_params), intent(in) :: wl_setup
    integer, intent(in) :: mpi_index, my_rank, num_walkers
    integer, intent(in) :: mpi_start_idx, mpi_start_idx_buffer, mpi_end_idx, mpi_end_idx_buffer
    integer, dimension(:,:), intent(in) :: window_indices, window_rank_index

    ! Input-Output
    real(real64), dimension(:), intent(inout) :: wl_logdos, wl_logdos_buffer, wl_logdos_write

    ! Internal
    integer :: i, j, ierror
    real(real64) :: beta_index, beta_original, beta_merge, beta_diff

    do i = 2, wl_setup%num_windows
      if (my_rank == window_rank_index(i,1)) then
        call MPI_Send(wl_logdos, wl_setup%bins, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
        call MPI_Send(mpi_start_idx, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
        call MPI_Send(mpi_end_idx, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
      end if
      if (my_rank == 0) then
        call MPI_Recv(wl_logdos_buffer, wl_setup%bins, MPI_DOUBLE_PRECISION, (i - 1)*num_walkers, 0, MPI_COMM_WORLD, &
                      MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_start_idx_buffer, 1, MPI_INT, (i - 1)*num_walkers, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        call MPI_Recv(mpi_end_idx_buffer, 1, MPI_INT, (i - 1)*num_walkers, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        scale_factor = 0.0_real64
        beta_diff = wl_logdos_write(2) - wl_logdos_write(1)
        do j = 0,  window_indices(i-1, 2) - window_indices(i, 1) - 1
          beta_original = wl_logdos_write(mpi_start_idx_buffer + j + 1) - wl_logdos_write(mpi_start_idx_buffer + j)
          beta_merge = wl_logdos_buffer(mpi_start_idx_buffer + j + 1) - wl_logdos_buffer(mpi_start_idx_buffer + j)
          if (ABS(beta_original - beta_merge) < beta_diff) then
            beta_diff = ABS(beta_original - beta_merge)
            beta_index = mpi_start_idx_buffer + j
          end if
        end do
        do j = beta_index, mpi_end_idx_buffer
          wl_logdos_write(j) = wl_logdos_buffer(j) + wl_logdos_write(beta_index) - wl_logdos_buffer(beta_index)
        end do
      end if
    end do
  end subroutine combine_dos

  subroutine save_data(setup, my_rank, bin_edges, wl_logdos, wl_hist, bin_energy, num_walkers &
                      rho_saved, radial_min, rho_of_E, rho_of_E_buffer)

    ! Input
    type(run_params), intent(in) :: setup
    integer, intent(in) :: my_rank, num_walkers
    integer, dimension(:) :: wl_hist
    real(real64) :: radial_min
    real(real64), dimension(:), intent(in) :: bin_edges, bin_energy, wl_logdos
    real(real64), dimension(:,:,:,:), intent(in) :: rho_of_E, rho_of_E_buffer

    ! Input-Output
    logical, intent(inout) :: rho_saved

    if (rho_saved .eqv. .False. .and. INT(radial_min) == 100) then
          
      rho_saved = .True.
      call MPI_REDUCE(rho_of_E, rho_of_E_buffer, setup%n_species*setup%n_species*setup%wc_range*wl_setup%bins, &
                      MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      if (my_rank == 0) then
        !========================================================
        print*, "Radial densities saved" ! Make this more visible
        !========================================================
        rho_of_E_buffer = rho_of_E_buffer/wl_setup%radial_samples/num_walkers
        call ncdf_radial_density_writer_across_energy("radial_densities/rho_of_E.dat", rho_of_E_buffer, shells, bin_energy, setup)
      end if
    end if

    if (my_rank == 0) then
      call ncdf_writer_1d("wl_dos_bins.dat", ierr, bin_edges)
      call ncdf_writer_1d("wl_dos.dat", ierr, wl_logdos_write)
    end if
  end subroutine save_data

  subroutine adjust_parameters(mpi_wl_hist, wl_f)
    real(real64), dimension(:), intent(inout) :: mpi_wl_hist
    real(real64), intent(inout) :: wl_f

    mpi_wl_hist = 0.0_real64
    wl_f = wl_f*1/2

  end subroutine adjust_parameters

  subroutine mpi_metadata(my_rank, mpi_index, num_walkers, rank_time, print_time)
    time_total = rank_time(mpi_index)

    call MPI_REDUCE(rank_time/num_walkers, print_time(:,1), num_windows, MPI_DOUBLE_PRECISION, &
    MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    rank_time = HUGE(time_total)
    rank_time(mpi_index) = time_total

    call MPI_REDUCE(rank_time(mpi_index), print_time(:,2), num_windows, MPI_DOUBLE_PRECISION, &
    MPI_MIN, 0, MPI_COMM_WORLD, ierror)

    rank_time = 0.0_real64
    rank_time(mpi_index) = time_total

    call MPI_REDUCE(rank_time(mpi_index), print_time(:,3), num_windows, MPI_DOUBLE_PRECISION, &
    MPI_MAX, 0, MPI_COMM_WORLD, ierror)

    if (my_rank == 0) then
      do i=1, num_windows
        write (6, '(a,i0,a,f8.2,a,f8.2,a,f8.2,a)') "MPI Window: ", i, " | Avg. time: ", print_time(i,1), &
        "s | Time min: ", print_time(i,2), "s Time max: " , print_time(i,3), "s"
      end do
      write (*, *)
    end if

  end subroutine mpi_metadata