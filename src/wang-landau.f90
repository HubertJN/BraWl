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

      ! Arrays for storing data
      type(run_params) :: setup
      type(wl_params) :: wl_setup
    
      ! Integers used in calculations
      integer :: i,j,k, ierr, accept, bins
      
      ! Temperature and temperature steps
      real(real64) :: temp, acceptance, beta, flatness
    
      ! wl variables and arrays
      real(real64) :: bin_width, energy_to_ry, wl_f, wl_f_prev, tolerance, flatness_tolerance
      real(real64) :: target_energy, min_val

      real(real64), allocatable :: bin_edges(:), wl_hist(:), wl_logdos(:)
      logical :: first_reset, hist_reset

      ! mpi variables
      real(real64) :: start, end, time_max, time_min
      integer :: mpi_bins, mpi_start_idx, mpi_end_idx, bin_overlap
      real(real64) :: mpi_width, scale_factor, scale_count
      real(real64), allocatable :: mpi_bin_edges(:), mpi_wl_hist(:), wl_logdos_buffer(:)

      ! Minimum value in array to be considered
      min_val = wl_setup%tolerance*1e-1_real64

      ! Get number of MPI processes
      call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierror)

      ! Allocate arrays
      bins = wl_setup%bins
      mpi_bins = bins/num_proc
      if (my_rank == num_proc-1) then
          mpi_bins = mpi_bins + MOD(bins, mpi_bins*(num_proc))
      end if

      bin_overlap = wl_setup%bin_overlap
      mpi_start_idx = max(my_rank*(bins/num_proc) + 1 - bin_overlap, 1)
      mpi_end_idx = min(my_rank*(bins/num_proc) + 1 + mpi_bins + bin_overlap, bins)
      mpi_bins = mpi_end_idx - mpi_start_idx + 1

      ! Adjust for spare bins after division
      if (my_rank < bins - (bins/num_proc)*num_proc) then
        mpi_start_idx = mpi_start_idx + my_rank
        mpi_end_idx = mpi_end_idx + my_rank
      end if

      allocate(bin_edges(bins+1))
      allocate(wl_hist(bins))
      allocate(wl_logdos(bins))
      if (my_rank == 0) then
        allocate(wl_logdos_buffer(bins))
      end if

      ! MPI arrays
      allocate(mpi_bin_edges(mpi_bins+1))
      allocate(mpi_wl_hist(mpi_bins))

      ! Set temperature
      temp = setup%T*k_b_in_Ry

      ! Conversion meV/atom to Rydberg
      energy_to_ry=setup%n_atoms/(eV_to_Ry*1000)

      mpi_width = (wl_setup%energy_max - wl_setup%energy_min)/real(num_proc)*energy_to_ry
      target_energy = wl_setup%energy_min*energy_to_ry + mpi_width*(my_rank + 0.5)

      j = 1
      bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)*energy_to_ry
      do i=1, bins+1
        bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i-1)*bin_width
        if (i > mpi_start_idx-1 .and. i < mpi_end_idx+2) then
            mpi_bin_edges(j) = bin_edges(i)
            j = j + 1
        end if
      end do

      wl_hist = 0.0_real64; wl_logdos = 0.0_real64; wl_f = wl_setup%wl_f; wl_f_prev = wl_f
      flatness = 0.0_real64; first_reset = .False.; hist_reset = .True.

      mpi_wl_hist = 0.0_real64
      !---------------------------------!

      ! Set up the lattice
      call initial_setup(setup, config)
  
      call lattice_shells(setup, shells, config)
    
      ! Are we swapping neighbours or on the whole lattice?
      if (setup%nbr_swap) then
        setup%mc_step => monte_carlo_step_nbr
      else
        setup%mc_step => monte_carlo_step_lattice
      end if

      if(my_rank == 0) then
        write(6,'(/,72("-"),/)')
        write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
        print*, "Number of atoms", setup%n_atoms
      end if
      call comms_wait()

      !---------!
      ! Burn in !
      !---------!

      beta = 1.0_real64/temp
      if (wl_setup%burn_in) then
          call wl_burn_in(setup, wl_setup, config, target_energy, MINVAL(mpi_bin_edges), MAXVAL(mpi_bin_edges))
      end if

      if(my_rank == 0) then
        write(*,*)
        write(6,'(27("-"),x,"Burn-in complete",x,27("-"),/)')
        write(*,*)
      end if

      call comms_wait()
      
      !--------------------!
      ! Target Temperature !
      !--------------------!
      i = 0
      j = 0
      tolerance = wl_setup%tolerance
      flatness_tolerance = wl_setup%flatness
      do while (wl_f > tolerance)
        if (hist_reset .eqv. .True.) then
          ! Start timer
          start = mpi_wtime()
          hist_reset = .False.
        end if
        acceptance = run_wl_sweeps(setup, wl_setup, config, bins, bin_edges, &
        mpi_bin_edges, mpi_bins, mpi_wl_hist, wl_logdos, wl_f)

        flatness = 100.0_real64*minval(mpi_wl_hist, MASK=(mpi_wl_hist > min_val))/(sum(mpi_wl_hist, &
        MASK=(mpi_wl_hist > min_val))/count(MASK=(mpi_wl_hist > min_val)))

        if (first_reset .eqv. .False. .and. minval(mpi_wl_hist, MASK=(mpi_wl_hist > min_val)) > 100) then ! First reset after system had time to explore
          first_reset = .True.
          mpi_wl_hist = 0.0_real64
        else
          !Check if we're flatness_tolerance% flat
          if (flatness > flatness_tolerance) then
            !Reset the histogram
            wl_f = wl_f * 0.5_real64
            mpi_wl_hist = 0.0_real64
            hist_reset = .True.
            !Reduce f
            wl_logdos = wl_logdos - minval(wl_logdos, MASK=(wl_logdos > min_val))
            wl_logdos = wl_logdos * merge(0, 1, wl_logdos < 0.0_real64)
          end if
        end if

        if (hist_reset .eqv. .True.) then
          ! End timer
          end = mpi_wtime()

          call MPI_REDUCE(end-start, time_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierror)
          call MPI_REDUCE(end-start, time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
          call comms_wait()

          if(my_rank == 0) then
              write(6,'(a,i0,a,f6.2,a,f12.10,a,f6.2,a,f6.2,a)', advance='no'), "Rank: ", my_rank, " Flatness: ", flatness, &
              "% W-L F: ", wl_f_prev, " Time min:", time_min, "s Time max:", time_max, "s"
              wl_f_prev = wl_f
              write(*,*)
          end if
        end if

        do i=1, num_proc-1 
          if (my_rank == i) then
            call MPI_Send(wl_logdos, bins, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
            call MPI_Send(mpi_start_idx, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, ierror)
            call MPI_Send(mpi_end_idx, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, ierror)
          end if
          if (my_rank == 0) then
            call MPI_Recv(wl_logdos_buffer, bins, MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_Recv(mpi_start_idx, 1, MPI_INT, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            call MPI_Recv(mpi_end_idx, 1, MPI_INT, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
            scale_factor = 0.0_real64
            scale_count = 0.0_real64
            do j=0, bin_overlap-1
              if (wl_logdos_buffer(mpi_start_idx+j) > min_val .and. wl_logdos(mpi_start_idx+j) > min_val) then
                scale_factor = scale_factor + EXP(wl_logdos(mpi_start_idx+j))-EXP(wl_logdos_buffer(mpi_start_idx+j))
                scale_count = scale_count + 1.0_real64
              end if
            end do
            scale_factor = scale_factor/scale_count
            do j=mpi_start_idx+bin_overlap, mpi_end_idx
              wl_logdos(j) = LOG(EXP(wl_logdos_buffer(j)) + scale_factor)
            end do
          end if
        end do
  
        if (my_rank == 0) then
          ! Write output files
          call ncdf_writer_1d("wl_dos_bins.dat", ierr, bin_edges)
          call ncdf_writer_1d("wl_dos.dat", ierr, wl_logdos)
          call ncdf_writer_1d("wl_hist.dat", ierr, wl_hist)
        end if
        
      end do
     

      
      if(my_rank == 0) then
          write(*, *)
          write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
      end if
  
  end subroutine wl_main    

  integer function bin_index(energy, bin_edges, bins) result(index)
      integer, intent(in) :: bins
      real(real64), intent(in) :: energy
      real(real64), dimension(:), intent(in) :: bin_edges
      real(real64) :: bin_range

      bin_range = bin_edges(bins+1) - bin_edges(1)
      index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  function run_wl_sweeps(setup, wl_setup, config, bins, bin_edges, &
    mpi_bin_edges, mpi_bins, mpi_wl_hist, wl_logdos, wl_f) result(acceptance)
      integer(int16), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      class(wl_params), intent(in) :: wl_setup
      integer, intent(in) :: bins, mpi_bins
      real(real64), dimension(:), intent(in) :: bin_edges, mpi_bin_edges
      real(real64), dimension(:), intent(inout) :: mpi_wl_hist, wl_logdos
      real(real64), intent(in) :: wl_f

      integer, dimension(4) :: rdm1, rdm2
      real(real64) :: e_swapped, e_unswapped, delta_e
      integer :: acceptance, i, ibin, jbin, mpi_ibin, mpi_jbin
      integer(int16) :: site1, site2

      ! Establish total energy before any moves
      e_unswapped = setup%full_energy(config)
      e_swapped = e_unswapped

      acceptance = 0.0_real64

      do i=1, wl_setup%mc_sweeps*setup%n_atoms
  
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

          mpi_ibin = bin_index(e_unswapped, mpi_bin_edges, mpi_bins)
          mpi_jbin = bin_index(e_swapped, mpi_bin_edges, mpi_bins)

          ! Only compute energy change if within limits where V is defined and within MPI region
          if (jbin > 0 .and. jbin < bins+1 .and. mpi_jbin > 0 .and. mpi_jbin < mpi_bins+1) then
                ! Add change in V into diff_energy
                delta_e = e_swapped - e_unswapped

                ! Accept or reject move
                if (genrand() .lt. exp(wl_logdos(ibin)-wl_logdos(jbin))) then
                    acceptance = acceptance + 1
                    e_unswapped = e_swapped
                else
                    call pair_swap(config, rdm1, rdm2)
                    jbin = ibin
                end if
              mpi_wl_hist(mpi_ibin) = mpi_wl_hist(mpi_ibin) + 1.0_real64
              wl_logdos(jbin) = wl_logdos(jbin) + wl_f
          else
              ! reject and reset
              call pair_swap(config, rdm1, rdm2)
          end if
      end do

  end function run_wl_sweeps

  subroutine wl_burn_in(setup, wl_setup, config, target_energy, min_e, max_e)
      integer(int16), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      class(wl_params), intent(in) :: wl_setup
      real(real64), intent(in) :: target_energy, min_e, max_e

      integer, dimension(4) :: rdm1, rdm2
      real(real64) :: e_swapped, e_unswapped, delta_e, beta
      integer :: i

      ! Establish total energy before any moves
      e_unswapped = setup%full_energy(config)

      do while (.True.)
        if (e_unswapped > min_e .and. e_unswapped < max_e) then
          exit
        end if
          ! Make one MC trial
          ! Generate random numbers
          rdm1 = setup%rdm_site()
          rdm2 = setup%rdm_site()

          call pair_swap(config, rdm1, rdm2)
  
          e_swapped = setup%full_energy(config)

          delta_e = e_swapped - e_unswapped

            ! Accept or reject move
          if (e_swapped > target_energy .and. delta_e < 0) then
            e_unswapped = e_swapped
        else if (e_swapped < target_energy .and. delta_e > 0) then
            e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
    end do
  end subroutine wl_burn_in

end module wang_landau
