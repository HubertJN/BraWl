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
    setup%T = setup%T*k_b_in_Ry
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
    call burn_in(...) ! Include radial density saving

    ! Perform initial pre-sampling
    call pre_sampling(...) ! Include radial density saving

    ! Have one sweep between radial densities
    ! Program timing such that it excludes radial density calculation

    ! Main Wang-Landau sampling loop
    do while(condition)

      ! Perform sweeps
      call sweeps(...)

      ! Calculate metrics from sweep
      call sweep_metrics(...)

      if (flatness_condition) then
        ! Average Density of States within MPI windows
        call dos_average(...)

        ! Combine Density of States across MPI windows
        call dos_combine(...)

        ! Save data with NetCDF 
        call save_data(...) ! Include radial density

        ! Adjust sampling parameters 
        call adjust_paramters(...)

        ! MPI metadata and MPI window optimisation
        call mpi_metadata(...)
        call mpi_window_optimise(...)
      end if

    end do

  end subroutine wl_main

  ! Subroutine that calculates energy range binning and populates bin_edges
  subroutine energy_binning(wl_setup, bin_edges)
    ! Subroutine input
    type(wl_params) :: wl_setup

    ! Subroutine output
    real(real64), dimension(:), intent(out) :: bin_edges

    ! Loop indices
    integer :: i

    ! Calculation variables
    real(real64) :: bin_width

    bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)
    do i = 1, wl_setup%bins + 1
      bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i - 1)*bin_width
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
  subroutine burn_in(setup, wl_setup, config, my_rank, mpi_index, window_rank_index, &
                     window_min_e, window_max_e)
    ! Subroutine input

    ! Target energy
    target_energy = (window_min_e + window_max_e)/2.0_real64

    ! Initiale configuration energy
    e_unswapped = congif%full_energy(config)

    ! Non-blocking MPI receive
    call MPI_IRECV(stop_burn_in, 1, MPI_LOGICAL, MPI_ANY_SOURCE, mpi_index, MPI_COMM_WORLD, request, ierr)

    do while(.True.)
      ! Check if MPI message received
      call MPI_TEST(request, flag, MPI_STATUS_IGNORE, ierr)

      ! Stop burn if other rank in window is burnt in
      ! or if burnt in send configuration to rest of window
      if (stop_burn_in .eqv. .True.) then
        call MPI_RECV(config, SIZE(config), MPI_SHORT, MPI_ANY_SOURCE, mpi_index, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) ! Check if you can put an array of accepted values in the source variable
        exit
      else if (e_unswapped > window_min_e .and. e_unswapped < window_max_e) then
        stop_burn_in = .True.
        call MPI_CANCEL(request, ierr)
        call MPI_REQUEST_FREE(request, ierr)
        do rank=window_rank_index(mpi_index, 1), window_rank_index(mpi_index, 2)
          call MPI_SEND(stop_burn, 1, MPI_LOGICAL, rank, mpi_index, MPI_COMM_WORLD, ierr)
          call MPI_SEND(config, SIZE(config), MPI_SHORT, rank, mpi_index, MPI_COMM_WORLD, ierr)
        end do
      end if

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
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do

  end subroutine burn_in