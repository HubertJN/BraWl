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
    integer :: ierror, num_proc

    ! Convert input variables to simulation variables
    eV_per_atom_to_Ry = setup%n_atoms/(eV_to_Ry*1000) ! Conversion meV/atom to Rydberg
    setup%T = setup%T*k_b_in_Ry
    wl_setup%energy_max = wl_setup%energy_max*eV_per_atom_to_Ry
    wl_setup%energy_min = wl_setup%energy_min*eV_per_atom_to_Ry

    ! Adjusting due to MPI processes saving data for the same bin
    wl_setup%radial_samples = wl_setup%radial_samples/num_walkers

    ! Get number of MPI processes
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierror)

    ! Number of MPI processes per MPI window
    num_walkers = num_proc/wl_setup%num_windows

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
  subroutine divide_range(wl_setup, window_intervals)
    ! Subroutine input
    type(wl_params) :: wl_setup

    ! Subroutine output
    integer, dimension(:, :), intent(out) :: window_intervals

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