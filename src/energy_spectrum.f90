!----------------------------------------------------------------------!
! Energy Spectrum Module                                               !
!                                                                      !
! This program runs MC for the given system and attempts to find       !
! lowest energy state.                                                 !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module energy_spectrum
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

  contains

  subroutine es_main(setup, es_setup, my_rank)
    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
    type(es_params) :: es_setup

    ! Integers used in calculations
    integer :: ierr, unique_energy, i, energy_min_loc, energy_max_loc

    ! Temperature and temperature steps
    real(real64) :: acceptance, step, energy_to_ry, min_energy

    real(real64), allocatable :: energy_spectrum(:), energy_spectrum_condensed(:), bin_edges(:)
    integer, allocatable :: energy_spectrum_sort(:)

    energy_to_ry = setup%n_atoms/(eV_to_Ry*1000)

    allocate (energy_spectrum(es_setup%unique_energy_count))
    allocate (bin_edges(es_setup%bins + 1))
    energy_spectrum = 10000.0_real64
    unique_energy = 1

    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells, config)

    min_energy = setup%full_energy(config)

    if (my_rank == 0) then
      write (6, '(/,72("-"),/)')
      write (6, '(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
      print *, "Number of atoms", setup%n_atoms
    end if

    acceptance = run_es_sweeps(setup, es_setup, config, min_energy)

    print *, "Energy Reached", min_energy/energy_to_ry, "meV"

    call comms_wait()
    if (my_rank == 0) then
      write (*, *)
      write (6, '(25("-"),x,"Simulation Complete!",x,25("-"))')
    end if

  end subroutine es_main

  function run_es_sweeps(setup, es_setup, config, min_energy) result(acceptance)
    integer(int16), dimension(:, :, :, :) :: config
    class(run_params), intent(in) :: setup
    class(es_params), intent(in) :: es_setup
    real(real64), intent(inout) :: min_energy

    integer, dimension(4) :: rdm1, rdm2
    real(real64) :: e_swapped, e_unswapped, eps, delta_e
    integer :: acceptance, i, cycle_loop, reject
    integer(int16) :: site1, site2

    ! Establish total energy before any moves
    e_unswapped = setup%full_energy(config)
    e_swapped = e_unswapped

    acceptance = 0.0_real64

    do i = 1, es_setup%mc_sweeps*setup%n_atoms
      ! Make one MC trial
      ! Generate random numbers
      rdm1 = setup%rdm_site()
      rdm2 = setup%rdm_site()

      ! Get what is on those sites
      site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
      site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

      ! Calculate energy if different species
      if (site1 /= site2) then
        call pair_swap(config, rdm1, rdm2)
        e_swapped = setup%full_energy(config)

        delta_e = e_swapped - e_unswapped

        ! Accept or reject move
        if (delta_e < 0) then
          e_unswapped = e_swapped
          if (e_swapped < min_energy) then
            min_energy = e_swapped
          end if
        else if (genrand() .lt. 0.001_real64) then ! to prevent getting stuck in local minimum
          e_unswapped = e_swapped
        else
          call pair_swap(config, rdm1, rdm2)
        end if
      end if
    end do

  end function run_es_sweeps
end module energy_spectrum

