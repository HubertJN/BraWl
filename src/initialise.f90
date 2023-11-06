!----------------------------------------------------!
! File with module containing functions and routines !
! to call                                            !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

module initialise

  use kinds
  use mpi_shared_data
  use c_functions
  use energetics
  use random_site
  use metropolis
  
  implicit none

  contains

  !----------------------------------------------------!
  ! Subroutine to allocate memory for rank 1 processor !
  !----------------------------------------------------!
  subroutine initialise_global_arrays(setup)
    type(run_params), intent(inout) :: setup
    ! Array for storing energy as a function of temperature
    allocate(av_energies_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(av_C_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(av_acceptance_of_T(setup%T_steps))
    ! Radial densities as a function of temperature
    allocate(av_rho_of_T(setup%n_species, setup%n_species, &
                     setup%radial_d, setup%T_steps))
    av_rho_of_T = 0.0_real64
  end subroutine initialise_global_arrays

  !---------------------------------------!
  ! Subroutine to set up function pointer !
  !---------------------------------------!
  subroutine initialise_function_pointers(setup)
    type(run_params), intent(inout) :: setup

    setup%full_energy => total_energy

    ! Setup functions used in calculations
    if(trim(setup%lattice) == 'simple_cubic') then
      if (setup%n_shells .eq. 1) then
        setup%nbr_energy => simple_cubic_energy_1shells
      else
        print*, 'Unsupported number of shells'
        stop
      end if
      setup%rdm_site => simple_cubic_random_site
      setup%rdm_nbr => simple_cubic_random_nbr
    else if(trim(setup%lattice) == 'bcc') then
      setup%lattice_vectors = reshape( (/-0.5, 0.5, 0.5,     &
                                          0.5,-0.5, 0.5,     &
                                          0.5, 0.5,-0.5  /), &
                                       (/ 3, 3 /))
      setup%basis_vectors   = (/ 0.0, 0.0, 0.0 /)
      if (setup%n_shells .eq. 1) then
        setup%nbr_energy => bcc_energy_1shells
      else if (setup%n_shells .eq. 2) then
        setup%nbr_energy => bcc_energy_2shells
      else if (setup%n_shells .eq. 3) then
        setup%nbr_energy => bcc_energy_3shells
      else if (setup%n_shells .eq. 4) then
        setup%nbr_energy => bcc_energy_4shells
      else if (setup%n_shells .eq. 5) then
        setup%nbr_energy => bcc_energy_5shells
      else
        print*, 'Unsupported number of shells'
        stop
      end if
      setup%rdm_site => bcc_random_site
      setup%rdm_nbr => bcc_random_nbr
!    else if(trim(setup%lattice) == 'fcc') then
!      if (setup%n_shells .eq. 1) then
!        setup%nbr_energy => fcc_energy_1shells
!      else if (setup%n_shells .eq. 2) then
!        setup%nbr_energy => fcc_energy_2shells
!      else if (setup%n_shells .eq. 3) then
!        setup%nbr_energy => fcc_energy_3shells
!      else if (setup%n_shells .eq. 4) then
!        setup%nbr_energy => fcc_energy_4shells
!      else if (setup%n_shells .eq. 5) then
!        setup%nbr_energy => fcc_energy_5shells
!      else if (setup%n_shells .eq. 6) then
!        setup%nbr_energy => fcc_energy_6shells
!      else
!        print*, 'Unsupported number of shells'
!        stop
!      end if
!      setup%rdm_site => fcc_random_site
!      setup%rdm_nbr => fcc_random_nbr
    ! Bomb if we ask for a lattice type we don't have!
    else
      print*, 'Lattice type not yet implemented!'
      stop
    end if
    ! Are we swapping neighbours or on the whole lattice?
    if (setup%nbr_swap) then
      setup%mc_step => monte_carlo_step_nbr
    else
      setup%mc_step => monte_carlo_step_lattice
    end if
  end subroutine initialise_function_pointers

  !----------------------------------------------------!
  ! Subroutine to allocate memory for rank 1 processor !
  !----------------------------------------------------!
  subroutine initialise_local_arrays(setup)
    type(run_params), intent(inout) :: setup
    ! Array for storing energy as a function of temperature
    allocate(energies_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(C_of_T(setup%T_steps))
    ! Array for storing energy as a function of temperature
    allocate(acceptance_of_T(setup%T_steps))
    ! Radial densities as a function of temperature
    allocate(rho_of_T(setup%n_species, setup%n_species, &
                     setup%radial_d, setup%T_steps))
    rho_of_T = 0.0_real64
    ! On-shell distances
    allocate(shells(setup%radial_d))
    shells = 0.0_real64
    ! Array for storing temperatures
    allocate(temperature(setup%T_steps))
    ! Allocate space for exchange parameters
    allocate(V_ex(setup%n_species, setup%n_species, setup%n_shells))
    ! Allocate array for storing configuration
    allocate(config(setup%n_basis, setup%n_1, setup%n_2, setup%n_3))
  end subroutine initialise_local_arrays

  !-----------------------------------!
  ! Subroutine to clean up at the end !
  !-----------------------------------!
  subroutine local_clean_up()
    deallocate(V_ex)
    deallocate(rho_of_T)
    deallocate(shells)
    deallocate(temperature)
    deallocate(config)
    deallocate(energies_of_T, C_of_T, acceptance_of_T)
  end subroutine local_clean_up

  !-----------------------------------!
  ! Subroutine to clean up at the end !
  !-----------------------------------!
  subroutine global_clean_up()
    deallocate(av_rho_of_T)
    deallocate(av_energies_of_T, av_C_of_T, av_acceptance_of_T)
  end subroutine global_clean_up

  !---------------------------------------------!
  ! Function to setup lattice. We put zeroes on !
  ! every 'site' not part of the lattice.       !
  !---------------------------------------------!
  subroutine initial_setup(setup, config)
    type(run_params) :: setup
    integer(int16), allocatable, dimension(:,:,:,:) :: config
    integer(int32), dimension(4) :: grid_dims
    integer :: i, j, k, l, n_sites
    integer(int16) :: m, n_species
    real(real64) :: rand
    integer(int32), dimension(setup%n_species) :: species_count, check

    check = 0

    n_sites = setup%n_1*setup%n_2*setup%n_3*setup%n_basis

    ! Make sure my arrays have been allocated and bomb if not
    if (.not. allocated(config)) then
      stop 'config not allocated in function initial_setup'
    end if

    ! Get the dimensions of the grid I am using
    grid_dims = shape(config)
    n_species = int(setup%n_species, kind=int16)

    do i=1, n_species-1
      species_count(i) = floor(real(n_sites)*setup%species_cs(i))
    end do
    species_count(n_species) = n_sites - sum(species_count(1:(n_species-1)))

    ! Setup the concentrations array
    call set_concentrations(setup)

    ! Set configuration to be zero
    config = 0_int32

    ! Set up the lattice
    ! Loop over lattice sites
    do l=1, grid_dims(4)
      do k=1, grid_dims(3)
        do j=1, grid_dims(2)
          do i=1, grid_dims(1)
            do while (config(i,j,k,l) .eq. 0_int32)
              ! Get a random number
              rand = genrand()
              ! Loop over species
              do m=1, n_species
                ! Decide which species to sit on that site
                if ((rand .ge. sum(setup%species_cs(0:(m-1)))) .and. &
                    (rand .le. sum(setup%species_cs(0:m)))) then
                  if (check(m) .lt. species_count(m)) then
                    config(i,j,k,l) = m
                    check(m) = check(m) + 1
                  end if
                end if
              end do ! Species
            end do ! While
          end do ! i
        end do ! j
      end do ! k
    end do ! l

  end subroutine initial_setup

  subroutine set_concentrations(parameters)
    type(run_params) :: parameters
    integer :: i

    parameters%species_cs = 0.0_real64

    ! Setup array of concentrations
    ! Sloppy coding, I know, but keeps input strictly to input file
    ! Currently handles up to six species
    if (parameters%n_species .eq. 2) then
      parameters%species_cs(1) = parameters%c_a
      parameters%species_cs(2) = 1.0_real64 - parameters%c_a
    else if (parameters%n_species .eq. 3) then
      parameters%species_cs(1) = parameters%c_a
      parameters%species_cs(2) = parameters%c_b
      parameters%species_cs(3) = 1.0_real64 - parameters%c_a&
                                            - parameters%c_b
    else if (parameters%n_species .eq. 4) then
      parameters%species_cs(1) = parameters%c_a
      parameters%species_cs(2) = parameters%c_b
      parameters%species_cs(3) = parameters%c_c
      parameters%species_cs(4) = 1.0_real64 - parameters%c_a&
                                            - parameters%c_b&
                                            - parameters%c_c
    else if (parameters%n_species .eq. 5) then
      parameters%species_cs(1) = parameters%c_a
      parameters%species_cs(2) = parameters%c_b
      parameters%species_cs(3) = parameters%c_c
      parameters%species_cs(4) = parameters%c_d
      parameters%species_cs(5) = 1.0_real64 - parameters%c_a&
                                            - parameters%c_b&
                                            - parameters%c_c&
                                            - parameters%c_d
    else if (parameters%n_species .eq. 6) then
      parameters%species_cs(1) = parameters%c_a
      parameters%species_cs(2) = parameters%c_b
      parameters%species_cs(3) = parameters%c_c
      parameters%species_cs(4) = parameters%c_d
      parameters%species_cs(5) = parameters%c_e
      parameters%species_cs(6) = 1.0_real64 - parameters%c_a&
                                            - parameters%c_b&
                                            - parameters%c_c&
                                            - parameters%c_d&
                                            - parameters%c_e
    else
      stop 'Unsupported number of species.'
    end if

    ! Error handling
    do i=1, parameters%n_species
      if ((parameters%species_cs(i) .lt. 0.0_real64) .or. &
          (parameters%species_cs(i) .gt. 1.0_real64)) then
        stop 'Bad concentrations specified'
      end if
    end do
  end subroutine set_concentrations
end module initialise
