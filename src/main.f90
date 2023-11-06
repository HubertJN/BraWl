!----------------------------------------------------!
! MC Routine for computing ground state of an alloy. !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

program main
  
  use initialise
  use shared_data
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use write_diagnostics
  use command_line
  use display

  implicit none

  ! How many times do we store the state for output
  integer(int32) :: N_store

  ! Integers used in calculations
  integer :: i,j, n_save, div_steps, ierr, accept
  
  ! Whether or not to seed the PRNG
  integer :: seedtime=1

  ! Temperature and temperature steps
  real(real64) :: beta, temp, sim_temp, current_energy, &
                  step_E, step_Esq, C, acceptance

  ! Runtime parameters type
  type(run_params) :: setup

  ! Which control file to read
  character(len=30) :: control

  ! Name of file for grid state at this temperature
  character(len=25) :: grid_file
  character(len=28) :: xyz_file

  !------------------------------!
  ! Parse command line arguments !
  !------------------------------!
  call parse_args()

  if(.not. get_arg('control', control)) then
    print*, 'Could not parse name of control file'
    print*, 'Defaulting to searching for "input.txt"'
    print*, ' '
    control = 'input.txt'
  end if

  ! Read the control file
  call read_control_file(control, setup)

  ! Initialise the arrays for storing data
  call initialise_arrays(setup)

  ! Read the exchange coefficients matrix
  call read_exchange(setup)

  ! Print it out as a sanity check for user
  call pretty_print_exchange(setup)

  ! Initialise the prng
  call f90_init_genrand(seedtime)

  ! Setup the concentrations array
  call set_concentrations(setup)

  ! Set up the lattice
  call initial_setup(setup)

  ! Check the average concentrations
  call print_particle_count(setup)
  
  ! make a directory for the grid states
  call execute_command_line('mkdir -p grids')

  N_store=0
  n_save=1000
  div_steps = setup%mc_steps/1000

  call lattice_shells(setup, global_config, shells, setup%radial_d)

  print*, '##########################'
  print*, '# Commencing Simulation! #'
  print*, '##########################', new_line('a')

  ! Loop over temperature steps
  do j=1, setup%T_steps

    order = 0.0_real64; step_E = 0.0_real64; step_Esq=0.0_real64
    acceptance = 0.0_real64
  
    ! Work out the temperature and corresponding beta
    temp = setup%T_start + real(j-1, real64)*setup%delta_T
    sim_temp = temp*t_conversion
    beta = 1.0_real64/sim_temp
  
    ! Store this in an array
    temperature(j) = temp
  
    !==================!
    ! Monte Carlo Loop !
    !==================!
    do i=1, setup%mc_steps
  
      ! Store if we need to
      if (modulo(i, n_save) == 0) then
        call store_state(order, setup)
      end if
  
      ! Make one MC move
      accept = setup%mc_step(beta)

      acceptance = acceptance + accept
  
      ! Write percentage progress to screen
      if (mod(i, div_steps) == 0) then
        current_energy = setup%full_energy()
        step_E   = step_E + current_energy
        step_Esq = step_Esq + current_energy**2
      end if
  
    end do

    ! Store the average energy at this temperature
    energies_of_T(j) = step_E/1000

    C = (step_Esq/1000 - (step_E/1000)**2)/(sim_temp*temp)
    print*, ''
    print*, 'C is ', C

    acceptance_of_T(j) = acceptance/real(setup%mc_steps)

    ! Store the specific heat capacity at this temperature
    C_of_T(j) = C
  
    ! Add this to the average state
    call average_state(order, setup, &
                       setup%mc_steps/n_save)

    ! Store it
    if (setup%lro) then
      do i=1, setup%n_species
        order_of_T(:,:,:,i,j) = order(i,:, :, :)
      end do
    end if

    write(grid_file, '(A16 I4.4 F2.1 A3)') 'grids/grid_at_T_', &
                                         int(temp), temp-int(temp), '.nc'
    write(xyz_file, '(A18 I4.4 F2.1 A4)') 'grids/config_at_T_', &
                                         int(temp), temp-int(temp),'.xyz'

  
    ! Write xyz file
    call xyz_writer(xyz_file, global_config, setup)

    ! Write grid to file
    call ncdf_grid_state_writer(grid_file , ierr, &
                           global_config, temp, setup)

    ! Compute the radial densities at the end of this temperature
    call radial_densities(setup, global_config, setup%radial_d, shells, rho_r_T, j)

    ! Write that we have completed a particular temperature
    write(*,'(a,a,f7.2,a,f5.1,a,a,f15.6)',advance='yes') &
     " Temperature ", temp, " 100% complete.", &
     " Energy is ", current_energy
  
  end do ! Loop over temperature

  print*, '##########################'
  print*, '#  Simulation Complete!  #'
  print*, '##########################', new_line('a')

  ! Print the state to the screen
  call display_grid(global_config)

  ! Check the average concentrations
  call print_particle_count(setup)

  ! Write energy diagnostics
  call diagnostics_writer('energy_diagnostics.dat', temperature, &
                          energies_of_T, C_of_T, acceptance_of_T)

  if (setup%lro) then
    ! Write order parameters to file
    call ncdf_order_writer('order_parameters.nc' , ierr, &
                           order_of_T, temperature, setup)
  end if

  ! Write the radial densities to file
  call ncdf_radial_density_writer('radial_density.nc', rho_r_T, &
                                  shells, temperature, energies_of_T, setup)

  ! Clean up
  call clean_up()

end program main
