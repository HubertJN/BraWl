!----------------------------------------------------!
! MC Routine for computing ground state of an alloy. !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

program main
  
  use initialise
  use comms
  use model_run
  use mpi_shared_data
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use write_diagnostics
  use command_line
  use display

  implicit none

  ! Random number seed
  integer(8) :: seed
  
  ! Whether or not to seed the PRNG
  integer :: seedtime=1

  ! Runtime parameters type
  type(run_params) :: setup

  ! Which control file to read
  character(len=30) :: control

  ! Start MPI
  call comms_initialise()

  !------------------------------!
  ! Parse command line arguments !
  !------------------------------!
  call parse_args()

  ! make a directory for the grid states, diagnostics, and radial_densities for each thread
  if(my_rank == 0) call execute_command_line('mkdir -p grids')
  if(my_rank == 0) call execute_command_line('mkdir -p diagnostics')
  if(my_rank == 0) call execute_command_line('mkdir -p radial_densities')

  if(.not. get_arg('control', control)) then
    if (my_rank == 0) then
      print*, 'Could not parse name of control file'
      print*, 'Defaulting to searching for "input.txt"'
      print*, ' '
    end if
    control = 'input.txt'
  end if

  if (my_rank .eq. 0) call print_parse()

  ! Read the control file
  call read_control_file(control, setup)

  if (my_rank .eq. 0) call echo_control_file(setup)

  ! Initialise some global arrays
  call initialise_global_arrays(setup)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Initialise some global arrays
  call initialise_local_arrays(setup)

  ! Read the exchange coefficients matrix
  call read_exchange(setup)

  ! Setup the concentrations array
  call set_concentrations(setup)

  ! Initialise the prng
  seed = f90_init_genrand(seedtime, int(my_rank, kind=C_INT))

  call comms_wait()
  print*, 'Thread ', my_rank, ' has seed ', seed
  call comms_wait()

  if (my_rank .eq. 0) then
    ! Print it out as a sanity check for user
    call pretty_print_exchange(setup)
  end if

  if (my_rank .eq. 0) then
    print*, '##########################'
    print*, '# Commencing Simulation! #'
    print*, '##########################'
  end if

  call run_model(setup, my_rank)

  call comms_reduce_results(setup)

  if (my_rank .eq. 0) then
    print*, '##########################'
    print*, '#  Simulation Complete!  #'
    print*, '##########################', new_line('a')

    ! Write energy diagnostics
    call diagnostics_writer('diagnostics/av_energy_diagnostics.dat', temperature, &
                            av_energies_of_T, av_C_of_T, av_acceptance_of_T)

    !Write the radial densities to file
    call ncdf_radial_density_writer('radial_densities/av_radial_density.nc', av_rho_of_T, &
                                  shells, temperature, av_energies_of_T, setup)

  end if

  ! Clean up
  call global_clean_up()

  ! Clean up
  call local_clean_up()

  ! Start MPI
  call comms_finalise()

end program main
