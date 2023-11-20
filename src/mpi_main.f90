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
  use io
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
  if(my_rank == 0) call write_info('s')

  if(my_rank == 0) then
    write(6,'(22("-"),x,"Parsing name of input file",x,22("-"),/)')
  end if
    if(.not. get_arg('control', control)) then
      if (my_rank == 0) then
        print*, 'Input file not specified with "control=<name>"'
        print*, 'Defaulting to searching for "input.txt"'
        print*, ' '
      end if
      control = 'input.txt'
    else
      if (my_rank == 0) then
        print*, 'Input file name is: ', control
        print*, ' '
      end if
    end if

  ! Read the control file
  call read_control_file(control, setup, my_rank)

  if(my_rank == 0) then
    call echo_control_file(setup)
    write(6,'(/,20("-"),x,"Parsed input file successfully",x,20("-"),/)')
  end if

  ! Initialise some global arrays
  call initialise_global_arrays(setup)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Initialise some global arrays
  call initialise_local_arrays(setup)

  if(my_rank == 0) then
    write(6,'(15("-"),x,"Reading atom-atom interaction parameters",x,15("-"),/)')
  end if

  ! Read the exchange coefficients matrix
  call read_exchange(setup)

  if (my_rank .eq. 0) then
    ! Print it out as a sanity check for user
    call pretty_print_exchange(setup)
  end if

  if(my_rank == 0) then
    write(6,'(72("-"),/)')
  end if

  ! Setup the concentrations array
  call set_concentrations(setup)

  if(my_rank == 0) then
    write(6,'(17("-"),x,"Initialising random number generators",x,16("-"),/)')
  end if

  ! Initialise the prng
  seed = f90_init_genrand(seedtime, int(my_rank, kind=C_INT))

  call comms_wait()
  print*, 'Thread ', my_rank, ' has seed ', seed
  call flush(6)

  call comms_wait()

  if(my_rank == 0) then
    write(6,'(/,72("-"),/)')
  end if

  if (my_rank .eq. 0) then
    write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"))')
  end if

  call run_model(setup, my_rank)

  call comms_reduce_results(setup)

  if (my_rank .eq. 0) then
    write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')


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
  call local_clean_up(setup)

  if(my_rank == 0) call write_info('f')

  ! Start MPI
  call comms_finalise()

end program main
