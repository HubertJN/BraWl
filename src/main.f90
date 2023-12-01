!----------------------------------------------------!
! MC Routine for computing ground state of an alloy. !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!
program main
  
  use initialise
  use comms
  use shared_data
  use metropolis
  use io
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use write_diagnostics
  use display

  implicit none

  ! Runtime parameters type
  type(run_params) :: setup

  ! Start MPI
  call comms_initialise()

  ! Print software info to the screen
  if(my_rank == 0) call write_info('s')

  ! Parse inputs
  call parse_inputs(setup, my_rank)

  ! Allocate space for atom-atom interaction parameters
  call initialise_interaction(setup)

  ! Read in atom-atom interaction
  call read_exchange(setup, my_rank)

  ! Initialise PNRG
  call initialise_pnrg(setup%seedtime)

  ! Initialise some function pointers
  call initialise_function_pointers(setup)

  ! Make directories for data
  call make_data_directories(my_rank)

  ! Initialise some global arrays
  call initialise_global_arrays(setup)

  ! Initialise some global arrays
  call initialise_local_arrays(setup)

  !---------------!
  ! Main Routines !
  !---------------!
  if (setup%mode == 301) then

    ! Metropolis with Kawasaki dynamics
    call metropolis_simulated_annealing(setup, my_rank)

  end if

  ! Clean up
  call global_clean_up()

  ! Clean up
  call clean_up_interaction(setup)

  ! Clean up
  call local_clean_up(setup)

  ! Print software info to the screen
  if(my_rank == 0) call write_info('f')

  ! Finalise MPI
  call comms_finalise()

end program main
