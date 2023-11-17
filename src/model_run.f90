!----------------------------------------------------!
! MC Routine for computing ground state of an alloy. !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

module model_run
  
  use initialise
  use mpi_shared_data
  use io
  use comms
  use kinds
  use c_functions
  use write_netcdf
  use write_xyz
  use write_diagnostics
  use command_line
  use display

  implicit none

  contains

  subroutine run_model(setup, my_rank)

    ! Rank of this processor
    integer, intent(in) :: my_rank

    ! Arrays for storing data
    type(run_params) :: setup
  
    ! How many times do we store the state for output
    integer(int32) :: N_store
  
    ! Integers used in calculations
    integer :: i,j, n_save, div_steps, ierr, accept
    
    ! Temperature and temperature steps
    real(real64) :: beta, temp, sim_temp, current_energy, &
                    step_E, step_Esq, C, acceptance
  
    ! Name of file for grid state and xyz file at this temperature
    character(len=34) :: grid_file
    character(len=36) :: xyz_file

    ! Name of file for writing diagnostics at the end
    character(len=43) :: diagnostics_file
  
    ! Name of file for writing radial densities at the end
    character(len=37) :: radial_file
  
    ! Setup the concentrations array
    call set_concentrations(setup)
  
    ! Set up the lattice
    call initial_setup(setup, config)

    call lattice_shells(setup, shells)
  
    N_store=0
    n_save=1000
    div_steps = setup%mc_steps/1000
  
    ! Loop over temperature steps
    do j=1, setup%T_steps
  
      step_E = 0.0_real64; step_Esq=0.0_real64
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
    
        ! Make one MC move
        accept = setup%mc_step(config, beta)
  
        acceptance = acceptance + accept

        ! Write percentage progress to screen
        if (mod(i, div_steps) .eq. 0) then
          current_energy = setup%full_energy(config)
          step_E   = step_E + current_energy
          step_Esq = step_Esq + current_energy**2
          if (my_rank .eq. 0) then
          write(6,'(a, a,f7.2,a,f5.1,a)',advance='no') achar(13), &
          " Temperature ", temp, "  ",(0.1*real(i)/real(div_steps)), "% complete on process 0."
          write(*,'(a, a,f7.2,a,f5.1,a)',advance='no') achar(13), &
           " Temperature ", temp, "  ",(0.1*real(i)/real(div_steps)), "% complete on process 0."
          end if
        end if
    
      end do

  
      ! Store the average energy at this temperature
      energies_of_T(j) = step_E/1000
  
      C = (step_Esq/1000 - (step_E/1000)**2)/(sim_temp*temp)

      acceptance_of_T(j) = acceptance/real(setup%mc_steps)
  
      ! Store the specific heat capacity at this temperature
      C_of_T(j) = C
    
      write(grid_file, '(A11 I3.3 A11 I4.4 F2.1 A3)') 'grids/proc_', my_rank, '_grid_at_T_', &
                                           int(temp), temp-int(temp), '.nc'
      write(xyz_file, '(A11 I3.3 A12 I4.4 F2.1 A4)') 'grids/proc_', my_rank, 'config_at_T_', &
                                           int(temp), temp-int(temp),'.xyz'
  
      ! Write xyz file
      call xyz_writer(xyz_file, config, setup)
  
      ! Write grid to file
      call ncdf_grid_state_writer(grid_file , ierr, &
                             config, temp, setup)
  
      ! Compute the radial densities at the end of this temperature
      call radial_densities(setup, config, setup%radial_d, shells, rho_of_T, j)
  
      if (my_rank ==0) then
      ! Write that we have completed a particular temperature
      write(*, '(a)', advance='no') new_line('a')
      end if
    
    end do ! Loop over temperature

    write(radial_file, '(A22 I3.3 A12)') 'radial_densities/proc_', my_rank, '_rho_of_T.nc'
    write(diagnostics_file, '(A17 I4.4 A22)') 'diagnostics/proc_', my_rank, &
                                         'energy_diagnostics.dat'
  
    
    ! Write energy diagnostics
    call diagnostics_writer(diagnostics_file, temperature, &
                            energies_of_T, C_of_T, acceptance_of_T)

    ! Write the radial densities to file
    call ncdf_radial_density_writer(radial_file, rho_of_T, &
                                  shells, temperature, energies_of_T, setup)

  
  end subroutine run_model
end module model_run
