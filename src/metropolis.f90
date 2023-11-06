!----------------------------------------------------!
! File with module containing functions and routines !
! to call                                            !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

module metropolis

  use kinds
  use mpi_shared_data
  use c_functions
  use energetics
  use random_site
  use analytics
  
  implicit none

  contains

  !---------------------------------------------------------------!
  ! Runs one MC step assuming pairs swapped across entire lattice !
  !---------------------------------------------------------------!
  function monte_carlo_step_lattice(setup, config, beta) result(accept)
    integer(int16), allocatable, dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(int16) :: site1, site2

    ! Generate random numbers
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_site()

    ! Get what is on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! Don't do anything if same species
    if (site1 == site2) then
      accept = 0
      return
    end if

    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    call pair_swap(config, rdm1, rdm2)

    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    delta_e = e_swapped - e_unswapped

    ! MC swap
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return
    else
      accept = 0
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_lattice

  !----------------------------------------------------------!
  ! Runs one MC step assuming only neighbours can be swapped !
  !----------------------------------------------------------!
  function monte_carlo_step_nbr(setup, config, beta) result(accept)
    integer(int16), allocatable, dimension(:,:,:,:) :: config
    class(run_params), intent(in) :: setup
    integer :: accept
    integer, dimension(4) :: rdm1, rdm2
    real(real64) , intent(in) :: beta
    real(real64) :: e_unswapped, e_swapped, delta_e
    integer(int16) :: site1, site2

    ! Generate random numbers
    rdm1 = setup%rdm_site()
    rdm2 = setup%rdm_nbr(rdm1)

    ! Get what is on those sites
    site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
    site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))

    ! Don't do anything if same species
    if (site1 == site2) then
      accept=1
      return
    end if

    e_unswapped = pair_energy(setup, config, rdm1, rdm2)

    call pair_swap(config, rdm1, rdm2)

    e_swapped = pair_energy(setup, config, rdm1, rdm2)   

    delta_e = e_swapped - e_unswapped

    ! MC swap
    if(delta_e .lt. 0.0_real64) then
      accept = 1
      return
    else if (genrand() .lt. exp(-beta*delta_E)) then
      accept = 1
      return
    else
      accept = 0
      call pair_swap(config, rdm1, rdm2)
    end if

  end function monte_carlo_step_nbr

end module metropolis
