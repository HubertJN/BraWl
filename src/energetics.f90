!-------------------------------------------------------!
! Module containing routines associated with energetics !
! of alloy.                                             !
!                                                       !
! C Woodgate, Warwick                              2020 !
!-------------------------------------------------------!

module energetics

  use kinds
  use mpi_shared_data
  use c_functions
  
  implicit none

  contains

  !-------------------------------------------------------------!
  ! Function to compute the full Hamiltonian for the BCC system !
  !-------------------------------------------------------------!
  function total_energy(setup,config) result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer :: i, j, k, l

    energy=0.0_real64

    do l=1, setup%n_3
      do k=1, setup%n_2
        do j=1, setup%n_1
          do i=1, setup%n_basis
            if (config(i,j,k,l) .eq. 0_int16) cycle
            energy = energy + setup%nbr_energy(config, j, k, l)
          end do
        end do
      end do
    end do
  end function total_energy

  !----------------------------------------------!
  ! Function to compute the bcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_shell1_energy(setup, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip1,im1,jp1,jm1,kp1,km1

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip1 = modulo(site_i, 2*setup%n_1) + 1
    im1 = modulo(site_i-2, 2*setup%n_1) + 1
    jp1 = modulo(site_j, 2*setup%n_2) + 1
    jm1 = modulo(site_j-2, 2*setup%n_2) + 1
    kp1 = modulo(site_k, 2*setup%n_3) + 1
    km1 = modulo(site_k-2, 2*setup%n_3) + 1
      
    allocate(nbrs(8))
    nbrs(1) = config(1, ip1, jp1, kp1)
    nbrs(2) = config(1, ip1, jp1, km1)
    nbrs(3) = config(1, ip1, jm1, kp1)
    nbrs(4) = config(1, ip1, jm1, km1)
    nbrs(5) = config(1, im1, jp1, kp1)
    nbrs(6) = config(1, im1, jp1, km1)
    nbrs(7) = config(1, im1, jm1, kp1)
    nbrs(8) = config(1, im1, jm1, km1)
    do i=1, 8
      energy = energy + V_ex(species, nbrs(i), 1)
    end do
    deallocate(nbrs)
  end function bcc_shell1_energy

  !----------------------------------------------!
  ! Function to compute the bcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_shell2_energy(setup, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip2,im2,jp2,jm2,kp2,km2

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip2 = modulo(site_i+1, 2*setup%n_1) + 1
    im2 = modulo(site_i-3, 2*setup%n_1) + 1
    jp2 = modulo(site_j+1, 2*setup%n_2) + 1
    jm2 = modulo(site_j-3, 2*setup%n_2) + 1
    kp2 = modulo(site_k+1, 2*setup%n_3) + 1
    km2 = modulo(site_k-3, 2*setup%n_3) + 1
      
    allocate(nbrs(6))
    nbrs(1) = config(1, ip2, site_j  , site_k  )
    nbrs(2) = config(1, im2, site_j  , site_k  )
    nbrs(3) = config(1, site_i  , jm2, site_k  )
    nbrs(4) = config(1, site_i  , jp2, site_k  )
    nbrs(5) = config(1, site_i  , site_j  , kp2)
    nbrs(6) = config(1, site_i  , site_j  , km2)
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i), 1)
    end do
    deallocate(nbrs)
  end function bcc_shell2_energy

  !----------------------------------------------!
  ! Function to compute the bcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_shell3_energy(setup, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, ip1,im1,jp1,jm1,kp1,km1
    integer :: ip2,im2,jp2,jm2,kp2,km2

    energy=0.0_real64
    
    ! Compute where my neighbours are
    ip1 = modulo(site_i, setup%n_1) + 1
    im1 = modulo(site_i-2, setup%n_1) + 1
    jp1 = modulo(site_j, setup%n_2) + 1
    jm1 = modulo(site_j-2, setup%n_2) + 1
    kp1 = modulo(site_k, setup%n_3) + 1
    km1 = modulo(site_k-2, setup%n_3) + 1
    ip2 = modulo(site_i+1, setup%n_1) + 1
    im2 = modulo(site_i-3, setup%n_1) + 1
    jp2 = modulo(site_j+1, setup%n_2) + 1
    jm2 = modulo(site_j-3, setup%n_2) + 1
    kp2 = modulo(site_k+1, setup%n_3) + 1
    km2 = modulo(site_k-3, setup%n_3) + 1
      
    allocate(nbrs(12))
    nbrs(1)  = config(1,   ip1,    jp1,    kp2)
    nbrs(2)  = config(1,   im1,    jm1,    km2)
    nbrs(3)  = config(1,   ip1,    jp2,    kp1)
    nbrs(4)  = config(1,   im1,    jm2,    km1)
    nbrs(5)  = config(1,   ip2,    jp1,    kp1)
    nbrs(6)  = config(1,   im2,    jm1,    km1)
    nbrs(7)  = config(1,site_i,    jp1,    km1)
    nbrs(8)  = config(1,site_i,    jm1,    kp1)
    nbrs(9)  = config(1,   ip1, site_j,    km1)
    nbrs(10) = config(1,   im1, site_j,    kp1)
    nbrs(11) = config(1,   ip1,    jm1, site_k)
    nbrs(12) = config(1,   im1,    jp1, site_k)
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i), 1)
    end do
    deallocate(nbrs)
  end function bcc_shell3_energy

  !----------------------------------------------!
  ! Function to compute the bcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_shell4_energy(setup, site_i, site_j, site_k, &
                             config, species)     &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt, upup, dndn, &
               fwfw, bwbw, ltlt, rtrt

    energy=0.0_real64
    
    up = modulo(  site_i, setup%n_1) + 1
    dn = modulo(site_i-2, setup%n_1) + 1
    lt = modulo(  site_j, setup%n_2) + 1
    rt = modulo(site_j-2, setup%n_2) + 1
    fw = modulo(  site_k, setup%n_3) + 1
    bw = modulo(site_k-2, setup%n_3) + 1
    upup = modulo(site_i+1, setup%n_1) + 1
    dndn = modulo(site_i-3, setup%n_1) + 1
    ltlt = modulo(site_j+1, setup%n_2) + 1
    rtrt = modulo(site_j-3, setup%n_2) + 1
    fwfw = modulo(site_k+1, setup%n_3) + 1
    bwbw = modulo(site_k-3, setup%n_3) + 1

    allocate(nbrs(24))
    nbrs(1)   = config(1,  dndn,   ltlt,     bw)
    nbrs(2)   = config(1,  dndn,     lt,   bwbw)
    nbrs(3)   = config(1,  dndn,     lt, site_k)
    nbrs(4)   = config(1,  dndn, site_j,     bw)
    nbrs(5)   = config(1,    dn,   ltlt,   bwbw)
    nbrs(6)   = config(1,    dn,   ltlt, site_k)
    nbrs(7)   = config(1,    dn,     lt,     fw)
    nbrs(8)   = config(1,    dn, site_k,   bwbw)
    nbrs(9)   = config(1,    dn,     rt,     bw)
    nbrs(10)  = config(1,    dn,     rt,     fw)
    nbrs(11)  = config(1,site_i,   ltlt,     bw)
    nbrs(12)  = config(1,site_i,     lt,   bwbw)
    nbrs(13)  = config(1,site_i,     rt,   fwfw)
    nbrs(14)  = config(1,site_i,   rtrt,     fw)
    nbrs(15)  = config(1,    up,     lt,     bw)
    nbrs(16)  = config(1,    up,     lt,     fw)
    nbrs(17)  = config(1,    up, site_j,   fwfw)
    nbrs(18)  = config(1,    up,     rt,     bw)
    nbrs(19)  = config(1,    up,   rtrt, site_k)
    nbrs(20)  = config(1,    up,   rtrt,   fwfw)
    nbrs(21)  = config(1,  upup, site_j,     fw)
    nbrs(22)  = config(1,  upup,     rt, site_k)
    nbrs(23)  = config(1,  upup,     rt,   fwfw)
    nbrs(24)  = config(1,  upup,   rtrt,     fw)
     
    ! Sum them
    do i=1, 24
      energy = energy + V_ex(species, nbrs(i),4)
    end do
    deallocate(nbrs)
  end function bcc_shell4_energy

  !----------------------------------------------!
  ! Function to compute the bcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_shell5_energy(setup, site_i, site_j, site_k, &
                             config,  species)     &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16), intent(in) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt

    energy=0.0_real64
    
    upup = modulo(site_i+1, setup%n_1) + 1
    dndn = modulo(site_i-3, setup%n_1) + 1
    ltlt = modulo(site_j+1, setup%n_2) + 1
    rtrt = modulo(site_j-3, setup%n_2) + 1
    fwfw = modulo(site_k+1, setup%n_3) + 1
    bwbw = modulo(site_k-3, setup%n_3) + 1

    allocate(nbrs(8))
    nbrs(1)  = config(1,  dndn,   ltlt,   bwbw)
    nbrs(2)  = config(1,  dndn, site_j, site_k)
    nbrs(3)  = config(1,site_i,   ltlt, site_k)
    nbrs(4)  = config(1,site_i, site_j,   bwbw)
    nbrs(5)  = config(1,site_i, site_j,   fwfw)
    nbrs(6)  = config(1,site_i,   rtrt, site_k)
    nbrs(7)  = config(1,  upup, site_j, site_k)
    nbrs(8)  = config(1,  upup,   rtrt,   fwfw)
     
    ! Sum them
    do i=1, 8
      energy = energy + V_ex(species, nbrs(i),5)
    end do
    deallocate(nbrs)
  end function bcc_shell5_energy

  !----------------------------------------------!
  ! Function to compute the fcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_energy_1shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy= bcc_shell1_energy(setup, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_1shells

  !----------------------------------------------!
  ! Function to compute the fcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_energy_2shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_2shells

  !----------------------------------------------!
  ! Function to compute the fcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_energy_3shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_3shells

  !----------------------------------------------!
  ! Function to compute the fcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_energy_4shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_4shells

  !----------------------------------------------!
  ! Function to compute the fcc neighbour energy !
  ! for a given site.                            !
  !----------------------------------------------!
  function bcc_energy_5shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy = bcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell3_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell4_energy(setup, site_i, site_j, site_k, config, species) &
           + bcc_shell5_energy(setup, site_i, site_j, site_k, config, species)
    
  end function bcc_energy_5shells

!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell1_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i, up, dn, fw, bw, lt, rt
!
!    energy=0.0_real64
!    
!    ! Compute where my neighbours are
!    up = modulo(site_i, 2*setup%n_x) + 1
!    dn = modulo(site_i-2, 2*setup%n_x) + 1
!    lt = modulo(site_j, 2*setup%n_y) + 1
!    rt = modulo(site_j-2, 2*setup%n_y) + 1
!    fw = modulo(site_k, 2*setup%n_z) + 1
!    bw = modulo(site_k-2, 2*setup%n_z) + 1
!      
!    allocate(nbrs(12))
!    nbrs(1)  = config(site_i, rt, fw)
!    nbrs(2)  = config(site_i, rt, bw)
!    nbrs(3)  = config(site_i, lt, fw)
!    nbrs(4)  = config(site_i, lt, bw)
!    nbrs(5)  = config(up, rt, site_k)
!    nbrs(6)  = config(up, lt, site_k)
!    nbrs(7)  = config(up, site_j, fw)
!    nbrs(8)  = config(up, site_j, bw)
!    nbrs(9)  = config(dn, rt, site_k)
!    nbrs(10) = config(dn, lt, site_k)
!    nbrs(11) = config(dn, site_j, fw)
!    nbrs(12) = config(dn, site_j, bw)
!    do i=1, 12
!      energy = energy + V_ex(species, nbrs(i), 1)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell1_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell2_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i
!    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt
!
!    energy=0.0_real64
!    
!    upup = modulo(site_i+1, 2*setup%n_x) + 1
!    dndn = modulo(site_i-3, 2*setup%n_x) + 1
!    ltlt = modulo(site_j+1, 2*setup%n_y) + 1
!    rtrt = modulo(site_j-3, 2*setup%n_y) + 1
!    fwfw = modulo(site_k+1, 2*setup%n_z) + 1
!    bwbw = modulo(site_k-3, 2*setup%n_z) + 1
!    allocate(nbrs(6))
!    nbrs(1)  = config(upup, site_j, site_k)
!    nbrs(2)  = config(dndn, site_j, site_k)
!    nbrs(3)  = config(site_i, ltlt, site_k)
!    nbrs(4)  = config(site_i, rtrt, site_k)
!    nbrs(5)  = config(site_i, site_j, fwfw)
!    nbrs(6)  = config(site_i, site_j, bwbw)
!    do i=1, 6
!      energy = energy + V_ex(species, nbrs(i), 2)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell2_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell3_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i, up, dn, fw, bw, lt, rt
!    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt
!
!    energy=0.0_real64
!    
!    up = modulo(site_i, 2*setup%n_x) + 1
!    dn = modulo(site_i-2, 2*setup%n_x) + 1
!    lt = modulo(site_j, 2*setup%n_y) + 1
!    rt = modulo(site_j-2, 2*setup%n_y) + 1
!    fw = modulo(site_k, 2*setup%n_z) + 1
!    bw = modulo(site_k-2, 2*setup%n_z) + 1
!      
!    upup = modulo(site_i+1, 2*setup%n_x) + 1
!    dndn = modulo(site_i-3, 2*setup%n_x) + 1
!    ltlt = modulo(site_j+1, 2*setup%n_y) + 1
!    rtrt = modulo(site_j-3, 2*setup%n_y) + 1
!    fwfw = modulo(site_k+1, 2*setup%n_z) + 1
!    bwbw = modulo(site_k-3, 2*setup%n_z) + 1
!
!    allocate(nbrs(24))
!    nbrs(1)   = config(dndn, lt, fw)
!    nbrs(2)   = config(dndn, lt, bw)
!    nbrs(3)   = config(dndn, rt, fw)
!    nbrs(4)   = config(dndn, rt, bw)
!    nbrs(5)   = config(upup, lt, fw)
!    nbrs(6)   = config(upup, lt, bw)
!    nbrs(7)   = config(upup, rt, fw)
!    nbrs(8)   = config(upup, rt, bw)
!    nbrs(9)   = config(up, ltlt, fw)
!    nbrs(10)  = config(dn, ltlt, fw)
!    nbrs(11)  = config(up, ltlt, bw)
!    nbrs(12)  = config(dn, ltlt, bw)
!    nbrs(13)  = config(up, rtrt, fw)
!    nbrs(14)  = config(dn, rtrt, fw)
!    nbrs(15)  = config(up, rtrt, bw)
!    nbrs(16)  = config(dn, rtrt, bw)
!    nbrs(17)  = config(up, lt, fwfw)
!    nbrs(18)  = config(dn, lt, fwfw)
!    nbrs(19)  = config(up, rt, fwfw)
!    nbrs(20)  = config(dn, rt, fwfw)
!    nbrs(21)  = config(up, lt, bwbw)
!    nbrs(22)  = config(dn, lt, bwbw)
!    nbrs(23)  = config(up, rt, bwbw)
!    nbrs(24)  = config(dn, rt, bwbw)
!    do i=1, 24
!      energy = energy + V_ex(species, nbrs(i), 3)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell3_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell4_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i
!    integer :: upup, dndn, fwfw, bwbw, ltlt, rtrt
!
!    energy=0.0_real64
!    
!    upup = modulo(site_i+1, 2*setup%n_x) + 1
!    dndn = modulo(site_i-3, 2*setup%n_x) + 1
!    ltlt = modulo(site_j+1, 2*setup%n_y) + 1
!    rtrt = modulo(site_j-3, 2*setup%n_y) + 1
!    fwfw = modulo(site_k+1, 2*setup%n_z) + 1
!    bwbw = modulo(site_k-3, 2*setup%n_z) + 1
!
!    allocate(nbrs(12))
!    nbrs(1)  = config(upup, site_j, bwbw)
!    nbrs(2)  = config(dndn, site_j, bwbw)
!    nbrs(3)  = config(site_i, ltlt, bwbw)
!    nbrs(4)  = config(site_i, rtrt, bwbw)
!    nbrs(5)  = config(upup, ltlt, site_k)
!    nbrs(6)  = config(dndn, ltlt, site_k)
!    nbrs(7)  = config(upup, rtrt, site_k)
!    nbrs(8)  = config(dndn, rtrt, site_k)
!    nbrs(9)  = config(upup, site_j, fwfw)
!    nbrs(10) = config(dndn, site_j, fwfw)
!    nbrs(11) = config(site_i, ltlt, fwfw)
!    nbrs(12) = config(site_i, rtrt, fwfw)
!    do i=1, 12
!      energy = energy + V_ex(species, nbrs(i), 4)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell4_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell5_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i, up, dn, fw, bw, lt, rt
!    integer :: upupup, dndndn, fwfwfw, bwbwbw, ltltlt, rtrtrt
!
!    energy=0.0_real64
!    
!    up = modulo(site_i, 2*setup%n_x) + 1
!    dn = modulo(site_i-2, 2*setup%n_x) + 1
!    lt = modulo(site_j, 2*setup%n_y) + 1
!    rt = modulo(site_j-2, 2*setup%n_y) + 1
!    fw = modulo(site_k, 2*setup%n_z) + 1
!    bw = modulo(site_k-2, 2*setup%n_z) + 1
!      
!    upupup = modulo(site_i+2, 2*setup%n_x) + 1
!    dndndn = modulo(site_i-4, 2*setup%n_x) + 1
!    ltltlt = modulo(site_j+2, 2*setup%n_y) + 1
!    rtrtrt = modulo(site_j-4, 2*setup%n_y) + 1
!    fwfwfw = modulo(site_k+2, 2*setup%n_z) + 1
!    bwbwbw = modulo(site_k-4, 2*setup%n_z) + 1
!
!    allocate(nbrs(24))
!    nbrs(1)   = config(up, site_j, bwbwbw)
!    nbrs(2)   = config(dn, site_j, bwbwbw)
!    nbrs(3)   = config(site_i, rt, bwbwbw)
!    nbrs(4)   = config(site_i, lt, bwbwbw)
!    nbrs(5)   = config(site_i, ltltlt, bw)
!    nbrs(6)   = config(upupup, site_j, bw)
!    nbrs(7)   = config(dndndn, site_j, bw)
!    nbrs(8)   = config(site_i, rtrtrt, bw)
!    nbrs(9)   = config(up, ltltlt, site_k)
!    nbrs(10)  = config(dn, ltltlt, site_k)
!    nbrs(11)  = config(upupup, lt, site_k)
!    nbrs(12)  = config(dndndn, lt, site_k)
!    nbrs(13)  = config(upupup, rt, site_k)
!    nbrs(14)  = config(dndndn, rt, site_k)
!    nbrs(15)  = config(up, rtrtrt, site_k)
!    nbrs(16)  = config(dn, rtrtrt, site_k)
!    nbrs(17)  = config(site_i, ltltlt, fw)
!    nbrs(18)  = config(upupup, site_j, fw)
!    nbrs(19)  = config(dndndn, site_j, fw)
!    nbrs(20)  = config(site_k, rtrtrt, fw)
!    nbrs(21)  = config(up, site_j, fwfwfw)
!    nbrs(22)  = config(dn, site_j, fwfwfw)
!    nbrs(23)  = config(site_i, rt, fwfwfw)
!    nbrs(24)  = config(site_i, lt, fwfwfw)
!    do i=1, 24
!      energy = energy + V_ex(species, nbrs(i), 5)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell5_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_shell6_energy(setup, site_i, site_j, site_k, &
!                             config, species)     &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16), intent(in) :: species
!    integer(int16), allocatable, dimension(:) :: nbrs
!    integer :: i, upup, dndn, fwfw, bwbw, ltlt, rtrt
!
!    energy=0.0_real64
!
!    upup = modulo(site_i+1, 2*setup%n_x) + 1
!    dndn = modulo(site_i-3, 2*setup%n_x) + 1
!    ltlt = modulo(site_j+1, 2*setup%n_y) + 1
!    rtrt = modulo(site_j-3, 2*setup%n_y) + 1
!    fwfw = modulo(site_k+1, 2*setup%n_z) + 1
!    bwbw = modulo(site_k-3, 2*setup%n_z) + 1
!
!    
!    allocate(nbrs(8))
!    nbrs(1)   = config(upup, ltlt, bwbw)
!    nbrs(2)   = config(dndn, ltlt, bwbw)
!    nbrs(3)   = config(upup, rtrt, bwbw)
!    nbrs(4)   = config(dndn, rtrt, bwbw)
!    nbrs(5)   = config(upup, ltlt, fwfw)
!    nbrs(6)   = config(dndn, ltlt, fwfw)
!    nbrs(7)   = config(upup, rtrt, fwfw)
!    nbrs(8)   = config(dndn, rtrt, fwfw)
!    do i=1, 8
!      energy = energy + V_ex(species, nbrs(i), 6)
!    end do
!    deallocate(nbrs)
!  end function fcc_shell6_energy
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_1shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy= fcc_shell1_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_1shells
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_2shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy = fcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell2_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_2shells
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_3shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy = fcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell3_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_3shells
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_4shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy = fcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell3_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell4_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_4shells
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_5shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy = fcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell3_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell4_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell5_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_5shells
!
!  !----------------------------------------------!
!  ! Function to compute the fcc neighbour energy !
!  ! for a given site.                            !
!  !----------------------------------------------!
!  function fcc_energy_6shells(setup, config, site_i, site_j, site_k) &
!           result(energy)
!    integer(int16), allocatable, dimension(:,:,:), intent(in) :: config
!    real(real64) :: energy
!    class(run_params), intent(in) :: setup
!    integer, intent(in) :: site_i, site_j, site_k
!    integer(int16) :: species
!
!    species = config(site_i, site_j, site_k)
!
!    energy = fcc_shell1_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell2_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell3_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell4_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell5_energy(setup, site_i, site_j, site_k, config, species) &
!           + fcc_shell6_energy(setup, site_i, site_j, site_k, config, species)
!    
!  end function fcc_energy_6shells

  !----------------------------------------------!
  ! Function to compute the neighbour energy for !
  ! a given site.                                !
  !----------------------------------------------!
  function simple_cubic_1shell_energy(setup, site_i, site_j, site_k, config, species) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species
    integer(int16), allocatable, dimension(:) :: nbrs
    integer :: i, up, dn, fw, bw, lt, rt

    energy=0.0_real64
    
    ! Compute where my neighbours are
    up = modulo(  site_i, setup%n_1) + 1
    dn = modulo(site_i-2, setup%n_1) + 1
    lt = modulo(  site_j, setup%n_2) + 1
    rt = modulo(site_j-2, setup%n_2) + 1
    fw = modulo(  site_k, setup%n_3) + 1
    bw = modulo(site_k-2, setup%n_3) + 1
      
    allocate(nbrs(6))

    ! Compute the energies of neighbours
    nbrs(1) = config(1,    up, site_j, site_k)
    nbrs(2) = config(1,    dn, site_j, site_k)
    nbrs(3) = config(1,site_i,     lt, site_k)
    nbrs(4) = config(1,site_i,     rt, site_k)
    nbrs(5) = config(1,site_i, site_j,     fw)
    nbrs(6) = config(1,site_i, site_j,     bw)
    
    ! Sum them
    do i=1, 6
      energy = energy + V_ex(species, nbrs(i),1)
    end do
 
    deallocate(nbrs)
  end function simple_cubic_1shell_energy

  !----------------------------------------------!
  ! Function to compute the neighbour energy for !
  ! a given site.                                !
  !----------------------------------------------!
  function simple_cubic_energy_1shells(setup, config, site_i, site_j, site_k) &
           result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    real(real64) :: energy
    class(run_params), intent(in) :: setup
    integer, intent(in) :: site_i, site_j, site_k
    integer(int16) :: species

    species = config(1,site_i, site_j, site_k)

    energy = simple_cubic_1shell_energy(setup, site_i, site_j, site_k, config, species)
    
  end function simple_cubic_energy_1shells

  !-----------------------------------------------!
  ! Function to compute energy of swapping a pair !
  !-----------------------------------------------!
  function pair_energy(setup, config, idx1, idx2)&
       result(energy)
    integer(int16), allocatable, dimension(:,:,:,:), intent(in) :: config
    type(run_params), intent(in) :: setup
    integer, dimension(4), intent(in) :: idx1, idx2
    real(real64) :: energy
    integer(int16) :: species1, species2

    species1 = config(idx1(1), idx1(2), idx1(3), idx1(4))
    species2 = config(idx2(1), idx2(2), idx2(3), idx2(4))

    energy = setup%nbr_energy(config, idx1(2), idx1(3), idx1(4)) &
           + setup%nbr_energy(config, idx2(2), idx2(3), idx2(4))
  end function pair_energy

end module energetics
