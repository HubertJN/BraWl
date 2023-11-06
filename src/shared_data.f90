!--------------------------------------------------!
! File with module containing all global variables !
!                                                  !
! C Woodgate, Warwick                         2020 !
!--------------------------------------------------!

module shared_data

  use kinds
  
  implicit none

  save

  ! Convert temperature input in Kelvin to simulation
  ! units by multiplying by this constant.
  ! Simulation units of energy are Ry.
  real(real64), parameter :: t_conversion &
                             =8.167333262e-5_real64/13.605693122_real64

  !---------------------------------!
  ! Type storing runtime parameters !
  !---------------------------------!
  type run_params
    ! Number of cubic cells in each direction
    integer :: n_x, n_y, n_z
    ! Number of species
    integer :: n_species
    ! Number of monte carlo steps
    integer :: mc_steps
    ! Lattice type
    character(len=20) :: lattice
    ! Exchange type
    character(len=50) :: exchange
    ! Concentrations
    real(real64) :: species_cs(0:6) = 0.0_real64
    ! Concentration of species a
    real(real64) :: c_a, c_b, c_c, c_d, c_e
    ! Chemical elements
    character(len=2) :: a,b,c,d,e
    ! Lattice parameter (for writing xyz file)
    real(real64) :: lattice_parameter
    ! Lattice vectors (for writing xyz file)
    real(real64), dimension(3,3) :: lattice_vectors
    ! Vector to second basis atom (for writing xyz file)
    real(real64), dimension(3) :: basis_vectors
    ! Inverse temperature
    real(real64) :: beta   
    ! If so, what is the start temperature
    real(real64) :: T_start
    ! How large do we want our steps in temperature to be?
    real(real64) :: delta_T
    ! Number of temperatures at which to run
    ! N.B. If this is one, we just run a simulation at one temperature
    integer :: T_steps
    ! Maximum neighbour distance to consider
    integer :: n_shells
    ! Distance of radial shells to consider
    integer :: radial_d
    ! Hamiltonian to use
    procedure(hamiltonian), pointer, pass :: full_energy => null()
    ! neighbour energy function to use
    procedure(neighbour), pointer, pass :: nbr_energy => null()
    ! Random site function to use
    procedure(rand_site), pointer, pass :: rdm_site => null()
    ! Do we want to swap neighbours
    logical :: nbr_swap
    ! Do we want to store long range order parameters
    logical :: lro
    ! Random neighbour function to use
    procedure(rand_neighbour), pointer, pass :: rdm_nbr => null()
    ! Monte Carlo step to call. (Neighbour swap or whole lattice swap)
    procedure(monte_carlo), pointer, pass :: mc_step => null()
  end type run_params

  !----------------------------------------!
  ! Interface for neighbour implementation !
  !----------------------------------------!
  interface
    function hamiltonian(setup)
      use kinds
      import :: run_params
      real(real64) :: hamiltonian
      class(run_params), intent(in) :: setup
    end function
    function neighbour(setup, site_i, site_j, site_k)
      use kinds
      import :: run_params
      real(real64) :: neighbour
      class(run_params), intent(in) :: setup
      integer, intent(in) :: site_i, site_j, site_k
    end function
    function rand_site(setup)
      use kinds
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4) :: rand_site
    end function
    function rand_neighbour(setup, site)
      use kinds
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4), intent(in) :: site
      integer, dimension(4) :: rand_neighbour
    end function
    function monte_carlo(setup, beta) result(accept)
      use kinds
      import :: run_params
      class(run_params), intent(in) :: setup
      integer :: accept
      real(real64), intent(in) :: beta
    end function monte_carlo
  end interface

  integer(int16), allocatable, dimension(:,:,:) :: global_config
  real(real64), allocatable, dimension(:) :: energies_of_T
  real(real64), allocatable, dimension(:) :: C_of_T
  real(real64), allocatable, dimension(:) :: acceptance_of_T
  real(real64), allocatable, dimension(:,:,:,:,:) :: order_of_T
  real(real64), allocatable, dimension(:,:,:,:) :: rho_of_T
  real(real64), allocatable, dimension(:) :: temperature
  real(real64), dimension(:,:,:,:), allocatable :: order
  real(real64), dimension(:,:,:), allocatable :: V_ex
  real(real64), allocatable, dimension(:,:,:,:) :: rho_r_T
  real(real64), allocatable, dimension(:) :: shells

  contains

  !----------------------------------!
  ! Subroutine to parse control file !
  !----------------------------------!
  subroutine read_control_file(filename, parameters)
    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(run_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios

    check = .false.

    ios=0; line=0

    open(15, file=filename, iostat=ios)

    if (ios .ne. 0) then
      stop 'Could not parse input file. Aborting...'
    end if

    write(*,'(a)', advance='no') new_line('a')
    print*, '###############################'
    print*, '#     Parsing input file      #'
    print*, '###############################'

    print*, '# Input file name: ', filename

    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('n_x')  
          read(buffer, *, iostat=ios) parameters%n_x
          print*, '# Read n_x = ', parameters%n_x
          check(1) = .true.
        case ('n_y')  
          read(buffer, *, iostat=ios) parameters%n_y
          print*, '# Read n_y = ', parameters%n_y
          check(2) = .true.
        case ('n_z')  
          read(buffer, *, iostat=ios) parameters%n_z
          print*, '# Read n_z = ', parameters%n_z
          check(3) = .true.
        case ('a')  
          read(buffer, *, iostat=ios) parameters%a
          print*, '# Read a = ', parameters%a
        case ('c_a')  
          read(buffer, *, iostat=ios) parameters%c_a
          print*, '# Read c_a = ', parameters%c_a
        case ('b')  
          read(buffer, *, iostat=ios) parameters%b
          print*, '# Read b = ', parameters%b
        case ('c_b')  
          read(buffer, *, iostat=ios) parameters%c_b
          print*, '# Read c_b = ', parameters%c_b
        case ('c')  
          read(buffer, *, iostat=ios) parameters%c
          print*, '# Read c = ', parameters%c
        case ('c_c')  
          read(buffer, *, iostat=ios) parameters%c_c
          print*, '# Read c_c = ', parameters%c_c
        case ('d')  
          read(buffer, *, iostat=ios) parameters%d
          print*, '# Read d = ', parameters%d
        case ('c_d')  
          read(buffer, *, iostat=ios) parameters%c_d
          print*, '# Read c_d = ', parameters%c_d
        case ('e')  
          read(buffer, *, iostat=ios) parameters%e
          print*, '# Read e = ', parameters%e
        case ('c_e')  
          read(buffer, *, iostat=ios) parameters%c_e
          print*, '# Read c_e = ', parameters%c_e
        case ('n_steps')  
          read(buffer, *, iostat=ios) parameters%mc_steps
          print*, '# Read n_steps = ', parameters%mc_steps
          check(4) = .true.
        case ('n_species')  
          read(buffer, *, iostat=ios) parameters%n_species
          print*, '# Read n_species = ', parameters%n_species
          check(5) = .true.
        case ('lattice_parameter')  
          read(buffer, *, iostat=ios) parameters%lattice_parameter
          print*, '# Read lattice = ', parameters%lattice_parameter
        case ('lattice')  
          read(buffer, *, iostat=ios) parameters%lattice
          print*, '# Read lattice = ', parameters%lattice
          check(6) = .true.
        case ('exchange')  
          read(buffer, *, iostat=ios) parameters%exchange
          print*, '# Read exchange = ', parameters%exchange
          check(7) = .true.
        case ('n_shells')  
          read(buffer, *, iostat=ios) parameters%n_shells
          print*, '# Read n_shells = ', parameters%n_shells
          check(8) = .true.
        case ('radial_d')  
          read(buffer, *, iostat=ios) parameters%radial_d
          print*, '# Read radial_d = ', parameters%radial_d
          check(8) = .true.
        case ('T_start')  
          read(buffer, *, iostat=ios) parameters%T_start
          print*, '# Read T_start = ', parameters%T_start
        case ('delta_T')  
          read(buffer, *, iostat=ios) parameters%delta_T
          print*, '# Read delta_T = ', parameters%delta_T
        case ('T_steps')  
          read(buffer, *, iostat=ios) parameters%T_steps
          print*, '# Read T_steps = ', parameters%T_steps
        case ('lro')  
          read(buffer, *, iostat=ios) parameters%lro
          print*, '# Read nbr_swap = ', parameters%lro
        case ('nbr_swap')  
          read(buffer, *, iostat=ios) parameters%nbr_swap
          print*, '# Read nbr_swap = ', parameters%nbr_swap
        case default
          print*, '# Skipping invalid label'
        end select
      end if
    end do

    do line=1, 8
      if(.not. check(line)) then
        print*, 'Missing parameter in input file', line
        stop
      end if
    end do

    print*, '# Finished parsing input file #'
    print*, '###############################', new_line('a')

  end subroutine read_control_file

  !-----------------------------------------------------!
  ! Subroutine to read in exchange parameters from file !
  !-----------------------------------------------------!
  subroutine read_exchange(setup)
    type(run_params) , intent(in) :: setup
    V_ex = 0.0_real64
    open(16, file=setup%exchange)
    read(16,*) V_ex   
    close(16)
  end subroutine read_exchange

end module shared_data
