!----------------------------------------------------------------------!
! io.f90                                                               !
!                                                                      !
! Module containing input/output routines.                             !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module io

  use kinds
  use mpi_shared_data
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Subroutine to write software version, date, and time               !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine write_info(point)

    character (len=1) point
    character (len=10) date,time

    ! Get date and time
    call date_and_time(date=date,time=time)

    write(6,'(/,72("="))')
    write(6,'(19x,"Bonte Warlo Version 1.0, 07.11.23")')
    write(6,'(72("-"))')

    if (point .eq. 's') then
      write(6,'(15x,"This run started at",1x,a," on",1x,a)')           &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
    else if(point .eq. 'f') then
      write(6,'(15x,"This run finished at",1x,a," on",1x,a)')          &
               time(1:2)//":"//time(3:4)//":"//time(5:6),              &
               date(7:8)//"."//date(5:6)//"."//date(1:4)
    endif

    write(6,'(72("="),/)' )

  end subroutine write_info

  !--------------------------------------------------------------------!
  ! Subroutine to parse control file                                   !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine read_control_file(filename, parameters)

    character(len=*), intent(in) :: filename
    logical, dimension(8) :: check
    type(run_params) :: parameters
    character(len=100) :: buffer, label
    integer :: line, pos, ios, ierr

    check = .false.

    ios=0; line=0

    open(15, file=filename, iostat=ios)

    do while (ios==0)

      read(15, "(A)", iostat=ios) buffer

      if(ios==0) then
        line=line+1

        pos = scan(buffer, '=')
        label=buffer(1:pos-1)
        buffer = buffer(pos+1:)

        select case (label)
        case ('n_1')  
          read(buffer, *, iostat=ios) parameters%n_1
          check(1) = .true.
        case ('n_2')  
          read(buffer, *, iostat=ios) parameters%n_2
          check(2) = .true.
        case ('n_3')  
          read(buffer, *, iostat=ios) parameters%n_3
          check(3) = .true.
        case ('n_basis')  
          read(buffer, *, iostat=ios) parameters%n_basis
          check(3) = .true.
        case ('a')  
          read(buffer, *, iostat=ios) parameters%a
        case ('c_a')  
          read(buffer, *, iostat=ios) parameters%c_a
        case ('b')  
          read(buffer, *, iostat=ios) parameters%b
        case ('c_b')  
          read(buffer, *, iostat=ios) parameters%c_b
        case ('c')  
          read(buffer, *, iostat=ios) parameters%c
        case ('c_c')  
          read(buffer, *, iostat=ios) parameters%c_c
        case ('d')  
          read(buffer, *, iostat=ios) parameters%d
        case ('c_d')  
          read(buffer, *, iostat=ios) parameters%c_d
        case ('e')  
          read(buffer, *, iostat=ios) parameters%e
        case ('c_e')  
          read(buffer, *, iostat=ios) parameters%c_e
        case ('n_steps')  
          read(buffer, *, iostat=ios) parameters%mc_steps
          check(4) = .true.
        case ('n_species')  
          read(buffer, *, iostat=ios) parameters%n_species
          check(5) = .true.
        case ('lattice_parameter')  
          read(buffer, *, iostat=ios) parameters%lattice_parameter
        case ('lattice')  
          read(buffer, *, iostat=ios) parameters%lattice
          check(6) = .true.
        case ('exchange')  
          read(buffer, *, iostat=ios) parameters%exchange
          check(7) = .true.
        case ('n_shells')  
          read(buffer, *, iostat=ios) parameters%n_shells
          check(8) = .true.
        case ('radial_d')  
          read(buffer, *, iostat=ios) parameters%radial_d
        case ('T_start')  
          read(buffer, *, iostat=ios) parameters%T_start
        case ('delta_T')  
          read(buffer, *, iostat=ios) parameters%delta_T
        case ('T_steps')  
          read(buffer, *, iostat=ios) parameters%T_steps
        case ('lro')  
          read(buffer, *, iostat=ios) parameters%lro
        case ('nbr_swap')  
          read(buffer, *, iostat=ios) parameters%nbr_swap
        case default
        end select
      end if
    end do

  end subroutine read_control_file

  subroutine print_parse()

    print*, '###############################'
    print*, '#     Parsing input file      #'

  end subroutine print_parse

  subroutine echo_control_file(parameters)
    type(run_params) :: parameters

    print*, '# Read n_1 = ', parameters%n_1
    print*, '# Read n_2 = ', parameters%n_2
    print*, '# Read n_3 = ', parameters%n_3
    print*, '# Read n_basis = ', parameters%n_basis
    print*, '# Read a = ', parameters%a
    print*, '# Read c_a = ', parameters%c_a
    print*, '# Read b = ', parameters%b
    print*, '# Read c_b = ', parameters%c_b
    if (parameters%n_species .gt. 2) then
      print*, '# Read c = ', parameters%c
      print*, '# Read c_c = ', parameters%c_c
    end if
    if (parameters%n_species .gt. 3) then
      print*, '# Read d = ', parameters%d
      print*, '# Read c_d = ', parameters%c_d
    end if
    if (parameters%n_species .gt. 4) then
      print*, '# Read e = ', parameters%e
      print*, '# Read c_e = ', parameters%c_e
    end if
    print*, '# Read n_steps = ', parameters%mc_steps
    print*, '# Read n_species = ', parameters%n_species
    print*, '# Read lattice = ', parameters%lattice_parameter
    print*, '# Read lattice = ', parameters%lattice
    print*, '# Read exchange = ', parameters%exchange
    print*, '# Read n_shells = ', parameters%n_shells
    print*, '# Read radial_d = ', parameters%radial_d
    print*, '# Read T_start = ', parameters%T_start
    print*, '# Read delta_T = ', parameters%delta_T
    print*, '# Read T_steps = ', parameters%T_steps
    print*, '# Read nbr_swap = ', parameters%lro
    print*, '# Read nbr_swap = ', parameters%nbr_swap
    print*, ' ', parameters%nbr_swap

  end subroutine echo_control_file

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

end module io
