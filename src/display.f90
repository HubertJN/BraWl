!----------------------------------------------------------------------!
! display.f90                                                          !
!                                                                      !
! Module for displaying some things to the screen.                     !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module display

  use kinds
  use mpi_shared_data
  
  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Routine to print the read in exchange coefficients to the screen   !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine pretty_print_exchange(setup, units)
    integer :: i,j,k
    character(len=*), optional :: units
    character(len=10) :: aunits
    type(run_params), intent(in) :: setup

    ! Can do mRy or meV, but default is meV
    if (present(units)) then
      aunits = units
    else
      aunits = 'meV'
    end if

    ! Case mRy
    if (trim(aunits) .eq. 'mRy') then
      write(*,"(x,'V_ij to be used in calculation (mRy):',/)") 
      do i=1, setup%interaction_range
        write(*,'(a, I2, a, a)') ' Interchange interaction on shell: ', i, new_line('a')
          write(*, '(A)', advance='no') ' '
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') '          ', setup%species_names(j)
        end do
          write(*, '(A)', advance='yes') ''
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') ' ', setup%species_names(j)
          do k=1, setup%n_species
            write(*, '(A, F8.3, A)', advance='no') '  ', 1000*V_ex(k,j,i), '  '
            if (k .eq. setup%n_species) write(*,'(a)') ''
          end do
          if (j .eq. setup%n_species) write(*,'(a)') ''
        end do
      end do
      write(*,'(a)') ''
    ! Case meV
    else if (trim(aunits) .eq. 'meV') then
      write(*,"(x,'V_ij to be used in calculation (meV):',/)") 
      do i=1, setup%interaction_range
        write(*,'(a, I2, a, a)') ' Interchange interaction on shell: ', i, new_line('a')
          write(*, '(A)', advance='no') ' '
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') '          ', setup%species_names(j)
        end do
          write(*, '(A)', advance='yes') ''
        do j=1, setup%n_species
          write(*, '(A, A)', advance='no') ' ', setup%species_names(j)
          do k=1, setup%n_species
            write(*, '(A, F8.3, A)', advance='no') '  ', 1000*13.606*V_ex(k,j,i), '  '
            if (k .eq. setup%n_species) write(*,'(a)') ''
          end do
          if (j .eq. setup%n_species) write(*,'(a)') ''
        end do
      end do
      write(*,'(a)') ''
    end if

  end subroutine pretty_print_exchange

  !--------------------------------------------------------------------!
  ! Routine to print the current state of the grid                     !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
  subroutine display_grid(grid, show_borders, title)

    integer(kind=int16), intent(in), dimension(:,:,:) :: grid
    logical, intent(in), optional :: show_borders
    character(len=*), optional :: title
    logical :: borders
    integer, dimension(3) :: sizes
    integer :: ix, iy, iz
    character(len=1) :: c
    character(len=4), parameter :: clrstr = char(27)//'[2j'

    borders = .true.

    if (present(show_borders)) borders = show_borders

    write(*,'(a)') clrstr

    sizes = shape(grid)
    if (present(title)) then
      write(*, '(a)') title
    end if

    do iz = sizes(3), 1, -1

      write(*, '(a, I2)') 'Layer for z = ', iz

      if (borders) write(*, '(a)') repeat('=', sizes(1)+2)

      do iy = sizes(2), 1, -1
        if (borders) write(*, '(a)', advance='no') '|'
        do ix = 1, sizes(1)
          if (grid(ix,iy,iz) .ne. 0) then
            write(c, '(I1)') grid(ix, iy, iz)
          else
            c = ' '
          end if
          write(*, '(a)', advance='no') c
        end do
        if (borders) write(*, '(a)', advance='no') '|'
        write(*, '(a)') ''
      end do
      if (borders) write(*, '(a)') repeat('=', sizes(1)+2)
    end do

  end subroutine display_grid

end module display
