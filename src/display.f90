!----------------------------------------------------!
! Module with various tools and utilities for code.  !
!                                                    !
! C Woodgate, Warwick                           2020 !
!----------------------------------------------------!

module display

  use kinds
  use mpi_shared_data
  
  implicit none

  contains

  !-------------------------------------------------------!
  ! Print the read in exchange coefficients to the screen !
  !-------------------------------------------------------!
  subroutine pretty_print_exchange(setup)
    integer :: i,j,k
    type(run_params), intent(in) :: setup

    write(*,'(a)') ''
    write(*,'(a, a)') 'V_ij to be used in calculation:'

    do i=1, setup%n_shells
      write(*,'(a, I2, a)') 'Interchange interaction on shell: ', i
      do j=1, setup%n_species
        do k=1, setup%n_species
          write(*, '(F9.6, A)', advance='no') V_ex(k,j,i), '  '
          if (k .eq. setup%n_species) write(*,'(a)') ''
        end do
      end do
    end do
    write(*,'(a)') ''
  end subroutine pretty_print_exchange

  !-------------------------------------------------!
  ! Subroutine to display current state of the grid !
  !-------------------------------------------------!
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
