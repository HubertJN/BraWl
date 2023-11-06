!=============================================================!
! Routines to write to .xyz format for visualisation of alloy !
! configurations.                                             !
!                                                             !
! C. D. Woodgate                                         2021 !
!=============================================================!
module write_xyz

  use kinds
  use mpi_shared_data
  use analytics

  implicit none

  contains

  subroutine xyz_writer(filename, configuration, setup)
    ! Simulation setup information
    type(run_params), intent(in) :: setup

    ! Data to write to file
    integer(int16), dimension(:,:,:,:), allocatable, intent(in) :: configuration

    ! Sizes of grid
    integer, dimension(4) :: grid_sizes

    ! Filename to which to write
    character(len=*), intent(in) :: filename

    integer :: n_particles, i, j, k, l

    real(real64), dimension(3) :: pos

    n_particles = total_particle_count(setup, configuration)

    grid_sizes = shape(configuration)

    open(unit=7, file=filename)

    write(7, *) n_particles

    write(7,*) ''

    do i=1, grid_sizes(1)
      do j=1, grid_sizes(2)
        do k=1, grid_sizes(3)
          do l=1, grid_sizes(4)
            pos = real(j)*setup%lattice_vectors(1,:) + &
                  real(k)*setup%lattice_vectors(2,:) + &
                  real(l)*setup%lattice_vectors(3,:) + &
                  real(i)*setup%basis_vectors
            if (configuration(i,j,k,l) .eq. 1_int16) then
              write(7,*) setup%a, pos(1),pos(2),pos(3)
            else if (configuration(i,j,k,l) .eq. 2_int16) then
              write(7,*) setup%b, pos(1),pos(2),pos(3)
            else if (configuration(i,j,k,l) .eq. 3_int16) then
              write(7,*) setup%c, pos(1),pos(2),pos(3)
            else if (configuration(i,j,k,l) .eq. 4_int16) then
              write(7,*) setup%d, pos(1),pos(2),pos(3)
            else if (configuration(i,j,k,l) .eq. 5_int16) then
              write(7,*) setup%e, pos(1),pos(2),pos(3)
            end if
          end do
        end do
      end do
    end do

    close(7)

  end subroutine xyz_writer
   
end module write_xyz
