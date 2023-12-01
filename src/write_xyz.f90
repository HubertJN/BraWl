!----------------------------------------------------------------------!
! write_xyz.f90                                                        !
!                                                                      !
! Module for writing xyz file                                          !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module write_xyz

  use kinds
  use mpi_shared_data
  use analytics

  implicit none

  contains

  !--------------------------------------------------------------------!
  ! Routine to write xyz file                                          !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
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
            if (configuration(i,j,k,l) .eq. 0) cycle
            pos = real(j)*setup%lattice_vectors(1,:) + &
                  real(k)*setup%lattice_vectors(2,:) + &
                  real(l)*setup%lattice_vectors(3,:) + &
                  real(i)*setup%basis_vectors
            write(7,*) setup%species_names(configuration(i,j,k,l)), pos(1),pos(2),pos(3)
          end do
        end do
      end do
    end do

    close(unit=7)

  end subroutine xyz_writer
   
end module write_xyz
