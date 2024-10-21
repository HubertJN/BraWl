!----------------------------------------------------------------------!
! Energy Spectrum Module                                               !
!                                                                      !
! This program run MC for the given system and records all visited     !
! energies                                                             !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module energy_spectrum
    use initialise
    use kinds
    use shared_data
    use c_functions
    use random_site
    use metropolis
    use mpi
  
    implicit none
  
    contains
  
    subroutine es_main(setup, es_setup, my_rank)
        ! Rank of this processor
        integer, intent(in) :: my_rank

        ! Arrays for storing data
        type(run_params) :: setup
        type(es_params) :: es_setup
      
        ! Integers used in calculations
        integer :: ierr, unique_energy, i, energy_min_loc, energy_max_loc
        
        ! Temperature and temperature steps
        real(real64) :: acceptance, step, energy_to_ry
      
        real(real64), allocatable :: energy_spectrum(:), energy_spectrum_condensed(:), bin_edges(:)
        integer, allocatable :: energy_spectrum_sort(:)
        
        energy_to_ry=setup%n_atoms/(eV_to_Ry*1000)

        allocate(energy_spectrum(es_setup%unique_energy_count))
        allocate(bin_edges(es_setup%bins+1))
        energy_spectrum = 0.0_real64
        unique_energy = 1
  
        ! Set up the lattice
        call initial_setup(setup, config)
    
        call lattice_shells(setup, shells, config)
      
        ! Are we swapping neighbours or on the whole lattice?
        if (setup%nbr_swap) then
          setup%mc_step => monte_carlo_step_nbr
        else
          setup%mc_step => monte_carlo_step_lattice
        end if
  
        if(my_rank == 0) then
          write(6,'(/,72("-"),/)')
          write(6,'(24("-"),x,"Commencing Simulation!",x,24("-"),/)')
          print*, "Number of atoms", setup%n_atoms
        end if
  
        acceptance = run_es_sweeps(setup, es_setup, config, energy_spectrum, unique_energy)

        print*, "Unique energies found: ", unique_energy-1
        
        allocate(energy_spectrum_condensed(unique_energy-1))
        allocate(energy_spectrum_sort(unique_energy-1))
        energy_spectrum_condensed = energy_spectrum(1:unique_energy-1)

        call merge_argsort(energy_spectrum_condensed, energy_spectrum_sort)
        energy_spectrum_sort = energy_spectrum_sort(SIZE(energy_spectrum_sort):1:-1 )
        energy_spectrum_condensed = energy_spectrum_condensed(energy_spectrum_sort)
        call ncdf_writer_1d("energy_spectrum.dat", ierr, energy_spectrum_condensed)

        energy_min_loc = minloc(abs(energy_spectrum_condensed-es_setup%energy_min*energy_to_ry), 1)
        energy_max_loc = minloc(abs(energy_spectrum_condensed-es_setup%energy_max*energy_to_ry), 1)

        energy_spectrum_condensed = energy_spectrum_condensed(energy_min_loc:energy_max_loc)

        !print*, energy_spectrum_condensed

        step = REAL(SIZE(energy_spectrum_condensed)/es_setup%bins)
        bin_edges(1) = MINVAL(energy_spectrum_condensed)
        do i=1, es_setup%bins
            bin_edges(i+1) = energy_spectrum_condensed(NINT(i*step))
        end do

        call ncdf_writer_1d("bin_edges.dat", ierr, bin_edges)

        write(*, *)
        write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
    
    end subroutine es_main    
  
    function run_es_sweeps(setup, es_setup, config, energy_spectrum, unique_energy) result(acceptance)
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        class(es_params), intent(in) :: es_setup
        real(real64), dimension(:), intent(inout) :: energy_spectrum
        integer, intent(inout) :: unique_energy
  
        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, eps, delta_e
        integer :: acceptance, i, cycle_loop, reject
        integer(int16) :: site1, site2
  
        eps = EPSILON(eps)
        cycle_loop = 0
        reject = 0
  
        ! Establish total energy before any moves
        e_unswapped = setup%full_energy(config)
        e_swapped = e_unswapped
  
        acceptance = 0.0_real64
  
        do i=1, es_setup%mc_sweeps*setup%n_atoms
            ! Make one MC trial
            ! Generate random numbers
            rdm1 = setup%rdm_site()
            rdm2 = setup%rdm_site()
  
            ! Get what is on those sites
            site1 = config(rdm1(1), rdm1(2), rdm1(3), rdm1(4))
            site2 = config(rdm2(1), rdm2(2), rdm2(3), rdm2(4))
  
            call pair_swap(config, rdm1, rdm2)
    
            ! Calculate energy if different species
            if (site1 /= site2) then
                e_swapped = setup%full_energy(config)
            else
                reject = reject - 1
            end if

            delta_e = e_swapped - e_unswapped

            ! Accept or reject move
            if (cycle_loop == 0 .and. delta_e < 0) then
                reject = 0
                e_unswapped = e_swapped
            else if (cycle_loop == 1 .and. delta_e > 0) then
                reject = 0
                e_unswapped = e_swapped
            else
                reject = reject + 1
                call pair_swap(config, rdm1, rdm2)
            end if
  
            if (MINVAL(ABS(energy_spectrum-e_swapped)) > eps) then
              energy_spectrum(unique_energy) = e_swapped
              unique_energy = unique_energy + 1
            end if

            if (unique_energy > es_setup%unique_energy_count) then
                exit
            end if

            if (reject > setup%n_atoms) then
                reject = 0
                cycle_loop = ABS(cycle_loop-1)
            end if
        end do
  
    end function run_es_sweeps

    subroutine merge_argsort(r,d)
              real(kind=8), intent(in), dimension(:) :: r
              integer, intent(out), dimension(size(r)) :: d

              integer, dimension(size(r)) :: il

              integer :: stepsize
              integer :: i,j,n,left,k,ksize

              n = size(r)

              do i=1,n
                  d(i)=i
              end do

              if ( n==1 ) return

              stepsize = 1
              do while (stepsize<n)
                  do left=1,n-stepsize,stepsize*2
                      i = left
                      j = left+stepsize
                      ksize = min(stepsize*2,n-left+1)
                      k=1

                      do while ( i<left+stepsize .and. j<left+ksize )
                          if ( r(d(i))>r(d(j)) ) then
                              il(k)=d(i)
                              i=i+1
                              k=k+1
                          else
                              il(k)=d(j)
                              j=j+1
                              k=k+1
                          endif
                      end do
                        if ( i<left+stepsize ) then
                          ! fill up remaining from left
                          il(k:ksize) = d(i:left+stepsize-1)
                      else
                          ! fill up remaining from right
                          il(k:ksize) = d(j:left+ksize-1)
                      endif
                      d(left:left+ksize-1) = il(1:ksize)
                  end do
                  stepsize=stepsize*2
              end do

              return
          end subroutine
  end module energy_spectrum
  
