!----------------------------------------------------------------------!
! Wang-Landau module                                                   !
!                                                                      !
! H. Naguszewski, Warwick                                         2024 !
!----------------------------------------------------------------------!

module wang_landau
  use initialise
  use kinds
  use shared_data
  use c_functions
  use random_site
  use metropolis
  use mpi

  implicit none

  contains

  subroutine wl_main(setup, wl_setup, my_rank)
      ! Rank of this processor
      integer, intent(in) :: my_rank
      integer :: ierror, num_proc

      ! Arrays for storing data
      type(run_params) :: setup
      type(tmmc_params) :: wl_setup
    
      ! Integers used in calculations
      integer :: i,j,k, ierr, accept, bins
      
      ! Temperature and temperature steps
      real(real64) :: temp, acceptance, beta, flatness
    
      ! wl variables and arrays
      real(real64) :: bin_width, energy_to_ry, target_energy, wl_f
      real(real64), allocatable :: bin_edges(:), wl_hist, wl_logdos
      logical :: first_reset

      ! Allocate arrays
      bins = wl_setup%bins

      allocate(bin_edges(bins+1))
      allocate(wl_hist(bins))
      allocate(wl_logdos(bins))

      ! Set temperature
      temp = setup%T*k_b_in_Ry

      ! Conversion meV/atom to Rydberg
      energy_to_ry=setup%n_atoms/(eV_to_Ry*1000)

      bin_width = (wl_setup%energy_max - wl_setup%energy_min)/real(wl_setup%bins)*energy_to_ry
      do i=1, bins+1
          bin_edges(i) = wl_setup%energy_min*energy_to_ry + (i-1)*bin_width
      end do

      wl_hist = 0.0_real64; wl_logdos = 0.0_real64; wl_f = 0.05_real64
      flatness = 0.0_real64; first_reset = .False.
      !---------------------------------!

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

      !---------!
      ! Burn in !
      !---------!
      beta = 1.0_real64/temp
      if (wl_setup%burn_in) then
          call wl_burn_in(setup, wl_setup, config, temp, target_energy, trans_matrix, bin_edges, bins)
      end if

      !--------------------!
      ! Target Temperature !
      !--------------------!
      do while (wl_f > wl_setup%tolerance)
          acceptance = run_wl_sweeps(setup, wl_setup, config, bins, bin_edges, bins_mpi, &
          bin_edges_mpi, wl_hist, wl_logdos, wl_f)

          flatness = 100.0_real64*minval(wl_hist)/(sum(wl_hist)/size(wl_hist))

          if (minval(wl_hist) > 10.0_real64) then
    
            if (first_reset == .False.) then
              first_reset = .True.
              wl_hist = 0.0_real64

            else
              !Check if we're 80% flat
              if flatness > 80.0_real64 then
  
                  !Reset the histogram
                  wl_hist = 0.0_real64
                  wl_f = wl_f * 0.5_real64
  
                  !Reduce f
                  wl_logdos = wl_logdos - minval(wl_logdos)
              end if
            end if

          if(my_rank == 0) then
              write(6,'(a,i0,a,f6.2,a,f6.2,a)',advance='yes') "Weight Update ", i, ": Accepted ", &
               (acceptance/(wl_setup%mc_sweeps*setup%n_atoms)*100.0), "% of Monte Carlo moves. Time taken: ", end-start, "s"
          end if
      end do
      
      if (my_rank == 0) then
          ! Write output files
          call ncdf_writer_1d("wl_dos_bins.dat", ierr, bin_edges)

          call ncdf_writer_1d("wl_dos.dat", ierr, wl_dos)
      end if
      
      if(my_rank == 0) then
          write(*, *)
          write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
      end if
  
  end subroutine wl_main    

  integer function bin_index(energy, bin_edges, bins) result(index)
      integer, intent(in) :: bins
      real(real64), intent(in) :: energy
      real(real64), dimension(:), intent(in) :: bin_edges
      real(real64) :: bin_range

      bin_range = bin_edges(bins+1) - bin_edges(1)
      index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
  end function bin_index

  function run_wl_sweeps(setup, wl_setup, config, bins, bin_edges, bins_mpi, &
      bin_edges_mpi, wl_hist, wl_logdos, wl_f) result(acceptance)
      integer(int16), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      class(tmmc_params), intent(in) :: wl_setup
      integer, intent(in) :: bins, bins_mpi
      real(real64), dimension(:), intent(in) :: bin_edges, bin_edges_mpi
      real(real64), dimension(:), intent(inout) :: wl_hist, wl_logdos
      real(real64) , intent(in) :: temp, wl_f

      integer, dimension(4) :: rdm1, rdm2
      real(real64) :: e_swapped, e_unswapped, delta_e
      integer :: acceptance, i, ibin, jbin, ibin_mpi, jbin_mpi

      ! Establish total energy before any moves
      e_unswapped = setup%full_energy(config)

      acceptance = 0.0_real64

      do i=1, wl_setup%mc_sweeps*setup%n_atoms
  
          ! Make one MC trial
          ! Generate random numbers
          rdm1 = setup%rdm_site()
          rdm2 = setup%rdm_site()

          call pair_swap(config, rdm1, rdm2)
  
          e_swapped = setup%full_energy(config)

          ibin = bin_index(e_unswapped, bin_edges, bins)
          jbin = bin_index(e_swapped, bin_edges, bins)

          ibin_mpi = bin_index(e_unswapped, bin_edges_mpi, bins_mpi)
          jbin_mpi = bin_index(e_swapped, bin_edges_mpi, bins_mpi)

          ! Only compute energy change if within limits where V is defined
          if (jbin > 0 .and. jbin < bins+1) then
              ! Only compute move if within mpi bin energy range
              if (jbin_mpi > 0 .and. jbin_mpi < bins_mpi+1) then
                  ! Add change in V into diff_energy
                  delta_e = e_swapped - e_unswapped

                  ! Accept or reject move
                  if (genrand() .lt. m.exp(wl_logdos[ibin]-wl_logdos[jbin])) then
                      acceptance = acceptance + 1
                      e_unswapped = e_swapped
                  else
                      call pair_swap(config, rdm1, rdm2)
                  end if
              else
                  call pair_swap(config, rdm1, rdm2)
              end if
              wl_hist(ibin) = wl_hist(ibin) + 1
              wl_logdos(jbin) = wl_logdos(jbin) + wl_f
          else
              ! reject and reset
              call pair_swap(config, rdm1, rdm2)
          end if
      end do

  end function run_wl_sweeps

  subroutine wl_burn_in(setup, wl_setup, config, temp, bin_edges, bins)
      integer(int16), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      class(tmmc_params), intent(in) :: wl_setup
      real(real64) , intent(in) :: temp

      integer, dimension(4) :: rdm1, rdm2
      real(real64) :: e_swapped, e_unswapped, delta_e, beta
      integer :: i, ibin, jbin

      ! Establish total energy before any moves
      e_unswapped = setup%full_energy(config)

      do i=1, wl_setup%burn_in_sweeps*setup%n_atoms
          ! Make one MC trial
          ! Generate random numbers
          rdm1 = setup%rdm_site()
          rdm2 = setup%rdm_site()

          call pair_swap(config, rdm1, rdm2)
  
          e_swapped = setup%full_energy(config)

          ibin = bin_index(e_unswapped, bin_edges, bins)
          jbin = bin_index(e_swapped, bin_edges, bins)

          delta_e = e_swapped - e_unswapped

          ! Accept or reject move
          if (delta_e < 0) then
              e_unswapped = e_swapped
          else 
            call pair_swap(config, rdm1, rdm2)
          end if
      end do
  end subroutine wl_burn_in

end module wang_landau
