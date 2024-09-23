!----------------------------------------------------------------------!
! tmmc.f90                                                             !
!                                                                      !
! Module containing routines implementing transition-matrix Monte      !
! Carlo (TMMC).                                                        !
!                                                                      !
! H. J. Naguszewski,  Warwick                                     2024 !
!----------------------------------------------------------------------!
module tmmc

    use initialise
    use kinds
    use shared_data
    use c_functions
    use random_site
    use metropolis
    use mpi

    implicit none

    contains

    !------------------------------------------------------------------!
    ! Main TMMC routine.                                               !
    !                                                                  !
    ! H. J. Naguszewski,  Warwick                                 2024 !
    !------------------------------------------------------------------!
    subroutine tmmc_main(setup, my_rank)

        ! Rank of this processor
        integer, intent(in) :: my_rank
        integer :: ierror, num_proc

        ! Arrays for storing data
        type(run_params) :: setup
        type(tmmc_params) :: tmmc_setup
      
        ! Integers used in calculations
        integer :: i,j,k, ierr, accept, bins
        
        ! Temperature and temperature steps
        real(real64) :: temp, acceptance, beta
      
        ! tmmc variables and arrays
        real(real64) :: bin_width, energy_to_ry, target_energy
        real(real64), allocatable :: bin_edges(:), probability_dist(:), bin_probability(:), bin_probability_buffer(:), &
        energy_bias(:)
        real(real64), allocatable :: trans_matrix(:,:), norm_trans_matrix(:,:), trans_matrix_buffer(:,:), energy_bias_all(:,:)
        ! MPI variables
        integer :: bins_mpi, mpi_start_idx, mpi_end_idx, bin_overlap
        real(real64) :: mpi_width, start, end, reduce_time, bias_time, tmmc_time
        real(real64), allocatable :: bin_edges_mpi(:), energy_bias_mpi(:)

        ! Get number of MPI processes
        call MPI_COMM_SIZE(MPI_COMM_WORLD, num_proc, ierror)
        !tmmc_setup%mc_sweeps = tmmc_setup%mc_sweeps/num_proc

        ! Allocate arrays
        bins = tmmc_setup%bins
        bins_mpi = bins/num_proc
        if (my_rank == num_proc-1) then
            bins_mpi = bins_mpi + MOD(bins, bins_mpi*(num_proc))
        end if

        bin_overlap = tmmc_setup%bin_overlap
        mpi_start_idx = max(my_rank*(bins/num_proc) + 1 - bin_overlap, 1)
        mpi_end_idx = min(my_rank*(bins/num_proc) + bins_mpi + bin_overlap, bins)
        bins_mpi = mpi_end_idx - mpi_start_idx + 1

        allocate(bin_edges(bins+1))
        allocate(probability_dist(bins))
        allocate(bin_probability(bins))
        allocate(trans_matrix(bins,bins))
        allocate(norm_trans_matrix(bins,bins))
        allocate(energy_bias(bins))
        if (my_rank == 0) then
            allocate(trans_matrix_buffer(bins,bins))
            allocate(bin_probability_buffer(bins))
            bin_probability_buffer = 0.0_real64
            allocate(energy_bias_all(bins,tmmc_setup%weight_update))
        end if

        ! MPI arrays
        allocate(bin_edges_mpi(bins_mpi+1))
        allocate(energy_bias_mpi(bins_mpi))

        ! Set temperature
        temp = setup%T*k_b_in_Ry

        ! Conversion meV/atom to Rydberg
        energy_to_ry=setup%n_atoms/(eV_to_Ry*1000)

        !---------------------------------!
        ! Initialise tmmc arrays and bins !
        !---------------------------------!
        mpi_width = (tmmc_setup%energy_max - tmmc_setup%energy_min)/real(num_proc)*energy_to_ry
        target_energy = tmmc_setup%energy_min*energy_to_ry + mpi_width*(my_rank + 0.5)

        j = 1
        bin_width = (tmmc_setup%energy_max - tmmc_setup%energy_min)/real(tmmc_setup%bins)*energy_to_ry
        do i=1, bins+1
            bin_edges(i) = tmmc_setup%energy_min*energy_to_ry + (i-1)*bin_width
            if (i > mpi_start_idx-1 .and. i < mpi_end_idx+2) then
                bin_edges_mpi(j) = bin_edges(i)
                j = j + 1
            end if
        end do

        energy_bias = 0.0_real64; probability_dist = 0.0_real64
        trans_matrix = 0.0_real64; norm_trans_matrix = 0.0_real64
        bin_probability = 0.0_real64

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
        if (tmmc_setup%burn_in) then
            call tmmc_burn_in(setup, tmmc_setup, config, temp, target_energy, trans_matrix, bin_edges, bins)
        end if

        !print*, my_rank, minval(bin_edges_mpi), maxval(bin_edges_mpi), target_energy, &
        !setup%full_energy(config)-minval(bin_edges_mpi), setup%full_energy(config)-maxval(bin_edges_mpi)

        !print*, my_rank, bin_index(setup%full_energy(config), bin_edges_mpi, bins_mpi)

        call comms_wait()

        !--------------------!
        ! Target Temperature !
        !--------------------!
        do i=1, tmmc_setup%weight_update
            start = MPI_Wtime()
            bin_probability = 0.0_real64
            acceptance = run_tmmc_sweeps(setup, tmmc_setup, config, temp, bins, &
            bin_edges, bins_mpi, bin_edges_mpi, energy_bias, trans_matrix, bin_probability)
            call comms_wait()
            tmmc_time = MPI_Wtime()

            call MPI_REDUCE(trans_matrix, trans_matrix_buffer, bins*bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
            call comms_wait()
            reduce_time = MPI_Wtime()

            if (my_rank == 0) then
                trans_matrix = trans_matrix_buffer
                call bias_from_tm(energy_bias, probability_dist, norm_trans_matrix,trans_matrix, &
                bins, bin_edges, bin_width, beta, 1)
            end if
            call comms_wait()
            bias_time = MPI_Wtime()
            call MPI_Bcast(energy_bias, bins, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
            call comms_wait()

            end = MPI_Wtime()
            if(my_rank == 0) then
                energy_bias_all(:,i) = energy_bias
                write(6,'(a,f6.2,a,f6.2,a,f6.2,a,f6.2,a)',advance='yes') "TMMC: ", tmmc_time-start, &
                 "s | R:", reduce_time-tmmc_time, "s | Bias:", bias_time-reduce_time, "s | B:", end-bias_time, "s"
                write(6,'(a,i0,a,f6.2,a,f6.2,a)',advance='yes') "Weight Update ", i, ": Accepted ", &
                 (acceptance/(tmmc_setup%mc_sweeps*setup%n_atoms)*100.0), "% of Monte Carlo moves. Time taken: ", end-start, "s"
            end if
        end do
        
        !call MPI_REDUCE(trans_matrix, trans_matrix_buffer, bins*bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        call MPI_REDUCE(bin_probability, bin_probability_buffer, bins, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

        if (my_rank == 0) then
            call bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, trans_matrix, &
            bins, bin_edges, bin_width, temp, 0)
            ! Normalize bins visited array
            bin_probability = bin_probability_buffer
            bin_probability = bin_probability/sum(bin_probability)

            ! Write output files
            call ncdf_writer_1d("dos_bins.dat", ierr, bin_edges)

            call ncdf_writer_1d("dos_probability.dat", ierr, probability_dist)

            call ncdf_writer_1d("bin_probability.dat", ierr, bin_probability)

            call ncdf_writer_2d("energy_bias_all.dat", ierr, energy_bias_all)

            call ncdf_writer_2d("transition_matrix.dat", ierr, trans_matrix)
        end if
        
        if(my_rank == 0) then
            write(*, *)
            write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
        end if
    
    end subroutine tmmc_main    

    !------------------------------------------------------------------!
    ! Routine for obtaining index of particular bin.                   !
    !                                                                  !
    ! H. J. Naguszewski,  Warwick                                 2024 !
    !------------------------------------------------------------------!
    function bin_index(energy, bin_edges, bins) result(index)

        integer, intent(in) :: bins
        real(real64), intent(in) :: energy
        real(real64), dimension(:), intent(in) :: bin_edges
        real(real64) :: bin_range

        bin_range = bin_edges(bins+1) - bin_edges(1)
        index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) &
                + 1
    end function bin_index

    !------------------------------------------------------------------!
    ! Routine for constructing the bias from the transition matrix.    !
    !                                                                  !
    ! H. J. Naguszewski,  Warwick                                 2024 !
    !------------------------------------------------------------------!
    subroutine bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, &
        trans_matrix, bins, bin_edges, bin_width, beta, bias_min)
        real(real64), dimension(:), intent(inout) :: energy_bias, probability_dist
        real(real64), dimension(:), intent(in) :: bin_edges
        real(real64), dimension(:,:), intent(inout) :: norm_trans_matrix
        real(real64), dimension(:,:), intent(in) :: trans_matrix
        real(real64), intent(in) :: bin_width, beta
        integer, intent(in) :: bins, bias_min

        integer, parameter :: lwmax=8000
        integer :: info, lwork
        real(real64), dimension(bins) :: wr, wi
        real(real64), dimension(lwmax) :: work
        real(real64), dimension(bins, bins) :: work_tm, vl, vr
        logical, dimension(bins) :: max_mask

        integer :: i, ii
        real(real32) :: Pnorm, bin_energy, min_bias, mincount
        external :: dgeev

        max_mask = .True.

        ! zero normalized transition matrix
        norm_trans_matrix = 0.0_real64

        ! Compute as appropriately normalised collection matrix
        do i=1, bins
            Pnorm = sum(trans_matrix(:,i))
            do ii=1, bins
                if (Pnorm > 0.0_real64) then
                    norm_trans_matrix(ii,i) = trans_matrix(ii,i)/Pnorm
                end if
            end do
        end do

        !-------------------------------------------------------------!
        ! Find the dominant eigenvector and store as probability_dist !
        !-------------------------------------------------------------!
        work_tm = norm_trans_matrix

        ! query work space
        lwork = -1
        call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins,   &
                   vr, bins, work, lwork, info)

        lwork = min(lwmax, int(work(1)))

        ! solve eigenproblem
        call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins,   &
                   vr, bins, work, lwork, info)

        ! check convergence
        IF( info.gt.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF
        
        probability_dist = abs(vr(:,maxloc(wr,1, mask=max_mask)))

        ! In case there are bins for which we have no data (yet)
        ! replace zeros with minimum non-zero probability
        if (bias_min .eq. 1) then
            mincount = minval(probability_dist,                       &
                              MASK=(probability_dist > 0.0_real64))

            do i=1, bins
                probability_dist(i) = max(probability_dist(i),mincount)
            end do
        end if

        probability_dist = probability_dist/sum(probability_dist)

        ! Construct energy_bias function needed for uniform energy 
        ! sampling
        do i=1, bins
            bin_energy = bin_edges(i) + 0.5*bin_width
            energy_bias(i) = log(probability_dist(i))/beta
        end do
        ! Shift energy_bias so that minimum value is zero
        min_bias = minval(energy_bias)
        energy_bias = energy_bias - min_bias
    end subroutine bias_from_tm

    !------------------------------------------------------------------!
    ! Routine to run TMMC sweeps.                                      !
    !                                                                  !
    ! H. J. Naguszewski,  Warwick                                 2024 !
    !------------------------------------------------------------------!
    function run_tmmc_sweeps(setup, tmmc_setup, config, temp, bin_edges, energy_bias, &
        trans_matrix, bin_probability) result(acceptance)
        Implicit None
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        class(tmmc_params), intent(in) :: tmmc_setup
        integer, intent(in) :: bins, bins_mpi
        real(real64), dimension(:), intent(in) :: energy_bias, bin_edges, bin_edges_mpi
        real(real64), dimension(:), intent(inout) :: bin_probability
        real(real64) , intent(in) :: temp
        real(real64), dimension(:,:), intent(inout) :: trans_matrix

        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, delta_e, beta, unswapped_bias, swapped_bias
        integer :: acceptance, i, ibin, jbin, ibin_mpi, jbin_mpi

        ! Store inverse temp
        beta = 1.0_real64/temp

        ! Establish total energy before any moves
        e_unswapped = setup%full_energy(config)

        acceptance = 0.0_real64

        do i=1, tmmc_setup%mc_sweeps*setup%n_atoms
    
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
                ! Probability of staying in ibin, ignoring energy_bias
                trans_matrix(ibin,ibin) = trans_matrix(ibin,ibin) &
                     + 1.0_real64 - min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                ! Probability of moving to jbin, ignoring energy_bias
                trans_matrix(jbin,ibin) = trans_matrix(jbin,ibin)     &
                + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                ! Only compute move if within mpi bin energy range
                if (jbin_mpi > 0 .and. jbin_mpi < bins_mpi+1) then
                    ! Add change in V into diff_energy
                    unswapped_bias = energy_bias(ibin)
                    swapped_bias = energy_bias(jbin)         
                    delta_e = (e_swapped + swapped_bias) - (e_unswapped + unswapped_bias)

                    ! Accept or reject move
                    if (genrand() .lt. exp(-beta*delta_e)) then
                        bin_probability(jbin) = bin_probability(jbin) + 1.0_real64
                        acceptance = acceptance + 1
                        e_unswapped = e_swapped
                    else
                        bin_probability(ibin) = bin_probability(ibin) + 1.0_real64
                        call pair_swap(config, rdm1, rdm2)
                    end if
                else
                    call pair_swap(config, rdm1, rdm2)
                end if
            else
                ! reject and reset
                call pair_swap(config, rdm1, rdm2)
            end if
        end do

    end function run_tmmc_sweeps

    subroutine tmmc_burn_in(setup, tmmc_setup, config, temp, target_energy, trans_matrix, bin_edges, bins)
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        class(tmmc_params), intent(in) :: tmmc_setup
        real(real64) , intent(in) :: temp, target_energy
        real(real64), dimension(:,:), intent(inout) :: trans_matrix
        real(real64), dimension(:), intent(in) :: bin_edges
        integer, intent(in) :: bins

        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, delta_e, beta
        integer :: i, ibin, jbin

        ! Store inverse temp
        beta = 1.0_real64/temp

        ! Establish total energy before any moves
        e_unswapped = setup%full_energy(config)

        do i=1, tmmc_setup%burn_in_sweeps*setup%n_atoms
            if (abs(e_unswapped - target_energy) < 1.0e-6_real64) then
                exit
            end if

            ! Make one MC trial
            ! Generate random numbers
            rdm1 = setup%rdm_site()
            rdm2 = setup%rdm_site()

            call pair_swap(config, rdm1, rdm2)
    
            e_swapped = setup%full_energy(config)

            ibin = bin_index(e_unswapped, bin_edges, bins)
            jbin = bin_index(e_swapped, bin_edges, bins)

            ! Only compute energy change if within limits where V is defined
            if (jbin > 0 .and. jbin < bins+1) then
                trans_matrix(ibin,ibin) = trans_matrix(ibin,ibin) &
                     + 1.0_real64 - min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                trans_matrix(jbin,ibin) = trans_matrix(jbin,ibin) + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))
            end if

            delta_e = e_swapped - e_unswapped

            ! Accept or reject move
            if (e_swapped > target_energy .and. delta_e < 0) then
                e_unswapped = e_swapped
            else if (e_swapped < target_energy .and. delta_e > 0) then
                e_unswapped = e_swapped
            else
              call pair_swap(config, rdm1, rdm2)
            end if
        end do
    end subroutine tmmc_burn_in

end module tmmc
