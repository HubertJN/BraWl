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

    subroutine tmmc_main(setup, tmmc_setup, my_rank) 
        ! Rank of this processor
        integer, intent(in) :: my_rank
        integer :: ierror

        ! Arrays for storing data
        type(run_params) :: setup
        type(tmmc_params) :: tmmc_setup
      
        ! Integers used in calculations
        integer :: i,j,k, ierr, accept, bias_min, bins
        
        ! Temperature and temperature steps
        real(real64) :: temp, acceptance, beta
      
        ! tmmc variables and arrays
        real(real64) :: bin_width, energy_to_ry, target_energy
        real(real64), allocatable :: bin_edges(:), probability_dist(:), bin_probability(:), energy_bias(:)
        real(real64), allocatable :: trans_matrix(:,:), norm_trans_matrix(:,:)
    
        ! Allocate arrays and find number of bins
        call tmmc_allocate_arrays(tmmc_setup, my_rank, bins, bin_edges, probability_dist, &
        bin_probability, trans_matrix, norm_trans_matrix, energy_bias)

        ! Set temperature
        temp = setup%T*k_b_in_Ry

        ! Conversion meV/atom to Rydberg
        energy_to_ry=setup%n_atoms/(eV_to_Ry*1000)

        !---------------------------------!
        ! Initialise tmmc arrays and bins !
        !---------------------------------!
        call tmmc_initialise_bins(tmmc_setup, energy_to_ry, my_rank, bins, bin_width, bin_edges, target_energy)
    
        energy_bias = 0.0_real64; probability_dist = 0.0_real64
        trans_matrix = 0.0_real64; norm_trans_matrix = 0.0_real64
        bin_probability = 0.0_real64
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
        if (tmmc_setup%burn_in) then
            call tmmc_burn_in(setup, tmmc_setup, config, temp, target_energy)
        end if

        !--------------------!
        ! Target Temperature !
        !--------------------!
        bias_min = 1 ! determines whether empty probability_dist entries are filled with lowest probability
        do i=1, tmmc_setup%weight_update
            if (i .eq. tmmc_setup%weight_update) then
                bias_min = 0
            end if
            acceptance = run_tmmc_sweeps(setup, tmmc_setup, config, temp, bins, &
            bin_edges, energy_bias, trans_matrix, bin_probability)
            
            call bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, trans_matrix, &
            bins, bin_edges, bin_width, temp, bias_min)

            if(my_rank == 0) then
                write(6,'(a,i0,a,f6.2,a)',advance='yes') "Weight Update ", i, ": Accepted ", &
                 (acceptance/(tmmc_setup%mc_sweeps*setup%n_atoms)*100.0), "% of Monte Carlo moves"
            end if
        end do

        
        ! Normalize bins visited array
        bin_probability = bin_probability/sum(bin_probability)

        ! Sync MPI threads
        call comms_wait()

        ! Merge transistion matrix and bins visited
        if (my_rank /= 0) then
            call MPI_Send(probability_dist, bins, mpi_real, 0, 0, MPI_COMM_WORLD, ierror)
        end if

        if (my_rank == 0) then
            call MPI_Recv(probability_dist, bins, mpi_real, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        end if

        ! Write output files
        !call ncdf_writer_1d("dos_bins.dat", ierr, bin_edges)

        !call ncdf_writer_1d("dos_probability.dat", ierr, probability_dist)

        !call ncdf_writer_1d("bin_probability.dat", ierr, bin_probability)

        !call ncdf_writer_2d("trans_matrix.dat", ierr, trans_matrix)
        
        if(my_rank == 0) then
            write(*, *)
            write(6,'(25("-"),x,"Simulation Complete!",x,25("-"))')
        end if
    
    end subroutine tmmc_main    

    integer function bin_index(energy, bin_edges, bins) result(index)
        integer, intent(in) :: bins
        real(real64), intent(in) :: energy
        real(real64), dimension(:), intent(in) :: bin_edges
        real(real64) :: bin_range

        bin_range = bin_edges(bins+1) - bin_edges(1)
        index = int(((energy - bin_edges(1))/(bin_range))*real(bins)) + 1
    end function bin_index

    subroutine bias_from_tm(energy_bias, probability_dist, norm_trans_matrix, &
        trans_matrix, bins, bin_edges, bin_width, temp, bias_min)
        real(real64), dimension(:), intent(inout) :: energy_bias, probability_dist
        real(real64), dimension(:), intent(in) :: bin_edges
        real(real64), dimension(:,:), intent(inout) :: norm_trans_matrix
        real(real64), dimension(:,:), intent(in) :: trans_matrix
        real(real64), intent(in) :: bin_width, temp
        integer, intent(in) :: bins, bias_min

        integer, parameter :: lwmax=5000
        integer :: info, lwork
        real(real64), dimension(bins) :: wr, wi
        real(real64), dimension(lwmax) :: work
        real(real64), dimension(bins, bins) :: work_tm, vl, vr
        logical, dimension(bins) :: max_mask

        integer :: i, ii
        real(real64) :: Pnorm, bin_energy, min_bias, mincount
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

        !--------------------------------------------------!
        ! Find the dominant eigenvector and store as probability_dist !
        !--------------------------------------------------!
        work_tm = norm_trans_matrix

        ! query work space
        lwork = -1
        call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins, vr, bins, work, lwork, info)
        lwork = min(lwmax, int(work(1)))

        ! solve eigenproblem
        call dgeev('V', 'V', bins, work_tm, bins, wr, wi, vl, bins, vr, bins, work, lwork, info)

        ! check convergence
        IF( info.gt.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF
        
        probability_dist = abs(vr(:,maxloc(wr,1, mask=max_mask)))

        ! In case there are bins for which we have no data (yet)
        ! replace zeros with minimum non-zero probability
        if (bias_min .eq. 1) then
            mincount = minval(probability_dist, MASK=(probability_dist > 0.0_real64))
            do i=1, bins
                probability_dist(i) = max(probability_dist(i),mincount)
            end do
        end if
        probability_dist = probability_dist/sum(probability_dist)
        !--------------------------------------------------!

        ! Construct energy_bias function needed for uniform energy sampling
        do i=1, bins
            bin_energy = bin_edges(i) + 0.5*bin_width
            energy_bias(i) = temp*log(probability_dist(i))
        end do
        ! Shift energy_bias so that minimum value is zero
        min_bias = minval(energy_bias)
        energy_bias = energy_bias - min_bias
    end subroutine bias_from_tm

    function run_tmmc_sweeps(setup, tmmc_setup, config, temp, bins, bin_edges, energy_bias, &
        trans_matrix, bin_probability) result(acceptance)
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        class(tmmc_params), intent(in) :: tmmc_setup
        integer, intent(in) :: bins
        real(real64), dimension(:), intent(in) :: energy_bias, bin_edges
        real(real64), dimension(:), intent(inout) :: bin_probability
        real(real64) , intent(in) :: temp
        real(real64), dimension(:,:), intent(inout) :: trans_matrix

        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, delta_e, beta, unswapped_bias, swapped_bias
        integer :: acceptance, i, ibin, jbin

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

            ! Only compute energy change if within limits where V is defined
            if (jbin > 0 .and. jbin < bins+1) then
                bin_probability(ibin) = bin_probability(ibin) + 1.0_real64

                trans_matrix(ibin,ibin) = trans_matrix(ibin,ibin) &
                     + 1.0_real64 - min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                ! Probability of moving to jbin, ignoring energy_bias
                trans_matrix(jbin,ibin) = trans_matrix(jbin,ibin) + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                ! Add change in V into diff_energy
                unswapped_bias = energy_bias(ibin)
                swapped_bias = energy_bias(jbin)         
                delta_e = (e_swapped + swapped_bias) - (e_unswapped + unswapped_bias)

                ! Accept or reject move
                if (genrand() .lt. exp(-beta*delta_e)) then
                    acceptance = acceptance + 1
                    e_unswapped = e_swapped
                else
                  call pair_swap(config, rdm1, rdm2)
                end if
            else
                ! reject and reset
                call pair_swap(config, rdm1, rdm2)
            end if
        end do

    end function run_tmmc_sweeps

    subroutine tmmc_burn_in(setup, tmmc_setup, config, temp, target_energy)
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        class(tmmc_params), intent(in) :: tmmc_setup
        real(real64) , intent(in) :: temp, target_energy

        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, delta_e, beta
        integer :: i

        ! Store inverse temp
        beta = 1.0_real64/temp

        ! Establish total energy before any moves
        e_unswapped = setup%full_energy(config)

        do i=1, tmmc_setup%burn_in_sweeps*setup%n_atoms
    
            ! Make one MC trial
            ! Generate random numbers
            rdm1 = setup%rdm_site()
            rdm2 = setup%rdm_site()

            call pair_swap(config, rdm1, rdm2)
    
            e_swapped = setup%full_energy(config)

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

    subroutine tmmc_allocate_arrays(tmmc_setup, my_rank, bins, bin_edges, probability_dist, &
        bin_probability, trans_matrix, norm_trans_matrix, energy_bias)
        class(tmmc_params), intent(in) :: tmmc_setup
        real(real64), allocatable, intent(inout) :: bin_edges(:), probability_dist(:), bin_probability(:), energy_bias(:)
        real(real64), allocatable, intent(inout) :: trans_matrix(:,:), norm_trans_matrix(:,:)
        integer, intent(in) :: my_rank
        integer, intent(inout) :: bins
        real(real64) :: rank_bins_real

        select case (tmmc_setup%use_mpi)
        case (.true.)
            rank_bins_real = tmmc_setup%bins/tmmc_setup%mpi_processes

            if (my_rank == 0) then
                rank_bins_real = rank_bins_real*(1+tmmc_setup%percent_overlap)
            else if (my_rank == tmmc_setup%mpi_processes-1) then
                rank_bins_real = rank_bins_real*(1+tmmc_setup%percent_overlap)
            else
                rank_bins_real = rank_bins_real*(1+tmmc_setup%percent_overlap*2)
            end if
            
            bins = int(rank_bins_real)
        case (.false.)
            bins = tmmc_setup%bins
        end select

        allocate(bin_edges(bins+1))
        allocate(probability_dist(bins))
        allocate(bin_probability(bins))
        allocate(trans_matrix(bins,bins))
        allocate(norm_trans_matrix(bins,bins))
        allocate(energy_bias(bins))
    end subroutine tmmc_allocate_arrays

    subroutine tmmc_initialise_bins(tmmc_setup, energy_to_ry, my_rank, bins, bin_width, bin_edges, target_energy)
        class(tmmc_params), intent(in) :: tmmc_setup
        integer, intent(in) :: my_rank, bins
        real(real64), intent(in) :: energy_to_ry
        real(real64), intent(inout) :: target_energy, bin_width
        real(real64), intent(inout), dimension(:) :: bin_edges
        real(real64) :: mpi_width
        integer :: i

        bin_width = (tmmc_setup%energy_max - tmmc_setup%energy_min)/real(tmmc_setup%bins)*energy_to_ry

        select case (tmmc_setup%use_mpi)
        case (.true.)
            mpi_width = (tmmc_setup%energy_max - tmmc_setup%energy_min)/real(tmmc_setup%mpi_processes)*energy_to_ry
            target_energy = tmmc_setup%energy_min*energy_to_ry + mpi_width*(my_rank + 0.5)

            select case (my_rank)
            case (0)
                do i=1, bins+1
                    bin_edges(i) = tmmc_setup%energy_min*energy_to_ry + (i-1)*bin_width
                end do
            case default
                do i=1, bins+1
                    bin_edges(i) = tmmc_setup%energy_min*energy_to_ry &
                    +mpi_width*(real(my_rank)-tmmc_setup%percent_overlap) + (i-1)*bin_width
                end do
            end select

        case (.false.)
            target_energy = (tmmc_setup%energy_min + (tmmc_setup%energy_max - tmmc_setup%energy_min)/2)*energy_to_ry

            do i=1, tmmc_setup%bins+1
                bin_edges(i) = tmmc_setup%energy_min*energy_to_ry + (i-1)*bin_width
            end do
        end select
    end subroutine tmmc_initialise_bins

    subroutine tmmc_merge()
        integer :: ierror

        real(real64), allocatable :: probability_dist_buffer(:), probability_dist(:)
        allocate(probability_dist_buffer(52))

        call MPI_Send(probability_dist, 52, mpi_real, 0, 0, MPI_COMM_WORLD, ierror)
        call MPI_Recv(probability_dist_buffer, 52, mpi_real, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)
        
    end subroutine tmmc_merge
end module tmmc
