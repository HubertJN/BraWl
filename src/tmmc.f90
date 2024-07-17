module tmmc
    use kinds
    use shared_data
    use c_functions

    implicit none

    contains

    subroutine tmmc_main
        !---------------------!
        ! setup sampling bins !
        !---------------------!
        integer, parameter :: bins=50
        integer :: i
        real(real64), dimension(2) :: energy_range
        real(real64), dimension((bins+1)) :: bin_edges
        real(real64) :: bin_width, bin_range, temp
        real(real64), dimension(bins) :: bias, histogram, gutter, statP
        real(real64), dimension(bins, bins) :: tm, norm_tm
    
        energy_range(1) = -6.1_real64
        energy_range(2) = 4.0_real64
    
        temp = 601
    
        bin_width = (energy_range(2) - energy_range(1))/real(bins)
    
        do i=1, bins+1
            bin_edges(i) = energy_range(1) + (i-1)*bin_width
        end do
    
        bin_range = bin_edges(bins) - bin_edges(1)
    
        bias = 0.0_real64; histogram = 0.0_real64; gutter = 0.0_real64; statP = 0.0_real64
        tm = 0.0_real64; norm_tm = 0.0_real64
    
        call bias_from_tm(bias, statP, norm_tm, tm, gutter, bins, bin_edges, bin_width, temp)
    
    end subroutine tmmc_main    

    real function bin_index(energy, bin_edges, bins)
        integer, intent(in) :: bins
        real(real64), intent(in) :: energy
        real(real64), dimension(:), intent(in) :: bin_edges
        real(real64) :: bin_range

        bin_range = bin_edges(bins) - bin_edges(1)
        bin_index = int(((energy - bin_edges(1))/(bin_range))*real(bins))
        return
    end function bin_index

    subroutine bias_from_tm(bias, statP, norm_tm, tm, gutter, bins, bin_edges, bin_width, temp)
        real(real64), dimension(:), intent(inout) :: bias, statP
        real(real64), dimension(:), intent(in) :: gutter, bin_edges
        real(real64), dimension(:,:), intent(inout) :: norm_tm
        real(real64), dimension(:,:), intent(in) :: tm
        real(real64), intent(in) :: bin_width, temp
        integer, intent(in) :: bins

        integer, parameter :: lwmax=1000
        integer :: info, lwork
        real(real64), dimension(bins) :: wr, wi
        real(real64), dimension(lwmax) :: work
        real(real64), dimension(bins, bins) :: work_tm, vl, vr

        integer :: i, ii
        real(real64) :: Pnorm, bin_energy, min_bias, mincount
        external :: dgeev

        ! zero normalized transition matrix
        norm_tm = 0.0_real64

        ! Compute as appropriately normalised collection matrix
        do i=1, bins
            Pnorm = sum(tm(:,i)) + gutter(i)
            do ii=1, bins
                if (Pnorm > 0.0_real64) then
                    norm_tm(ii,i) = tm(ii,i)/Pnorm
                end if
            end do
        end do

        !--------------------------------------------------!
        ! Find the dominant eigenvector and store as statP !
        !--------------------------------------------------!
        work_tm = norm_tm

        ! query work space
        lwork = -1
        call dgeev('N', 'V', bins, work_tm, bins, wr, wi, vl, bins, vr, bins, work, lwork, info)
        lwork = min(lwmax, int(work(1)))

        ! solve eigenproblem
        call dgeev('N', 'V', bins, work_tm, bins, wr, wi, vl, bins, vr, bins, work, lwork, info)

        ! check convergence
        IF( info.gt.0 ) THEN
            WRITE(*,*)'The algorithm failed to compute eigenvalues.'
            STOP
        END IF
        
        statP = abs(vr(:,maxloc(wr,1)))

        ! In case there are bins for which we have no data (yet)
        ! replace zeros with minimum non-zero probability
        mincount = minval(statP, MASK=(statP > 0.0_real64))
        do i=1, bins
            statP(i) = max(statP(i),mincount)
        end do
        statP = statP/sum(statP)
        !--------------------------------------------------!

        ! Construct bias function needed for uniform energy sampling
        do i=1, bins
            bin_energy = bin_edges(i) + 0.5*bin_width
            bias(i) = temp*log(statP(i))
        end do
        ! Shift bias so that minimum value is zero
        min_bias = minval(bias)
        bias = bias - min_bias
    end subroutine bias_from_tm

    function run_tmmc_sweeps(setup, config, temp, bins, bin_edges, bias, trans_matrix, gutter) result(acceptance)
        integer(int16), dimension(:,:,:,:) :: config
        class(run_params), intent(in) :: setup
        real(real64), dimension(:), intent(in) :: bias, bin_edges
        real(real64) , intent(in) :: temp
        integer, intent(in) :: bins

        real(real64), dimension(:), intent(inout) :: gutter
        real(real64), dimension(:,:), intent(inout) :: trans_matrix

        integer, dimension(4) :: rdm1, rdm2
        real(real64) :: e_swapped, e_unswapped, delta_e, beta, unswapped_bias, swapped_bias
        integer(int16) :: site1, site2
        integer :: acceptance, i, ibin, jbin

        ! Store inverse temp
        beta = 1.0_real64/temp

        ! Establish total energy before any moves
        e_unswapped = setup%full_energy(config)

        do i=1, setup%mc_steps
    
            ! Make one MC trial
            ! Generate random numbers
            rdm1 = setup%rdm_site()
            rdm2 = setup%rdm_site()

            call pair_swap(config, rdm1, rdm2)
    
            e_swapped = setup%full_energy(config)

            ibin = bin_index(e_unswapped, bin_edges, bins)
            jbin = bin_index(e_swapped, bin_edges, bins)

            trans_matrix(ibin,ibin) = trans_matrix(ibin,ibin) + 1.0_real64 - min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

            ! Only compute energy change if within limits where V is defined
            if (jbin > 0 .and. jbin < bins+1) then
    
                ! Probability of moving to jbin, ignoring bias
                trans_matrix(jbin,ibin) = trans_matrix(jbin,ibin) + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))

                ! Add change in V into diff_energy
                unswapped_bias = bias(ibin)
                swapped_bias = bias(jbin)         
                delta_e = (e_swapped + swapped_bias) - (e_unswapped + unswapped_bias)

                ! Accept or reject move
                if (genrand() .lt. exp(-beta*delta_e)) then
                    acceptance = acceptance + 1
                    e_unswapped = e_swapped
                else
                  call pair_swap(config, rdm1, rdm2)
                end if
            else
                ! Add the transition probability for leaving range of interest into gutter
                gutter(ibin) = gutter(ibin) + min(1.0_real64, exp(-beta*(e_swapped - e_unswapped)))
    
                ! reject and reset
                call pair_swap(config, rdm1, rdm2)
            end if
        end do

    end function run_tmmc_sweeps

end module tmmc