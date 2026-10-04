!*****************************************************************************************
!>
!  Random feasible QPs of sizes 2 to 500, checked by their KKT conditions
!  (stationarity, primal and dual feasibility, complementarity) and against
!  the known solution, with tolerances scaled by `epsilon` and the data.

    program test_kkt_random

    use daqp_module
    use daqp_test_utils

    implicit none

    integer, parameter :: ns(9) = [2, 3, 5, 10, 20, 50, 100, 200, 500]
    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:)
    real(wp) :: kappas(3), stat, prim, dual, comp, tol, scale
    type(daqp_type) :: qp
    integer :: in, im, ik, n, m, ms, nact
    integer(ip) :: istat
    character(len=128) :: label

    call rng_seed(31415)
    if (precision(1.0_wp) > 10) then
        kappas = [1.0_wp, 1.0e2_wp, 1.0e5_wp]
    else
        kappas = [1.0_wp, 1.0e1_wp, 1.0e2_wp]
    end if

    do in = 1, size(ns)
        n = ns(in)
        do im = 0, 2
            if (n == 500 .and. im > 1) cycle
            m = max(1, im*n + n/2)   ! about n/2, 1.5n, 2.5n
            ms = min(m, n)/3
            nact = min(m, n)/2
            do ik = 1, size(kappas)
                if (n == 500 .and. ik > 1) cycle
                call generate_qp(n, m, ms, nact, kappas(ik), xref, H, f, A, bu, bl)
                if (allocated(x)) deallocate(x, lam)
                allocate(x(n), lam(m))
                call qp%setup(H, f, A, bu, bl, istat)
                write(label,'("n=",I0,", m=",I0,", kappa=",ES8.1)') n, m, kappas(ik)
                call check(istat == daqp_success, 'setup: '//trim(label))
                call qp%solve(x, lam, istat)
                call check(istat == daqp_optimal, 'exit flag: '//trim(label))
                call kkt_residuals(ms, bu, bl, x, lam, stat, prim, dual, comp, H=H, f=f, A=A)
                scale = maxval(abs(xref)) + maxval(abs(f)) + maxval(abs(lam))
                tol = kkt_tolerance(kappas(ik), scale, qp%primal_tol)
                call check(stat <= tol*(1.0_wp + maxval(abs(H))), 'stationarity: '//trim(label))
                call check(prim <= tol, 'primal feasibility: '//trim(label))
                call check(dual <= tol, 'dual feasibility: '//trim(label))
                call check(comp <= tol*(1.0_wp + scale), 'complementarity: '//trim(label))
                call check(maxval(abs(x-xref)) <= tol*kappas(ik), 'solution: '//trim(label))
                if (stat > tol*(1.0_wp + maxval(abs(H))) .or. prim > tol .or. &
                    maxval(abs(x-xref)) > tol*kappas(ik)) &
                    write(*,'(A,5ES10.2)') trim(label)//': ', stat, prim, dual, comp, maxval(abs(x-xref))
            end do
        end do
    end do

    call report('test_kkt_random')

    end program test_kkt_random
!*****************************************************************************************
