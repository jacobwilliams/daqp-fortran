!*****************************************************************************************
!>
!  Degenerate problems: duplicated rows, linearly dependent active rows, more
!  active constraints than variables at the solution, weakly active
!  constraints, ties, and simple bounds given as simple bounds or as general rows.

    program test_degenerate

    use daqp_module
    use daqp_test_utils

    implicit none

    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), &
                             x2(:), lam2(:), A2(:,:), xs(:), ax(:), mult(:)
    integer, allocatable :: dup(:)
    type(daqp_type) :: qp
    integer :: k, i, j, n, m, ms, nact, nd
    integer(ip) :: istat
    real(wp) :: tol, stat, prim, dual, comp

    call rng_seed(2718)
    tol = max(1.0e-6_wp, 1.0e3_wp*sqrt(epsilon(1.0_wp)))*10.0_wp

    ! --- duplicated rows give the same solution
    do k = 1, 10
        n = 15; m = 45; ms = 5; nact = 10
        call generate_qp(n, m, ms, nact, 1.0e2_wp, xref, H, f, A, bu, bl)
        nd = 8
        dup = [(rand_int(1, m-ms), i=1,nd)]
        A2 = reshape([transpose(A), transpose(A(dup,:))], [n, m-ms+nd])
        A2 = transpose(A2)
        allocate(x(n), lam(m+nd))
        call qp%setup(H, f, A2, [bu, bu(ms+dup)], [bl, bl(ms+dup)], istat)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal, 'duplicated rows: exit flag')
        call check(maxval(abs(x-xref)) < tol, 'duplicated rows: solution')
        call kkt_residuals(ms, [bu, bu(ms+dup)], [bl, bl(ms+dup)], x, lam, stat, prim, dual, comp, &
                           H=H, f=f, A=A2)
        call check(stat < tol .and. prim < tol .and. dual < tol .and. comp < tol, &
                   'duplicated rows: KKT')
        deallocate(x, lam)
    end do

    ! --- a vertex with more active rows than variables: some weakly active,
    !     some linearly dependent (sums of two others)
    do k = 1, 10
        n = 8 + k
        m = 3*n
        if (allocated(A)) deallocate(A, bu, bl)
        allocate(A(m,n), bu(m), bl(m), mult(m), x(n), lam(m))
        xs = [(randn(), i=1,n)]
        do i = 1, m
            if (i > 2*n .and. mod(i,2) == 0) then
                A(i,:) = A(i-2*n,:) + A(i-2*n+1,:)
            else
                A(i,:) = [(randn(), j=1,n)]
            end if
        end do
        mult = 0.0_wp
        do i = 1, n/2
            mult(i) = 0.1_wp + urand()
        end do
        H = reshape([((merge(1.0_wp, 0.0_wp, i == j), i=1,n), j=1,n)], [n,n])
        f = -xs - matmul(transpose(A), mult)
        ax = matmul(A, xs)
        do i = 1, m
            if (i <= 3*n/2 .or. i > 2*n) then
                bu(i) = ax(i)         ! active (weakly, if the multiplier is 0)
            else
                bu(i) = ax(i) + 0.5_wp
            end if
            bl(i) = -daqp_inf
        end do
        call qp%setup(H, f, A, bu, bl, istat)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal, 'over-determined vertex: exit flag')
        call check(maxval(abs(x-xs)) < tol, 'over-determined vertex: solution')
        call kkt_residuals(0, bu, bl, x, lam, stat, prim, dual, comp, H=H, f=f, A=A)
        call check(stat < tol .and. prim < tol .and. dual < tol .and. comp < tol, &
                   'over-determined vertex: KKT')
        deallocate(A, bu, bl, mult, x, lam)
    end do

    ! --- ties: the same violated constraint three times (equal violations)
    allocate(x(2), lam(3))
    call qp%setup(reshape([1.0_wp,0.0_wp,0.0_wp,1.0_wp],[2,2]), [-1.0_wp,-1.0_wp], &
                  reshape([1.0_wp,1.0_wp,1.0_wp, 1.0_wp,1.0_wp,1.0_wp],[3,2]), &
                  [1.0_wp,1.0_wp,1.0_wp], [(-daqp_inf, i=1,3)], istat)
    call qp%solve(x, lam, istat)
    call check(istat == daqp_optimal, 'ties: exit flag')
    call check(maxval(abs(x-0.5_wp)) < tol, 'ties: solution')
    call check(abs(sum(lam)-0.5_wp) < tol, 'ties: multipliers')
    deallocate(x, lam)

    ! --- simple bounds as simple bounds or as general rows
    do k = 1, 5
        n = 12; m = 30; ms = 12; nact = 8
        call generate_qp(n, m, ms, nact, 1.0e2_wp, xref, H, f, A, bu, bl)
        allocate(x(n), lam(m), x2(n), lam2(m))
        call qp%setup(H, f, A, bu, bl, istat)
        call qp%solve(x, lam, istat)
        A2 = reshape([((merge(1.0_wp, 0.0_wp, i == j), i=1,n), j=1,n)], [n,n])
        A2 = reshape([transpose(A2), transpose(A)], [n, m])
        A2 = transpose(A2)
        call qp%setup(H, f, A2, bu, bl, istat)
        call check(istat == daqp_success, 'bounds as rows: setup')
        call qp%solve(x2, lam2, istat)
        call check(istat == daqp_optimal, 'bounds as rows: exit flag')
        call check(maxval(abs(x-x2)) < tol .and. maxval(abs(x-xref)) < tol, 'bounds as rows: solution')
        call check(maxval(abs(lam-lam2)) < tol*(1.0_wp+maxval(abs(lam))), 'bounds as rows: multipliers')
        deallocate(x, lam, x2, lam2)
    end do

    call report('test_degenerate')

    end program test_degenerate
!*****************************************************************************************
