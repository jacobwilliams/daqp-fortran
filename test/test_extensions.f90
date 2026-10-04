!*****************************************************************************************
!>
!  The extensions of DAQP, through the object-oriented interface (ports of
!  upstream's tests in `interfaces/daqp-julia/test/core_tests.jl` where they
!  exist): branch and bound, hierarchical QPs, AVIs, the elimination of
!  equalities, individual soft weights, a prefactorized Hessian, primal and
!  dual starts, the time limit, and redundant constraints.

    program test_extensions

    use daqp_module
    use daqp_test_utils

    implicit none

    real(wp) :: tol

    tol = max(1.0e-4_wp, 1.0e3_wp*sqrt(epsilon(1.0_wp)))

    call test_bnb()
    call test_hierarchical()
    call test_avi()
    call test_eq_elimination()
    call test_soft_weights()
    call test_factored()
    call test_starts()
    call test_time_limit()
    call test_minrep()
    call test_unconstrained_shortcut()

    call report('test_extensions')

    contains

    function eye(k) result(e)
        integer, intent(in) :: k
        real(wp) :: e(k,k)
        integer :: j
        e = 0.0_wp
        do j = 1, k
            e(j,j) = 1.0_wp
        end do
    end function eye

!*****************************************************************************************

    subroutine test_bnb()
        real(wp), allocatable :: H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), x1(:)
        integer(ip), allocatable :: sense(:)
        type(daqp_type) :: qp
        integer(ip) :: istat, iter1, i, k, nb
        real(wp) :: x3(3), l5(5), x8(8), l8(8), f1, tb

        tb = max(1.0e-5_wp, 1.0e3_wp*epsilon(1.0_wp))

        call rng_seed(11)
        nb = 10
        do k = 1, 3
            call generate_miqp(20, 60, 20, int(nb), H, f, A, bu, bl, sense)
            if (allocated(x)) deallocate(x, lam)
            allocate(x(20), lam(60))
            call qp%setup(H, f, A, bu, bl, istat, sense=sense)
            call check(istat == daqp_success, 'bnb: setup')
            call qp%solve(x, lam, istat)
            call check(istat == daqp_optimal, 'bnb: exit flag')
            call check(all(abs(x(1:nb)) < tb .or. abs(x(1:nb)-1.0_wp) < tb), &
                       'bnb: binary feasible')
            call check(sum(x(1:nb)) <= real(nb/2,wp) + tb, 'bnb: cardinality')
            ! the second solve is warm started at the root: fewer iterations
            iter1 = qp%iter
            f1 = qp%fval
            x1 = x
            call qp%solve(x, lam, istat)
            call check(istat == daqp_optimal .and. abs(qp%fval-f1) < tb*(1.0_wp+abs(f1)), &
                       'bnb: repeated solve')
            call check(qp%iter < iter1, 'bnb: the root is warm started')
            ! an optimal incumbent
            call qp%setup(H, f, A, bu, bl, istat, sense=sense, primal_start=x1)
            call qp%solve(x, lam, istat)
            call check(istat == daqp_optimal .and. maxval(abs(x-x1)) < tb, 'bnb: incumbent')
        end do

        ! upstream's small example: x = (0,1,1)
        call qp%setup(reshape([1.0_wp,0.5_wp,0.0_wp, 0.5_wp,1.0_wp,0.5_wp, 0.0_wp,0.5_wp,1.0_wp],[3,3]), &
                      [1.0_wp,0.0_wp,0.0_wp], reshape([1.0_wp,1.0_wp, 2.0_wp,1.0_wp, 3.0_wp,0.0_wp],[2,3]), &
                      [1.0_wp,1.0_wp,1.0_wp,daqp_inf,daqp_inf], [0.0_wp,0.0_wp,0.0_wp,4.0_wp,1.0_wp], &
                      istat, sense=[daqp_binary,daqp_binary,daqp_binary,0_ip,0_ip])
        call qp%solve(x3, l5, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x3-[0.0_wp,1.0_wp,1.0_wp])) < tol, &
                   'bnb: small example')

        ! a binary constraint at a zero-dual endpoint does not create branches
        call qp%setup(eye(8), [(0.0_wp, i=1,8)], bupper=[(1.0_wp, i=1,8)], blower=[(0.0_wp, i=1,8)], &
                      istat=istat, sense=[(daqp_binary, i=1,8)])
        call qp%solve(x8, l8, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x8)) < tol .and. qp%outer_iter == 1, &
                   'bnb: zero-dual endpoint')

        ! binary constraints need a Hessian
        call qp%setup(f=[1.0_wp], bupper=[1.0_wp], blower=[0.0_wp], istat=istat, sense=[daqp_binary])
        call check(istat == daqp_invalid_input, 'bnb: an LP is rejected')
    end subroutine test_bnb

!*****************************************************************************************

    subroutine test_hierarchical()
        type(daqp_type) :: qp
        integer(ip) :: istat, i, j, n, mh, ml, nl, m
        real(wp) :: x3(3), l6(6), x2(2), l3(3), l4(4), te
        real(wp), allocatable :: L(:,:), H(:,:), f(:), A(:,:), bu(:), bl(:), x0(:), x(:), lam(:), x1(:)

        te = max(1.0e-8_wp, 1.0e3_wp*epsilon(1.0_wp))

        ! upstream's example (no H: minimal norm): x = (1, 0.5, -1)
        call qp%setup(A=reshape([1.0_wp,1.0_wp,3.0_wp, 1.0_wp,-1.0_wp,1.0_wp, 1.0_wp,0.0_wp,-1.0_wp],[3,3]), &
                      bupper=[1.0_wp,1.0_wp,1.0_wp,1.0_wp,0.5_wp,20.0_wp], &
                      blower=[-1.0_wp,-1.0_wp,-1.0_wp,-daqp_inf,0.5_wp,10.0_wp], istat=istat, &
                      break_points=[3_ip,4_ip,5_ip,6_ip])
        call check(istat == daqp_success, 'hierarchical: setup')
        call qp%solve(x3, l6, istat)
        call check(istat > 0 .and. maxval(abs(x3-[1.0_wp,0.5_wp,-1.0_wp])) < tol, 'hierarchical: example')

        ! degenerate: conflicting equalities x1 = 4 and x1 = 8 in a soft level
        call qp%setup(A=reshape([1.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp],[3,2]), &
                      bupper=[4.0_wp,8.0_wp,1.0_wp], blower=[4.0_wp,8.0_wp,1.0_wp], istat=istat, &
                      break_points=[0_ip,2_ip,3_ip])
        call qp%solve(x2, l3, istat)
        call check(istat > 0 .and. maxval(abs(x2-[6.0_wp,1.0_wp])) < tol, 'hierarchical: degenerate')

        ! linearly dependent equalities in a soft level are resolved by slacks
        call qp%setup(eye(2), [0.0_wp,0.0_wp], reshape([1.0_wp,2.0_wp,1.0_wp,2.0_wp],[2,2]), &
                      [10.0_wp,10.0_wp,1.0_wp,4.0_wp], [-10.0_wp,-10.0_wp,1.0_wp,4.0_wp], istat, &
                      break_points=[2_ip,4_ip])
        call qp%solve(x2, l4, istat)
        call check(istat > 0 .and. maxval(abs(x2-0.75_wp)) < tol, 'hierarchical: dependent equalities')
        call qp%update(istat)
        call qp%solve(x2, l4, istat)
        call check(istat > 0 .and. maxval(abs(x2-0.75_wp)) < tol, 'hierarchical: after an update')

        ! a solve shifts the bounds of the soft levels: an update is needed before the next one
        call rng_seed(1)
        n = 10; mh = 3; ml = 3; nl = 2
        m = mh + nl*ml
        allocate(L(n,n), A(m,n), x(n), lam(m))
        do j = 1, n
            do i = 1, n
                L(i,j) = randn()
            end do
        end do
        H = matmul(transpose(L), L)/real(n,wp) + eye(n)
        f = [(randn(), i=1,n)]
        do j = 1, n
            do i = 1, m
                A(i,j) = randn()
            end do
        end do
        x0 = [(randn(), i=1,n)]
        bu = matmul(A, x0) + [(0.5_wp*urand(), i=1,m)]
        bl = matmul(A, x0) - [(0.5_wp*urand(), i=1,m)]
        x0 = [(randn(), i=1,n)]
        bu(mh+1:) = matmul(A(mh+1:,:), x0) + 0.05_wp
        bl(mh+1:) = bu(mh+1:) - 0.1_wp
        call qp%setup(H, f, A, bu, bl, istat, break_points=[(mh + i*ml, i=0,nl)])
        call qp%solve(x, lam, istat)
        call check(istat > 0, 'hierarchical: solve')
        x1 = x
        call qp%solve(x, lam, istat)
        call check(istat == daqp_unsupported, 'hierarchical: a second solve needs an update')
        call qp%update(istat)
        call qp%solve(x, lam, istat)
        call check(istat > 0 .and. maxval(abs(x-x1)) < te*(1.0_wp+maxval(abs(x1))), &
                   'hierarchical: re-solve')

        ! a hierarchy with a singular Hessian is unsupported
        H = 0.0_wp
        H(1,1) = 1.0_wp
        call qp%setup(H, f, A, bu, bl, istat, break_points=[(mh + i*ml, i=0,nl)])
        call check(istat == daqp_unsupported, 'hierarchical: singular H')
    end subroutine test_hierarchical

!*****************************************************************************************

    subroutine test_avi()
        real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), b(:), x(:), lam(:)
        type(daqp_type) :: qp
        type(daqp_result) :: ra
        integer(ip) :: istat, k, i
        real(wp) :: tol_avi

        tol_avi = max(1.0e-4_wp, 1.0e3_wp*sqrt(epsilon(1.0_wp)))*10.0_wp
        call rng_seed(21)
        do k = 1, 5
            call generate_avi(30, 100, xref, H, f, A, b)
            if (allocated(x)) deallocate(x, lam)
            allocate(x(30), lam(100))
            call qp%setup(H, f, A, b, [(-daqp_inf, i=1,100)], istat, is_avi=.true.)
            call check(istat == daqp_success, 'avi: setup')
            call qp%solve(x, lam, istat)
            call check(istat > 0, 'avi: exit flag')
            call check(maxval(abs(x-xref)) < tol_avi*(1.0_wp+maxval(abs(xref))), 'avi: solution')
            call daqp_avi(30, 100, 0_ip, H, f, b, [(-daqp_inf, i=1,100)], x, ra, A=A)
            call check(ra%exitflag > 0 .and. maxval(abs(x-xref)) < tol_avi*(1.0_wp+maxval(abs(xref))), &
                       'avi: daqp_avi')
        end do
    end subroutine test_avi

!*****************************************************************************************

    subroutine test_eq_elimination()
        type(daqp_type) :: qp, full
        real(wp), allocatable :: H(:,:), A(:,:), f(:), bu(:), bl(:), x(:), lam(:), x2(:), lam2(:), &
                                 xref(:), L(:,:)
        integer(ip), allocatable :: sense(:)
        type(daqp_result) :: r1, r2
        integer(ip) :: istat, i, j, n, neq, nineq, ifield, ifield2
        real(wp) :: te

        te = max(1.0e-8_wp, 1.0e3_wp*epsilon(1.0_wp))

        ! the elimination of equalities, forced: the same as without it
        n = 10; neq = 6
        H = eye(n)
        A = reshape([eye(n), [(1.0_wp, i=1,n)]], [n, n+1])
        A = transpose(A(:, [(i, i=1,neq), n+1]))
        bu = [(10.0_wp, i=1,n), (0.0_wp, i=1,neq), 10.0_wp]
        bl = [(-10.0_wp, i=1,n), (0.0_wp, i=1,neq), -10.0_wp]
        sense = [(0_ip, i=1,n), (daqp_equality, i=1,neq), 0_ip]
        allocate(x(n), lam(n+neq+1), x2(n), lam2(n+neq+1))
        qp%eq_reduction = daqp_eq_reduction_on
        call qp%setup(H, [(0.0_wp, i=1,n)], A, bu, bl, istat, sense=sense)
        call check(istat == daqp_success, 'eq: setup')
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x)) < tol, 'eq: solution')
        bu(n+neq+1) = 9.0_wp
        bl(n+neq+1) = -9.0_wp
        call qp%update(istat, bupper=bu, blower=bl)
        call check(istat == daqp_success, 'eq: update of a reduced problem')

        ! random problems with many equalities: the same as without the elimination
        call rng_seed(2024)
        n = 30; neq = 20; nineq = 40
        allocate(L(n,n))
        do j = 1, n
            do i = 1, n
                L(i,j) = randn()
            end do
        end do
        H = matmul(transpose(L), L)/real(n,wp) + eye(n)
        f = [(randn(), i=1,n)]
        deallocate(A)
        allocate(A(neq+nineq,n))
        do j = 1, n
            do i = 1, neq+nineq
                A(i,j) = randn()
            end do
        end do
        xref = [(randn(), i=1,n)]
        bu = [xref + [(urand(), i=1,n)], matmul(A(1:neq,:), xref), &
              matmul(A(neq+1:,:), xref) + [(0.3_wp*urand(), i=1,nineq)]]
        bl = [xref - [(urand(), i=1,n)], matmul(A(1:neq,:), xref), [(-daqp_inf, i=1,nineq)]]
        sense = [(0_ip, i=1,n), (daqp_equality, i=1,neq), (0_ip, i=1,nineq)]
        deallocate(x, lam, x2, lam2)
        allocate(x(n), lam(n+neq+nineq), x2(n), lam2(n+neq+nineq))
        call daqp_quadprog(n, n+neq+nineq, n, bu, bl, x2, r1, lam2, H=H, f=f, A=A, sense=sense)
        call check(r1%exitflag == daqp_optimal, 'eq: quadprog (auto)')
        full%eq_reduction = daqp_eq_reduction_off
        call full%setup(H, f, A, bu, bl, istat, sense=sense)
        call full%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-x2)) < te, 'eq: auto = off')
        call check(maxval(abs(lam-lam2)) < 1.0e2_wp*te*(1.0_wp+maxval(abs(lam))), 'eq: multipliers')
        call qp%setup(H, f, A, bu, bl, istat, sense=sense)
        call qp%solve(x, lam, istat)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-x2)) < te, 'eq: solved twice')
        call check(abs(qp%fval-r1%fval) < te*(1.0_wp+abs(r1%fval)), 'eq: objective')

        ! consecutive updates with a forced reduction
        n = 40; neq = 16
        do ifield = 1, 3
            do ifield2 = 1, 3, 2
                H = eye(n)
                deallocate(x, lam)
                allocate(x(n), lam(n))
                bu = [(0.0_wp, i=1,neq), (1.0_wp, i=1,n-neq)]
                bl = [(0.0_wp, i=1,neq), (-daqp_inf, i=1,n-neq)]
                sense = [(daqp_equality, i=1,neq), (0_ip, i=1,n-neq)]
                call qp%setup(H, [(-2.0_wp, i=1,n)], eye(n), bu, bl, istat, sense=sense)
                call qp%solve(x, lam, istat)
                do j = 1, 2
                    select case (merge(ifield, ifield2, j == 1))
                    case (1)
                        call qp%update(istat, H=H)
                    case (2)
                        call qp%update(istat, f=[(-2.0_wp, i=1,n)])
                    case (3)
                        call qp%update(istat, A=eye(n))
                    end select
                    call check(istat == daqp_success, 'eq: update')
                end do
                call qp%solve(x, lam, istat)
                call check(istat == daqp_optimal .and. maxval(abs(x-[(0.0_wp, i=1,neq), (1.0_wp, i=1,n-neq)])) &
                           < te .and. abs(qp%fval+36.0_wp) < 1.0e2_wp*te, 'eq: consecutive updates')
            end do
        end do
        call daqp_quadprog(2, 2, 2, [1.0_wp,1.0_wp], [1.0_wp,1.0_wp], x(1:2), r2, H=eye(2), f=[0.0_wp,0.0_wp])
        call check(r2%exitflag == daqp_optimal .and. maxval(abs(x(1:2)-1.0_wp)) < tol, 'eq: determined x')
    end subroutine test_eq_elimination

!*****************************************************************************************

    subroutine test_soft_weights()
        type(daqp_type) :: qp
        integer(ip) :: istat
        real(wp) :: x1(1), l1(1), ts

        ts = max(1.0e-9_wp, 1.0e3_wp*epsilon(1.0_wp))
        ! min 0.5x^2 - 10x, x <= 0 (soft): w*s + s^2/(2 rho)
        qp%rho_soft = 0.5_wp
        qp%w_soft = 2.0_wp
        call qp%setup(eye(1), [-10.0_wp], bupper=[0.0_wp], blower=[-daqp_inf], istat=istat, &
                      sense=[daqp_soft])
        call qp%solve(x1, l1, istat)
        call check(istat == daqp_soft_optimal .and. abs(x1(1)-8.0_wp/3.0_wp) < ts .and. &
                   abs(l1(1)-22.0_wp/3.0_wp) < ts, 'soft weights: uniform')
        qp%w_soft = 3.0_wp
        call qp%solve(x1, l1, istat)
        call check(abs(x1(1)-7.0_wp/3.0_wp) < ts, 'soft weights: w_soft = 3')
        qp%w_soft = 2.0_wp
        call qp%solve(x1, l1, istat)
        call check(abs(x1(1)-8.0_wp/3.0_wp) < ts, 'soft weights: w_soft = 2 again')
        ! individual weights (the first allocation invalidates the factorization)
        call qp%set_soft_weights(istat, rho_u=[0.1_wp])
        call check(istat == daqp_success, 'soft weights: set')
        call qp%solve(x1, l1, istat)
        call check(istat == daqp_soft_optimal .and. abs(x1(1)-8.0_wp/11.0_wp) < ts .and. &
                   abs(l1(1)-102.0_wp/11.0_wp) < ts, 'soft weights: individual')
        call qp%set_soft_weights(istat, rho_u=[0.2_wp])
        call qp%solve(x1, l1, istat)
        call check(abs(x1(1)-4.0_wp/3.0_wp) < ts .and. abs(l1(1)-26.0_wp/3.0_wp) < ts, &
                   'soft weights: individual, changed')
        call qp%set_soft_weights(istat, rho_u=[0.1_wp, 0.2_wp])
        call check(istat == daqp_invalid_input, 'soft weights: wrong size')
    end subroutine test_soft_weights

!*****************************************************************************************

    subroutine test_factored()
        real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), R(:)
        type(daqp_type) :: qp
        integer(ip) :: istat, i

        call rng_seed(31)
        call generate_qp(30, 90, 15, 20, 1.0e2_wp, xref, H, f, A, bu, bl)
        allocate(x(30), lam(90))
        R = cholesky_packed(H)
        call qp%setup(f=f, A=A, bupper=bu, blower=bl, istat=istat, R=R)
        call check(istat == daqp_success, 'factored: setup')
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'factored: solution')
        ! a factor certifies positive definiteness: even a positive eps_prox is not used
        qp%eps_prox = 1.0e-3_wp
        call qp%setup(f=f, A=A, bupper=bu, blower=bl, istat=istat, R=R)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'factored: positive eps_prox')
        call qp%set_defaults()
        ! diagonal
        H = 0.0_wp
        do i = 1, 30
            H(i,i) = 0.5_wp + urand()
        end do
        call qp%setup(H, f, A, bu, bl, istat)
        call qp%solve(xref, lam, istat)
        call qp%setup(f=f, A=A, bupper=bu, blower=bl, istat=istat, R=cholesky_packed(H))
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'factored: diagonal')
    end subroutine test_factored

!*****************************************************************************************

    subroutine test_starts()
        real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), lam0(:)
        type(daqp_type) :: qp
        integer(ip) :: istat
        real(wp) :: x2(2), l3(3)

        call rng_seed(41)
        call generate_qp(100, 500, 50, 80, 1.0e2_wp, xref, H, f, A, bu, bl)
        allocate(x(100), lam(500))
        ! primal start: one iteration
        call qp%setup(H, f, A, bu, bl, istat, primal_start=xref)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'primal start: solution')
        ! (in single precision, the generated solution is not accurate enough to find
        ! every active constraint)
        if (precision(1.0_wp) > 10) call check(qp%iter == 1, 'primal start: iterations')
        lam0 = lam
        ! dual start: one iteration
        call qp%setup(H, f, A, bu, bl, istat, dual_start=lam0)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'dual start: solution')
        if (precision(1.0_wp) > 10) call check(qp%iter == 1, 'dual start: iterations')
        ! recovery from a degenerate starting point
        call qp%setup(eye(2), [0.0_wp,0.0_wp], reshape([1.0_wp,1.0_wp],[1,2]), [1.0_wp,1.0_wp,2.0_wp], &
                      [-daqp_inf,-daqp_inf,-daqp_inf], istat, primal_start=[1.0_wp,1.0_wp])
        call qp%solve(x2, l3, istat)
        call check(istat == daqp_optimal .and. maxval(abs(x2)) < tol, 'primal start: degenerate')
    end subroutine test_starts

!*****************************************************************************************

    subroutine test_time_limit()
        real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), Q(:,:), &
                                 target(:), center(:), width(:)
        integer(ip), allocatable :: sense(:)
        type(daqp_type) :: qp
        integer(ip) :: istat, i, j, nt, nbt, mt

        call rng_seed(51)
        call generate_qp(100, 500, 50, 80, 1.0e2_wp, xref, H, f, A, bu, bl)
        allocate(x(100), lam(500))
        qp%time_limit = 1.0e-9_wp
        call qp%setup(H, f, A, bu, bl, istat)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_time_limit, 'time limit: reached')
        qp%time_limit = 100.0_wp
        call qp%solve(x, lam, istat, cold=.true.)
        call check(istat == daqp_optimal .and. maxval(abs(x-xref)) < tol, 'time limit: generous')
        ! branch and bound: the limit is also checked across the tree
        nt = 30; nbt = 14; mt = 6
        allocate(Q(nt,nt))
        do j = 1, nt
            do i = 1, nt
                Q(i,j) = randn()
            end do
        end do
        H = matmul(transpose(Q), Q)/real(nt,wp) + 0.2_wp*eye(nt)
        target = [(0.15_wp + 0.7_wp*urand(), i=1,nt)]
        f = -matmul(H, target)
        deallocate(A)
        allocate(A(mt,nt))
        do i = 1, mt
            A(i,1:nbt) = [(0.2_wp + urand(), j=1,nbt)]
            A(i,nbt+1:) = [(0.1_wp*randn(), j=nbt+1,nt)]
        end do
        center = matmul(A, target)
        width = [(0.15_wp + 0.15_wp*urand(), i=1,mt)]
        bu = [(1.0_wp, i=1,nbt), (2.0_wp, i=nbt+1,nt), center + width]
        bl = [(0.0_wp, i=1,nbt), (-2.0_wp, i=nbt+1,nt), center - width]
        sense = [(daqp_binary, i=1,nbt), (0_ip, i=nbt+1,nt+mt)]
        deallocate(x, lam)
        allocate(x(nt), lam(nt+mt))
        qp%time_limit = 1.0e-9_wp
        call qp%setup(H, f, A, bu, bl, istat, sense=sense)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_time_limit .and. qp%outer_iter <= 32, 'time limit: branch and bound')
    end subroutine test_time_limit

!*****************************************************************************************

    subroutine test_minrep()
        integer(ip) :: red(5)
        ! the box |x| <= 1 (as general constraints) and x1 + x2 <= 5 (redundant)
        call daqp_minrep(reshape([1.0_wp,-1.0_wp,0.0_wp,0.0_wp,1.0_wp, 0.0_wp,0.0_wp,1.0_wp,-1.0_wp,1.0_wp], &
                                 [5,2]), [1.0_wp,1.0_wp,1.0_wp,1.0_wp,5.0_wp], 0_ip, red)
        call check(all(red == [0_ip,0_ip,0_ip,0_ip,1_ip]), 'minrep')
        call check(daqp_first_violating([0.5_wp,0.5_wp], reshape([1.0_wp,1.0_wp],[1,2]), &
                                        [1.0_wp,1.0_wp,0.9_wp], [0.0_wp,0.0_wp,0.0_wp], 2_ip, &
                                        1.0e-9_wp) == 3, 'first violating')
    end subroutine test_minrep

!*****************************************************************************************

    subroutine test_unconstrained_shortcut()
        real(wp), allocatable :: A(:,:), b(:)
        real(wp) :: x(10)
        type(daqp_result) :: r
        integer :: i, j
        call rng_seed(61)
        allocate(A(10000,10), b(10000))
        do j = 1, 10
            do i = 1, 10000
                A(i,j) = randn()
            end do
        end do
        do i = 1, 10000
            b(i) = urand()
        end do
        call daqp_quadprog(10, 10000, 0, b, [(-daqp_inf, i=1,10000)], x, r, H=eye(10), &
                           f=[(0.0_wp, i=1,10)], A=A)
        call check(r%exitflag == daqp_optimal .and. maxval(abs(x)) < tol .and. r%iter == 1, &
                   'unconstrained shortcut')
    end subroutine test_unconstrained_shortcut

    end program test_extensions
!*****************************************************************************************
