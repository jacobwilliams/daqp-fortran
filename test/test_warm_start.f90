!*****************************************************************************************
!>
!  Warm and hot starts: `update` with new data gives the same as a fresh
!  setup (in fewer iterations), a warm start from the optimal working set
!  finishes in at most one iteration, a cold start repeats a fresh solve, and a
!  copied object solves independently.

    program test_warm_start

    use daqp_module
    use daqp_test_utils

    implicit none

    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), &
                             x2(:), lam2(:), f0(:), df(:), df2(:), db(:), bu0(:), bl0(:), H2(:,:), A2(:,:)
    integer(ip), allocatable :: active(:)
    logical, allocatable :: at_lower(:)
    type(daqp_type) :: qp, qp2, fresh
    integer :: k, i, n, m, ms, nact, iter_hot, iter_cold, iter0
    integer(ip) :: istat
    real(wp) :: tol

    call rng_seed(161803)
    tol = max(1.0e-8_wp, 1.0e3_wp*sqrt(epsilon(1.0_wp)))
    n = 30; m = 90; ms = 15; nact = 20
    call generate_qp(n, m, ms, nact, 1.0e2_wp, xref, H, f, A, bu, bl)
    allocate(x(n), lam(m), x2(n), lam2(m))

    ! --- warm start from the optimal working set: at most one iteration
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x, lam, istat)
    iter0 = qp%iter
    call qp%get_working_set(active, at_lower)
    call check(size(active) == nact, 'working set size')
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x2, lam2, istat, active=active, at_lower=at_lower)
    call check(istat == daqp_optimal, 'warm start: exit flag')
    call check(qp%iter <= 1, 'warm start: iterations')
    call check(maxval(abs(x-x2)) < tol .and. maxval(abs(lam-lam2)) < tol*(1.0_wp+maxval(abs(lam))), &
               'warm start: solution')

    ! --- cold start: the same as a fresh solve
    call qp%solve(x2, lam2, istat, cold=.true.)
    call check(qp%iter == iter0, 'cold start: iterations')
    call check(maxval(abs(x-x2)) < tol, 'cold start: solution')

    ! --- an MPC-like sequence: update f and the bounds; hot starts agree with
    !     fresh setups, in fewer iterations
    f0 = f; bu0 = bu; bl0 = bl
    df = [(randn(), i=1,n)]
    ! the bounds move with a translation of the feasible set (which stays feasible)
    df2 = [(0.1_wp*randn(), i=1,n)]
    db = [df2(1:ms), matmul(A, df2)]
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x, lam, istat)
    iter_hot = 0
    iter_cold = 0
    do k = 1, 20
        f = f0 + 0.05_wp*real(k,wp)*df
        bu = bu0 + 0.01_wp*real(k,wp)*db
        bl = bl0 + 0.01_wp*real(k,wp)*db
        call qp%update(istat, f=f, bupper=bu, blower=bl)
        call check(istat == daqp_success, 'update: status')
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal, 'hot start: exit flag')
        iter_hot = iter_hot + qp%iter
        call fresh%setup(H, f, A, bu, bl, istat)
        call fresh%solve(x2, lam2, istat)
        call check(istat == daqp_optimal, 'fresh: exit flag')
        iter_cold = iter_cold + fresh%iter
        call check(maxval(abs(x-x2)) < tol, 'hot start: solution')
        call check(maxval(abs(lam-lam2)) < tol*(1.0_wp+maxval(abs(lam2))), 'hot start: multipliers')
    end do
    call check(iter_hot < iter_cold, 'hot starts take fewer iterations')
    write(*,'(A,I0,A,I0)') 'iterations (hot, cold): ', iter_hot, ', ', iter_cold

    ! --- update H and A: the same as a fresh setup
    call generate_qp(n, m, ms, nact, 1.0e2_wp, xref, H2, f, A2, bu, bl)
    call qp%update(istat, H=H2, f=f, A=A2, bupper=bu, blower=bl)
    call check(istat == daqp_success, 'update H, A: status')
    call qp%solve(x, lam, istat)
    call check(istat == daqp_optimal, 'update H, A: exit flag')
    call check(maxval(abs(x-xref)) < tol*1.0e2_wp, 'update H, A: solution')

    ! --- a copied object solves independently
    call qp%setup(H, f0, A, bu0, bl0, istat)
    call qp%solve(x, lam, istat)
    qp2 = qp
    call qp2%update(istat, f=2.0_wp*f0)
    call qp2%solve(x2, lam2, istat)
    call check(istat == daqp_optimal, 'copy: exit flag')
    call qp%solve(x2, lam2, istat)
    call check(maxval(abs(x-x2)) < tol, 'copy: the original is unchanged')
    call qp2%destroy()
    call qp%solve(x2, lam2, istat)
    call check(istat == daqp_optimal .and. maxval(abs(x-x2)) < tol, 'copy: destroying the copy')

    ! --- invalid working sets
    call qp%set_working_set([1_ip, 1_ip], istat)
    call check(istat == daqp_invalid_input, 'duplicate in working set')
    call qp%set_working_set([0_ip], istat)
    call check(istat == daqp_invalid_input, 'index out of range')

    ! --- not set up
    call fresh%destroy()
    call fresh%solve(x, lam, istat)
    call check(istat == daqp_not_setup, 'not set up')

    call report('test_warm_start')

    end program test_warm_start
!*****************************************************************************************
