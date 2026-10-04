!*****************************************************************************************
!>
!  Ports of upstream's core tests (`interfaces/daqp-julia/test/core_tests.jl`):
!  random QPs and LPs with a known solution, the iteration limit, one-sided
!  bounds, a zero-row equality, a trivially infeasible problem, an LP that
!  used to cycle, warm starts, the semi-proximal method, and soft constraints.

    program test_core_upstream

    use daqp_module
    use daqp_core, only: daqp_quadprog
    use daqp_test_utils

    implicit none

    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:), &
                             x1(:), bl_inf(:)
    integer(ip), allocatable :: sense(:), active(:)
    logical, allocatable :: at_lower(:)
    type(daqp_type) :: qp
    type(daqp_result) :: res
    type(daqp_settings) :: s
    integer :: k, j, n, m, ms, nact
    integer(ip) :: istat
    real(wp) :: tol, kappa

    call rng_seed(1234)
    tol = max(1.0e-4_wp, sqrt(epsilon(1.0_wp))*1.0e3_wp)

    ! --- quadprog: random QPs (n=100, m=500, ms=50, 80 active, cond 1e2)
    n = 100; m = 500; ms = 50; nact = 80; kappa = 1.0e2_wp
    allocate(x(n), lam(m))
    do k = 1, 20
        call generate_qp(n, m, ms, nact, kappa, xref, H, f, A, bu, bl)
        call daqp_quadprog(n, m, ms, bu, bl, x, res, lam, H=H, f=f, A=A)
        call check(res%exitflag == daqp_optimal, 'quadprog exit flag')
        call check(norm2(x-xref) < tol, 'quadprog solution')
    end do
    ! iteration limit
    s%iter_limit = 1
    call daqp_quadprog(n, m, ms, bu, bl, x, res, lam, H=H, f=f, A=A, settings=s)
    call check(res%exitflag == daqp_iteration_limit, 'iteration limit')

    ! --- one-sided bounds give the same as blower = -inf
    do k = 1, 5
        call generate_qp(n, m, ms, nact, kappa, xref, H, f, A, bu, bl)
        bl_inf = [(-daqp_inf, j=1, m)]
        call daqp_quadprog(n, m, ms, bu, bl_inf, xref, res, H=H, f=f, A=A)
        call qp%setup(H, f, A, bu, bl_inf, istat)
        call qp%solve(x, lam, istat)
        call check(istat == daqp_optimal, 'one-sided exit flag')
        call check(norm2(x-xref) < tol, 'one-sided solution')
    end do

    ! --- the object interface (setup + solve), and stationarity
    call generate_qp(n, m, ms, nact, kappa, xref, H, f, A, bu, bl)
    call qp%setup(H, f, A, bu, bl, istat)
    call check(istat == daqp_success, 'setup status')
    call qp%solve(x, lam, istat)
    call check(istat == daqp_optimal, 'model exit flag')
    call check(norm2(x-xref) < tol, 'model solution')
    call check(norm2(matmul(H,x) + f + [lam(1:ms), [(0.0_wp, j=ms+1,n)]] + &
                     matmul(transpose(A), lam(ms+1:m))) < tol, 'model stationarity')

    ! --- warm start from the optimal working set: one iteration
    call qp%get_working_set(active, at_lower)
    call generate_qp(n, m, ms, nact, kappa, xref, H, f, A, bu, bl) ! (a different QP...)
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x, lam, istat)
    call qp%get_working_set(active, at_lower)
    call qp%setup(H, f, A, bu, bl, istat)     ! (...then this one again, warm started)
    call qp%solve(x, lam, istat, active=active, at_lower=at_lower)
    call check(istat == daqp_optimal, 'warm start exit flag')
    call check(norm2(x-xref) < tol, 'warm start solution')
    call check(qp%iter == 1, 'warm start iterations')

    ! --- LPs (via the proximal loop)
    do k = 1, 20
        call generate_lp(n, m, ms, xref, f, A, bu, bl)
        call daqp_quadprog(n, m, ms, bu, bl, x, res, lam, f=f, A=A)
        call check(res%exitflag == daqp_optimal, 'linprog exit flag')
        call check(abs(dot_product(f, xref-x)) < tol*(1.0_wp+abs(dot_product(f, xref))), 'linprog objective')
    end do

    ! --- an LP that used to cycle (a boundary step must carry its blocking
    !     constraint into the next inner solve)
    call lp_cycle()

    ! --- zero-row equality: infeasible
    call qp%setup(reshape([1.0_wp],[1,1]), [0.0_wp], reshape([0.0_wp],[1,1]), &
                  [1.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], istat, sense=[0_ip, daqp_equality])
    call check(istat == daqp_infeasible, 'zero-row equality')

    ! --- trivially infeasible (zero rows of A with infeasible bounds)
    call trivial_infeasible()

    ! --- small infeasible problem: x1 <= 1 and x1 >= 2
    call qp%setup(reshape([1.0_wp,0.0_wp,0.0_wp,1.0_wp],[2,2]), [0.0_wp,0.0_wp], &
                  reshape([1.0_wp,-1.0_wp,0.0_wp,0.0_wp],[2,2]), &
                  [1.0_wp,-2.0_wp], [-daqp_inf,-daqp_inf], istat)
    deallocate(x, lam); allocate(x(2), lam(2))
    call qp%solve(x, lam, istat)
    call check(istat == daqp_infeasible, 'small infeasible problem')

    ! --- unconstrained optimum (feasible): nothing active
    call qp%setup(reshape([1.0_wp,0.0_wp,0.0_wp,1.0_wp],[2,2]), [0.0_wp,0.0_wp], &
                  reshape([1.0_wp,0.0_wp,0.0_wp,1.0_wp],[2,2]), [1.0_wp,1.0_wp], &
                  [-daqp_inf,-daqp_inf], istat)
    call qp%solve(x, lam, istat)
    call qp%get_working_set(active, at_lower)
    call check(istat == daqp_optimal .and. size(active) == 0, 'unconstrained')

    ! --- degenerate starting point: min ||x||^2, x1+x2 <= 2, x <= 1, warm started
    !     at x = (1,1) (all three active)
    call qp%setup(reshape([1.0_wp,0.0_wp,0.0_wp,1.0_wp],[2,2]), [0.0_wp,0.0_wp], &
                  reshape([1.0_wp,1.0_wp],[1,2]), [1.0_wp,1.0_wp,2.0_wp], &
                  [(-daqp_inf, k=1,3)], istat)
    deallocate(lam); allocate(lam(3))
    call qp%solve(x, lam, istat, active=[1_ip,2_ip,3_ip])
    call check(istat == daqp_optimal .and. norm2(x) < tol, 'degenerate starting point')

    ! --- the semi-proximal method
    call semi_proximal()

    ! --- soft constraints: min 0.5x^2 - 10x, x <= 0 (soft), rho = 0.5, w = 2
    call soft_weights()

    call report('test_core_upstream')

    contains

    subroutine lp_cycle()
        real(wp) :: fc(4), Ac(5,4), buc(7), blc(7), xc(4), lc(7), xcref(4)
        type(daqp_result) :: r
        fc = [0.39565845534355815_wp, 1.0547224186098394_wp, -0.6079736130218228_wp, &
              0.08664354239083986_wp]
        Ac(1,:) = [0.02567471103491733_wp, 1.7939107275243051_wp, -1.5930403662318176_wp, &
                   -1.6268049010944303_wp]
        Ac(2,:) = [0.17430767539690895_wp, 0.9393983598218338_wp, 1.3853828760377995_wp, &
                   -0.12156671583176686_wp]
        Ac(3,:) = [0.1501841861988604_wp, 1.3229055909658318_wp, 0.06340819575037879_wp, &
                   0.00444070515780678_wp]
        Ac(4,:) = [-0.3381833190415828_wp, -0.8034139813242723_wp, 1.1295816385119204_wp, &
                   -0.13300267414467912_wp]
        Ac(5,:) = [0.8222305102124798_wp, -0.4876066796950338_wp, 1.6040189329364452_wp, &
                   1.4786982612616595_wp]
        buc = [0.6219967993947606_wp, -1.7148261687971498_wp, -4.035856723121925_wp, &
               0.42802173327763904_wp, -2.1046708743722897_wp, 2.1244179526639546_wp, &
               3.3321342110568404_wp]
        blc = [-0.5562536036016433_wp, -2.5667350684115884_wp, -5.200480703737678_wp, &
               -0.5024572486929929_wp, -2.1713208424238117_wp, 1.4858484259752665_wp, &
               1.8907523366788457_wp]
        xcref = [0.30757892308706225_wp, -1.7148261687971498_wp, 0.7856067762605017_wp, &
                 0.27583126130300034_wp]
        call daqp_quadprog(4, 7, 2, buc, blc, xc, r, lc, f=fc, A=Ac)
        call check(r%exitflag == daqp_optimal, 'cycling LP exit flag')
        call check(maxval(abs(xc-xcref)) < tol, 'cycling LP solution')
        call check(r%iter <= 20, 'cycling LP iterations')
    end subroutine lp_cycle

    subroutine trivial_infeasible()
        real(wp) :: Ht(4,4), ft(4), At(16,4), bt(16)
        integer(ip) :: is
        type(daqp_type) :: q
        Ht(1,:) = [6.837677669279314_wp, 1.3993262799977795_wp, 1.9781574256330445_wp, 0.7988389688453156_wp]
        Ht(2,:) = [1.3993262799977795_wp, 4.91607513347457_wp, 0.8347008717503388_wp, 0.964319980996552_wp]
        Ht(3,:) = [1.9781574256330445_wp, 0.8347008717503388_wp, 5.6186371819867755_wp, 0.23421356059787485_wp]
        Ht(4,:) = [0.7988389688453156_wp, 0.964319980996552_wp, 0.23421356059787485_wp, 5.512564828518534_wp]
        ft = [-0.9018168388545096_wp, 1.3888380439021342_wp, -3.2050167583822065_wp, 6.2604158413126205_wp]
        At(1,:) = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
        At(2,:) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
        At(3,:) = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
        At(4,:) = [-1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp]
        At(5,:) = [1.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
        At(6,:) = [0.5_wp, 1.0_wp, 0.5_wp, 1.0_wp]
        At(7,:) = [-1.0_wp, 0.0_wp, -1.0_wp, 0.0_wp]
        At(8,:) = [-0.5_wp, -1.0_wp, -0.5_wp, -1.0_wp]
        At(9,:) = [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
        At(10,:) = [-1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
        At(11,:) = [0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp]
        At(12,:) = [0.0_wp, -1.0_wp, 0.0_wp, 0.0_wp]
        At(13,:) = [0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp]
        At(14,:) = [0.0_wp, 0.0_wp, -1.0_wp, 0.0_wp]
        At(15,:) = [0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp]
        At(16,:) = [0.0_wp, 0.0_wp, 0.0_wp, -1.0_wp]
        bt = [2.2693082025353517_wp, 1.3445735938597536_wp, -0.2693082025353519_wp, &
              0.6554264061402464_wp, 2.2330893356345_wp, 1.172286796929877_wp, &
              -0.23308933563449985_wp, 0.8277132030701232_wp, 1.0_wp, 0.5_wp, 1.0_wp, &
              0.5_wp, 0.5_wp, 2.0_wp, 0.5_wp, 2.0_wp]
        call q%setup(Ht, ft, At, bt, [(-daqp_inf, k=1,16)], is)
        call check(is == daqp_infeasible, 'trivially infeasible')
    end subroutine trivial_infeasible

    subroutine semi_proximal()
        real(wp) :: x2(2), l2(2), x3(3), l3(3)
        real(wp) :: Hd(3,3), xd(3), fd(3)
        integer(ip) :: is
        type(daqp_type) :: q
        ! rank-1 Hessian: the x2 direction is singular
        q%eps_prox = -1.0e-3_wp
        call q%setup(reshape([1.0_wp,0.0_wp,0.0_wp,0.0_wp],[2,2]), [1.0_wp,1.0_wp], &
                     bupper=[2.0_wp,2.0_wp], blower=[-2.0_wp,-2.0_wp], istat=is)
        call q%solve(x2, l2, is)
        call check(is == daqp_optimal, 'semi-proximal exit flag')
        call check(abs(x2(1)+1.0_wp) < tol .and. abs(x2(2)+2.0_wp) < tol, 'semi-proximal solution')
        ! dense singular Hessian (ones): a full shift
        Hd = 1.0_wp
        xd = [1.0_wp, -1.0_wp, 0.5_wp]
        fd = -matmul(Hd,xd) - [0.5_wp, 0.75_wp, 1.0_wp]
        call q%set_defaults()
        q%eps_prox = -1.0e-8_wp
        call q%setup(Hd, fd, bupper=xd, blower=xd-1.0_wp, istat=is)
        call q%solve(x3, l3, is)
        call check(is == daqp_optimal, 'dense singular exit flag')
        call check(maxval(abs(x3-xd)) < tol, 'dense singular solution')
        ! singular quadratic with no linear term
        call q%set_defaults()
        call q%setup(reshape([1.0_wp,0.0_wp,0.0_wp,0.0_wp],[2,2]), &
                     bupper=[2.0_wp,2.0_wp], blower=[-2.0_wp,1.0_wp], istat=is)
        call q%solve(x2, l2, is)
        call check(is == daqp_optimal, 'singular without f exit flag')
        call check(maxval(abs(x2-[0.0_wp,1.0_wp])) < tol, 'singular without f solution')
        ! zero Hessian: all at the lower bound
        q%eps_prox = -1.0e-2_wp
        Hd = 0.0_wp
        call q%setup(Hd, [1.0_wp,1.0_wp,1.0_wp], bupper=[5.0_wp,5.0_wp,5.0_wp], &
                     blower=[-5.0_wp,-5.0_wp,-5.0_wp], istat=is)
        call q%solve(x3, l3, is)
        call check(is > 0, 'zero Hessian exit flag')
        call check(norm2(x3+5.0_wp) < tol, 'zero Hessian solution')
        ! an ill-conditioned positive definite Hessian is not regularized
        call q%set_defaults()
        call q%setup(reshape([0.5_wp*(1.0_wp+1.0e-7_wp), 0.5_wp*(1.0e-7_wp-1.0_wp), &
                              0.5_wp*(1.0e-7_wp-1.0_wp), 0.5_wp*(1.0_wp+1.0e-7_wp)],[2,2]), &
                     [0.0_wp,0.0_wp], bupper=[1.0_wp,1.0_wp], blower=[-1.0_wp,-1.0_wp], istat=is)
        call q%solve(x2, l2, is)
        call check(is == daqp_optimal, 'ill-conditioned PD exit flag')
    end subroutine semi_proximal

    subroutine soft_weights()
        real(wp) :: x1s(1), l1s(1), tsoft
        integer(ip) :: is
        type(daqp_type) :: q
        tsoft = max(1.0e-8_wp, 1.0e3_wp*epsilon(1.0_wp))
        q%rho_soft = 0.5_wp
        q%w_soft = 2.0_wp
        call q%setup(reshape([1.0_wp],[1,1]), [-10.0_wp], bupper=[0.0_wp], blower=[-daqp_inf], &
                     istat=is, sense=[daqp_soft])
        call q%solve(x1s, l1s, is)
        call check(is == daqp_soft_optimal, 'soft exit flag')
        call check(abs(x1s(1)-8.0_wp/3.0_wp) < tsoft, 'soft solution')
        call check(abs(l1s(1)-22.0_wp/3.0_wp) < tsoft, 'soft multiplier')
        q%w_soft = 3.0_wp
        call q%solve(x1s, l1s, is)
        call check(is == daqp_soft_optimal .and. abs(x1s(1)-7.0_wp/3.0_wp) < tsoft, &
                   'soft solution (w = 3)')
    end subroutine soft_weights

    end program test_core_upstream
!*****************************************************************************************
