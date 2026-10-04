!*****************************************************************************************
!>
!  Infeasible, unbounded, nonconvex, and semidefinite problems; LPs (`H = 0`
!  and no `H`); equality, immutable, and soft constraints; invalid input and
!  the out-of-memory status.

    program test_infeasible

    use daqp_module
    use daqp_test_utils

    implicit none

    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:)
    type(daqp_type) :: qp
    integer :: i, n, m, ms
    integer(ip) :: istat
    real(wp) :: tol, x1(1), l1(1), x2(2), l2(2), l3(3)

    call rng_seed(4242)
    tol = max(1.0e-6_wp, 1.0e3_wp*sqrt(epsilon(1.0_wp)))

    ! --- inconsistent bounds (lower > upper): detected by setup
    call qp%setup(eye(2), [1.0_wp,1.0_wp], bupper=[1.0_wp,1.0_wp], blower=[0.0_wp,2.0_wp], istat=istat)
    call check(istat == daqp_infeasible, 'inconsistent bounds')

    ! --- infeasible general constraints: x1+x2 >= 3, x1 <= 1, x2 <= 1
    call qp%setup(eye(2), [0.0_wp,0.0_wp], reshape([1.0_wp,1.0_wp],[1,2]), &
                  [1.0_wp,1.0_wp,daqp_inf], [-daqp_inf,-daqp_inf,3.0_wp], istat)
    call check(istat == daqp_success, 'infeasible: setup')
    call qp%solve(x2, l3, istat)
    call check(istat == daqp_infeasible, 'infeasible: exit flag')

    ! --- random infeasible problem: contradicting rows
    n = 10; m = 30; ms = 5
    call generate_qp(n, m, ms, 5, 1.0e2_wp, xref, H, f, A, bu, bl)
    A(2,:) = -A(1,:)
    bl(ms+2) = -bl(ms+1) + 1.0_wp
    bu(ms+2) = daqp_inf
    bu(ms+1) = daqp_inf
    allocate(x(n), lam(m))
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x, lam, istat)
    call check(istat == daqp_infeasible, 'contradicting rows')

    ! --- nonconvex: indefinite H (diagonal and dense)
    call qp%setup(reshape([1.0_wp,0.0_wp,0.0_wp,-1.0_wp],[2,2]), [0.0_wp,0.0_wp], &
                  bupper=[1.0_wp,1.0_wp], blower=[-1.0_wp,-1.0_wp], istat=istat)
    call check(istat == daqp_nonconvex, 'indefinite diagonal H')
    call qp%setup(reshape([1.0_wp,2.0_wp,2.0_wp,1.0_wp],[2,2]), [0.0_wp,0.0_wp], &
                  bupper=[1.0_wp,1.0_wp], blower=[-1.0_wp,-1.0_wp], istat=istat)
    ! (not in single precision: the regularization, which starts at sqrt(zero_tol), can
    ! then make H + eps*I positive definite, and the proximal loop finds a stationary point)
    if (precision(1.0_wp) > 10) call check(istat == daqp_nonconvex, 'indefinite dense H')

    ! --- LP without H: min -x1 - x2, x1 + 2 x2 <= 4, 3 x1 + x2 <= 6, x >= 0
    call qp%setup(f=[-1.0_wp,-1.0_wp], A=reshape([1.0_wp,3.0_wp,2.0_wp,1.0_wp],[2,2]), &
                  bupper=[daqp_inf,daqp_inf,4.0_wp,6.0_wp], blower=[0.0_wp,0.0_wp,-daqp_inf,-daqp_inf], &
                  istat=istat)
    call check(istat == daqp_success, 'LP: setup')
    deallocate(lam); allocate(lam(4))
    call qp%solve(x2, lam, istat)
    call check(istat == daqp_optimal, 'LP: exit flag')
    call check(maxval(abs(x2-[1.6_wp,1.2_wp])) < tol, 'LP: solution')
    call check(abs(qp%fval+2.8_wp) < tol, 'LP: objective')

    ! --- the same LP with H = 0 (proximal loop)
    call qp%setup(zeros(2), [-1.0_wp,-1.0_wp], reshape([1.0_wp,3.0_wp,2.0_wp,1.0_wp],[2,2]), &
                  [daqp_inf,daqp_inf,4.0_wp,6.0_wp], [0.0_wp,0.0_wp,-daqp_inf,-daqp_inf], istat)
    call qp%solve(x2, lam, istat)
    call check(istat == daqp_optimal, 'LP (H = 0): exit flag')
    call check(maxval(abs(x2-[1.6_wp,1.2_wp])) < tol, 'LP (H = 0): solution')

    ! --- unbounded LP: min -x1, x2 <= 1
    call qp%setup(f=[-1.0_wp,0.0_wp], bupper=[daqp_inf,1.0_wp], blower=[-daqp_inf,-daqp_inf], &
                  istat=istat)
    call qp%solve(x2, l2, istat)
    call check(istat == daqp_unbounded, 'unbounded LP')

    ! --- semidefinite H (rank 1, dense): min 0.5 (x1+x2)^2 - x1, x in [-1,1]^2
    call qp%setup(reshape([1.0_wp,1.0_wp,1.0_wp,1.0_wp],[2,2]), [-1.0_wp,0.0_wp], &
                  bupper=[1.0_wp,1.0_wp], blower=[-1.0_wp,-1.0_wp], istat=istat)
    call qp%solve(x2, l2, istat)
    call check(istat == daqp_optimal, 'semidefinite: exit flag')
    ! (the optimal value is -1, at x = (1,-1))
    call check(abs(qp%fval - objective()) < tol, 'semidefinite: objective')
    call check(abs(qp%fval - (-1.0_wp)) < tol, 'semidefinite: optimal value')
    call check(maxval(abs(x2-[1.0_wp,-1.0_wp])) < tol, 'semidefinite: solution')

    ! --- equality constraints: by flag and by equal bounds
    call qp%setup(eye(2), [-1.0_wp,-1.0_wp], reshape([1.0_wp,-1.0_wp],[1,2]), &
                  [daqp_inf,daqp_inf,0.5_wp], [-daqp_inf,-daqp_inf,0.5_wp], istat)
    call qp%solve(x2, l3, istat)
    call check(istat == daqp_optimal .and. maxval(abs(x2-[1.25_wp,0.75_wp])) < tol, &
               'equality (equal bounds)')
    call qp%setup(eye(2), [-1.0_wp,-1.0_wp], reshape([1.0_wp,-1.0_wp],[1,2]), &
                  [daqp_inf,daqp_inf,0.5_wp], [-daqp_inf,-daqp_inf,0.5_wp], istat, &
                  sense=[0_ip,0_ip,daqp_equality])
    call qp%solve(x2, l3, istat)
    call check(istat == daqp_optimal .and. maxval(abs(x2-[1.25_wp,0.75_wp])) < tol, &
               'equality (flag)')
    call check(abs(l3(3)+0.25_wp) < tol, 'equality: multiplier')  ! (Hx + f + A'lam = 0)
    ! inconsistent equalities: x1 - x2 = 0.5 and 2 x1 - 2 x2 = 3
    call qp%setup(eye(2), [-1.0_wp,-1.0_wp], reshape([1.0_wp,2.0_wp,-1.0_wp,-2.0_wp],[2,2]), &
                  [daqp_inf,daqp_inf,0.5_wp,3.0_wp], [-daqp_inf,-daqp_inf,0.5_wp,3.0_wp], istat)
    call check(istat == daqp_overdetermined .or. istat == daqp_infeasible, &
               'inconsistent equalities')

    ! --- soft constraints: min 0.5x^2 - 10x, x <= 0 (soft): s^2/(2 rho) + w s
    qp%rho_soft = 0.5_wp
    call qp%setup(eye(1), [-10.0_wp], bupper=[0.0_wp], blower=[-daqp_inf], istat=istat, &
                  sense=[daqp_soft])
    call qp%solve(x1, l1, istat)
    call check(istat == daqp_soft_optimal, 'soft (quadratic): exit flag')
    call check(abs(x1(1)-10.0_wp/3.0_wp) < tol, 'soft (quadratic): solution')
    call check(abs(qp%soft_slack-10.0_wp/3.0_wp) < tol, 'soft (quadratic): slack')
    qp%w_soft = 2.0_wp
    call qp%solve(x1, l1, istat)
    call check(abs(x1(1)-8.0_wp/3.0_wp) < tol, 'soft (linear + quadratic): solution')
    qp%w_soft = 20.0_wp  ! a large linear weight: the constraint holds (exact penalty)
    call qp%solve(x1, l1, istat)
    call check(istat == daqp_optimal .and. abs(x1(1)) < tol, 'soft (exact penalty)')
    call qp%set_defaults()

    ! --- invalid input
    call qp%setup(eye(2), [1.0_wp], bupper=[1.0_wp], blower=[0.0_wp], istat=istat)
    call check(istat == daqp_invalid_input, 'invalid: size of f')
    call qp%setup(eye(2), [1.0_wp,1.0_wp], bupper=[1.0_wp], blower=[0.0_wp,0.0_wp], istat=istat)
    call check(istat == daqp_invalid_input, 'invalid: size of blower')
    call qp%setup(eye(2), [1.0_wp,1.0_wp], bupper=[1.0_wp,1.0_wp,1.0_wp], &
                  blower=[0.0_wp,0.0_wp,0.0_wp], istat=istat)
    call check(istat == daqp_invalid_input, 'invalid: more simple bounds than variables')
    call qp%setup(eye(2), [1.0_wp,1.0_wp], bupper=[1.0_wp,1.0_wp], blower=[0.0_wp,0.0_wp], &
                  istat=istat, ms=1)
    call check(istat == daqp_invalid_input, 'invalid: ms')

    ! --- out of memory (the packed triangles would not fit the default integers)
    deallocate(f, A)
    n = 70000
    allocate(f(n), A(1,n))
    f = 1.0_wp
    A = 1.0_wp
    call qp%setup(f=f, A=A, bupper=[1.0_wp], blower=[-1.0_wp], istat=istat)
    call check(istat == daqp_out_of_memory, 'out of memory')

    call report('test_infeasible')

    contains

    pure function eye(k) result(e)
        integer, intent(in) :: k
        real(wp) :: e(k,k)
        integer :: j
        e = 0.0_wp
        do j = 1, k
            e(j,j) = 1.0_wp
        end do
    end function eye

    pure function zeros(k) result(z)
        integer, intent(in) :: k
        real(wp) :: z(k,k)
        z = 0.0_wp
    end function zeros

    real(wp) function objective()
        objective = 0.5_wp*(x2(1)+x2(2))**2 - x2(1)
    end function objective

    end program test_infeasible
!*****************************************************************************************
