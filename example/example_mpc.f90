!*****************************************************************************************
!>
!  Repeated solves with hot starts, as in model predictive control: a double
!  integrator is steered to the origin with a bounded input, solving a QP
!  over a horizon at each step. Only the linear term (which depends on the
!  current state) changes between the QPs, so `update` keeps the factorization
!  of `H`, and each solve starts from the previous working set.

    program example_mpc

    use daqp_module, wp => daqp_wp

    implicit none

    integer, parameter :: nh = 20              !! horizon (number of inputs)
    real(wp), parameter :: dt = 0.1_wp         !! time step
    real(wp), parameter :: umax = 1.0_wp       !! input bound
    real(wp), parameter :: r = 0.1_wp          !! weight of the input

    type(daqp_type) :: qp
    real(wp) :: Phi(2,2), Gam(2,nh), H(nh,nh), Fx(nh,2), state(2), f(nh), u(nh), lam(nh), &
                bupper(nh), blower(nh), Ak(2,2), Bk(2)
    integer :: i, k, istat, total_iter

    ! prediction: state(k) = Phi^k state(0) + sum_j Phi^(k-1-j) B u(j)
    Ak = reshape([1.0_wp, 0.0_wp, dt, 1.0_wp], [2,2])
    Bk = [0.5_wp*dt**2, dt]
    ! cost: sum_k |state(k)|^2 + r |u(k)|^2  =  0.5 u'Hu + (Fx state(0))'u + const
    H = 0.0_wp
    Fx = 0.0_wp
    Phi = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2,2])
    Gam = 0.0_wp
    do k = 1, nh
        ! the state after step k: Phi state(0) + Gam u, with Gam(:,j) = Ak^(k-j) Bk (j <= k)
        Gam = 0.0_wp
        do i = 1, k
            Gam(:,i) = matmul(mpow(Ak, k-i), Bk)
        end do
        Phi = matmul(Ak, Phi)
        H = H + 2.0_wp*matmul(transpose(Gam), Gam)
        Fx = Fx + 2.0_wp*matmul(transpose(Gam), Phi)
    end do
    do i = 1, nh
        H(i,i) = H(i,i) + 2.0_wp*r
    end do
    bupper = umax
    blower = -umax

    state = [5.0_wp, 0.0_wp]
    f = matmul(Fx, state)
    call qp%setup(H, f, bupper=bupper, blower=blower, istat=istat)
    if (istat /= daqp_success) error stop 'setup failed'

    total_iter = 0
    write(*,'(A)') ' step   position   velocity      input  iterations'
    do k = 1, 40
        call qp%solve(u, lam, istat)   ! hot start from the previous working set
        if (istat <= 0) error stop 'solve failed'
        total_iter = total_iter + qp%iter
        write(*,'(I5,3F11.5,I12)') k, state, u(1), qp%iter
        state = matmul(Ak, state) + Bk*u(1)  ! apply the first input
        f = matmul(Fx, state)
        call qp%update(istat, f=f)            ! new linear term only
    end do
    write(*,'(A,I0)') 'total iterations: ', total_iter

    contains

    function mpow(Am, p) result(P2)
        real(wp), intent(in) :: Am(2,2)
        integer, intent(in) :: p
        real(wp) :: P2(2,2)
        integer :: j
        P2 = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2,2])
        do j = 1, p
            P2 = matmul(Am, P2)
        end do
    end function mpow

    end program example_mpc
!*****************************************************************************************
