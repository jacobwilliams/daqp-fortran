!*****************************************************************************************
!>
!  A small QP, set up and solved with the object-oriented interface:
!
!      min  0.5*(x1^2 + x2^2) - x1 - 2 x2
!      s.t. 0 <= x1 <= 0.5,  x2 >= 0          (simple bounds)
!           x1 + x2 <= 1                       (a general constraint)
!
!  The solution is x = (0, 1).

    program example_simple

    use daqp_module, wp => daqp_wp

    implicit none

    real(wp), parameter :: H(2,2) = reshape([1.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [2,2])
    real(wp), parameter :: f(2) = [-1.0_wp, -2.0_wp]
    real(wp), parameter :: A(1,2) = reshape([1.0_wp, 1.0_wp], [1,2])
    real(wp), parameter :: bupper(3) = [0.5_wp, daqp_inf, 1.0_wp] ! the first 2: bounds on x
    real(wp), parameter :: blower(3) = [0.0_wp, 0.0_wp, -daqp_inf]

    type(daqp_type) :: qp
    real(wp) :: x(2), lam(3)
    integer :: istat

    call qp%setup(H, f, A, bupper, blower, istat)
    if (istat /= daqp_success) error stop 'setup failed'

    call qp%solve(x, lam, istat)
    if (istat <= 0) error stop 'solve failed'

    write(*,'(A,2F10.6)') 'x      = ', x
    write(*,'(A,3F10.6)') 'lambda = ', lam   ! > 0 at an upper bound, < 0 at a lower bound
    write(*,'(A,F10.6)')  'f      = ', qp%fval
    write(*,'(A,I0)')     'iter   = ', qp%iter

    end program example_simple
!*****************************************************************************************
