!*****************************************************************************************
!>
!  A mixed-integer QP, solved by branch and bound: choose at most two of four
!  items (binary variables) to get as close as possible to target values,
!
!      min  0.5*||x - t||^2   s.t.  x(i) in {0, 1},  x1 + x2 + x3 + x4 <= 2
!
!  with `t = (0.9, 0.2, 0.8, 0.6)`: the solution is x = (1, 0, 1, 0).

    program example_bnb

    use daqp_module, wp => daqp_wp

    implicit none

    real(wp), parameter :: t(4) = [0.9_wp, 0.2_wp, 0.8_wp, 0.6_wp]

    type(daqp_type) :: qp
    real(wp) :: H(4,4), x(4), lam(5)
    integer :: i, istat

    H = 0.0_wp
    do i = 1, 4
        H(i,i) = 1.0_wp
    end do
    ! simple bounds 0 <= x <= 1 (binary), then the cardinality constraint
    call qp%setup(H, -t, reshape([1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp], [1,4]), &
                  bupper=[1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp], &
                  blower=[0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, -daqp_inf], istat=istat, &
                  sense=[daqp_binary, daqp_binary, daqp_binary, daqp_binary, 0])
    if (istat /= daqp_success) error stop 'setup failed'

    call qp%solve(x, lam, istat)
    if (istat <= 0) error stop 'solve failed'

    write(*,'(A,4F6.2)') 'x     = ', x
    write(*,'(A,I0)')    'nodes = ', qp%outer_iter
    write(*,'(A,I0)')    'iter  = ', qp%iter

    end program example_bnb
!*****************************************************************************************
