!*****************************************************************************************
!>
!  The real kind of the build (`-DREAL32`, `-DREAL64`, `-DREAL128`): a QP
!  solved to a tolerance scaled by `epsilon`. (CI runs every test in each kind.)

    program test_kinds

    use daqp_module
    use daqp_test_utils

    implicit none

    real(wp), allocatable :: xref(:), H(:,:), f(:), A(:,:), bu(:), bl(:), x(:), lam(:)
    type(daqp_type) :: qp
    integer(ip) :: istat
    real(wp) :: tol

    write(*,'(A,I0,A,I0)') 'real kind: ', daqp_wp, ', precision: ', precision(1.0_wp)
    call check(daqp_wp == wp, 'kind of the package and of the tests')

    call rng_seed(1)
    call generate_qp(20, 60, 10, 10, 1.0e1_wp, xref, H, f, A, bu, bl)
    allocate(x(20), lam(60))
    call qp%setup(H, f, A, bu, bl, istat)
    call qp%solve(x, lam, istat)
    tol = max(1.0e3_wp*epsilon(1.0_wp), 10.0_wp*qp%primal_tol)*1.0e2_wp
    call check(istat == daqp_optimal, 'exit flag')
    call check(maxval(abs(x-xref)) < tol*(1.0_wp+maxval(abs(xref))), 'solution')
    ! the defaults are upstream's values, floored at multiples of epsilon
    call check(qp%sing_tol >= 3.7e-11_wp .and. qp%sing_tol >= 1.0e3_wp*epsilon(1.0_wp), 'sing_tol')
    call check(qp%primal_tol >= 1.0e-6_wp, 'primal_tol')

    call report('test_kinds')

    end program test_kinds
!*****************************************************************************************
