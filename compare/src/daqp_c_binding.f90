!*****************************************************************************************
!>
!  Interfaces to the upstream DAQP C library (v0.10.3), for the comparison.
!
!  The derived types mirror `DAQPProblem`, `DAQPSettings`, and `DAQPResult`
!  of `upstream/daqp/include/types.h` and `api.h`, field by field (order and
!  type): a mismatch here produces plausible-looking garbage.

    module daqp_c_binding

    use, intrinsic :: iso_c_binding

    implicit none

    private

    integer(c_int), parameter, public :: daqp_c_eq_reduction_off = -1 !! `DAQP_EQ_REDUCTION_OFF`

    type, bind(c), public :: daqp_c_problem
        !! `DAQPProblem`
        integer(c_int) :: n = 0
        integer(c_int) :: m = 0
        integer(c_int) :: ms = 0
        type(c_ptr) :: H = c_null_ptr        !! `n x n`, row major
        type(c_ptr) :: f = c_null_ptr
        type(c_ptr) :: A = c_null_ptr        !! `(m-ms) x n`, row major
        type(c_ptr) :: bupper = c_null_ptr
        type(c_ptr) :: blower = c_null_ptr
        type(c_ptr) :: sense = c_null_ptr
        type(c_ptr) :: break_points = c_null_ptr
        integer(c_int) :: nh = 0
        integer(c_int) :: problem_type = 0
    end type daqp_c_problem

    type, bind(c), public :: daqp_c_settings
        !! `DAQPSettings`
        real(c_double) :: primal_tol
        real(c_double) :: dual_tol
        real(c_double) :: zero_tol
        real(c_double) :: pivot_tol
        real(c_double) :: progress_tol
        integer(c_int) :: cycle_tol
        integer(c_int) :: iter_limit
        real(c_double) :: fval_bound
        real(c_double) :: eps_prox
        real(c_double) :: eta_prox
        real(c_double) :: rho_soft
        real(c_double) :: rel_subopt
        real(c_double) :: abs_subopt
        real(c_double) :: sing_tol
        real(c_double) :: refactor_tol
        real(c_double) :: time_limit
        real(c_double) :: w_soft
        integer(c_int) :: eq_reduction
    end type daqp_c_settings

    type, bind(c), public :: daqp_c_result
        !! `DAQPResult`
        type(c_ptr) :: x = c_null_ptr
        type(c_ptr) :: lam = c_null_ptr
        real(c_double) :: fval = 0
        real(c_double) :: soft_slack = 0
        integer(c_int) :: exitflag = 0
        integer(c_int) :: iter = 0
        integer(c_int) :: nodes = 0
        real(c_double) :: solve_time = 0
        real(c_double) :: setup_time = 0
    end type daqp_c_result

    interface
        subroutine daqp_default_settings(settings) bind(c, name='daqp_default_settings')
            import :: daqp_c_settings
            implicit none
            type(daqp_c_settings), intent(inout) :: settings
        end subroutine daqp_default_settings
        subroutine daqp_c_quadprog(res, qp, settings) bind(c, name='daqp_quadprog')
            import :: daqp_c_result, daqp_c_problem, daqp_c_settings
            implicit none
            type(daqp_c_result), intent(inout) :: res
            type(daqp_c_problem), intent(inout) :: qp
            type(daqp_c_settings), intent(inout) :: settings
        end subroutine daqp_c_quadprog
        integer(c_int) function setup_daqp(qp, work, setup_time) bind(c, name='setup_daqp')
            import :: daqp_c_problem, c_ptr, c_int
            implicit none
            type(daqp_c_problem), intent(inout) :: qp
            type(c_ptr), value :: work
            type(c_ptr), value :: setup_time
        end function setup_daqp
        subroutine daqp_c_solve(res, work) bind(c, name='daqp_solve')
            import :: daqp_c_result, c_ptr
            implicit none
            type(daqp_c_result), intent(inout) :: res
            type(c_ptr), value :: work
        end subroutine daqp_c_solve
        integer(c_int) function daqp_c_update_ldp(mask, work, qp) bind(c, name='daqp_update_ldp')
            import :: daqp_c_problem, c_ptr, c_int
            implicit none
            integer(c_int), value :: mask
            type(c_ptr), value :: work
            type(daqp_c_problem), intent(inout) :: qp
        end function daqp_c_update_ldp
        integer(c_int) function setup_daqp_main(qp, work, setup_time, init_mask) bind(c, name='setup_daqp_main')
            import :: daqp_c_problem, c_ptr, c_int
            implicit none
            type(daqp_c_problem), intent(inout) :: qp
            type(c_ptr), value :: work
            type(c_ptr), value :: setup_time
            integer(c_int), value :: init_mask
        end function setup_daqp_main
        integer(c_int) function daqp_c_set_soft_weights(work, rho_l, rho_u, w_l, w_u) &
                bind(c, name='daqp_set_soft_weights')
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value :: work, rho_l, rho_u, w_l, w_u
        end function daqp_c_set_soft_weights
        subroutine daqp_c_primal_init_active(qp, x) bind(c, name='daqp_primal_init_active')
            import :: daqp_c_problem, c_double
            implicit none
            type(daqp_c_problem), intent(inout) :: qp
            real(c_double), intent(in) :: x(*)
        end subroutine daqp_c_primal_init_active
        subroutine daqp_c_dual_init_active(qp, lam) bind(c, name='daqp_dual_init_active')
            import :: daqp_c_problem, c_double
            implicit none
            type(daqp_c_problem), intent(inout) :: qp
            real(c_double), intent(in) :: lam(*)
        end subroutine daqp_c_dual_init_active
        subroutine daqp_c_set_primal_start(work, x) bind(c, name='daqp_set_primal_start')
            import :: c_ptr, c_double
            implicit none
            type(c_ptr), value :: work
            real(c_double), intent(in) :: x(*)
        end subroutine daqp_c_set_primal_start
        subroutine daqp_c_minrep(is_redundant, A, b, n, m, ms) bind(c, name='daqp_minrep')
            import :: c_int, c_double
            implicit none
            integer(c_int), intent(out) :: is_redundant(*)
            real(c_double), intent(inout) :: A(*)
            real(c_double), intent(inout) :: b(*)
            integer(c_int), value :: n, m, ms
        end subroutine daqp_c_minrep
        ! helpers (c_shim.c)
        type(c_ptr) function cmp_ws_new(settings) bind(c, name='cmp_ws_new')
            import :: daqp_c_settings, c_ptr
            implicit none
            type(daqp_c_settings), intent(inout) :: settings
        end function cmp_ws_new
        subroutine cmp_ws_free(work) bind(c, name='cmp_ws_free')
            import :: c_ptr
            implicit none
            type(c_ptr), value :: work
        end subroutine cmp_ws_free
        integer(c_int) function cmp_ws_n_active(work) bind(c, name='cmp_ws_n_active')
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value :: work
        end function cmp_ws_n_active
        subroutine cmp_ws_working_set(work, ws, lower) bind(c, name='cmp_ws_working_set')
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value :: work
            integer(c_int), intent(out) :: ws(*)
            integer(c_int), intent(out) :: lower(*)
        end subroutine cmp_ws_working_set
        integer(c_int) function cmp_ws_set_working_set(work, na, ids, lower) &
                bind(c, name='cmp_ws_set_working_set')
            import :: c_ptr, c_int
            implicit none
            type(c_ptr), value :: work
            integer(c_int), value :: na
            integer(c_int), intent(in) :: ids(*)
            integer(c_int), intent(in) :: lower(*)
        end function cmp_ws_set_working_set
    end interface

    public :: daqp_default_settings, daqp_c_quadprog, setup_daqp, daqp_c_solve, daqp_c_update_ldp
    public :: setup_daqp_main, daqp_c_set_soft_weights, daqp_c_primal_init_active, daqp_c_dual_init_active
    public :: daqp_c_set_primal_start, daqp_c_minrep
    public :: cmp_ws_new, cmp_ws_free, cmp_ws_n_active, cmp_ws_working_set, cmp_ws_set_working_set

    end module daqp_c_binding
!*****************************************************************************************
