!*****************************************************************************************
!> author: Jacob Williams
!
!  The one-to-one port of the DAQP solver (dual active-set method for dense
!  convex quadratic programs).
!
!  This module is a modern Fortran translation of the C code of
!  [DAQP](https://github.com/darnstrom/daqp) v0.10.3
!  (Copyright (c) 2022 Daniel Arnström, MIT licence).
!  Changed from the original: translated to Fortran, 1-based indexing,
!  column-major storage of the constraint rows, and without the
!  branch-and-bound, hierarchical, AVI, equality-elimination,
!  per-constraint soft weight, and timing parts of upstream.
!
!  The QP is
!
!      min  0.5 x'Hx + f'x
!      s.t. blower(1:ms) <= x(1:ms)  <= bupper(1:ms)
!           blower(ms+1:) <= A x     <= bupper(ms+1:)
!
!  It is turned into the least-distance problem (LDP)
!  `min 0.5 ||u||^2 s.t. dlower <= M u <= dupper`, with `H = R'R`,
!  `u = R x + R'\f` and `M = A R^{-1}`, which is solved by a dual active-set
!  method that updates the `LDL'` factors of the working set's Gram matrix.
!
!### Storage
!
!  * `R` (holding \( R^{-1} \), or \( R \) while its inverse is deferred) is
!    the upper triangle packed by rows: element `(i,j)`, `i<=j`, is
!    `R(ridx(i,j,n))`, and a row is contiguous.
!  * `M(:,k)` is the k-th general constraint row of the LDP (so a row is
!    contiguous), and `At(:,k)` the k-th row of `A`.
!  * `Hc(j,i) = H(i,j)`: the Hessian in the memory order of the C code.
!  * `L` is the unit lower triangle packed by rows (with a placeholder
!    for the diagonal): element `(i,j)`, `j<=i`, is `L(lidx(i,j))`.
!
!  Working set positions and constraint indices are 1-based; `0` marks an
!  empty index (upstream's `-1`).

    module daqp_core

    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    private

    ! exit flags (upstream's values)
    integer(ip), parameter, public :: daqp_exit_optimal_inexact        = 4  !! optimal, found with the noise floor; violates a constraint by more than `primal_tol`
    integer(ip), parameter, public :: daqp_exit_soft_optimal           = 2  !! optimal, with a soft constraint violated
    integer(ip), parameter, public :: daqp_exit_optimal                = 1  !! optimal
    integer(ip), parameter, public :: daqp_exit_infeasible             = -1 !! primal infeasible
    integer(ip), parameter, public :: daqp_exit_cycle                  = -2 !! cycling detected
    integer(ip), parameter, public :: daqp_exit_unbounded              = -3 !! unbounded (LP)
    integer(ip), parameter, public :: daqp_exit_iterlimit              = -4 !! iteration limit reached
    integer(ip), parameter, public :: daqp_exit_nonconvex              = -5 !! the Hessian could not be factored
    integer(ip), parameter, public :: daqp_exit_overdetermined_initial = -6 !! inconsistent equalities in the initial working set
    integer(ip), parameter, public :: daqp_exit_invalid_input          = -101 !! invalid input (Fortran interface)
    integer(ip), parameter, public :: daqp_exit_out_of_memory          = -102 !! an allocation failed (Fortran interface)
    integer(ip), parameter, public :: daqp_exit_not_setup              = -103 !! the solver was not set up (Fortran interface)

    ! constraint flags (bits of `sense`)
    integer(ip), parameter, public :: daqp_active        = 1   !! the constraint is in the working set
    integer(ip), parameter, public :: daqp_lower         = 2   !! the active bound is the lower one
    integer(ip), parameter, public :: daqp_immutable     = 4   !! the constraint never leaves (or enters) the working set
    integer(ip), parameter, public :: daqp_soft          = 8   !! the constraint may be violated, at a penalty
    integer(ip), parameter, public :: daqp_binary        = 16  !! binary constraint (not supported by the port)
    integer(ip), parameter, public :: daqp_slack_fixed   = 32  !! the slack of the soft constraint is zero
    integer(ip), parameter, public :: daqp_set_aside     = 64  !! temporarily set aside (proximal step)
    integer(ip), parameter, public :: daqp_auto_equality = 128 !! active and immutable set by the detection of equal bounds
    integer(ip), parameter, public :: daqp_equality      = daqp_active + daqp_immutable !! an equality constraint

    ! update masks
    integer(ip), parameter, public :: daqp_update_rinv          = 1
    integer(ip), parameter, public :: daqp_update_m             = 2
    integer(ip), parameter, public :: daqp_update_v             = 4
    integer(ip), parameter, public :: daqp_update_d             = 8
    integer(ip), parameter, public :: daqp_update_sense         = 16
    integer(ip), parameter, public :: daqp_update_hierarchy     = 32
    integer(ip), parameter, public :: daqp_update_unconstrained = 64

    ! workspace state masks
    integer(ip), parameter :: state_pending = daqp_update_rinv + daqp_update_m + daqp_update_v + &
                                              daqp_update_d + daqp_update_sense + daqp_update_hierarchy
    integer(ip), parameter :: state_unconstrained     = 256
    integer(ip), parameter :: state_rinv_normalized   = 512
    integer(ip), parameter :: state_ill_conditioned   = 2048
    integer(ip), parameter :: state_cholesky_pending  = 4096
    integer(ip), parameter :: state_noise_floor       = 8192

    integer(ip), parameter :: empty_ind = 0                  !! an empty index (upstream's -1)
    integer(ip), parameter :: unconstrained_optimal = -2     !! return value of [[check_unconstrained]]

    real(wp), parameter, public :: daqp_inf = 1.0e30_wp      !! "infinite" bound

    ! default settings
    ! (the tolerances are floored at a multiple of epsilon, which only changes
    ! them in single precision: in double and quadruple precision, they are
    ! upstream's values)
    real(wp), parameter, public :: daqp_default_prim_tol     = max(1.0e-6_wp, 100.0_wp*epsilon(1.0_wp))
    real(wp), parameter, public :: daqp_default_dual_tol     = max(1.0e-12_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter, public :: daqp_default_zero_tol     = max(1.0e-11_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter, public :: daqp_default_prog_tol     = max(1.0e-14_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter, public :: daqp_default_pivot_tol    = 1.0e-8_wp
    integer(ip), parameter, public :: daqp_default_cycle_tol = 10
    real(wp), parameter, public :: daqp_default_eta          = -1.0_wp
    integer(ip), parameter, public :: daqp_default_iter_limit = 10000
    real(wp), parameter, public :: daqp_default_rho_soft     = 1.0e-6_wp
    real(wp), parameter, public :: daqp_default_w_soft       = 0.0_wp
    real(wp), parameter, public :: daqp_default_sing_tol     = max(3.7e-11_wp, 1000.0_wp*epsilon(1.0_wp))
    real(wp), parameter, public :: daqp_default_refactor_tol = 1.0e-9_wp
    real(wp), parameter, public :: daqp_default_eps_prox     = -max(1.0e-6_wp, 100.0_wp*epsilon(1.0_wp))

    ! internal constants
    real(wp), parameter :: auto_eta_cap        = max(1.0e-6_wp, 1000.0_wp*epsilon(1.0_wp)) ! (upstream's 1e-6, floored in single precision)
    integer(ip), parameter :: refine_min_iter  = 5
    real(wp), parameter :: refine_cond         = 1.0e6_wp
    real(wp), parameter :: refine_gain         = 1.0e3_wp
    real(wp), parameter :: refine_pivot        = 1.0e-6_wp
    real(wp), parameter :: add_noise_gain      = 10.0_wp
    real(wp), parameter :: prox_eps_max        = 1.0e-3_wp
    real(wp), parameter :: hessian_cond_eps    = 0.1_wp
    real(wp), parameter :: cond_defer_margin   = 100.0_wp
    integer(ip), parameter :: prox_face_start  = 16
    integer(ip), parameter :: prox_face_steps  = 3

    ! what `R` holds
    integer(ip), parameter :: rinv_none  = 0 !! no Hessian (LP): `R = I`
    integer(ip), parameter :: rinv_dense = 1 !! dense packed `R`
    integer(ip), parameter :: rinv_diag  = 2 !! diagonal Hessian: `R(1:n)` holds the diagonal of `R^{-1}`

    type, public :: daqp_settings
        !! Solver settings (upstream's `DAQPSettings`, with its defaults).
        real(wp)    :: primal_tol   = daqp_default_prim_tol     !! tolerance for primal feasibility
        real(wp)    :: dual_tol     = daqp_default_dual_tol     !! tolerance for dual feasibility
        real(wp)    :: zero_tol     = daqp_default_zero_tol     !! values below are regarded as zero
        real(wp)    :: pivot_tol    = daqp_default_pivot_tol    !! pivots of the `LDL'` below are reordered
        real(wp)    :: progress_tol = daqp_default_prog_tol     !! minimum objective progress (cycle guard)
        integer(ip) :: cycle_tol    = daqp_default_cycle_tol    !! iterations without progress before cycling is assumed
        integer(ip) :: iter_limit   = daqp_default_iter_limit   !! maximum number of iterations
        real(wp)    :: fval_bound   = daqp_inf             !! upper bound of the objective (infeasible above)
        real(wp)    :: eps_prox     = daqp_default_eps_prox     !! proximal regularization (negative: automatic, only if needed)
        real(wp)    :: eta_prox     = daqp_default_eta          !! tolerance of the proximal outer loop (negative: automatic)
        real(wp)    :: rho_soft     = daqp_default_rho_soft     !! reciprocal quadratic weight of the soft constraints
        real(wp)    :: rel_subopt   = 0.0_wp               !! (branch and bound only; unused)
        real(wp)    :: abs_subopt   = 0.0_wp               !! (branch and bound only; unused)
        real(wp)    :: sing_tol     = daqp_default_sing_tol     !! pivots of the `LDL'` below mark a singular working set
        real(wp)    :: refactor_tol = daqp_default_refactor_tol !! pivots below trigger a refactorization at a solution
        real(wp)    :: w_soft       = daqp_default_w_soft       !! linear weight of the soft constraints
    end type daqp_settings

    type, public :: daqp_result
        !! Result of a solve (upstream's `DAQPResult`, without the timing).
        real(wp)    :: fval       = 0.0_wp  !! objective function value
        real(wp)    :: soft_slack = 0.0_wp  !! largest violation of a soft constraint
        integer(ip) :: exitflag   = daqp_exit_not_setup !! exit flag
        integer(ip) :: iter       = 0       !! number of iterations
        integer(ip) :: nodes      = 0       !! number of outer (proximal) iterations
    end type daqp_result

    type, public :: daqp_workspace
        !! The workspace of the solver (upstream's `DAQPWorkspace`).
        !! All arrays are allocated in [[daqp_setup]], so a solve allocates nothing.
        integer(ip) :: n  = 0  !! number of variables
        integer(ip) :: m  = 0  !! number of constraints (including the simple bounds)
        integer(ip) :: ms = 0  !! number of simple bounds (the first `ms` constraints)
        integer(ip) :: ns = 0  !! number of soft constraints
        type(daqp_settings) :: settings !! the settings
        ! the problem
        logical :: has_H = .false.      !! a Hessian was given (else an LP)
        logical :: has_f = .false.      !! a linear term was given
        logical :: has_sense = .false.  !! constraint flags were given
        real(wp), allocatable :: Hc(:,:)      !! `Hc(j,i) = H(i,j)`
        real(wp), allocatable :: f(:)         !! linear term
        real(wp), allocatable :: At(:,:)      !! `At(:,k)` is the k-th row of `A` (`n x (m-ms)`)
        real(wp), allocatable :: bupper(:)    !! upper bounds
        real(wp), allocatable :: blower(:)    !! lower bounds
        integer(ip), allocatable :: sense_in(:) !! the given constraint flags
        ! the LDP
        integer(ip) :: rmode = rinv_none      !! what `R` holds
        real(wp), allocatable :: R(:)         !! packed upper triangular `R^{-1}` (see the module doc)
        real(wp), allocatable :: Mr(:,:)      !! `Mr(:,k)`: the k-th general constraint of the LDP
        real(wp), allocatable :: dupper(:)    !! upper bounds of the LDP
        real(wp), allocatable :: dlower(:)    !! lower bounds of the LDP
        real(wp), allocatable :: scaling(:)   !! normalization of the constraints
        real(wp), allocatable :: Mu(:)        !! `M'u` of the latest feasibility scan
        logical :: has_v = .false.            !! `v` is used
        real(wp), allocatable :: v(:)         !! `v = R'\f`
        integer(ip), allocatable :: sense(:)  !! constraint flags
        ! iterates
        real(wp), allocatable :: x(:)         !! the primal iterate (also `u` of the LDP)
        real(wp), allocatable :: xold(:)      !! the previous primal iterate (proximal loop)
        real(wp), allocatable :: lam(:)       !! dual iterate
        real(wp), allocatable :: lam_star(:)  !! constrained stationary point
        real(wp) :: fval = 0.0_wp             !! `||u||^2` (plus soft penalties)
        ! LDL' factors of the working set
        real(wp), allocatable :: L(:)         !! packed unit lower triangle
        real(wp), allocatable :: D(:)         !! pivots
        real(wp), allocatable :: xldl(:)      !! work array of the forward substitution
        real(wp), allocatable :: zldl(:)      !! work array (`xldl/D`)
        integer(ip) :: reuse_ind = 0          !! number of rows of the forward substitution that can be reused
        integer(ip), allocatable :: WS(:)     !! the working set (constraint indices)
        integer(ip) :: n_active = 0           !! number of active constraints
        integer(ip) :: iterations = 0         !! number of iterations of the latest solve
        integer(ip) :: sing_ind = empty_ind   !! position of the constraint that made the working set singular
        logical, allocatable :: prox_mask(:)  !! directions with a proximal regularization
        integer(ip) :: n_prox = 0             !! number of regularized directions
        real(wp) :: soft_slack = 0.0_wp       !! largest violation of a soft constraint
        integer(ip) :: nh = 1                 !! number of outer iterations (proximal loop)
        integer(ip) :: state = 0              !! state mask
        logical :: is_setup = .false.         !! set up successfully
    end type daqp_workspace

    public :: daqp_setup, daqp_solve, daqp_update_ldp, daqp_destroy, daqp_quadprog
    public :: daqp_ldp, daqp_ldp2qp_solution, daqp_reset_workspace
    public :: daqp_activate_constraints, daqp_deactivate_constraints
    public :: daqp_set_working_set, daqp_refresh_soft_weights

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Position of element `(i,j)`, `i<=j`, of an `n x n` upper triangle packed by rows.

    pure integer(ip) function ridx(i,j,n)

    integer(ip), intent(in) :: i !! row
    integer(ip), intent(in) :: j !! column
    integer(ip), intent(in) :: n !! dimension

    ridx = ((i-1)*(2*n-i))/2 + j

    end function ridx
!*****************************************************************************************

!*****************************************************************************************
!>
!  Position of element `(i,j)`, `j<=i`, of a lower triangle packed by rows.

    pure integer(ip) function lidx(i,j)

    integer(ip), intent(in) :: i !! row
    integer(ip), intent(in) :: j !! column

    lidx = ((i-1)*i)/2 + j

    end function lidx
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `flag` is set in `s`.

    pure logical function has(s,flag)

    integer(ip), intent(in) :: s    !! flags
    integer(ip), intent(in) :: flag !! flag to test

    has = iand(s,flag) /= 0

    end function has
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dot product with four accumulators (upstream's `dot_row`, which keeps its
!  order of summation under fast math).

    pure real(wp) function dot_row(n,a,b)

    integer(ip), intent(in) :: n    !! length
    real(wp), intent(in) :: a(n)    !! first vector
    real(wp), intent(in) :: b(n)    !! second vector

    real(wp) :: s0, s1, s2, s3
    integer(ip) :: i

    s0 = 0.0_wp; s1 = 0.0_wp; s2 = 0.0_wp; s3 = 0.0_wp
    i = 1
    do while (i+3 <= n)
        s0 = s0 + a(i)*b(i)
        s1 = s1 + a(i+1)*b(i+1)
        s2 = s2 + a(i+2)*b(i+2)
        s3 = s3 + a(i+3)*b(i+3)
        i = i + 4
    end do
    do while (i <= n)
        s0 = s0 + a(i)*b(i)
        i = i + 1
    end do
    dot_row = (s0+s1)+(s2+s3)

    end function dot_row
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dot product, summed in order.

    pure real(wp) function dot_seq(n,a,b)

    integer(ip), intent(in) :: n    !! length
    real(wp), intent(in) :: a(n)    !! first vector
    real(wp), intent(in) :: b(n)    !! second vector

    integer(ip) :: i

    dot_seq = 0.0_wp
    do i = 1, n
        dot_seq = dot_seq + a(i)*b(i)
    end do

    end function dot_seq
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange `x` and `xold` (no copy).

    subroutine swap_x(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    real(wp), allocatable :: tmp(:)

    call move_alloc(work%x, tmp)
    call move_alloc(work%xold, work%x)
    call move_alloc(tmp, work%xold)

    end subroutine swap_x
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange `lam` and `lam_star` (no copy).

    subroutine swap_lam(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    real(wp), allocatable :: tmp(:)

    call move_alloc(work%lam, tmp)
    call move_alloc(work%lam_star, work%lam)
    call move_alloc(tmp, work%lam_star)

    end subroutine swap_lam
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reset the working set (upstream's `reset_daqp_workspace`).

    subroutine daqp_reset_workspace(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    work%sing_ind  = empty_ind
    work%n_active  = 0
    work%reuse_ind = 0

    end subroutine daqp_reset_workspace
!*****************************************************************************************

!*****************************************************************************************
!>
!  Element `j` of constraint row `id` of the LDP (for a simple bound, only if
!  `R` is dense, and `j>=id`).

    pure real(wp) function row_elem(work,id,j)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id !! constraint
    integer(ip), intent(in) :: j  !! column

    if (id <= work%ms) then
        row_elem = work%R(ridx(id,j,work%n))
    else
        row_elem = work%Mr(j,id-work%ms)
    end if

    end function row_elem
!*****************************************************************************************

!*****************************************************************************************
!>
!  [[dot_row]] of constraint rows `id1` and `id2` of the LDP, over columns `j:n`.

    pure real(wp) function row_dot(work,id1,id2,j)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id1 !! first constraint
    integer(ip), intent(in) :: id2 !! second constraint
    integer(ip), intent(in) :: j   !! first column

    integer(ip) :: n, ms, len, k1, k2

    n = work%n
    ms = work%ms
    len = n-j+1
    if (id1 <= ms) then
        k1 = ridx(id1,j,n)
        if (id2 <= ms) then
            k2 = ridx(id2,j,n)
            row_dot = dot_row(len, work%R(k1:k1+len-1), work%R(k2:k2+len-1))
        else
            row_dot = dot_row(len, work%R(k1:k1+len-1), work%Mr(j:n,id2-ms))
        end if
    else
        if (id2 <= ms) then
            k2 = ridx(id2,j,n)
            row_dot = dot_row(len, work%Mr(j:n,id1-ms), work%R(k2:k2+len-1))
        else
            row_dot = dot_row(len, work%Mr(j:n,id1-ms), work%Mr(j:n,id2-ms))
        end if
    end if

    end function row_dot
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add constraint `add_ind` to the `LDL'` factors of the working set
!  (upstream's `daqp_update_LDL_add`). A nonzero `rho` marks a free soft
!  slack, and is added to the diagonal.

    subroutine update_LDL_add(work,add_ind,rho)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: add_ind !! constraint to add
    real(wp), intent(in) :: rho        !! reciprocal soft weight (or 0)

    integer(ip) :: i, j, id, na, off, rs, start_col, ns_active
    logical :: mi_null, mk_null
    real(wp) :: s, tmp

    work%sing_ind = empty_ind
    na = work%n_active
    off = (na*(na+1))/2  ! the new row starts after this
    ns_active = 0

    ! di <-- Mi' Mi
    mi_null = add_ind <= work%ms .and. work%rmode /= rinv_dense
    if (add_ind <= work%ms) then
        start_col = add_ind
    else
        start_col = 1
    end if
    if (mi_null) then
        s = 1.0_wp
    else
        s = row_dot(work, add_ind, add_ind, start_col)
    end if

    ! a nonzero rho marks a free soft slack and contributes to the diagonal
    if (rho /= 0.0_wp) then
        s = s + rho
        ns_active = ns_active + 1
    end if
    work%D(na+1) = s

    if (na == 0) return

    ! l <-- Mk*m
    do i = 1, na
        id = work%WS(i)
        if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) &
            ns_active = ns_active + 1
        if (id <= work%ms) then
            mk_null = work%rmode /= rinv_dense
            j = max(start_col, id)
        else
            mk_null = .false.
            j = start_col
        end if
        if (mk_null) then
            if (mi_null) then
                s = 0.0_wp
            else
                s = row_elem(work, add_ind, j)
            end if
        else if (mi_null) then
            s = row_elem(work, id, j)
        else
            s = row_dot(work, id, add_ind, j)
        end if
        work%L(off+i) = s
    end do

    ! forward substitution: l <-- L\(Mk*m)
    do i = 1, na
        s = work%L(off+i)
        rs = ((i-1)*i)/2
        do j = 1, i-1
            s = s - work%L(rs+j)*work%L(off+j)
        end do
        work%L(off+i) = s
    end do

    ! scale: l_i <-- l_i/d_i, and d_new -= l'Dl
    s = work%D(na+1)
    do i = 1, na
        tmp = work%L(off+i)
        work%L(off+i) = work%L(off+i)/work%D(i)
        s = s - tmp*work%L(off+i)
    end do
    work%D(na+1) = s

    ! check for singularity
    if (work%D(na+1) < work%settings%sing_tol .or. na >= work%n + ns_active) work%sing_ind = na+1

    end subroutine update_LDL_add
!*****************************************************************************************

!*****************************************************************************************
!>
!  Remove the constraint at position `rm_ind` from the `LDL'` factors
!  (upstream's `daqp_update_LDL_remove`: algorithm C1 of Gill et al., 1974).

    subroutine update_LDL_remove(work,rm_ind)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: rm_ind  !! position in the working set

    integer(ip) :: i, j, jj, rr, na, nc, n_update, idx, w0
    real(wp) :: p, beta, dbar, alpha

    na = work%n_active
    if (na == rm_ind) return
    n_update = na - rm_ind
    w0 = rm_ind - 1  ! w(k) = zldl(w0+k) (zldl is obsolete here)

    ! remove column rm_ind (and move the rows below it up)
    do i = rm_ind+1, na
        nc = 0
        do j = 1, i-1
            if (j /= rm_ind) then
                nc = nc + 1
                work%L(lidx(i-1,nc)) = work%L(lidx(i,j))
            else
                work%zldl(w0+i-rm_ind) = work%L(lidx(i,j))
            end if
        end do
    end do

    ! low-rank update of the L2 block
    alpha = work%D(rm_ind)
    do jj = 1, n_update
        i = rm_ind + jj  ! (old) row to update
        p = work%zldl(w0+jj)
        dbar = work%D(i) + alpha*p*p
        work%D(i-1) = dbar
        beta = p*alpha/dbar
        alpha = work%D(i)*alpha/dbar
        do rr = jj+1, n_update
            idx = lidx(rm_ind+rr-1, rm_ind+jj-1)
            work%zldl(w0+rr) = work%zldl(w0+rr) - p*work%L(idx)
            work%L(idx) = work%L(idx) + beta*work%zldl(w0+rr)
        end do
    end do

    end subroutine update_LDL_remove
!*****************************************************************************************

!*****************************************************************************************
!>
!  Reciprocal quadratic weight of a soft constraint.

    pure real(wp) function soft_rho(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    soft_rho = work%settings%rho_soft

    end function soft_rho
!*****************************************************************************************

!*****************************************************************************************
!>
!  Linear weight of a soft constraint.

    pure real(wp) function soft_w(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    soft_w = work%settings%w_soft

    end function soft_w
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the soft constraints have a linear penalty.

    pure logical function has_l1(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    has_l1 = work%settings%w_soft /= 0.0_wp

    end function has_l1
!*****************************************************************************************

!*****************************************************************************************
!>
!  Signed violation of a free soft constraint for multiplier `lam`.

    pure real(wp) function soft_residual(work,id,lam)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id !! constraint
    real(wp), intent(in) :: lam   !! multiplier

    real(wp) :: w

    w = soft_w(work)
    if (w == 0.0_wp) then
        if (lam == 0.0_wp) then
            soft_residual = 0.0_wp
        else
            soft_residual = soft_rho(work)*lam
        end if
    else
        if (has(work%sense(id),daqp_lower)) then
            soft_residual = soft_rho(work)*(lam + w)
        else
            soft_residual = soft_rho(work)*(lam - w)
        end if
    end if

    end function soft_residual
!*****************************************************************************************

!*****************************************************************************************
!>
!  Contribution to the objective from the slack of a soft constraint.

    pure real(wp) function soft_penalty(work,id,lam)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id !! constraint
    real(wp), intent(in) :: lam   !! multiplier

    real(wp) :: w

    w = soft_w(work)
    if (w == 0.0_wp) then
        soft_penalty = soft_rho(work)*lam*lam
    else if (has(work%sense(id),daqp_slack_fixed)) then
        soft_penalty = 0.0_wp
    else
        soft_penalty = soft_rho(work)*(lam*lam-w*w)
    end if

    end function soft_penalty
!*****************************************************************************************

!*****************************************************************************************
!>
!  Slack of the soft constraint that is active at working set position `i`.

    pure real(wp) function soft_slack(work,i)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: i !! working set position

    integer(ip) :: id

    id = work%WS(i)
    if (has(work%sense(id),daqp_slack_fixed)) then
        soft_slack = 0.0_wp
    else
        soft_slack = soft_residual(work, id, work%lam_star(i))
    end if

    end function soft_slack
!*****************************************************************************************

!*****************************************************************************************
!>
!  Largest violation of a soft constraint, in the units of the original problem.

    pure real(wp) function max_soft_slack(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, id
    real(wp) :: s

    max_soft_slack = 0.0_wp
    do i = 1, work%n_active
        id = work%WS(i)
        if (.not. has(work%sense(id),daqp_soft)) cycle
        s = soft_slack(work, i)
        if (s < 0.0_wp) s = -s
        s = s/work%scaling(id)
        if (s > max_soft_slack) max_soft_slack = s
    end do

    end function max_soft_slack
!*****************************************************************************************

!*****************************************************************************************
!>
!  Use the noise floor when adding constraints (cycling that persists after a
!  refactorization). Returns false if the cycling is to be reported.

    logical function set_noise_floor(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (work%n_prox > 0 .or. has(work%state,state_noise_floor)) then
        set_noise_floor = .false.
    else
        work%state = ior(work%state, state_noise_floor)
        set_noise_floor = .true.
    end if

    end function set_noise_floor
!*****************************************************************************************

!*****************************************************************************************
!>
!  The (negative) rounding level of `u = -M'lam`, below which violations are not
!  added (0 if the noise floor is not used).

    pure real(wp) function noise_floor(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i
    real(wp) :: s

    noise_floor = 0.0_wp
    if (.not. has(work%state,state_noise_floor)) return
    s = 0.0_wp
    do i = 1, work%n_active
        s = s + abs(work%lam_star(i))
    end do
    noise_floor = (-add_noise_gain*epsilon(1.0_wp))*s

    end function noise_floor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Remove the constraint at working set position `rm_ind`.

    recursive subroutine remove_constraint(work,rm_ind)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: rm_ind !! working set position

    integer(ip) :: i, id

    id = work%WS(rm_ind)
    work%sense(id) = iand(work%sense(id), not(daqp_active))
    call update_LDL_remove(work, rm_ind)
    work%n_active = work%n_active - 1

    do i = rm_ind, work%n_active
        work%WS(i) = work%WS(i+1)
        work%lam(i) = work%lam(i+1)
    end do
    ! only work before the removed constraint can be reused
    if (rm_ind-1 < work%reuse_ind) work%reuse_ind = rm_ind-1

    ! check if the removal led to singularity (can happen due to numerics)
    if (work%n_active > 0) then
        if (work%D(work%n_active) < work%settings%sing_tol) then
            work%sing_ind = work%n_active
            return
        end if
    end if
    call pivot_last(work) ! pivot for improved numerics

    end subroutine remove_constraint
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add a constraint, keeping the slack state that is marked in `sense`.

    recursive subroutine add_constraint_keep_slack(work,add_ind,lam)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: add_ind !! constraint
    real(wp), intent(in) :: lam        !! its multiplier

    real(wp) :: rho

    work%sense(add_ind) = ior(work%sense(add_ind), daqp_active)
    rho = 0.0_wp
    if (has(work%sense(add_ind),daqp_soft) .and. .not. has(work%sense(add_ind),daqp_slack_fixed)) &
        rho = soft_rho(work)
    call update_LDL_add(work, add_ind, rho)
    work%n_active = work%n_active + 1
    work%WS(work%n_active) = add_ind
    work%lam(work%n_active) = lam

    call pivot_last(work) ! pivot for improved numerics

    end subroutine add_constraint_keep_slack
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add a constraint to the working set (upstream's `daqp_add_constraint`).

    recursive subroutine add_constraint(work,add_ind,lam_in)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: add_ind !! constraint
    real(wp), intent(in) :: lam_in     !! its multiplier

    real(wp) :: lam, w, magnitude, y
    logical :: lower

    lam = lam_in
    ! mark whether the slack is zero, given the multiplier
    if (has(work%sense(add_ind),daqp_soft)) then
        work%sense(add_ind) = iand(work%sense(add_ind), not(daqp_immutable))
        w = soft_w(work)
        lower = has(work%sense(add_ind),daqp_lower)
        if (w > 0.0_wp) then
            y = lam; if (lower) y = -lam
            magnitude = 0.0_wp; if (y >= w) magnitude = w
            lam = magnitude; if (lower) lam = -magnitude
        end if
        y = lam; if (lower) y = -lam
        if (w > 0.0_wp .and. y < w) then
            work%sense(add_ind) = ior(work%sense(add_ind), daqp_slack_fixed)
        else
            work%sense(add_ind) = iand(work%sense(add_ind), not(daqp_slack_fixed))
        end if
    end if
    call add_constraint_keep_slack(work, add_ind, lam)

    end subroutine add_constraint
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute `u = -M_W'lam_star` and `fval = ||u||^2` (plus the soft penalties).

    subroutine compute_primal_and_fval(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, id, n
    real(wp) :: fval, li
    logical :: l1

    n = work%n
    fval = 0.0_wp
    l1 = has_l1(work)
    work%x = 0.0_wp
    do i = 1, work%n_active
        id = work%WS(i)
        li = work%lam_star(i)
        if (id <= work%ms) then ! simple constraint
            if (work%rmode == rinv_dense) then
                k = ridx(id,id,n)
                do j = id, n
                    work%x(j) = work%x(j) - work%R(k)*li
                    k = k + 1
                end do
            else
                work%x(id) = work%x(id) - li
            end if
        else ! general constraint
            k = id - work%ms
            do j = 1, n
                work%x(j) = work%x(j) - work%Mr(j,k)*li
            end do
        end if
        if (has(work%sense(id),daqp_soft)) then
            if (l1) then
                fval = fval + soft_penalty(work, id, li)
            else
                fval = fval + work%settings%rho_soft*li*li
            end if
        end if
    end do
    do j = 1, n
        fval = fval + work%x(j)*work%x(j)
    end do
    work%fval = fval

    end subroutine compute_primal_and_fval
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute `Mu = M'u` for all the general constraints (four rows at a time).

    subroutine compute_Mu(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: rows, n, row, k
    real(wp) :: s0, s1, s2, s3, uk

    rows = work%m - work%ms
    n = work%n
    row = 1
    do while (row+3 <= rows)
        s0 = 0.0_wp; s1 = 0.0_wp; s2 = 0.0_wp; s3 = 0.0_wp
        do k = 1, n
            uk = work%x(k)
            s0 = s0 + work%Mr(k,row)*uk
            s1 = s1 + work%Mr(k,row+1)*uk
            s2 = s2 + work%Mr(k,row+2)*uk
            s3 = s3 + work%Mr(k,row+3)*uk
        end do
        work%Mu(row)   = s0
        work%Mu(row+1) = s1
        work%Mu(row+2) = s2
        work%Mu(row+3) = s3
        row = row + 4
    end do
    do while (row <= rows)
        work%Mu(row) = dot_seq(n, work%Mr(:,row), work%x)
        row = row + 1
    end do

    end subroutine compute_Mu
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add the most violated constraint to the working set. Returns false if no
!  constraint is violated (primal feasible).

    logical function add_infeasible(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: j, n, add_ind
    real(wp) :: ep, min_val, bound, Mu, min_cand, noise
    logical :: isupper

    n = work%n
    ep = -work%settings%primal_tol
    min_val = 0.0_wp
    isupper = .false.
    add_ind = empty_ind
    noise = noise_floor(work) ! 0 unless cycling persisted

    ! simple bounds
    do j = 1, work%ms
        ! never activate immutable or already active constraints
        if (iand(work%sense(j), daqp_active+daqp_immutable) /= 0) cycle
        if (work%rmode /= rinv_dense) then
            Mu = work%x(j)
        else
            Mu = dot_seq(n-j+1, work%R(ridx(j,j,n):), work%x(j:n))
        end if
        bound = ep*work%scaling(j)
        if (bound > noise) bound = noise
        min_cand = work%dupper(j) - Mu
        if (min_cand < min_val .and. min_cand < bound) then
            add_ind = j; isupper = .true.
            min_val = min_cand
        else
            min_cand = Mu - work%dlower(j)
            if (min_cand < min_val .and. min_cand < bound) then
                add_ind = j; isupper = .false.
                min_val = min_cand
            end if
        end if
    end do

    ! general two-sided constraints
    call compute_Mu(work)
    do j = work%ms+1, work%m
        if (iand(work%sense(j), daqp_active+daqp_immutable) /= 0) cycle
        Mu = work%Mu(j-work%ms)
        bound = ep*work%scaling(j)
        if (bound > noise) bound = noise
        min_cand = work%dupper(j) - Mu
        if (min_cand < min_val .and. min_cand < bound) then
            add_ind = j; isupper = .true.
            min_val = min_cand
        else
            min_cand = Mu - work%dlower(j)
            if (min_cand < min_val .and. min_cand < bound) then
                add_ind = j; isupper = .false.
                min_val = min_cand
            end if
        end if
    end do

    ! no constraint is infeasible
    if (add_ind == empty_ind) then
        add_infeasible = .false.
        return
    end if

    ! otherwise add the infeasible constraint to the working set
    if (isupper) then
        work%sense(add_ind) = iand(work%sense(add_ind), not(daqp_lower))
    else
        work%sense(add_ind) = ior(work%sense(add_ind), daqp_lower)
    end if
    call swap_lam(work) ! lam = lam_star
    if (isupper) then
        call add_constraint(work, add_ind, -min_val)
    else
        call add_constraint(work, add_ind, min_val)
    end if
    add_infeasible = .true.

    end function add_infeasible
!*****************************************************************************************

!*****************************************************************************************
!>
!  Take the step `lam <- lam + alpha*(lam_star-lam)` (`lam + alpha*lam_star`
!  if the working set is singular), stopping at the first multiplier that
!  reaches zero (which removes its constraint) or that passes `w` (which
!  switches the state of its slack). Returns false if the full step can be taken.

    logical function remove_blocking(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, ind, rm_ind, ns_active
    logical :: singular, l1, lower, fixed, release
    real(wp) :: alpha, alpha_cand, y, ystar, p, target, rm_target, tol, w, lam, rho, lo

    singular = work%sing_ind /= empty_ind
    l1 = has_l1(work)
    alpha = daqp_inf
    rm_ind = empty_ind
    rm_target = 0.0_wp
    ! blocking beyond dual_tol, or beyond zero_tol for a singular direction
    if (singular) then
        tol = work%settings%zero_tol
    else
        tol = work%settings%dual_tol
    end if

    do i = 1, work%n_active
        ind = work%WS(i)
        if (has(work%sense(ind),daqp_immutable)) cycle
        ! fold the sign of the multiplier (dual feasibility <=> y >= 0)
        lower = has(work%sense(ind),daqp_lower)
        ystar = work%lam_star(i); if (lower) ystar = -ystar

        if (.not. l1 .or. .not. has(work%sense(ind),daqp_soft)) then
            ! blocked when the multiplier reaches zero
            if (ystar >= -tol) cycle
            target = 0.0_wp
        else
            ! the multiplier is confined to [0,w] while the slack is zero and
            ! to [w,inf) otherwise; the state switches when it leaves its range
            w = soft_w(work)
            fixed = has(work%sense(ind),daqp_slack_fixed)
            target = w; if (fixed) target = 0.0_wp
            lo = target; if (singular) lo = 0.0_wp
            if (ystar >= lo - tol) then
                lo = w; if (singular) lo = 0.0_wp
                if (.not. fixed .or. w == 0.0_wp .or. ystar <= lo + tol) cycle
                target = w ! a zero slack is released
            end if
        end if

        y = work%lam(i); if (lower) y = -y
        if (singular) then
            p = ystar
        else
            p = ystar - y
        end if
        alpha_cand = (target-y)/p
        if (target /= 0.0_wp .and. alpha_cand < 0.0_wp) alpha_cand = 0.0_wp
        if (alpha_cand < alpha) then
            alpha = alpha_cand
            rm_ind = i
            rm_target = target
        end if
    end do

    if (rm_ind == empty_ind) then ! either dual feasible or primal infeasible
        remove_blocking = .false.
        return
    end if
    remove_blocking = .true.

    ! a zero-length transition cannot make progress when the CSP is singular
    if (singular .and. alpha <= 0.0_wp) rm_target = 0.0_wp

    ! update lambda
    if (singular) then
        do i = 1, work%n_active
            work%lam(i) = work%lam(i) + alpha*work%lam_star(i)
        end do
    else
        do i = 1, work%n_active
            work%lam(i) = work%lam(i) + alpha*(work%lam_star(i)-work%lam(i))
        end do
    end if

    work%sing_ind = empty_ind
    ind = work%WS(rm_ind)
    if (rm_target == 0.0_wp) then ! the constraint leaves the working set
        call remove_constraint(work, rm_ind)
        return
    end if

    ! the slack switches state, which only adds or removes rho on the diagonal
    lam = rm_target; if (has(work%sense(ind),daqp_lower)) lam = -rm_target
    release = has(work%sense(ind),daqp_slack_fixed)
    if (release) then
        work%sense(ind) = iand(work%sense(ind), not(daqp_slack_fixed))
    else
        work%sense(ind) = ior(work%sense(ind), daqp_slack_fixed)
    end if

    ! nothing in the factorization depends on the diagonal of the last row,
    ! so a slack there can switch without forming its row of M*M' again
    if (rm_ind == work%n_active .and. .not. singular) then
        rho = soft_rho(work)
        if (release) then
            work%D(rm_ind) = work%D(rm_ind) + rho
        else
            work%D(rm_ind) = work%D(rm_ind) - rho
        end if
        work%lam(rm_ind) = lam
        ! the shift in this row's right-hand side changed
        if (work%reuse_ind > rm_ind-1) work%reuse_ind = rm_ind-1
        ns_active = 0
        do i = 1, work%n_active
            if (has(work%sense(work%WS(i)),daqp_soft) .and. &
                .not. has(work%sense(work%WS(i)),daqp_slack_fixed)) ns_active = ns_active + 1
        end do
        if (work%D(rm_ind) < work%settings%sing_tol .or. rm_ind-1 >= work%n + ns_active) then
            work%sing_ind = rm_ind
        else
            call pivot_last(work) ! the new diagonal may be a worse pivot
        end if
    else
        call remove_constraint(work, rm_ind)
        if (work%sing_ind == empty_ind) call add_constraint_keep_slack(work, ind, lam)
    end if

    end function remove_blocking
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the constrained stationary point `lam_star` (solve
!  `M_W M_W' lam_star = -d_W` with the `LDL'` factors).

    subroutine compute_CSP(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, id, na, rs
    real(wp) :: s
    logical :: l1

    na = work%n_active
    l1 = has_l1(work)
    ! forward substitution (xldl <-- L\d)
    do i = work%reuse_ind+1, na
        id = work%WS(i)
        if (has(work%sense(id),daqp_lower)) then
            s = -work%dlower(id)
        else
            s = -work%dupper(id)
        end if
        ! linear weight of a nonzero slack
        if (l1) then
            if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) &
                s = s - soft_residual(work, id, 0.0_wp)
        end if
        rs = ((i-1)*i)/2
        do j = 1, i-1
            s = s - work%L(rs+j)*work%xldl(j)
        end do
        work%xldl(i) = s
    end do
    ! scale with D
    do i = work%reuse_ind+1, na
        work%zldl(i) = work%xldl(i)/work%D(i)
    end do
    ! backward substitution (lam_star <-- L'\z)
    do i = na, 1, -1
        s = work%zldl(i)
        do j = na, i+1, -1
            s = s - work%lam_star(j)*work%L(lidx(j,i))
        end do
        work%lam_star(i) = s
    end do
    work%reuse_ind = na ! save the forward substitution

    end subroutine compute_CSP
!*****************************************************************************************

!*****************************************************************************************
!>
!  The direction of a singular working set (stored in `lam_star`).

    subroutine compute_singular_direction(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, s, id
    logical :: lower, flip
    real(wp) :: w, y

    s = work%sing_ind
    ! backward substitution (p <-- L'\(-l))
    do i = s-1, 1, -1
        work%lam_star(i) = -work%L(lidx(s,i))
        do j = s-1, i+1, -1
            work%lam_star(i) = work%lam_star(i) - work%lam_star(j)*work%L(lidx(j,i))
        end do
    end do
    work%lam_star(s) = 1.0_wp

    ! orient the direction such that it is a descent direction
    id = work%WS(s)
    lower = has(work%sense(id),daqp_lower)
    flip = lower
    if (has(work%sense(id),daqp_soft) .and. has(work%sense(id),daqp_slack_fixed)) then
        w = soft_w(work)
        y = work%lam(s); if (lower) y = -y
        if (w > 0.0_wp .and. y >= w-work%settings%dual_tol) flip = .not. lower
    end if
    if (flip) then
        do i = 1, s
            work%lam_star(i) = -work%lam_star(i)
        end do
    end if

    end subroutine compute_singular_direction
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `-d_W'lam` exceeds `fval = ||u||^2` (its value at a solution) by
!  more than `fval` and the active residuals allowed by `primal_tol`.

    pure logical function inconsistent_dual(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, id
    real(wp) :: dl, tol, l

    inconsistent_dual = .false.
    dl = 0.0_wp
    tol = 0.0_wp
    do i = 1, work%n_active
        id = work%WS(i)
        if (has(work%sense(id),daqp_soft)) return
        l = work%lam_star(i)
        if (has(work%sense(id),daqp_lower)) then
            dl = dl + l*work%dlower(id)
        else
            dl = dl + l*work%dupper(id)
        end if
        tol = tol + abs(l)*work%scaling(id)
    end do
    inconsistent_dual = -dl > 2.0_wp*work%fval + work%settings%primal_tol*tol

    end function inconsistent_dual
!*****************************************************************************************

!*****************************************************************************************
!>
!  Swap the last two constraints of the working set if the next to last pivot is small.

    recursive subroutine pivot_last(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: rm_ind, na, ind_old
    real(wp) :: lam_old

    na = work%n_active
    if (na <= 1) return
    if (work%sing_ind /= empty_ind) return
    rm_ind = na - 1
    if (work%D(rm_ind) < work%settings%pivot_tol .and. work%D(rm_ind) < work%D(na)) then
        ind_old = work%WS(rm_ind)
        ! binaries never swap order (since this order is exploited)
        if (has(work%sense(ind_old),daqp_binary) .and. has(work%sense(work%WS(na)),daqp_binary)) return
        lam_old = work%lam(rm_ind)
        call remove_constraint(work, rm_ind) ! pivot_last might be recursively called here
        if (work%sing_ind /= empty_ind) return ! abort if D becomes singular
        ! reordering only: the slack state of the constraint is unchanged
        call add_constraint_keep_slack(work, ind_old, lam_old)
    end if

    end subroutine pivot_last
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add constraint `id` to the working set, with a multiplier that reproduces
!  the slack state that is marked in `sense`.

    subroutine activate_constraint(work,id)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: id !! constraint

    real(wp) :: lam, w

    lam = 1.0_wp
    w = 0.0_wp
    if (has(work%sense(id),daqp_soft)) w = soft_w(work)
    if (w > 0.0_wp) then
        if (has(work%sense(id),daqp_slack_fixed)) then
            lam = 0.9_wp*w
        else
            lam = w + 1.0_wp
        end if
    end if
    if (has(work%sense(id),daqp_lower)) lam = -lam
    call add_constraint(work, id, lam)

    end subroutine activate_constraint
!*****************************************************************************************

!*****************************************************************************************
!>
!  Remove the last constraint of a singular working set. Equalities become
!  mutable (re-added if violated). Returns the constraint.

    integer(ip) function drop_singular_last(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: id

    id = work%WS(work%n_active)
    work%n_active = work%n_active - 1
    work%sense(id) = iand(work%sense(id), not(daqp_active))
    if (.not. has(work%sense(id),daqp_binary)) work%sense(id) = iand(work%sense(id), not(daqp_immutable))
    work%sing_ind = empty_ind
    if (work%reuse_ind > work%n_active) work%reuse_ind = work%n_active
    drop_singular_last = id

    end function drop_singular_last
!*****************************************************************************************

!*****************************************************************************************
!>
!  Activate the constraints that are marked active in `sense`. Equalities are
!  activated before inequalities. Returns 1, or
!  `daqp_exit_overdetermined_initial` if the equalities are inconsistent.

    integer(ip) function daqp_activate_constraints(work) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, id, first_mutable
    real(wp) :: dependency_residual, dependency_scale, bound, term

    first_mutable = work%m + 1
    exitflag = 1
    do i = 1, work%m
        if (.not. has(work%sense(i),daqp_active)) cycle
        if (.not. has(work%sense(i),daqp_immutable)) then
            if (i < first_mutable) first_mutable = i
            cycle
        end if
        call activate_constraint(work, i)
        if (work%sing_ind /= empty_ind) then
            ! the new equality is linearly dependent on the active equalities
            dependency_residual = 0.0_wp
            dependency_scale = 1.0_wp
            call compute_singular_direction(work)
            do j = 1, work%n_active
                id = work%WS(j)
                if (has(work%sense(id),daqp_lower)) then
                    bound = work%dlower(id)
                else
                    bound = work%dupper(id)
                end if
                term = work%lam_star(j)*bound
                dependency_residual = dependency_residual + term
                dependency_scale = dependency_scale + abs(term)
            end do
            ! the dependency might only be numerical => keep as a mutable constraint
            id = drop_singular_last(work)
            if (dependency_residual > work%settings%primal_tol*dependency_scale .or. &
                dependency_residual < -work%settings%primal_tol*dependency_scale) &
                exitflag = daqp_exit_overdetermined_initial
        end if
    end do

    ! activate the active inequalities
    do i = first_mutable, work%m
        if (.not. has(work%sense(i),daqp_active) .or. has(work%sense(i),daqp_immutable)) cycle
        call activate_constraint(work, i)
        if (work%sing_ind /= empty_ind) then
            ! drop the dependent constraint and leave the remaining mutable ones inactive
            id = drop_singular_last(work)
            do j = i+1, work%m
                if (.not. has(work%sense(j),daqp_immutable)) &
                    work%sense(j) = iand(work%sense(j), not(daqp_active))
            end do
            return
        end if
    end do

    end function daqp_activate_constraints
!*****************************************************************************************

!*****************************************************************************************
!>
!  Deactivate all the active constraints that are mutable (i.e., not equalities).

    subroutine daqp_deactivate_constraints(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, id

    do i = 1, work%n_active
        id = work%WS(i)
        if (has(work%sense(id),daqp_immutable)) cycle
        work%sense(id) = iand(work%sense(id), not(daqp_active))
    end do
    call daqp_reset_workspace(work)

    end subroutine daqp_deactivate_constraints
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve `L D L' dlam = r` (`r` and `dlam` in `xldl`, `zldl` is used as scratch).

    subroutine solve_working_set(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, na, rs
    real(wp) :: s

    na = work%n_active
    ! forward substitution L*y = r
    do i = 1, na
        s = work%xldl(i)
        rs = ((i-1)*i)/2
        do j = 1, i-1
            s = s - work%L(rs+j)*work%xldl(j)
        end do
        work%xldl(i) = s
    end do
    ! scale by D^{-1}
    do i = 1, na
        work%zldl(i) = work%xldl(i)/work%D(i)
    end do
    ! backward substitution L'*dlam = z
    do i = na, 1, -1
        s = work%zldl(i)
        do j = na, i+1, -1
            s = s - work%xldl(j)*work%L(lidx(j,i))
        end do
        work%xldl(i) = s
    end do

    end subroutine solve_working_set
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y <-- y - M_W'*dlam` for the working set `W`.
!  (`y` and `dlam` may be components of `work` that this routine does not
!  otherwise reference.)

    subroutine sub_working_set_rows(work,dlam,y)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(in) :: dlam(:)  !! multipliers of the working set
    real(wp), intent(inout) :: y(:)  !! vector of size `n`

    integer(ip) :: i, j, k, id, n
    real(wp) :: dl

    n = work%n
    do i = 1, work%n_active
        dl = dlam(i)
        id = work%WS(i)
        if (id <= work%ms) then
            if (work%rmode == rinv_dense) then
                k = ridx(id,id,n)
                do j = id, n
                    y(j) = y(j) - work%R(k)*dl
                    k = k + 1
                end do
            else
                y(id) = y(id) - dl
            end if
        else
            k = id - work%ms
            do j = 1, n
                y(j) = y(j) - work%Mr(j,k)*dl
            end do
        end if
    end do

    end subroutine sub_working_set_rows
!*****************************************************************************************

!*****************************************************************************************
!>
!  One step of iterative refinement for the active constraints: solve
!  `LDL' dlam = r` with the residuals `r = M_W u - d_W`, and update
!  `lam_star += dlam`, `u -= M_W' dlam`.

    subroutine refine_active(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, id, n, na
    real(wp) :: Mu, d, fval

    n = work%n
    na = work%n_active
    ! xldl and zldl are used as scratch, which invalidates the cached CSP
    work%reuse_ind = 0

    ! r(i) = M_i*u - d_i in xldl
    do i = 1, na
        id = work%WS(i)
        if (id <= work%ms) then
            if (work%rmode == rinv_dense) then
                Mu = 0.0_wp
                k = ridx(id,id,n)
                do j = id, n
                    Mu = Mu + work%R(k)*work%x(j)
                    k = k + 1
                end do
            else
                Mu = work%x(id)
            end if
        else
            Mu = dot_seq(n, work%Mr(:,id-work%ms), work%x)
        end if
        if (has(work%sense(id),daqp_lower)) then
            d = work%dlower(id)
        else
            d = work%dupper(id)
        end if
        work%xldl(i) = Mu - d
        ! a nonzero soft slack adds a diagonal term to the CSP system
        if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) &
            work%xldl(i) = work%xldl(i) - soft_residual(work, id, work%lam_star(i))
    end do

    call solve_working_set(work) ! xldl = dlam

    do i = 1, na
        work%lam_star(i) = work%lam_star(i) + work%xldl(i)
    end do
    call sub_working_set_rows(work, work%xldl, work%x)

    ! recompute fval, since both u and lam_star changed
    fval = 0.0_wp
    do i = 1, na
        id = work%WS(i)
        if (has(work%sense(id),daqp_soft)) fval = fval + soft_penalty(work, id, work%lam_star(i))
    end do
    do j = 1, n
        fval = fval + work%x(j)*work%x(j)
    end do
    work%fval = fval

    end subroutine refine_active
!*****************************************************************************************

!*****************************************************************************************
!>
!  Residual of active constraint `id` at `x`: `A_id x - b_id` (`b` the active side).

    pure real(wp) function active_residual(work,id)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id !! constraint

    real(wp) :: val

    if (id <= work%ms) then
        val = work%x(id)
    else
        val = dot_seq(work%n, work%At(:,id-work%ms), work%x)
    end if
    if (has(work%sense(id),daqp_lower)) then
        active_residual = val - work%blower(id)
    else
        active_residual = val - work%bupper(id)
    end if

    end function active_residual
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y <-- Rinv*y` (as [[daqp_ldp2qp_solution]], without `v`).

    subroutine apply_Rinv(work,y)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(inout) :: y(:) !! vector of size `n`

    integer(ip) :: i, j, k, n

    n = work%n
    if (work%rmode == rinv_dense) then
        k = 1
        do i = 1, n
            y(i) = y(i)*work%R(k)
            k = k + 1
            do j = i+1, n
                y(i) = y(i) + work%R(k)*y(j)
                k = k + 1
            end do
        end do
        do i = 1, work%ms
            y(i) = y(i)/work%scaling(i)
        end do
    else if (work%rmode == rinv_diag) then
        do i = 1, n
            y(i) = y(i)*work%R(i)
        end do
    end if

    end subroutine apply_Rinv
!*****************************************************************************************

!*****************************************************************************************
!>
!  One step of iterative refinement of `x` on the active constraints, done
!  directly in `x` to avoid the cancellation in `x = Rinv*(u-v)`.

    subroutine refine_primal(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, id, n, na
    real(wp) :: dfval, scale2, min_D, err, val

    n = work%n
    na = work%n_active
    dfval = 0.0_wp
    if (na == 0 .or. work%sing_ind /= empty_ind) return

    ! skip if x is accurate
    if (.not. has(work%state,state_ill_conditioned) .and. work%n_prox == 0) then
        scale2 = work%fval
        min_D = 1.0_wp
        do i = 1, na
            if (work%D(i) < min_D) min_D = work%D(i)
        end do
        if (work%has_v) then
            do i = 1, n
                if (work%v(i)*work%v(i) > scale2) scale2 = work%v(i)*work%v(i)
            end do
        end if
        err = refine_gain*epsilon(1.0_wp)
        if (err*err*scale2 <= work%settings%primal_tol*work%settings%primal_tol*min_D*min_D) return
    end if
    do i = 1, na
        if (has(work%sense(work%WS(i)),daqp_soft)) return
    end do

    ! xldl and zldl are used as scratch
    work%reuse_ind = 0

    ! r = S*(A_W x - b_W)
    do i = 1, na
        id = work%WS(i)
        val = active_residual(work, id)
        dfval = dfval + work%lam_star(i)*val
        work%xldl(i) = val*work%scaling(id)
    end do

    call solve_working_set(work) ! r = dlam

    ! du = -M_W'*dlam
    do j = 1, n
        work%zldl(j) = 0.0_wp
    end do
    call sub_working_set_rows(work, work%xldl, work%zldl)
    do j = 1, n
        dfval = dfval + 0.5_wp*work%zldl(j)*work%zldl(j)
    end do

    call apply_Rinv(work, work%zldl) ! dx
    do i = 1, n
        work%x(i) = work%x(i) + work%zldl(i)
    end do

    do i = 1, na
        work%lam_star(i) = work%lam_star(i) + work%xldl(i)*work%scaling(work%WS(i))
    end do
    work%fval = work%fval + 2.0_wp*dfval ! fval is twice the objective function value

    end subroutine refine_primal
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the LDP with the dual active-set method (upstream's `daqp_ldp`).
!  Returns the exit flag.

    integer(ip) function daqp_ldp(work) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: iter, i, cycle_counter, refine_adds
    logical :: tried_repair, to_guard
    real(wp) :: best_fval, fval_bound, min_D

    exitflag = daqp_exit_iterlimit
    work%soft_slack = 0.0_wp ! only set when a solution is found
    tried_repair = .false.
    cycle_counter = 0
    refine_adds = 0  ! refinements without progress that added a constraint
    best_fval = -1.0_wp
    fval_bound = 2.0_wp*work%settings%fval_bound ! the internal objective is twice the nominal
    work%state = iand(work%state, not(state_noise_floor))

    ! clean up a singular working set on entry
    if (work%sing_ind /= empty_ind .and. work%sing_ind == work%n_active) then
        if (.not. has(work%sense(work%WS(work%sing_ind)),daqp_immutable)) then
            i = work%WS(work%n_active)
            work%n_active = work%n_active - 1
            work%sense(i) = iand(work%sense(i), not(daqp_active))
            work%sing_ind = empty_ind
            if (work%reuse_ind > work%n_active) work%reuse_ind = work%n_active
        end if
    end if

    do iter = 1, work%settings%iter_limit - 1

        if (work%sing_ind == empty_ind) then
            call compute_CSP(work)
            ! check dual feasibility of the CSP
            if (.not. remove_blocking(work)) then ! lam_star >= 0 (dual feasible)
                call compute_primal_and_fval(work)
                ! fval termination criterion
                if (work%fval > fval_bound) then
                    exitflag = daqp_exit_infeasible
                    exit
                end if
                to_guard = .false.
                ! try to add an infeasible constraint
                if (.not. add_infeasible(work)) then ! primal feasible: optimum found

                    min_D = work%D(1)
                    do i = 2, work%n_active
                        if (work%D(i) < min_D) min_D = work%D(i)
                    end do

                    ! if the LDL is truly ill-conditioned, refactor for a better pivot ordering
                    if (work%n_active > 2 .and. .not. tried_repair .and. &
                        min_D < work%settings%refactor_tol) then
                        tried_repair = .true.
                        ! correct LOWER/UPPER (important for equality constraints)
                        do i = 1, work%n_active
                            if (work%lam(i) >= 0.0_wp) then
                                work%sense(work%WS(i)) = iand(work%sense(work%WS(i)), not(daqp_lower))
                            else
                                work%sense(work%WS(i)) = ior(work%sense(work%WS(i)), daqp_lower)
                            end if
                        end do
                        call daqp_reset_workspace(work)
                        i = daqp_activate_constraints(work)
                        cycle ! try again with a new LDL factorization
                    end if

                    ! if the LDL is near-singular, apply one step of iterative
                    ! refinement before declaring optimal (at most two that add
                    ! a constraint without progress)
                    if (work%n_active > 0 .and. min_D < refine_pivot .and. refine_adds < 2) then
                        call refine_active(work)
                        ! a constraint added after the refinement goes through the cycle guard
                        if (add_infeasible(work)) then
                            refine_adds = refine_adds + 1
                            to_guard = .true.
                        end if
                    end if

                    if (.not. to_guard) then
                        ! check for an inconsistent dual
                        if (inconsistent_dual(work)) then
                            exitflag = daqp_exit_infeasible
                            exit
                        end if
                        ! softening was needed if a soft constraint ended up violated
                        work%soft_slack = max_soft_slack(work)
                        if (work%soft_slack > work%settings%primal_tol) then
                            exitflag = daqp_exit_soft_optimal
                        else
                            exitflag = daqp_exit_optimal
                        end if
                        exit
                    end if
                end if

                ! cycle guard
                if (work%fval - best_fval < work%settings%progress_tol) then
                    cycle_counter = cycle_counter + 1
                    if (cycle_counter-1 > work%settings%cycle_tol) then
                        if (tried_repair) then
                            if (.not. set_noise_floor(work)) then
                                exitflag = daqp_exit_cycle
                                exit
                            end if
                            cycle_counter = 0
                            best_fval = -1.0_wp
                        else ! cycling -> try to reorder and refactorize the LDL
                            tried_repair = .true.
                            call daqp_reset_workspace(work)
                            i = daqp_activate_constraints(work)
                            cycle_counter = 0
                            best_fval = -1.0_wp
                        end if
                    end if
                else ! progress was made
                    best_fval = work%fval
                    cycle_counter = 0
                    refine_adds = 0
                end if
            end if
        else ! singular case
            call compute_singular_direction(work)
            if (.not. remove_blocking(work)) then
                exitflag = daqp_exit_infeasible
                exit
            end if
        end if
    end do

    work%iterations = iter

    end function daqp_ldp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Compute the QP solution `x = Rinv*(u-v)` from the LDP solution, and scale
!  the multipliers back.

    subroutine daqp_ldp2qp_solution(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i

    if (work%has_v) then
        do i = 1, work%n
            work%x(i) = work%x(i) - work%v(i)
        end do
    end if
    if (work%rmode /= rinv_none) call apply_Rinv(work, work%x)
    do i = 1, work%n_active
        work%lam_star(i) = work%lam_star(i)*work%scaling(work%WS(i))
    end do

    end subroutine daqp_ldp2qp_solution
!*****************************************************************************************

!*****************************************************************************************
!>
!  The proximal regularization, scaled by the Hessian.

    pure real(wp) function prox_reg_scaled(work,hessian_scale) result(eps)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(in) :: hessian_scale    !! largest absolute diagonal of `H`

    real(wp) :: fl

    eps = abs(work%settings%eps_prox) ! negative eps_prox selects the automatic mode
    fl = sqrt(work%settings%zero_tol)*hessian_scale
    if (eps > 0.0_wp .and. eps < fl) eps = fl

    end function prox_reg_scaled
!*****************************************************************************************

!*****************************************************************************************
!>
!  Largest absolute diagonal element of `H`.

    pure real(wp) function hessian_scale_of(work) result(hs)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i

    hs = 0.0_wp
    do i = 1, work%n
        if (abs(work%Hc(i,i)) > hs) hs = abs(work%Hc(i,i))
    end do

    end function hessian_scale_of
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form the packed Cholesky factor `R` of `H` (with reciprocal diagonal).
!  Returns false if a pivot is not positive.

    logical function form_R(work,regularize_all,eps,min_pivot,max_pivot)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: regularize_all    !! regularize all the directions
    real(wp), intent(in) :: eps              !! the regularization
    real(wp), intent(out) :: min_pivot       !! smallest pivot (of the unregularized directions)
    real(wp), intent(out) :: max_pivot       !! largest pivot

    integer(ip) :: i, j, k, n, di, kik, kij
    real(wp) :: pivot, inv_diag, s, s0, s1, s2, s3, c

    n = work%n
    min_pivot = daqp_inf
    max_pivot = 0.0_wp
    form_R = .false.

    ! pack (symmetrized) H
    k = 1
    do i = 1, n
        work%R(k) = work%Hc(i,i)
        if (regularize_all) then
            work%R(k) = work%R(k) + eps
        else if (work%n_prox > 0) then
            if (work%prox_mask(i)) work%R(k) = work%R(k) + eps
        end if
        k = k + 1
        do j = i+1, n
            work%R(k) = 0.5_wp*(work%Hc(j,i) + work%Hc(i,j))
            k = k + 1
        end do
    end do

    do i = 1, n
        di = ridx(i,i,n)
        pivot = work%R(di)
        kik = i  ! position of R(k,i)
        do k = 1, i-1
            pivot = pivot - work%R(kik)*work%R(kik)
            kik = kik + n - k
        end do
        if (pivot <= work%settings%zero_tol) return
        if (pivot < min_pivot) then
            if (regularize_all .or. work%n_prox == 0) then
                min_pivot = pivot
            else if (.not. work%prox_mask(i)) then
                min_pivot = pivot
            end if
        end if
        if (pivot > max_pivot) max_pivot = pivot
        inv_diag = 1.0_wp/sqrt(pivot)
        ! four entries of row i per sweep over k, since R(k,j:j+3) are contiguous
        ! (the same order of summation; the sweep is too short to pay off for i < 3)
        j = i + 1
        if (i >= 3) then
            do while (j+3 <= n)
                s0 = work%R(di+j-i); s1 = work%R(di+j-i+1)
                s2 = work%R(di+j-i+2); s3 = work%R(di+j-i+3)
                kik = i
                kij = j
                do k = 1, i-1
                    c = work%R(kik)
                    s0 = s0 - c*work%R(kij);   s1 = s1 - c*work%R(kij+1)
                    s2 = s2 - c*work%R(kij+2); s3 = s3 - c*work%R(kij+3)
                    kik = kik + n - k
                    kij = kij + n - k
                end do
                work%R(di+j-i) = s0*inv_diag;   work%R(di+j-i+1) = s1*inv_diag
                work%R(di+j-i+2) = s2*inv_diag; work%R(di+j-i+3) = s3*inv_diag
                j = j + 4
            end do
        end if
        do while (j <= n)
            s = work%R(di+j-i)
            kik = i
            kij = j
            do k = 1, i-1
                s = s - work%R(kik)*work%R(kij)
                kik = kik + n - k
                kij = kij + n - k
            end do
            work%R(di+j-i) = s*inv_diag
            j = j + 1
        end do
        work%R(di) = inv_diag
    end do
    form_R = .true.

    end function form_R
!*****************************************************************************************

!*****************************************************************************************
!>
!  Invert the packed factor `R` in place.

    subroutine invert_R(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n, dk, di
    real(wp) :: rkk, rki

    n = work%n
    do k = 1, n
        dk = ridx(k,k,n)
        rkk = work%R(dk)
        do j = k+1, n
            work%R(dk+j-k) = work%R(dk+j-k)*(-rkk)
        end do
        do i = k+1, n
            di = ridx(i,i,n)
            work%R(dk+i-k) = work%R(dk+i-k)*work%R(di)
            rki = work%R(dk+i-k)
            do j = i+1, n
                work%R(dk+j-k) = work%R(dk+j-k) - work%R(di+j-i)*rki
            end do
        end do
    end do

    end subroutine invert_R
!*****************************************************************************************

!*****************************************************************************************
!>
!  Complete the inverse of the factor, and check the conditioning of `H`.
!  Returns false if `H` needs to be refactored with a regularization.

    logical function finish_Rinv(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n
    real(wp) :: hinv_max, hmax, hii, s2

    n = work%n
    work%state = iand(work%state, not(state_cholesky_pending))
    call invert_R(work)
    ! cond(H) >= max (H^-1)_ii * max H_ii
    hinv_max = 0.0_wp
    hmax = 0.0_wp
    k = 1
    do i = 1, n
        hii = work%Hc(i,i)
        s2 = 0.0_wp
        do j = i, n
            s2 = s2 + work%R(k)*work%R(k)
            k = k + 1
        end do
        if (s2 > hinv_max) hinv_max = s2
        if (hii > hmax) hmax = hii
    end do
    ! regularize an ill-conditioned Hessian, or mark it for refinement
    if (work%n_prox == 0 .and. real(n,wp)*epsilon(1.0_wp)*hinv_max*hmax > hessian_cond_eps) then
        finish_Rinv = .false.
        return
    end if
    if (hinv_max*hmax > refine_cond) work%state = ior(work%state, state_ill_conditioned)
    finish_Rinv = .true.

    end function finish_Rinv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Factor `H` (upstream's `daqp_update_R`): diagonal, dense, or regularized
!  (proximal) when it is singular. With `defer_inverse`, the inverse of a
!  well-conditioned factor is deferred until a constrained solve needs it.

    integer(ip) function update_R(work,defer_inverse) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: defer_inverse !! defer the inverse of the factor

    integer(ip) :: i, j, n, regularization_tries
    real(wp) :: eps, zero_tol, factor_tol, hessian_scale, acceptance_tol, &
                dmin, dmax, hi, min_pivot, max_pivot, d, cond, eps_mach, ptol
    logical :: regularize_all, force_prox, is_diagonal, ok, zero_row

    n = work%n
    eps = work%settings%eps_prox
    zero_tol = work%settings%zero_tol
    factor_tol = zero_tol
    hessian_scale = 0.0_wp
    regularization_tries = 0
    eps_mach = epsilon(1.0_wp)
    flag = 1

    force_prox = work%settings%eps_prox > 0.0_wp
    regularize_all = force_prox

    ! reset the semi-proximal mask for this factorization
    work%prox_mask = .false.
    work%n_prox = 0
    work%state = iand(work%state, not(state_rinv_normalized + state_ill_conditioned + &
                                      state_cholesky_pending))

    if (.not. work%has_H) then ! LP: all directions need proximal regularization
        if (work%has_f) work%n_prox = n
        work%scaling(1:work%ms) = 1.0_wp
        return
    end if

    ! check if diagonal
    is_diagonal = .true.
    do i = 1, n
        if (abs(work%Hc(i,i)) > hessian_scale) hessian_scale = abs(work%Hc(i,i))
        do j = i+1, n
            if (work%Hc(j,i) > zero_tol .or. work%Hc(j,i) < -zero_tol) then
                is_diagonal = .false.
                exit
            end if
        end do
        if (.not. is_diagonal) exit
    end do

    if (force_prox) then
        if (.not. is_diagonal) hessian_scale = hessian_scale_of(work)
        eps = prox_reg_scaled(work, hessian_scale)
        if (eps <= 0.0_wp) then
            flag = daqp_exit_nonconvex
            return
        end if
        work%n_prox = n
        work%prox_mask = .true.
    end if

    ! diagonal case
    if (is_diagonal) then
        if (hessian_scale > 0.0_wp) factor_tol = zero_tol*hessian_scale
        eps = prox_reg_scaled(work, hessian_scale)
        acceptance_tol = min(factor_tol, zero_tol)
        work%rmode = rinv_diag
        dmin = daqp_inf
        dmax = 0.0_wp
        do i = 1, n
            hi = work%Hc(i,i)
            if (force_prox .or. hi <= factor_tol) then
                if (.not. force_prox) then
                    work%prox_mask(i) = .true.
                    work%n_prox = work%n_prox + 1
                end if
                hi = hi + eps
            end if
            if (hi <= acceptance_tol) then
                flag = daqp_exit_nonconvex
                return
            end if
            hi = sqrt(hi)
            work%R(i) = 1.0_wp/hi
            if (i <= work%ms) work%scaling(i) = hi
            if (hi < dmin) dmin = hi
            if (hi > dmax) dmax = hi
        end do
        if (dmax*dmax > refine_cond*dmin*dmin) work%state = ior(work%state, state_ill_conditioned)
        return
    end if

    ! not diagonal
    work%rmode = rinv_dense
    if (.not. regularize_all) then
        ! zero rows of H are decoupled => regularize only them (semi-proximal)
        hessian_scale = hessian_scale_of(work)
        do i = 1, n
            if (work%Hc(i,i) /= 0.0_wp) cycle
            zero_row = .true.
            do j = 1, n
                if (work%Hc(j,i) /= 0.0_wp .or. work%Hc(i,j) /= 0.0_wp) then
                    zero_row = .false.
                    exit
                end if
            end do
            if (.not. zero_row) cycle
            work%prox_mask(i) = .true.
            work%n_prox = work%n_prox + 1
        end do
        if (work%n_prox > 0) then
            eps = prox_reg_scaled(work, hessian_scale)
            if (eps <= 0.0_wp) then
                flag = daqp_exit_nonconvex
                return
            end if
        end if
    end if

    do ! form R (retried with a larger regularization)
        ok = form_R(work, regularize_all, eps, min_pivot, max_pivot)
        if (ok) then
            if ((regularize_all .or. work%n_prox > 0) .and. .not. force_prox) then
                ptol = sqrt(zero_tol)
            else
                ptol = zero_tol
            end if
            if (min_pivot <= ptol*max_pivot) ok = .false.
        end if
        if (ok) then
            ! defer the inverse unless cond(H) might be close to the limits
            ! for which finish_Rinv regularizes H
            if (defer_inverse) then
                dmin = daqp_inf
                dmax = 0.0_wp
                do i = 1, n
                    d = work%R(ridx(i,i,n))
                    if (d < dmin) dmin = d
                    if (d > dmax) dmax = d
                end do
                cond = cond_defer_margin*(dmax*dmax)/(dmin*dmin)
                if (ieee_is_finite(dmin) .and. ieee_is_finite(dmax) .and. dmin > 0.0_wp) then
                    if (ieee_is_finite(cond)) then
                        if (real(n,wp)*eps_mach*cond <= hessian_cond_eps) then
                            work%state = ior(work%state, state_cholesky_pending)
                            return
                        end if
                    end if
                end if
            end if
            if (finish_Rinv(work)) return
        end if

        ! regularize the Hessian
        if (regularize_all) then
            if (eps <= 0.0_wp .or. regularization_tries >= 16) then
                flag = daqp_exit_nonconvex
                return
            end if
            regularization_tries = regularization_tries + 1
            eps = eps*2.0_wp
        else
            hessian_scale = hessian_scale_of(work)
            eps = prox_reg_scaled(work, hessian_scale)
            if (eps <= 0.0_wp) then
                flag = daqp_exit_nonconvex
                return
            end if
            regularize_all = .true.
            work%n_prox = n
            work%prox_mask = .true.
        end if
    end do

    end function update_R
!*****************************************************************************************

!*****************************************************************************************
!>
!  The proximal regularization that the factor of `H` was formed with.

    real(wp) function get_proximal_regularization(work) result(eps)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i
    real(wp) :: recovered, rinv

    eps = 0.0_wp
    if (work%n_prox == 0 .or. .not. work%has_H) return

    if (work%rmode == rinv_diag .and. work%n_prox < work%n) then
        i = findloc(work%prox_mask, .true., dim=1)
        eps = 1.0_wp/(work%R(i)*work%R(i)) - work%Hc(i,i)
        return
    end if
    if (work%rmode == rinv_diag) then
        eps = prox_reg_scaled(work, hessian_scale_of(work))
        return
    end if

    ! semi-proximal: recover eps from a regularized row of Rinv (e_i/sqrt(eps))
    if (work%n_prox < work%n) then
        i = findloc(work%prox_mask, .true., dim=1)
        if (i <= work%ms .and. has(work%state,state_rinv_normalized)) then
            rinv = 1.0_wp/work%scaling(i)
        else
            rinv = work%R(ridx(i,i,work%n))
        end if
        eps = 1.0_wp/(rinv*rinv)
        return
    end if

    ! handle the eps-shift correctly for simple bounds
    rinv = work%R(1)
    if (work%ms > 0) rinv = rinv/work%scaling(1)
    recovered = 1.0_wp/(rinv*rinv) - work%Hc(1,1)

    eps = prox_reg_scaled(work, hessian_scale_of(work))
    if (eps <= 0.0_wp) then
        eps = 0.0_wp
        return
    end if
    do while (1.5_wp*eps < recovered)
        eps = eps*2.0_wp
    end do

    end function get_proximal_regularization
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y <-- R'^{-1}... `: transform a linear term in place, `v = Rinv'*f`
!  (upstream's `daqp_update_v`, with `f` given in `v`).

    subroutine transform_v(work,v)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(inout) :: v(:) !! `f` on input, `v` on output (size `n`)

    integer(ip) :: i, j, n, stop_id, dj
    real(wp) :: fj

    n = work%n
    if (work%rmode /= rinv_dense) then ! Rinv = I (or diagonal)
        if (work%rmode == rinv_diag) then
            do i = 1, n
                v(i) = v(i)*work%R(i)
            end do
        end if
        return
    end if
    stop_id = 0
    if (has(work%state,state_rinv_normalized)) stop_id = work%ms
    do j = n, 1, -1
        dj = ridx(j,j,n)
        if (j > stop_id) then
            fj = v(j)
        else ! take the scaling in Rinv into account
            fj = v(j)/work%scaling(j)
        end if
        do i = n, j+1, -1
            v(i) = v(i) + work%R(dj+i-j)*fj
        end do
        v(j) = work%R(dj)*fj
    end do

    end subroutine transform_v
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form `v = Rinv'*f` (upstream's `daqp_update_v(qp->f, work)`).

    subroutine update_v(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (.not. work%has_v .or. .not. work%has_f) return
    work%v = work%f
    call transform_v(work, work%v)

    end subroutine update_v
!*****************************************************************************************

!*****************************************************************************************
!>
!  `Mr(:,k0:k0+nb-1) <-- Rinv'*Mr(:,k0:k0+nb-1)` in place, for `nb <= 4` rows of
!  `A` (upstream's `daqp_rinv_product_block`): the columns of `Rinv` are taken
!  in pairs, in descending order, and each element is summed from its diagonal
!  term downwards, as upstream does.

    subroutine rinv_product_block(work,k0,nb)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: k0 !! first row of A (column of `Mr`)
    integer(ip), intent(in) :: nb !! number of rows (1 to 4)

    integer(ip) :: n, c, j, b, idx, k1, k2, k3
    real(wp) :: s00, s01, s10, s11, s20, s21, s30, s31, r0, r1, rd, x0, x1, x2, x3

    n = work%n
    if (nb == 4) then
        k1 = k0 + 1; k2 = k0 + 2; k3 = k0 + 3
        c = n - 1
        do while (c >= 1) ! columns c and c+1
            rd = work%R(ridx(c+1,c+1,n))
            s01 = work%Mr(c+1,k0)*rd; s11 = work%Mr(c+1,k1)*rd
            s21 = work%Mr(c+1,k2)*rd; s31 = work%Mr(c+1,k3)*rd
            s00 = 0.0_wp; s10 = 0.0_wp; s20 = 0.0_wp; s30 = 0.0_wp
            idx = ridx(c,c,n)
            do j = c, 1, -1
                r0 = work%R(idx)
                r1 = work%R(idx+1)
                x0 = work%Mr(j,k0); x1 = work%Mr(j,k1)
                x2 = work%Mr(j,k2); x3 = work%Mr(j,k3)
                s00 = s00 + x0*r0; s01 = s01 + x0*r1
                s10 = s10 + x1*r0; s11 = s11 + x1*r1
                s20 = s20 + x2*r0; s21 = s21 + x2*r1
                s30 = s30 + x3*r0; s31 = s31 + x3*r1
                idx = idx - (n-j+1)  ! R(j-1,c)
            end do
            work%Mr(c,k0) = s00; work%Mr(c+1,k0) = s01
            work%Mr(c,k1) = s10; work%Mr(c+1,k1) = s11
            work%Mr(c,k2) = s20; work%Mr(c+1,k2) = s21
            work%Mr(c,k3) = s30; work%Mr(c+1,k3) = s31
            c = c - 2
        end do
    else ! the remaining rows, one at a time (the same order of summation)
        do b = k0, k0+nb-1
            c = n - 1
            do while (c >= 1)
                s01 = work%Mr(c+1,b)*work%R(ridx(c+1,c+1,n))
                s00 = 0.0_wp
                idx = ridx(c,c,n)
                do j = c, 1, -1
                    s00 = s00 + work%Mr(j,b)*work%R(idx)
                    s01 = s01 + work%Mr(j,b)*work%R(idx+1)
                    idx = idx - (n-j+1)
                end do
                work%Mr(c,b) = s00
                work%Mr(c+1,b) = s01
                c = c - 2
            end do
        end do
    end if
    if (mod(n,2) == 1) then ! the first column is left over, and has a single term
        rd = work%R(1)
        do b = k0, k0+nb-1
            work%Mr(1,b) = work%Mr(1,b)*rd
        end do
    end if

    end subroutine rinv_product_block
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form `M = A*Rinv` (rows normalized).

    integer(ip) function update_M(work) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n, ns, mA

    n = work%n
    mA = work%m - work%ms
    ! the rows of Rinv of the simple bounds are scaled if Rinv is normalized
    ns = 0
    if (has(work%state,state_rinv_normalized)) ns = work%ms
    select case (work%rmode)
    case (rinv_dense)
        do k = 1, mA
            do j = 1, ns ! undo the scaling in Rinv
                work%Mr(j,k) = work%At(j,k)/work%scaling(j)
            end do
            do j = ns+1, n
                work%Mr(j,k) = work%At(j,k)
            end do
        end do
        ! Mr(:,k) <-- Rinv'*Mr(:,k), in place, four rows of A at a time
        k = 1
        do while (k+3 <= mA)
            call rinv_product_block(work, k, 4)
            k = k + 4
        end do
        if (k <= mA) call rinv_product_block(work, k, mA-k+1)
    case (rinv_diag)
        do k = 1, mA
            do i = 1, n
                work%Mr(i,k) = work%At(i,k)*work%R(i)
            end do
        end do
    case default ! copy A to M
        do k = 1, mA
            work%Mr(:,k) = work%At(:,k)
        end do
    end select

    call daqp_reset_workspace(work) ! internal factorizations need to be redone
    flag = normalize_M(work)

    end function update_M
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form `d = b + M*v` (scaled).

    subroutine update_d(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n
    real(wp) :: s

    n = work%n
    work%reuse_ind = 0 ! the right-hand side changed => cannot reuse intermediate results
    do i = 1, work%m
        work%dupper(i) = work%bupper(i)*work%scaling(i)
        work%dlower(i) = work%blower(i)*work%scaling(i)
    end do

    if (.not. work%has_v) return
    ! simple bounds
    if (work%rmode == rinv_dense) then
        k = 1
        do i = 1, work%ms
            s = 0.0_wp
            do j = i, n
                s = s + work%R(k)*work%v(j)
                k = k + 1
            end do
            work%dupper(i) = work%dupper(i) + s
            work%dlower(i) = work%dlower(i) + s
        end do
    else
        do i = 1, work%ms
            work%dupper(i) = work%dupper(i) + work%v(i)
            work%dlower(i) = work%dlower(i) + work%v(i)
        end do
    end if
    ! general bounds
    do i = work%ms+1, work%m
        k = i - work%ms
        s = 0.0_wp
        do j = 1, n
            s = s + work%Mr(j,k)*work%v(j)
        end do
        work%dupper(i) = work%dupper(i) + s
        work%dlower(i) = work%dlower(i) + s
    end do

    end subroutine update_d
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check the bounds for trivial infeasibility, and detect equality constraints
!  (equal bounds). Returns 1 if the working set has to be activated again,
!  0, or `daqp_exit_infeasible`.

    integer(ip) function check_bounds(work) result(do_activate)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i
    real(wp) :: diff

    do_activate = 0
    do i = 1, work%m
        if (has(work%sense(i),daqp_immutable) .and. .not. has(work%sense(i),daqp_auto_equality)) cycle
        diff = work%bupper(i) - work%blower(i)
        if (diff < -work%settings%primal_tol) then ! trivial infeasibility
            do_activate = daqp_exit_infeasible
            return
        else if (diff < work%settings%zero_tol .and. .not. has(work%sense(i),daqp_soft)) then
            ! unmarked equality constraint (blower == bupper)
            if (.not. has(work%sense(i),daqp_auto_equality) .or. .not. has(work%sense(i),daqp_active)) &
                do_activate = 1
            work%sense(i) = ior(work%sense(i), daqp_active + daqp_immutable + daqp_auto_equality)
        else if (has(work%sense(i),daqp_auto_equality)) then
            work%sense(i) = iand(work%sense(i), not(daqp_active + daqp_immutable + daqp_auto_equality))
            do_activate = 1
        end if
    end do

    end function check_bounds
!*****************************************************************************************

!*****************************************************************************************
!>
!  Normalize the rows of `Rinv` of the simple bounds.

    subroutine normalize_Rinv(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n
    real(wp) :: scaling_i

    if (has(work%state,state_rinv_normalized)) return
    work%state = ior(work%state, state_rinv_normalized)
    if (work%rmode /= rinv_dense) return
    n = work%n
    do i = 1, work%ms
        k = ridx(i,i,n)
        scaling_i = 0.0_wp
        do j = 0, n-i
            scaling_i = scaling_i + work%R(k+j)*work%R(k+j)
        end do
        scaling_i = 1.0_wp/sqrt(scaling_i)
        work%scaling(i) = scaling_i ! needed to retrieve the solution
        do j = 0, n-i
            work%R(k+j) = work%R(k+j)*scaling_i
        end do
    end do

    end subroutine normalize_Rinv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Normalize the general constraints of the LDP. Returns 0 or `daqp_exit_infeasible`.

    integer(ip) function normalize_M(work) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, n
    real(wp) :: scaling_i, zero_tol

    flag = 0
    n = work%n
    zero_tol = work%settings%zero_tol
    do i = work%ms+1, work%m
        k = i - work%ms
        scaling_i = 0.0_wp
        do j = 1, n
            scaling_i = scaling_i + work%Mr(j,k)*work%Mr(j,k)
        end do
        if (scaling_i < zero_tol) then
            ! keep downstream transformations well-defined for constraints
            ! that are omitted from the normalized LDP
            work%scaling(i) = 1.0_wp
            if (work%bupper(i) < -zero_tol .or. work%blower(i) > zero_tol) then
                if (iand(work%sense(i), daqp_immutable+daqp_active) /= daqp_immutable .and. &
                    .not. has(work%sense(i),daqp_soft)) then
                    flag = daqp_exit_infeasible
                    return
                end if
            end if
            work%sense(i) = daqp_immutable ! ignore a zero-row constraint
            cycle
        end if
        scaling_i = 1.0_wp/sqrt(scaling_i)
        work%scaling(i) = scaling_i
        do j = 1, n
            work%Mr(j,k) = work%Mr(j,k)*scaling_i
        end do
    end do

    end function normalize_M
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve `R'v = f`, `R x = -v` with the factor `R` (reciprocal diagonal), for
!  the unconstrained optimum while the inverse is deferred.

    subroutine unconstrained_cholesky(work,vs,feasible)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(inout) :: vs(:) !! scratch for `v` (size `n`)
    logical, intent(inout) :: feasible !! set to false if `x` is not finite

    integer(ip) :: i, j, k, n, off
    real(wp) :: vi, s

    n = work%n
    if (work%has_f) then
        vs(1:n) = work%f
    else
        vs(1:n) = 0.0_wp
    end if
    k = 1
    do i = 1, n ! (row-wise, to traverse R contiguously)
        vs(i) = vs(i)*work%R(k)
        vi = vs(i)
        k = k + 1
        do j = i+1, n
            vs(j) = vs(j) - work%R(k)*vi
            k = k + 1
        end do
    end do
    do i = n, 1, -1
        s = -vs(i)
        off = ridx(i,i,n)
        do j = i+1, n
            s = s - work%R(off+j-i)*work%x(j)
        end do
        work%x(i) = s*work%R(off)
        if (.not. ieee_is_finite(work%x(i))) feasible = .false.
    end do

    end subroutine unconstrained_cholesky
!*****************************************************************************************

!*****************************************************************************************
!>
!  Check whether the unconstrained optimum is feasible (and hence the solution).
!
!  Returns 0 if the unconstrained optimum was not computed, 1 if it was, but is
!  not optimal (`d` is then formed, unscaled), and `unconstrained_optimal` if it is
!  the solution (in `x`).

    integer(ip) function check_unconstrained(work,mask) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: mask  !! update mask

    integer(ip) :: i, j, k, n
    real(wp) :: s, primal_tol
    logical :: feasible

    flag = 0
    if (.not. has(mask,daqp_update_unconstrained)) return
    if (iand(mask, daqp_update_rinv+daqp_update_m+daqp_update_v+daqp_update_d) == 0) return
    if (work%n_prox > 0) return ! not a standard QP
    do i = 1, work%m ! no equalities
        if (iand(work%sense(i), daqp_active+daqp_immutable) /= 0) return
    end do

    n = work%n
    primal_tol = work%settings%primal_tol
    feasible = .true.

    ! compute x_unc, temporarily in x
    call swap_x(work)

    if (has(work%state,state_cholesky_pending)) then
        ! form v only when checking the candidate: R'v = f, then R x = -v
        if (work%has_v) then
            call unconstrained_cholesky(work, work%v, feasible)
        else
            call unconstrained_cholesky(work, work%xldl, feasible)
        end if
    else if (work%has_v) then
        select case (work%rmode)
        case (rinv_dense)
            k = 1
            do i = 1, n
                s = 0.0_wp
                do j = i, n
                    s = s + work%R(k)*work%v(j)
                    k = k + 1
                end do
                work%x(i) = -s
            end do
            if (has(work%state,state_rinv_normalized)) then
                do i = 1, work%ms
                    work%x(i) = work%x(i)/work%scaling(i)
                end do
            end if
        case (rinv_diag)
            do i = 1, n
                work%x(i) = -work%R(i)*work%v(i)
            end do
        case default
            do i = 1, n
                work%x(i) = -work%v(i)
            end do
        end select
    else
        work%x = 0.0_wp ! no linear term: the unconstrained optimum is x = 0
    end if

    ! check the simple bounds
    do i = 1, work%ms
        work%dupper(i) = work%bupper(i) - work%x(i)
        work%dlower(i) = work%blower(i) - work%x(i)
        if (work%dupper(i) < -primal_tol .or. work%dlower(i) > primal_tol) feasible = .false.
    end do
    ! check the general constraints
    do i = work%ms+1, work%m
        s = dot_seq(n, work%At(:,i-work%ms), work%x)
        work%dupper(i) = work%bupper(i) - s
        work%dlower(i) = work%blower(i) - s
        if (work%dupper(i) < -primal_tol .or. work%dlower(i) > primal_tol) feasible = .false.
    end do
    if (feasible) then
        call daqp_reset_workspace(work)
        work%state = ior(work%state, state_unconstrained)
        flag = unconstrained_optimal
        return
    end if
    ! switch back, so that a warm start is preserved
    call swap_x(work)
    flag = 1

    end function check_unconstrained
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form the LDP of the problem in the workspace, as marked by `mask_in`
!  (upstream's `daqp_update_ldp`). Returns 0, or a negative exit flag.

    integer(ip) function daqp_update_ldp(work,mask_in) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: mask_in !! the parts to update (`daqp_update_*`)

    integer(ip) :: mask, unconstrained_flag, i
    logical :: do_activate

    do_activate = .false.
    unconstrained_flag = 0

    ! also form what an earlier update left pending
    mask = ior(mask_in, iand(work%state, state_pending))
    work%state = ior(iand(work%state, state_rinv_normalized + state_ill_conditioned + &
                                      state_cholesky_pending), &
                     iand(mask, state_pending))

    ! update the constraint flags
    if (has(mask,daqp_update_sense)) then
        if (.not. work%has_sense) then ! all constraints are inequalities
            work%sense = 0
        else
            work%sense = work%sense_in
            do_activate = .true.
        end if
    end if

    ! check the bounds early
    if (iand(mask, daqp_update_m+daqp_update_v+daqp_update_d+daqp_update_sense) /= 0) then
        flag = check_bounds(work)
        if (flag < 0) return
        if (flag == 1) do_activate = .true.
    end if

    ! form R first; the inverse is deferred until after the candidate check
    if (has(mask,daqp_update_rinv)) then
        flag = update_R(work, .true.)
        if (flag < 0) return
    end if

    ! update v (if Rinv still holds R, the unconstrained check forms v itself)
    if (.not. has(work%state,state_cholesky_pending) .and. &
        iand(mask, daqp_update_rinv+daqp_update_v) /= 0) call update_v(work)

    unconstrained_flag = check_unconstrained(work, mask)
    if (unconstrained_flag == unconstrained_optimal) then
        ! Rinv (or R), v, and sense are formed, but not M and d, which depend on them
        work%state = iand(work%state, not(daqp_update_rinv + daqp_update_v + daqp_update_sense))
        work%state = ior(work%state, daqp_update_d)
        if (has(mask,daqp_update_rinv)) work%state = ior(work%state, daqp_update_m)
        flag = 0
        return
    end if

    ! a constrained solve needs Rinv (and M formed from it)
    if (has(work%state,state_cholesky_pending)) then
        work%state = ior(work%state, daqp_update_m)
        mask = ior(mask, daqp_update_m)
        if (.not. finish_Rinv(work)) then
            ! refactor with regularization => v and d from the check are stale
            work%state = ior(work%state, daqp_update_rinv)
            mask = ior(mask, daqp_update_rinv)
            flag = update_R(work, .false.)
            if (flag < 0) return
            unconstrained_flag = 0
            call update_v(work)
        else if (unconstrained_flag == 0 .and. iand(mask, daqp_update_rinv+daqp_update_v) /= 0) then
            call update_v(work) ! not formed by the check
        end if
    end if

    ! update M
    if (iand(mask, daqp_update_rinv+daqp_update_m) /= 0) then
        flag = update_M(work)
        if (flag < 0) return
        do_activate = .true. ! update_M cleared the working set
    end if

    call normalize_Rinv(work)

    ! update d
    if (iand(mask, daqp_update_rinv+daqp_update_m+daqp_update_v+daqp_update_d) /= 0) then
        if (unconstrained_flag == 1) then ! d is already computed: normalize it
            do i = 1, work%m
                work%dupper(i) = work%dupper(i)*work%scaling(i)
                work%dlower(i) = work%dlower(i)*work%scaling(i)
            end do
            work%reuse_ind = 0
        else
            call update_d(work)
        end if
    end if

    flag = 0
    ! an empty working set can be one that a reset left out
    if (do_activate .or. work%n_active == 0) then
        call daqp_reset_workspace(work)
        flag = daqp_activate_constraints(work)
    end if
    if (flag < 0) return
    work%state = iand(work%state, not(state_pending)) ! everything has been formed
    flag = 0

    end function daqp_update_ldp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the QP solution in `x` violates a hard constraint by more than `primal_tol`.

    pure logical function violates_hard(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i
    real(wp) :: val, tol

    violates_hard = .false.
    tol = work%settings%primal_tol
    do i = 1, work%m
        if (has(work%sense(i),daqp_soft)) cycle
        if (has(work%sense(i),daqp_immutable) .and. .not. has(work%sense(i),daqp_active)) cycle
        if (i <= work%ms) then
            val = work%x(i)
        else
            val = dot_seq(work%n, work%At(:,i-work%ms), work%x)
        end if
        if (val > work%bupper(i)+tol .or. val < work%blower(i)-tol) then
            violates_hard = .true.
            return
        end if
    end do

    end function violates_hard
!*****************************************************************************************

!*****************************************************************************************
!>
!  Change `eps` of the semi-proximal directions to `eps_new` (only their
!  columns of `Rinv` and `M` are scaled; the working set has to be refactored
!  afterwards).

    subroutine prox_rescale(work,eps,eps_new)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in) :: eps     !! current regularization
    real(wp), intent(in) :: eps_new !! new regularization

    integer(ip) :: i, j, k, n, ms
    real(wp) :: r, norm2, sc

    n = work%n
    ms = work%ms
    do i = 1, n
        if (.not. work%prox_mask(i)) cycle
        r = sqrt((work%Hc(i,i)+eps)/(work%Hc(i,i)+eps_new))
        if (work%rmode /= rinv_dense) then
            work%R(i) = work%R(i)*r
            if (i <= ms) work%scaling(i) = work%scaling(i)/r
        else if (i <= ms .and. has(work%state,state_rinv_normalized)) then
            work%scaling(i) = work%scaling(i)/r
        else
            work%R(ridx(i,i,n)) = work%R(ridx(i,i,n))*r
        end if
    end do
    do i = ms+1, work%m
        k = i - ms
        norm2 = 1.0_wp
        do j = 1, n
            if (.not. work%prox_mask(j) .or. work%Mr(j,k) == 0.0_wp) cycle
            r = sqrt((work%Hc(j,j)+eps)/(work%Hc(j,j)+eps_new))
            norm2 = norm2 + (r*r-1.0_wp)*work%Mr(j,k)*work%Mr(j,k)
            work%Mr(j,k) = work%Mr(j,k)*r
        end do
        if (norm2 == 1.0_wp) cycle
        sc = 1.0_wp/sqrt(norm2)
        do j = 1, n
            work%Mr(j,k) = work%Mr(j,k)*sc
        end do
        work%scaling(i) = work%scaling(i)*sc
    end do

    end subroutine prox_rescale
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the certificate of infeasibility from [[daqp_ldp]] (the dependency
!  `lam_star` of a singular working set) is valid and implies a violation above
!  `primal_tol`.

    pure logical function prox_is_infeasible(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, id
    logical :: lower, wrong_sign
    real(wp) :: q, gap, norm, wrong

    prox_is_infeasible = .false.
    if (work%sing_ind == empty_ind) return
    gap = 0.0_wp
    norm = 0.0_wp
    wrong = 0.0_wp
    do i = 1, work%n_active
        id = work%WS(i)
        lower = has(work%sense(id),daqp_lower)
        q = work%lam_star(i)*work%scaling(id)
        if (lower) then
            gap = gap - q*work%blower(id)
        else
            gap = gap - q*work%bupper(id)
        end if
        if (lower) then
            wrong_sign = q > 0.0_wp
        else
            wrong_sign = q < 0.0_wp
        end if
        if (.not. has(work%sense(id),daqp_immutable) .and. wrong_sign) then
            if (work%bupper(id) < daqp_inf .and. work%blower(id) > -daqp_inf) then
                gap = gap - abs(q)*(work%bupper(id)-work%blower(id))
            else
                wrong = wrong + abs(q)
            end if
        end if
        norm = norm + abs(q)
    end do
    prox_is_infeasible = wrong <= 1.0e-8_wp*norm .and. gap > work%settings%primal_tol*norm

    end function prox_is_infeasible
!*****************************************************************************************

!*****************************************************************************************
!>
!  Project the latest proximal step onto the active face (the steepest descent
!  direction on the face, in the metric of `H+E`), stored as `x-xold`.
!  Returns false if no projection was done.

    logical function prox_project_step(work,eps)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in) :: eps !! the regularization

    integer(ip) :: i, j, k, id, n, na
    real(wp) :: s

    n = work%n
    na = work%n_active
    prox_project_step = .false.
    ! near a vertex, the inner solver resolves the remaining face cheaper
    if (na >= 9*n/10 .or. work%sing_ind /= empty_ind .or. .not. work%has_H) return
    do i = 1, na
        if (has(work%sense(work%WS(i)),daqp_soft)) return
    end do
    ! r = Rinv'*E*(x-xold), in xold
    do i = 1, n
        if (work%prox_mask(i)) then
            work%xold(i) = eps*(work%x(i)-work%xold(i))
        else
            work%xold(i) = 0.0_wp
        end if
    end do
    call transform_v(work, work%xold)
    ! r <-- r - M_W'*(M_W*M_W')^{-1}*M_W*r
    do i = 1, na
        id = work%WS(i)
        if (id > work%ms) then
            s = dot_seq(n, work%Mr(:,id-work%ms), work%xold)
        else if (work%rmode == rinv_dense) then
            s = 0.0_wp
            k = ridx(id,id,n)
            do j = id, n
                s = s + work%R(k)*work%xold(j)
                k = k + 1
            end do
        else
            s = work%xold(id)
        end if
        work%xldl(i) = s
    end do
    call solve_working_set(work)
    call sub_working_set_rows(work, work%xldl, work%xold)
    call apply_Rinv(work, work%xold)
    do i = 1, n
        work%xold(i) = work%x(i) - work%xold(i)
    end do
    work%reuse_ind = 0
    prox_project_step = .true.

    end function prox_project_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  Step length `-g'd/d'Hd` that minimizes the objective along `d = x-xold`
!  (`d` in `xldl`). Returns `daqp_inf` if `d` has no curvature, -1 if `d` is not
!  a descent direction.

    real(wp) function prox_curvature_step(work) result(step)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, n
    real(wp) :: gd, dhd, dd, hmax, hii, s

    n = work%n
    gd = 0.0_wp; dhd = 0.0_wp; dd = 0.0_wp; hmax = 0.0_wp
    do i = 1, n
        work%xldl(i) = work%x(i) - work%xold(i)
    end do
    if (work%has_f) then
        do i = 1, n
            gd = gd + work%f(i)*work%xldl(i)
        end do
    end if
    if (work%has_H) then
        do i = 1, n
            hii = abs(work%Hc(i,i))
            s = 0.0_wp
            do j = 1, n
                s = s + work%Hc(j,i)*work%xldl(j)
            end do
            work%zldl(i) = s
            if (hii > hmax) hmax = hii
        end do
        do i = 1, n
            gd = gd + work%x(i)*work%zldl(i)
            dhd = dhd + work%xldl(i)*work%zldl(i)
            dd = dd + work%xldl(i)*work%xldl(i)
        end do
    end if
    if (gd >= 0.0_wp) then
        step = -1.0_wp
    else if (dhd > work%settings%zero_tol*hmax*dd) then
        step = -gd/dhd
    else
        step = daqp_inf
    end if

    end function prox_curvature_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  First inactive constraint that blocks `x + s*(x-xold)` for `s < step`
!  (`step` is shortened to the blocking step).

    integer(ip) function prox_blocking_constraint(work,step,lower) result(ind)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(inout) :: step  !! step length
    logical, intent(inout) :: lower  !! the blocking bound is a lower one

    integer(ip) :: i, j, k, n
    real(wp) :: ad, ax, sb

    ind = empty_ind
    n = work%n
    do i = 1, work%m
        if (iand(work%sense(i), daqp_active+daqp_immutable+daqp_set_aside) /= 0) cycle
        if (i <= work%ms) then
            ax = work%x(i)
            ad = ax - work%xold(i)
        else
            k = i - work%ms
            ad = 0.0_wp
            ax = 0.0_wp
            do j = 1, n
                ax = ax + work%At(j,k)*work%x(j)
                ad = ad + work%At(j,k)*(work%x(j)-work%xold(j))
            end do
        end if
        if (ad > 0.0_wp .and. work%bupper(i) < daqp_inf) then
            sb = (work%bupper(i)-ax)/ad
        else if (ad < 0.0_wp .and. work%blower(i) > -daqp_inf) then
            sb = (work%blower(i)-ax)/ad
        else
            cycle
        end if
        if (sb < step) then
            step = sb
            lower = ad < 0.0_wp
            ind = i
        end if
    end do

    end function prox_blocking_constraint
!*****************************************************************************************

!*****************************************************************************************
!>
!  Step along the latest proximal step `d = x-xold`: move `x` to the minimizer
!  along `d`, or to the first blocking constraint, which is added to the working
!  set (dependent ones are set aside).
!
!  Returns 1 if `x` was moved, 0 otherwise, and `daqp_exit_unbounded` for an LP
!  with an unblocked descent direction.

    integer(ip) function prox_step(work,s_prev,eps,projected_in) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(inout) :: s_prev   !! exact step length of the latest step (-1: none)
    real(wp), intent(in) :: eps         !! the regularization
    logical, intent(in) :: projected_in !! the step is projected onto the active face

    integer(ip) :: i, k, ind, id, n_projections
    logical :: lower, moved, skipped, first, null_direction, projected, leave
    real(wp) :: s, s_exact, ad

    lower = .false.
    moved = .false.
    skipped = .false.
    first = .true.
    projected = projected_in
    n_projections = 0
    do
        s = prox_curvature_step(work)
        if (s < 0.0_wp) exit
        null_direction = s >= daqp_inf
        if (first) then ! lagged (Barzilai-Borwein) step length
            s_exact = s
            if (s_exact < daqp_inf .and. s_prev >= 0.0_wp) then
                if (s_prev < 2.0_wp*s_exact) then
                    s = s_prev
                else
                    s = 2.0_wp*s_exact
                end if
            end if
            if (s_exact < daqp_inf) then
                s_prev = s_exact
            else
                s_prev = -1.0_wp
            end if
            first = .false.
        end if
        ind = prox_blocking_constraint(work, s, lower)
        ! roundoff in a projected direction is amplified by a long step:
        ! keep the inner iterate if the step would leave the active face
        if (projected .and. s > 0.0_wp .and. s < daqp_inf) then
            leave = .false.
            do i = 1, work%n_active
                id = work%WS(i)
                if (id <= work%ms) then
                    ad = work%xldl(id)
                else
                    ad = dot_seq(work%n, work%At(:,id-work%ms), work%xldl)
                end if
                if (abs(s*ad) > work%settings%primal_tol) then
                    leave = .true.
                    exit
                end if
            end do
            if (leave) exit
        end if
        if (ind == empty_ind) then
            if (s < daqp_inf) then ! the minimizer along d
                do k = 1, work%n
                    work%x(k) = work%x(k) + s*work%xldl(k)
                end do
                moved = .true.
            else if (.not. work%has_H .and. .not. moved .and. .not. skipped) then
                flag = daqp_exit_unbounded
                return
            end if
            exit
        end if
        s_prev = -1.0_wp ! blocked: the working set changes
        ! advance to the blocking constraint and activate it
        if (s >= 0.0_wp) then
            do k = 1, work%n
                work%x(k) = work%x(k) + s*work%xldl(k)
            end do
        end if
        moved = .true.
        if (lower) then
            work%sense(ind) = ior(work%sense(ind), daqp_lower)
            call add_constraint(work, ind, -1.0_wp)
        else
            work%sense(ind) = iand(work%sense(ind), not(daqp_lower))
            call add_constraint(work, ind, 1.0_wp)
        end if
        if (work%sing_ind == empty_ind) then
            ! along a null direction of H, continue on the new face
            if (null_direction .and. work%settings%eps_prox < 0.0_wp .and. &
                work%nh >= prox_face_start) then
                n_projections = n_projections + 1
                if (n_projections-1 < prox_face_steps) then
                    if (prox_project_step(work, eps)) then
                        projected = .true.
                        cycle
                    end if
                end if
            end if
            exit
        end if
        ! linearly dependent on the active constraints: set it aside
        id = drop_singular_last(work)
        work%sense(id) = ior(work%sense(id), daqp_set_aside)
        skipped = .true.
    end do
    if (skipped) then
        do i = 1, work%m
            work%sense(i) = iand(work%sense(i), not(daqp_set_aside))
        end do
    end if
    flag = 0
    if (moved) flag = 1

    end function prox_step
!*****************************************************************************************

!*****************************************************************************************
!>
!  The outer proximal-point (or semi-proximal) loop, for a semidefinite `H`
!  or an LP (upstream's `daqp_prox`). Returns the exit flag.

    integer(ip) function daqp_prox(work) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, total_iter, nx, step_flag
    logical :: center_relaxed, is_lp, all_pd, adaptive, rescaled, converged, projected
    real(wp) :: s_prev, max_diff, tol_stat, eta, eps, eps_max, hmax, eps0, prox_norm

    total_iter = 0
    s_prev = -1.0_wp
    center_relaxed = .false.
    exitflag = daqp_exit_iterlimit ! if no iteration can be taken
    eta = work%settings%eta_prox

    work%nh = 0 ! counts the outer iterations
    nx = work%n
    is_lp = work%rmode == rinv_none
    if (is_lp) then
        eps = 1.0_wp
    else
        eps = get_proximal_regularization(work)
    end if

    ! for a positive definite H, the inner QP equals the original problem
    all_pd = .not. is_lp .and. work%n_prox == 0

    ! eps of semi-proximal directions can be changed cheaply (prox_rescale);
    ! a failed inner problem (often due to a small eps) is resolved with eps_max
    adaptive = .not. is_lp .and. .not. all_pd .and. work%n_prox < nx
    eps_max = eps
    rescaled = .false.
    if (adaptive) then
        hmax = hessian_scale_of(work)
        ! reset an eps that an earlier solve has raised
        eps0 = abs(work%settings%eps_prox)
        if (eps0 < sqrt(work%settings%zero_tol)*hmax) eps0 = sqrt(work%settings%zero_tol)*hmax
        if (eps > 1.01_wp*eps0) then
            call prox_rescale(work, eps, eps0)
            eps = eps0
            rescaled = .true.
        end if
        if (prox_eps_max*hmax > eps) then
            eps_max = prox_eps_max*hmax
        else
            eps_max = eps
        end if
    end if

    ! a negative eta selects an automatic tolerance
    if (.not. all_pd .and. eta < 0.0_wp) then
        eta = auto_eta_cap
        if (work%settings%dual_tol /= daqp_default_dual_tol .and. 0.1_wp*work%settings%dual_tol < eta) &
            eta = 0.1_wp*work%settings%dual_tol
    end if

    do while (total_iter < work%settings%iter_limit)

        ! perturb the problem: form v = R'\(f - eps_mask*x_old)
        if (is_lp) then
            if (total_iter > 0) then
                if (work%iterations == 1) then
                    eps = eps*10.0_wp
                else
                    eps = eps*0.9_wp
                end if
            end if
            if (eps > 1.0e3_wp) eps = 1.0e3_wp
            do i = 1, nx
                work%v(i) = work%f(i)*eps - work%x(i)
            end do
        else
            if (work%n_prox == nx) then ! full shift
                if (work%has_f) then
                    do i = 1, nx
                        work%v(i) = work%f(i) - eps*work%x(i)
                    end do
                else
                    do i = 1, nx
                        work%v(i) = -eps*work%x(i)
                    end do
                end if
            else ! regularize only the singular directions
                do i = 1, nx
                    if (work%prox_mask(i)) then
                        if (work%has_f) then
                            work%v(i) = work%f(i) - eps*work%x(i)
                        else
                            work%v(i) = -eps*work%x(i)
                        end if
                    else
                        if (work%has_f) then
                            work%v(i) = work%f(i) - 0.0_wp*work%x(i)
                        else
                            work%v(i) = -0.0_wp*work%x(i)
                        end if
                    end if
                end do
            end if
            call transform_v(work, work%v)
        end if

        call update_d(work)
        if (rescaled) then ! the working set is factored anew after a rescaling
            call daqp_reset_workspace(work)
            i = daqp_activate_constraints(work)
            rescaled = .false.
        end if

        call swap_x(work) ! xold <-- x

        ! solve the (regularized) least-distance problem
        work%nh = work%nh + 1
        exitflag = daqp_ldp(work)

        total_iter = total_iter + work%iterations
        if (adaptive .and. eps < eps_max .and. total_iter < work%settings%iter_limit) then
            if (exitflag == daqp_exit_cycle .or. &
                (exitflag == daqp_exit_infeasible .and. .not. prox_is_infeasible(work))) then
                call prox_rescale(work, eps, eps_max)
                eps = eps_max
                rescaled = .true.
                s_prev = -1.0_wp
                work%x = work%xold ! the center
                cycle
            end if
        end if
        if (exitflag < 0) exit ! inner solver failed
        call daqp_ldp2qp_solution(work)

        if (eps == 0.0_wp) exit ! no regularization -> single outer step

        ! H positive definite: the first solve gives the exact solution
        if (all_pd) then
            exitflag = daqp_exit_optimal
            exit
        end if

        ! convergence check: fixed point ||x - x_old||_inf < tol_stat
        if (is_lp) then
            tol_stat = eta*eps
        else
            tol_stat = eta/eps
        end if
        converged = .true.
        do i = 1, nx
            max_diff = work%x(i) - work%xold(i)
            if (max_diff > tol_stat .or. max_diff < -tol_stat) then
                converged = .false.
                exit
            end if
        end do
        if (converged) then
            if (center_relaxed .and. total_iter < work%settings%iter_limit) then
                center_relaxed = .false.
                cycle ! confirm convergence from the feasible iterate
            end if
            exitflag = daqp_exit_optimal
            exit
        end if

        ! unchanged working set => accelerate by moving the center along the step;
        ! after many outer steps, also accelerate small working-set changes,
        ! projecting the step onto the new face first
        center_relaxed = .false.
        projected = .false.
        if (work%iterations /= 1) then ! the working set has changed
            s_prev = -1.0_wp
            if (.not. is_lp .and. work%settings%eps_prox < 0.0_wp .and. &
                work%nh >= prox_face_start .and. work%iterations <= 3) then
                projected = prox_project_step(work, eps)
            end if
        end if
        if ((work%iterations == 1 .or. projected) .and. work%n_active < nx .and. &
            total_iter < work%settings%iter_limit) then
            step_flag = prox_step(work, s_prev, eps, projected)
            if (step_flag == daqp_exit_unbounded) then
                exitflag = daqp_exit_unbounded
                exit
            end if
            center_relaxed = step_flag == 1
        end if
    end do

    ! finalize
    if (total_iter >= work%settings%iter_limit) exitflag = daqp_exit_iterlimit
    ! refine x (skipped for short solves)
    if (exitflag > 0 .and. total_iter > refine_min_iter) call refine_primal(work)
    if (is_lp) then
        do i = 1, work%n_active
            work%lam_star(i) = work%lam_star(i)/eps ! rescale the dual variables
        end do
    else
        ! correct the regularized objective
        prox_norm = 0.0_wp
        do i = 1, nx
            if (work%prox_mask(i)) prox_norm = prox_norm + work%x(i)*work%x(i)
        end do
        work%fval = work%fval + eps*prox_norm
    end if
    work%iterations = total_iter

    end function daqp_prox
!*****************************************************************************************

!*****************************************************************************************
!>
!  Deallocate the workspace.

    subroutine daqp_destroy(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    type(daqp_settings) :: settings

    settings = work%settings
    work = daqp_workspace()  ! deallocates every component
    work%settings = settings

    end subroutine daqp_destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate `a(n)`, unless it already has that size. A nonzero `istat` is kept.

    subroutine resize1(a,n,istat)

    real(wp), allocatable, intent(inout) :: a(:) !! array
    integer(ip), intent(in) :: n                 !! size
    integer(ip), intent(inout) :: istat          !! status (nonzero if an allocation failed)

    integer :: stat

    if (allocated(a)) then
        if (size(a) == n) return
        deallocate(a)
    end if
    allocate(a(n), stat=stat)
    if (stat /= 0) istat = stat

    end subroutine resize1
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate `a(n1,n2)`, unless it already has that shape. A nonzero `istat` is kept.

    subroutine resize2(a,n1,n2,istat)

    real(wp), allocatable, intent(inout) :: a(:,:) !! array
    integer(ip), intent(in) :: n1                  !! number of rows
    integer(ip), intent(in) :: n2                  !! number of columns
    integer(ip), intent(inout) :: istat            !! status (nonzero if an allocation failed)

    integer :: stat

    if (allocated(a)) then
        if (size(a,1) == n1 .and. size(a,2) == n2) return
        deallocate(a)
    end if
    allocate(a(n1,n2), stat=stat)
    if (stat /= 0) istat = stat

    end subroutine resize2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate `a(n)`, unless it already has that size. A nonzero `istat` is kept.

    subroutine resize1i(a,n,istat)

    integer(ip), allocatable, intent(inout) :: a(:) !! array
    integer(ip), intent(in) :: n                    !! size
    integer(ip), intent(inout) :: istat             !! status (nonzero if an allocation failed)

    integer :: stat

    if (allocated(a)) then
        if (size(a) == n) return
        deallocate(a)
    end if
    allocate(a(n), stat=stat)
    if (stat /= 0) istat = stat

    end subroutine resize1i
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate `a(n)`, unless it already has that size. A nonzero `istat` is kept.

    subroutine resize1l(a,n,istat)

    logical, allocatable, intent(inout) :: a(:) !! array
    integer(ip), intent(in) :: n                !! size
    integer(ip), intent(inout) :: istat         !! status (nonzero if an allocation failed)

    integer :: stat

    if (allocated(a)) then
        if (size(a) == n) return
        deallocate(a)
    end if
    allocate(a(n), stat=stat)
    if (stat /= 0) istat = stat

    end subroutine resize1l
!*****************************************************************************************

!*****************************************************************************************
!>
!  `B = A'` (blocked, to limit the strided accesses).

    subroutine transpose_into(A,B)

    real(wp), intent(in) :: A(:,:)  !! `(p,q)`
    real(wp), intent(out) :: B(:,:) !! `(q,p)`

    integer(ip), parameter :: nb = 32
    integer(ip) :: i, j, i0, j0

    do j0 = 1, size(A,2), nb
        do i0 = 1, size(A,1), nb
            do i = i0, min(i0+nb-1, int(size(A,1),ip))
                do j = j0, min(j0+nb-1, int(size(A,2),ip))
                    B(j,i) = A(i,j)
                end do
            end do
        end do
    end do

    end subroutine transpose_into
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up the workspace for a QP, and form its LDP (upstream's `setup_daqp`).
!  The settings are taken from `work%settings`, which the caller sets first.
!
!  The data is given in the Fortran layout: `H(n,n)`, `A(m-ms,n)`.
!  Returns 1, or a negative exit flag.

    integer(ip) function daqp_setup(work,n,m,ms,bupper,blower,H,f,A,sense,init_mask,At) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: n   !! number of variables
    integer(ip), intent(in) :: m   !! number of constraints (including the simple bounds)
    integer(ip), intent(in) :: ms  !! number of simple bounds
    real(wp), intent(in) :: bupper(m) !! upper bounds
    real(wp), intent(in) :: blower(m) !! lower bounds
    real(wp), intent(in), optional :: H(n,n)    !! Hessian (absent: an LP)
    real(wp), intent(in), optional :: f(n)      !! linear term
    real(wp), intent(in), optional :: A(m-ms,n) !! constraint matrix (absent if `m == ms`)
    integer(ip), intent(in), optional :: sense(m) !! constraint flags
    integer(ip), intent(in), optional :: init_mask !! extra update mask (`daqp_update_unconstrained`)
    real(wp), intent(in), optional :: At(n,m-ms) !! the transpose of `A` (instead of `A`: a row of `A` per
                                                 !! column, the internal layout, which saves a transposed copy)

    integer(ip) :: istat, nw, mask

    work%is_setup = .false.
    work%n = n
    work%m = m
    work%ms = ms
    work%has_H = present(H)
    work%has_f = present(f)
    work%has_sense = present(sense)
    work%ns = 0
    if (present(sense)) work%ns = count(iand(sense, daqp_soft) /= 0)

    ! allocate (arrays of the right size are kept, for repeated setups)
    nw = n + work%ns  ! to account for soft constraints
    ! the packed triangles are indexed by default integers
    if ((int(nw,int64)+1_int64)*(int(nw,int64)+2_int64)/2_int64 > int(huge(1_ip),int64) .or. &
        int(n,int64)*int(max(n,m-ms),int64) > int(huge(1_ip),int64)) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if
    istat = 0
    ! the problem
    if (work%has_H) then
        call resize2(work%Hc, n, n, istat)
    else
        call resize2(work%Hc, 0, 0, istat)
    end if
    call resize1(work%f, n, istat)
    call resize2(work%At, n, m-ms, istat)
    call resize1(work%bupper, m, istat)
    call resize1(work%blower, m, istat)
    call resize1i(work%sense_in, m, istat)
    ! the iterates
    call resize1(work%lam, nw+1, istat)
    call resize1(work%lam_star, nw+1, istat)
    call resize1i(work%WS, nw+1, istat)
    call resize1(work%D, nw+1, istat)
    call resize1(work%xldl, nw+1, istat)
    call resize1(work%zldl, nw+1, istat)
    call resize1(work%L, ((nw+1)*(nw+2))/2, istat)
    call resize1(work%x, n, istat)
    call resize1(work%xold, n, istat)
    call resize1l(work%prox_mask, n, istat)
    ! the LDP
    call resize1(work%scaling, m, istat)
    call resize2(work%Mr, n, m-ms, istat)
    call resize1(work%Mu, m-ms, istat)
    call resize1(work%dupper, m, istat)
    call resize1(work%dlower, m, istat)
    call resize1i(work%sense, m, istat)
    call resize1(work%v, n, istat)
    if (work%has_H) then
        call resize1(work%R, (n*(n+1))/2, istat)
    else
        call resize1(work%R, 0, istat)
    end if
    if (istat /= 0) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if

    if (work%has_H) work%Hc = transpose(H)
    if (work%has_f) then
        work%f = f
    else
        work%f = 0.0_wp
    end if
    if (present(At)) then
        work%At = At
    else if (present(A)) then
        call transpose_into(A, work%At)
    else
        work%At = 0.0_wp
    end if
    work%bupper = bupper
    work%blower = blower
    work%sense_in = 0
    if (work%has_sense) work%sense_in = sense

    ! (the work arrays are written before they are read, as upstream's malloc'ed ones)
    work%x = 0.0_wp  ! an uninitialized iterate is 0
    work%xold = 0.0_wp
    work%D(1) = 0.0_wp
    work%prox_mask = .false.
    work%n_prox = 0
    work%state = 0
    work%soft_slack = 0.0_wp
    work%nh = 1
    work%fval = 0.0_wp
    call daqp_reset_workspace(work)

    work%scaling = 1.0_wp
    work%v = 0.0_wp
    work%sense = 0
    work%has_v = work%has_f
    if (work%has_H) then
        work%rmode = rinv_dense
    else
        work%rmode = rinv_none
    end if

    ! always update M, d and sense
    mask = daqp_update_m + daqp_update_d + daqp_update_sense
    if (present(init_mask)) mask = ior(mask, init_mask)
    if (work%has_H) mask = ior(mask, daqp_update_rinv)
    if (work%has_f) mask = ior(mask, daqp_update_v)

    ! for an LP, mark all directions as needing proximal regularization
    if (.not. work%has_H .and. work%has_f) work%n_prox = n

    istat = daqp_update_ldp(work, mask)
    if (istat < 0) then
        call daqp_destroy(work)
        flag = istat
        return
    end if

    ! a singular quadratic needs v for the proximal linear term even without f
    if (work%n_prox > 0 .and. .not. work%has_v .and. work%has_H) then
        work%has_v = .true.
        work%v = 0.0_wp
    end if

    work%is_setup = .true.
    flag = 1

    end function daqp_setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the working set for a warm start: the constraints `active` (at their
!  lower bound where `at_lower`). Equality (immutable) constraints are kept.
!  Returns 1, or `daqp_exit_overdetermined_initial`.

    integer(ip) function daqp_set_working_set(work,active,at_lower) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: active(:)  !! indices of the active constraints
    logical, intent(in) :: at_lower(:)    !! whether each is active at its lower bound

    integer(ip) :: i, id

    do i = 1, work%m
        if (has(work%sense(i),daqp_immutable)) cycle
        work%sense(i) = iand(work%sense(i), not(daqp_active))
    end do
    do i = 1, size(active)
        id = active(i)
        if (has(work%sense(id),daqp_immutable)) cycle
        work%sense(id) = ior(work%sense(id), daqp_active)
        if (at_lower(i)) then
            work%sense(id) = ior(work%sense(id), daqp_lower)
        else
            work%sense(id) = iand(work%sense(id), not(daqp_lower))
        end if
    end do
    call daqp_reset_workspace(work)
    flag = daqp_activate_constraints(work)

    end function daqp_set_working_set
!*****************************************************************************************

!*****************************************************************************************
!>
!  Refactor the working set after a change of `rho_soft` or `w_soft` on a
!  live workspace (upstream's `daqp_refresh_soft_weights`).

    subroutine daqp_refresh_soft_weights(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, flag

    ! only an active soft constraint makes the factorization stale
    do i = 1, work%n_active
        if (has(work%sense(work%WS(i)),daqp_soft)) then
            call daqp_reset_workspace(work)
            flag = daqp_activate_constraints(work)
            return
        end if
    end do

    end subroutine daqp_refresh_soft_weights
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the QP from the current working set (upstream's `daqp_solve`).

    subroutine daqp_solve(work,x,res,lam)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(out) :: x(:)               !! the solution (size `n`)
    type(daqp_result), intent(out) :: res       !! exit flag, objective, iterations
    real(wp), intent(out), optional :: lam(:)   !! multipliers (size `m`)

    integer(ip) :: i, id
    real(wp) :: l
    logical :: wrong_sign

    if (.not. work%is_setup) then
        res%exitflag = daqp_exit_not_setup
        x = 0.0_wp
        if (present(lam)) lam = 0.0_wp
        return
    end if

    work%nh = 1
    if (.not. has(work%state,state_unconstrained)) then
        if (work%n_prox == 0) then
            res%exitflag = daqp_ldp(work)
            if (res%exitflag > 0) then
                call daqp_ldp2qp_solution(work) ! retrieve the QP solution
                call refine_primal(work)        ! refine x (if it might be inaccurate)
                ! constraints were only added above the rounding level
                if (has(work%state,state_noise_floor)) then
                    if (violates_hard(work)) res%exitflag = daqp_exit_optimal_inexact
                end if
            end if
        else
            res%exitflag = daqp_prox(work)
        end if
    else ! unconstrained optimum
        work%iterations = 1
        work%fval = 0.0_wp
        work%soft_slack = 0.0_wp
        res%exitflag = daqp_exit_optimal
    end if

    ! package the result
    x = work%x
    if (present(lam)) then
        lam = 0.0_wp
        do i = 1, work%n_active
            id = work%WS(i)
            l = work%lam_star(i)
            ! report a multiplier of the wrong sign (within dual_tol) as zero
            if (has(work%sense(id),daqp_lower)) then
                wrong_sign = l > 0.0_wp
            else
                wrong_sign = l < 0.0_wp
            end if
            if (.not. has(work%sense(id),daqp_immutable) .and. &
                .not. has(work%sense(id),daqp_soft) .and. wrong_sign) then
                lam(id) = 0.0_wp
            else
                lam(id) = l
            end if
        end do
    end if

    ! shift back the function value
    if (work%has_v .and. work%rmode /= rinv_none) then ! QP
        res%fval = work%fval
        do i = 1, work%n
            res%fval = res%fval - work%v(i)*work%v(i)
        end do
        res%fval = res%fval*0.5_wp
    else if (work%has_f) then ! LP
        res%fval = 0.0_wp
        do i = 1, work%n
            res%fval = res%fval + work%f(i)*work%x(i)
        end do
    else ! no linear term (upstream leaves this undefined)
        res%fval = 0.5_wp*work%fval
    end if

    res%soft_slack = work%soft_slack
    res%iter = work%iterations
    res%nodes = work%nh

    end subroutine daqp_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up and solve a QP in one call (upstream's `daqp_quadprog`, with the
!  check of the unconstrained optimum, and without the equality elimination).

    subroutine daqp_quadprog(n,m,ms,bupper,blower,x,res,lam,H,f,A,sense,settings)

    integer(ip), intent(in) :: n   !! number of variables
    integer(ip), intent(in) :: m   !! number of constraints (including the simple bounds)
    integer(ip), intent(in) :: ms  !! number of simple bounds
    real(wp), intent(in) :: bupper(m) !! upper bounds
    real(wp), intent(in) :: blower(m) !! lower bounds
    real(wp), intent(out) :: x(n)     !! the solution
    type(daqp_result), intent(out) :: res     !! exit flag, objective, iterations
    real(wp), intent(out), optional :: lam(m) !! multipliers
    real(wp), intent(in), optional :: H(n,n)    !! Hessian (absent: an LP)
    real(wp), intent(in), optional :: f(n)      !! linear term
    real(wp), intent(in), optional :: A(m-ms,n) !! constraint matrix
    integer(ip), intent(in), optional :: sense(m) !! constraint flags
    type(daqp_settings), intent(in), optional :: settings !! settings (default: upstream's)

    type(daqp_workspace) :: work
    integer(ip) :: flag

    if (present(settings)) work%settings = settings
    flag = daqp_setup(work, n, m, ms, bupper, blower, H, f, A, sense, &
                      init_mask=daqp_update_unconstrained)
    if (flag < 0) then
        res%exitflag = flag
        x = 0.0_wp
        if (present(lam)) lam = 0.0_wp
        return
    end if
    call daqp_solve(work, x, res, lam)
    call daqp_destroy(work)

    end subroutine daqp_quadprog
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_core
!*****************************************************************************************
