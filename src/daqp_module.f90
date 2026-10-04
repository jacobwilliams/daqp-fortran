!*****************************************************************************************
!> author: Jacob Williams
!
!  Object-oriented interface to the Fortran port of
!  [DAQP](https://github.com/darnstrom/daqp), a dual active-set solver for
!  dense convex quadratic programs:
!
!      min  0.5 x'Hx + f'x
!      s.t. blower(1:ms)  <= x(1:ms) <= bupper(1:ms)
!           blower(ms+1:) <= A x     <= bupper(ms+1:)
!
!  `H` must be symmetric positive definite, or positive semidefinite (then
!  a proximal-point outer loop is used); with no `H`, the problem is an LP.
!
!### Multipliers
!
!  The multipliers satisfy \( Hx + f + A_{all}^T \lambda = 0 \), with
!  \( A_{all} = [I_{ms}\; 0; A] \): \( \lambda_i \ge 0 \) when the upper bound
!  of row `i` is active, \( \lambda_i \le 0 \) when its lower bound is active,
!  and 0 for an inactive row. [[daqp_lower_positive]] converts them to the
!  convention \( \lambda \ge 0 \) at a lower bound (\( Hx + f = A_{all}^T \lambda \)).
!
!### Example
!
!```fortran
!    type(daqp_type) :: qp
!    call qp%setup(H, f, A, bupper, blower, istat)
!    call qp%solve(x, lam, istat)
!```

    module daqp_module

    use daqp_kinds
    use daqp_core

    implicit none

    private

    integer, parameter :: wp = daqp_wp
    integer, parameter :: ip = daqp_ip

    ! re-export
    public :: daqp_wp, daqp_ip, daqp_inf
    public :: daqp_settings, daqp_result, daqp_workspace
    public :: daqp_active, daqp_lower, daqp_immutable, daqp_soft, daqp_equality

    ! statuses (0 or positive = success)
    integer(ip), parameter, public :: daqp_success               = 0   !! setup or update succeeded
    integer(ip), parameter, public :: daqp_optimal               = daqp_exit_optimal         !! solved
    integer(ip), parameter, public :: daqp_soft_optimal          = daqp_exit_soft_optimal    !! solved, a soft constraint is violated
    integer(ip), parameter, public :: daqp_optimal_inexact       = daqp_exit_optimal_inexact !! solved after cycling, violates a constraint slightly
    integer(ip), parameter, public :: daqp_infeasible            = daqp_exit_infeasible      !! the problem is infeasible
    integer(ip), parameter, public :: daqp_cycling               = daqp_exit_cycle           !! cycling
    integer(ip), parameter, public :: daqp_unbounded             = daqp_exit_unbounded       !! the problem is unbounded (LP)
    integer(ip), parameter, public :: daqp_iteration_limit       = daqp_exit_iterlimit       !! the iteration limit was reached
    integer(ip), parameter, public :: daqp_nonconvex             = daqp_exit_nonconvex       !! `H` could not be factored
    integer(ip), parameter, public :: daqp_overdetermined        = daqp_exit_overdetermined_initial !! inconsistent equalities in the initial working set
    integer(ip), parameter, public :: daqp_invalid_input         = daqp_exit_invalid_input   !! invalid input
    integer(ip), parameter, public :: daqp_out_of_memory         = daqp_exit_out_of_memory   !! an allocation failed
    integer(ip), parameter, public :: daqp_not_setup             = daqp_exit_not_setup       !! `setup` was not called (or failed)

    type, public :: daqp_type
        !! A DAQP solver for one QP (copyable).

        ! options (upstream's settings, with its defaults)
        real(wp)    :: primal_tol   = daqp_default_prim_tol  !! tolerance for primal feasibility
        real(wp)    :: dual_tol     = daqp_default_dual_tol  !! tolerance for dual feasibility
        real(wp)    :: zero_tol     = daqp_default_zero_tol  !! values below are regarded as zero
        real(wp)    :: pivot_tol    = daqp_default_pivot_tol !! pivots of the `LDL'` below are reordered
        real(wp)    :: progress_tol = daqp_default_prog_tol  !! minimum objective progress (cycle guard)
        integer(ip) :: cycle_tol    = daqp_default_cycle_tol !! iterations without progress before cycling is assumed
        integer(ip) :: iter_limit   = daqp_default_iter_limit !! maximum number of iterations
        real(wp)    :: fval_bound   = daqp_inf   !! the problem is regarded as infeasible if the objective exceeds this
        real(wp)    :: eps_prox     = daqp_default_eps_prox  !! proximal regularization (negative: automatic, only if `H` is singular)
        real(wp)    :: eta_prox     = daqp_default_eta       !! tolerance of the proximal loop (negative: automatic)
        real(wp)    :: rho_soft     = daqp_default_rho_soft  !! soft constraints: the penalty is `s**2/(2*rho_soft)` per unit slack `s` (normalized)
        real(wp)    :: w_soft       = daqp_default_w_soft    !! soft constraints: linear penalty `w_soft*s` (0: purely quadratic)
        real(wp)    :: sing_tol     = daqp_default_sing_tol  !! pivots below mark a singular working set
        real(wp)    :: refactor_tol = daqp_default_refactor_tol !! pivots below trigger a refactorization at a solution

        ! results of the latest solve
        integer(ip) :: status     = daqp_not_setup !! status of the latest call
        integer(ip) :: iter       = 0      !! iterations
        integer(ip) :: outer_iter = 0      !! outer (proximal) iterations
        real(wp)    :: fval       = 0.0_wp !! objective function value
        real(wp)    :: soft_slack = 0.0_wp !! largest violation of a soft constraint

        type(daqp_workspace), private :: work !! the solver's workspace

    contains

        procedure, public :: setup
        procedure, public :: solve
        procedure, public :: update
        procedure, public :: set_working_set
        procedure, public :: get_working_set
        procedure, public :: info
        procedure, public :: set_defaults
        procedure, public :: destroy
        procedure, public :: n_variables
        procedure, public :: n_constraints
        procedure, private :: copy_settings

    end type daqp_type

    public :: daqp_lower_positive

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Restore the default options.

    subroutine set_defaults(me)

    class(daqp_type), intent(inout) :: me

    type(daqp_settings) :: s

    me%primal_tol   = s%primal_tol
    me%dual_tol     = s%dual_tol
    me%zero_tol     = s%zero_tol
    me%pivot_tol    = s%pivot_tol
    me%progress_tol = s%progress_tol
    me%cycle_tol    = s%cycle_tol
    me%iter_limit   = s%iter_limit
    me%fval_bound   = s%fval_bound
    me%eps_prox     = s%eps_prox
    me%eta_prox     = s%eta_prox
    me%rho_soft     = s%rho_soft
    me%w_soft       = s%w_soft
    me%sing_tol     = s%sing_tol
    me%refactor_tol = s%refactor_tol

    end subroutine set_defaults
!*****************************************************************************************

!*****************************************************************************************
!>
!  Copy the options into the workspace.

    subroutine copy_settings(me)

    class(daqp_type), intent(inout) :: me

    me%work%settings%primal_tol   = me%primal_tol
    me%work%settings%dual_tol     = me%dual_tol
    me%work%settings%zero_tol     = me%zero_tol
    me%work%settings%pivot_tol    = me%pivot_tol
    me%work%settings%progress_tol = me%progress_tol
    me%work%settings%cycle_tol    = me%cycle_tol
    me%work%settings%iter_limit   = me%iter_limit
    me%work%settings%fval_bound   = me%fval_bound
    me%work%settings%eps_prox     = me%eps_prox
    me%work%settings%eta_prox     = me%eta_prox
    me%work%settings%rho_soft     = me%rho_soft
    me%work%settings%w_soft       = me%w_soft
    me%work%settings%sing_tol     = me%sing_tol
    me%work%settings%refactor_tol = me%refactor_tol

    end subroutine copy_settings
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up the solver for a QP: factor `H` and form the LDP. Call once per QP
!  (or after a change of the dimensions).
!
!  The number of simple bounds is `ms = size(bupper) - size(A,1)` (`size(bupper)`
!  without `A`). Omitted arguments must be followed by keyword arguments.

    subroutine setup(me, H, f, A, bupper, blower, istat, ms, sense)

    class(daqp_type), intent(inout) :: me
    real(wp), intent(in), optional, contiguous :: H(:,:)  !! Hessian `(n,n)`, symmetric (absent: an LP)
    real(wp), intent(in), optional, contiguous :: f(:)    !! linear term `(n)` (absent: zero)
    real(wp), intent(in), optional, contiguous :: A(:,:)  !! general constraints `(m-ms,n)` (absent: only simple bounds)
    real(wp), intent(in), contiguous :: bupper(:)         !! upper bounds `(m)` (`daqp_inf` for none)
    real(wp), intent(in), contiguous :: blower(:)         !! lower bounds `(m)` (`-daqp_inf` for none)
    integer(ip), intent(out) :: istat         !! `daqp_success` or a negative status
    integer(ip), intent(in), optional :: ms   !! number of simple bounds (checked if present)
    integer(ip), intent(in), optional, contiguous :: sense(:) !! constraint flags `(m)` (`daqp_equality`, `daqp_soft`, ...)

    integer(ip) :: n, m, nms, ma

    call me%destroy()
    istat = daqp_invalid_input

    ! dimensions
    m = size(bupper)
    if (present(H)) then
        n = size(H,1)
    else if (present(f)) then
        n = size(f)
    else if (present(A)) then
        n = size(A,2)
    else
        n = 0
    end if
    ma = 0
    if (present(A)) ma = size(A,1)
    nms = m - ma
    if (n < 1 .or. size(blower) /= m .or. nms < 0 .or. nms > n) then
        me%status = istat
        return
    end if
    if (present(ms)) then
        if (ms /= nms) then
            me%status = istat; return
        end if
    end if
    if (present(H)) then
        if (size(H,2) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(f)) then
        if (size(f) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(A)) then
        if (size(A,2) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(sense)) then
        if (size(sense) /= m) then
            me%status = istat; return
        end if
        if (any(iand(sense, 16_ip) /= 0)) then ! binary constraints are not supported
            me%status = istat; return
        end if
    end if

    call me%copy_settings()
    istat = daqp_setup(me%work, n, m, nms, bupper, blower, H, f, A, sense)
    if (istat > 0) istat = daqp_success
    me%status = istat

    end subroutine setup
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the QP. By default, the solve starts from the working set of the
!  previous solve (a hot start, which suits a sequence of similar QPs), or
!  from the empty working set (plus the equalities) after `setup`.
!
!  With `active` (and `at_lower`), it starts from that working set instead;
!  with `cold = .true.`, from the empty working set.

    subroutine solve(me, x, lam, istat, active, at_lower, cold)

    class(daqp_type), intent(inout) :: me
    real(wp), intent(out) :: x(:)       !! solution `(n)`
    real(wp), intent(out) :: lam(:)     !! multipliers `(m)` (see the module documentation for the sign)
    integer(ip), intent(out) :: istat   !! status (positive: solved)
    integer(ip), intent(in), optional :: active(:)  !! warm start: indices of the active constraints
    logical, intent(in), optional :: at_lower(:)    !! warm start: whether each is at its lower bound (default: upper)
    logical, intent(in), optional :: cold           !! start from the empty working set

    type(daqp_result) :: res
    logical :: soft_changed

    x = 0.0_wp
    lam = 0.0_wp
    if (.not. me%work%is_setup) then
        istat = daqp_not_setup
        me%status = istat
        return
    end if
    if (size(x) /= me%work%n .or. size(lam) /= me%work%m) then
        istat = daqp_invalid_input
        me%status = istat
        return
    end if

    soft_changed = me%work%settings%rho_soft /= me%rho_soft .or. &
                   me%work%settings%w_soft /= me%w_soft
    call me%copy_settings()
    if (soft_changed) call daqp_refresh_soft_weights(me%work)
    if (present(active)) then
        call me%set_working_set(active, istat, at_lower)
        if (istat < 0) return
    else if (present(cold)) then
        if (cold) then
            call daqp_deactivate_constraints(me%work)
            istat = daqp_activate_constraints(me%work)
            if (istat < 0) then
                me%status = istat
                return
            end if
        end if
    end if

    call daqp_solve(me%work, x, res, lam)
    istat          = res%exitflag
    me%status      = istat
    me%iter        = res%iter
    me%outer_iter  = res%nodes
    me%fval        = res%fval
    me%soft_slack  = res%soft_slack

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replace the data of the QP (same dimensions), and recompute only what
!  depends on it: a new `H` refactors; new `f` or bounds keep the
!  factorizations and the working set (for a hot start).

    subroutine update(me, istat, H, f, A, bupper, blower)

    class(daqp_type), intent(inout) :: me
    integer(ip), intent(out) :: istat          !! `daqp_success` or a negative status
    real(wp), intent(in), optional :: H(:,:)   !! new Hessian `(n,n)` (only if one was given to `setup`)
    real(wp), intent(in), optional :: f(:)     !! new linear term `(n)` (only if one was given to `setup`)
    real(wp), intent(in), optional :: A(:,:)   !! new general constraints `(m-ms,n)`
    real(wp), intent(in), optional :: bupper(:) !! new upper bounds `(m)`
    real(wp), intent(in), optional :: blower(:) !! new lower bounds `(m)`

    integer(ip) :: mask, n, m, ms

    if (.not. me%work%is_setup) then
        istat = daqp_not_setup
        me%status = istat
        return
    end if
    istat = daqp_invalid_input
    n = me%work%n
    m = me%work%m
    ms = me%work%ms
    mask = 0
    if (present(H)) then
        if (.not. me%work%has_H .or. size(H,1) /= n .or. size(H,2) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(f)) then
        if (.not. me%work%has_f .or. size(f) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(A)) then
        if (size(A,1) /= m-ms .or. size(A,2) /= n) then
            me%status = istat; return
        end if
    end if
    if (present(bupper)) then
        if (size(bupper) /= m) then
            me%status = istat; return
        end if
    end if
    if (present(blower)) then
        if (size(blower) /= m) then
            me%status = istat; return
        end if
    end if

    if (present(H)) then
        me%work%Hc = transpose(H)
        mask = ior(mask, daqp_update_rinv)
    end if
    if (present(f)) then
        me%work%f = f
        mask = ior(mask, daqp_update_v)
    end if
    if (present(A)) then
        me%work%At = transpose(A)
        mask = ior(mask, daqp_update_m)
    end if
    if (present(bupper)) then
        me%work%bupper = bupper
        mask = ior(mask, daqp_update_d)
    end if
    if (present(blower)) then
        me%work%blower = blower
        mask = ior(mask, daqp_update_d)
    end if

    call me%copy_settings()
    istat = daqp_update_ldp(me%work, mask)
    if (istat > 0) istat = daqp_success
    me%status = istat

    end subroutine update
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the working set for the next solve (a warm start). The indices must
!  be distinct and in `1..m`; equality constraints are always active.
!  Constraints that are linearly dependent on the others are dropped.

    subroutine set_working_set(me, active, istat, at_lower)

    class(daqp_type), intent(inout) :: me
    integer(ip), intent(in) :: active(:)          !! indices of the active constraints
    integer(ip), intent(out) :: istat             !! `daqp_success` or a negative status
    logical, intent(in), optional :: at_lower(:)  !! whether each is at its lower bound (default: upper)

    integer(ip) :: i, j, m
    logical, allocatable :: lower(:)

    if (.not. me%work%is_setup) then
        istat = daqp_not_setup
        me%status = istat
        return
    end if
    istat = daqp_invalid_input
    m = me%work%m
    if (present(at_lower)) then
        if (size(at_lower) /= size(active)) then
            me%status = istat; return
        end if
    end if
    do i = 1, size(active)
        if (active(i) < 1 .or. active(i) > m) then
            me%status = istat; return
        end if
        do j = 1, i-1
            if (active(j) == active(i)) then
                me%status = istat; return
            end if
        end do
    end do

    allocate(lower(size(active)), stat=i)
    if (i /= 0) then
        istat = daqp_out_of_memory
        me%status = istat
        return
    end if
    lower = .false.
    if (present(at_lower)) lower = at_lower
    istat = daqp_set_working_set(me%work, active, lower)
    if (istat > 0) istat = daqp_success
    me%status = istat

    end subroutine set_working_set
!*****************************************************************************************

!*****************************************************************************************
!>
!  The working set at the latest solution (for a later warm start).

    subroutine get_working_set(me, active, at_lower)

    class(daqp_type), intent(in) :: me
    integer(ip), allocatable, intent(out) :: active(:)  !! indices of the active constraints
    logical, allocatable, intent(out) :: at_lower(:)    !! whether each is at its lower bound

    integer(ip) :: i, na

    na = me%work%n_active
    if (.not. me%work%is_setup) na = 0
    allocate(active(na), at_lower(na))
    do i = 1, na
        active(i) = me%work%WS(i)
        at_lower(i) = iand(me%work%sense(active(i)), daqp_lower) /= 0
    end do

    end subroutine get_working_set
!*****************************************************************************************

!*****************************************************************************************
!>
!  Information on the latest solve.

    subroutine info(me, status, iter, outer_iter, fval, soft_slack, n_active)

    class(daqp_type), intent(in) :: me
    integer(ip), intent(out), optional :: status     !! status of the latest call
    integer(ip), intent(out), optional :: iter       !! iterations
    integer(ip), intent(out), optional :: outer_iter !! outer (proximal) iterations
    real(wp), intent(out), optional :: fval          !! objective function value
    real(wp), intent(out), optional :: soft_slack    !! largest violation of a soft constraint
    integer(ip), intent(out), optional :: n_active   !! number of active constraints

    if (present(status))     status     = me%status
    if (present(iter))       iter       = me%iter
    if (present(outer_iter)) outer_iter = me%outer_iter
    if (present(fval))       fval       = me%fval
    if (present(soft_slack)) soft_slack = me%soft_slack
    if (present(n_active))   n_active   = me%work%n_active

    end subroutine info
!*****************************************************************************************

!*****************************************************************************************
!>
!  Number of variables (0 if not set up).

    pure integer(ip) function n_variables(me)

    class(daqp_type), intent(in) :: me

    n_variables = me%work%n

    end function n_variables
!*****************************************************************************************

!*****************************************************************************************
!>
!  Number of constraints, including the simple bounds (0 if not set up).

    pure integer(ip) function n_constraints(me)

    class(daqp_type), intent(in) :: me

    n_constraints = me%work%m

    end function n_constraints
!*****************************************************************************************

!*****************************************************************************************
!>
!  Free the memory (the options are kept).

    subroutine destroy(me)

    class(daqp_type), intent(inout) :: me

    call daqp_destroy(me%work)
    me%status = daqp_not_setup
    me%iter = 0
    me%outer_iter = 0
    me%fval = 0.0_wp
    me%soft_slack = 0.0_wp

    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  Convert multipliers from DAQP's convention (\( \lambda \le 0 \) at a lower
!  bound, \( Hx + f + A^T\lambda = 0 \)) to the convention \( \lambda \ge 0 \)
!  at a lower bound (\( Hx + f = A^T\lambda \)).

    elemental real(wp) function daqp_lower_positive(lam)

    real(wp), intent(in) :: lam !! multiplier in DAQP's convention

    daqp_lower_positive = -lam

    end function daqp_lower_positive
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_module
!*****************************************************************************************
