!*****************************************************************************************
!> author: Jacob Williams
!
!  Constants, derived types, and small helpers of the Fortran port of
!  [DAQP](https://github.com/darnstrom/daqp) v0.10.3
!  (Copyright (c) 2022 Daniel Arnström, MIT licence): upstream's
!  `constants.h` and `types.h`, translated to Fortran (changed from the
!  original: 1-based indexing, allocatable components instead of pointers).
!
!### Storage
!
!  Matrices are kept in upstream's memory order: a C row-major matrix `X[r][c]`
!  is the Fortran array `X(c,r)`. So `Hc(j,i)` is `H(i,j)`, `At(:,k)` is the
!  k-th row of `A`, and `Mr(:,k)` the k-th general constraint row of the LDP.
!  Triangular factors are packed by rows (see [[ridx]] and [[lidx]]).
!
!  Constraint indices and working set positions are 1-based; `empty_ind` (0)
!  marks an empty index (upstream's -1).

    module daqp_types

    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    public
    private :: wp, ip, int64

    ! exit flags (upstream's values)
    integer(ip), parameter :: daqp_exit_optimal_inexact        = 4  !! optimal, found with the noise floor; violates a constraint by more than `primal_tol`
    integer(ip), parameter :: daqp_exit_no_freedom             = 3  !! (hierarchical QP) a level failed: the solution of the levels before
    integer(ip), parameter :: daqp_exit_soft_optimal           = 2  !! optimal, with a soft constraint violated
    integer(ip), parameter :: daqp_exit_optimal                = 1  !! optimal
    integer(ip), parameter :: daqp_exit_infeasible             = -1 !! primal infeasible
    integer(ip), parameter :: daqp_exit_cycle                  = -2 !! cycling detected
    integer(ip), parameter :: daqp_exit_unbounded              = -3 !! unbounded (LP)
    integer(ip), parameter :: daqp_exit_iterlimit              = -4 !! iteration limit reached
    integer(ip), parameter :: daqp_exit_nonconvex              = -5 !! the Hessian could not be factored
    integer(ip), parameter :: daqp_exit_overdetermined_initial = -6 !! inconsistent equalities in the initial working set
    integer(ip), parameter :: daqp_exit_timelimit              = -7 !! time limit reached
    integer(ip), parameter :: daqp_exit_unsupported            = -8 !! unsupported combination (e.g., hierarchy with a proximal term)
    integer(ip), parameter :: daqp_exit_invalid_input          = -101 !! invalid input (Fortran interface)
    integer(ip), parameter :: daqp_exit_out_of_memory          = -102 !! an allocation failed (Fortran interface)
    integer(ip), parameter :: daqp_exit_not_setup              = -103 !! the solver was not set up (Fortran interface)

    ! constraint flags (bits of `sense`)
    integer(ip), parameter :: daqp_active        = 1   !! the constraint is in the working set
    integer(ip), parameter :: daqp_lower         = 2   !! the active bound is the lower one
    integer(ip), parameter :: daqp_immutable     = 4   !! the constraint never leaves (or enters) the working set
    integer(ip), parameter :: daqp_soft          = 8   !! the constraint may be violated, at a penalty
    integer(ip), parameter :: daqp_binary        = 16  !! the constraint is active at one of its bounds (branch and bound)
    integer(ip), parameter :: daqp_slack_fixed   = 32  !! the slack of the soft constraint is zero
    integer(ip), parameter :: daqp_set_aside     = 64  !! temporarily set aside (proximal step)
    integer(ip), parameter :: daqp_auto_equality = 128 !! active and immutable set by the detection of equal bounds
    integer(ip), parameter :: daqp_equality      = daqp_active + daqp_immutable !! an equality constraint

    ! update masks
    integer(ip), parameter :: daqp_update_rinv          = 1
    integer(ip), parameter :: daqp_update_m             = 2
    integer(ip), parameter :: daqp_update_v             = 4
    integer(ip), parameter :: daqp_update_d             = 8
    integer(ip), parameter :: daqp_update_sense         = 16
    integer(ip), parameter :: daqp_update_hierarchy     = 32
    integer(ip), parameter :: daqp_update_unconstrained = 64
    integer(ip), parameter :: daqp_update_eliminate     = 128 !! lets `eq_reduction = auto` eliminate the equalities

    ! workspace state masks
    integer(ip), parameter :: state_pending = daqp_update_rinv + daqp_update_m + daqp_update_v + &
                                              daqp_update_d + daqp_update_sense + daqp_update_hierarchy
    integer(ip), parameter :: state_unconstrained     = 256
    integer(ip), parameter :: state_rinv_normalized   = 512
    integer(ip), parameter :: state_incumbent         = 1024
    integer(ip), parameter :: state_ill_conditioned   = 2048
    integer(ip), parameter :: state_cholesky_pending  = 4096
    integer(ip), parameter :: state_noise_floor       = 8192

    ! equality-reduction policy (`settings%eq_reduction`)
    integer(ip), parameter :: daqp_eq_reduction_off  = -1 !! never eliminate
    integer(ip), parameter :: daqp_eq_reduction_auto = 0  !! eliminate for a solve from scratch, if worthwhile
    integer(ip), parameter :: daqp_eq_reduction_on   = 1  !! always eliminate (also for warm-started workspaces)

    ! problem types
    integer(ip), parameter :: daqp_problem_qp       = 0 !! a QP (or LP)
    integer(ip), parameter :: daqp_problem_avi      = 1 !! an affine variational inequality (nonsymmetric `H`)
    integer(ip), parameter :: daqp_problem_factored = 2 !! a QP whose `H` is given by its Cholesky factor

    integer(ip), parameter :: empty_ind = 0                  !! an empty index (upstream's -1)
    integer(ip), parameter :: unconstrained_optimal = -2     !! return value of `check_unconstrained`
    integer(ip), parameter :: eq_not_reduced = 1             !! return value of `daqp_eq_update`

    real(wp), parameter :: daqp_inf = 1.0e30_wp      !! "infinite" bound

    ! default settings (upstream's; the tolerances are floored at a multiple of
    ! epsilon, which only changes them in single precision)
    real(wp), parameter :: daqp_default_prim_tol     = max(1.0e-6_wp, 100.0_wp*epsilon(1.0_wp))
    real(wp), parameter :: daqp_default_dual_tol     = max(1.0e-12_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter :: daqp_default_zero_tol     = max(1.0e-11_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter :: daqp_default_prog_tol     = max(1.0e-14_wp, 10.0_wp*epsilon(1.0_wp))
    real(wp), parameter :: daqp_default_pivot_tol    = 1.0e-8_wp
    integer(ip), parameter :: daqp_default_cycle_tol = 10
    real(wp), parameter :: daqp_default_eta          = -1.0_wp
    integer(ip), parameter :: daqp_default_iter_limit = 10000
    real(wp), parameter :: daqp_default_rho_soft     = max(1.0e-6_wp, 1.0e4_wp*epsilon(1.0_wp)) ! (above sing_tol)
    real(wp), parameter :: daqp_default_w_soft       = 0.0_wp
    real(wp), parameter :: daqp_default_sing_tol     = max(3.7e-11_wp, 1000.0_wp*epsilon(1.0_wp))
    real(wp), parameter :: daqp_default_refactor_tol = 1.0e-9_wp
    real(wp), parameter :: daqp_default_eps_prox     = -max(1.0e-6_wp, 100.0_wp*epsilon(1.0_wp))

    ! internal constants
    real(wp), parameter :: auto_eta_cap        = max(1.0e-6_wp, 1000.0_wp*epsilon(1.0_wp)) ! (upstream's 1e-6, floored in single precision)
    integer(ip), parameter :: refine_min_iter  = 5
    real(wp), parameter :: refine_cond         = 1.0e6_wp
    real(wp), parameter :: refine_gain         = 1.0e3_wp
    real(wp), parameter :: refine_pivot        = 1.0e-6_wp
    real(wp), parameter :: add_noise_gain      = 10.0_wp
    real(wp), parameter :: prox_eps_max        = 1.0e-3_wp
    real(wp), parameter :: hessian_cond_max    = 1.0e8_wp
    real(wp), parameter :: hessian_cond_eps    = 0.1_wp
    real(wp), parameter :: cond_defer_margin   = 100.0_wp
    integer(ip), parameter :: prox_face_start  = 16
    integer(ip), parameter :: prox_face_steps  = 3
    real(wp), parameter :: avi_pivot_trigger   = 0.05_wp
    real(wp), parameter :: avi_retry_rho_reduction = 16.0_wp
    integer(ip), parameter :: eq_path_ldp = 0  !! reduced problem: identity Hessian, no linear term
    integer(ip), parameter :: eq_path_qp  = 1  !! reduced problem: reduced Hessian W'HW
    integer(ip), parameter :: eq_path_lp  = 2  !! reduced problem: no Hessian
    integer(ip), parameter :: eq_min_count = 5
    integer(ip), parameter :: eq_min_ratio = 10
    integer(ip), parameter :: eq_min_dim = 12
    integer(ip), parameter :: eq_diag_min_ratio = 4
    integer(ip), parameter :: bnb_lower_bit = 16 !! the bit that marks a lower bound in a BnB constraint id

    ! what `R` holds
    integer(ip), parameter :: rinv_none  = 0 !! no Hessian (LP): `R = I`
    integer(ip), parameter :: rinv_dense = 1 !! dense packed `R`
    integer(ip), parameter :: rinv_diag  = 2 !! diagonal Hessian: `R(1:n)` holds the diagonal of `R^{-1}`

    type :: daqp_settings
        !! Solver settings (upstream's `DAQPSettings`, with its defaults).
        real(wp)    :: primal_tol   = daqp_default_prim_tol     !! tolerance for primal feasibility
        real(wp)    :: dual_tol     = daqp_default_dual_tol     !! tolerance for dual feasibility
        real(wp)    :: zero_tol     = daqp_default_zero_tol     !! values below are regarded as zero
        real(wp)    :: pivot_tol    = daqp_default_pivot_tol    !! pivots of the `LDL'` below are reordered
        real(wp)    :: progress_tol = daqp_default_prog_tol     !! minimum objective progress (cycle guard)
        integer(ip) :: cycle_tol    = daqp_default_cycle_tol    !! iterations without progress before cycling is assumed
        integer(ip) :: iter_limit   = daqp_default_iter_limit   !! maximum number of iterations
        real(wp)    :: fval_bound   = daqp_inf                  !! upper bound of the objective (infeasible above)
        real(wp)    :: eps_prox     = daqp_default_eps_prox     !! proximal regularization (negative: automatic, only if needed)
        real(wp)    :: eta_prox     = daqp_default_eta          !! tolerance of the proximal outer loop (negative: automatic)
        real(wp)    :: rho_soft     = daqp_default_rho_soft     !! reciprocal quadratic weight of the soft constraints
        real(wp)    :: rel_subopt   = 0.0_wp                    !! relative suboptimality tolerance (branch and bound)
        real(wp)    :: abs_subopt   = 0.0_wp                    !! absolute suboptimality tolerance (branch and bound)
        real(wp)    :: sing_tol     = daqp_default_sing_tol     !! pivots of the `LDL'` below mark a singular working set
        real(wp)    :: refactor_tol = daqp_default_refactor_tol !! pivots below trigger a refactorization at a solution
        real(wp)    :: time_limit   = 0.0_wp                    !! time limit of a solve in seconds (0: none)
        real(wp)    :: w_soft       = daqp_default_w_soft       !! linear weight of the soft constraints
        integer(ip) :: eq_reduction = daqp_eq_reduction_auto    !! equality-reduction policy (`daqp_eq_reduction_*`)
    end type daqp_settings

    type :: daqp_result
        !! Result of a solve (upstream's `DAQPResult`).
        real(wp)    :: fval       = 0.0_wp  !! objective function value
        real(wp)    :: soft_slack = 0.0_wp  !! largest violation of a soft constraint
        integer(ip) :: exitflag   = daqp_exit_not_setup !! exit flag
        integer(ip) :: iter       = 0       !! number of iterations
        integer(ip) :: nodes      = 0       !! nodes (branch and bound), or outer iterations (proximal, AVI)
        real(wp)    :: solve_time = 0.0_wp  !! time of the solve [s]
        real(wp)    :: setup_time = 0.0_wp  !! time of the setup [s]
    end type daqp_result

    type :: daqp_problem
        !! A problem (upstream's `DAQPProblem`), in upstream's memory order.
        integer(ip) :: n  = 0  !! number of variables
        integer(ip) :: m  = 0  !! number of constraints (including the simple bounds)
        integer(ip) :: ms = 0  !! number of simple bounds
        logical :: has_H = .false.      !! a Hessian is given (else an LP)
        logical :: has_f = .false.      !! a linear term is given
        logical :: has_sense = .false.  !! constraint flags are given
        integer(ip) :: problem_type = daqp_problem_qp !! `daqp_problem_*`
        real(wp), allocatable :: Hc(:,:)      !! `Hc(j,i) = H(i,j)` (not for `daqp_problem_factored`)
        real(wp), allocatable :: Rf(:)        !! the Cholesky factor of `H`, packed by rows (`daqp_problem_factored`)
        real(wp), allocatable :: f(:)         !! linear term
        real(wp), allocatable :: At(:,:)      !! `At(:,k)` is the k-th row of `A`
        real(wp), allocatable :: bupper(:)    !! upper bounds
        real(wp), allocatable :: blower(:)    !! lower bounds
        integer(ip), allocatable :: sense(:)  !! constraint flags
        integer(ip) :: nh = 1                 !! number of levels (hierarchical QP)
        integer(ip), allocatable :: break_points(:) !! the last constraint of each level (hierarchical QP)
    end type daqp_problem

    type :: daqp_ldp_data
        !! The parts of the workspace that describe the LDP that is solved
        !! (upstream's `DAQPLDPData`): an equality elimination swaps those of its
        !! reduced problem into the workspace.
        type(daqp_problem) :: qp
        integer(ip) :: n = 0, m = 0, ms = 0
        real(wp), allocatable :: Mr(:,:), dupper(:), dlower(:), R(:), v(:), scaling(:), Mu(:)
        integer(ip) :: rmode = rinv_none
        logical :: has_v = .false.
        integer(ip), allocatable :: sense(:)
        logical :: has_weights = .false.
        real(wp), allocatable :: rho_ls(:), rho_us(:), w_ls(:), w_us(:)
        integer(ip) :: state = 0
        integer(ip) :: n_prox = 0
        integer(ip), allocatable :: bin_ids(:)
        integer(ip) :: nb = 0
    end type daqp_ldp_data

    type :: daqp_node
        !! A node of the branch-and-bound tree (upstream's `DAQPNode`).
        integer(ip) :: bin_id = 0   !! binary constraint that is fixed (with the lower-bound bit)
        integer(ip) :: depth = 0    !! depth in the tree (-1: the root)
        integer(ip) :: ws_start = 0 !! start of the node's warm start in `tree_ws` (0-based)
        integer(ip) :: ws_end = 0   !! end of the node's warm start in `tree_ws` (0-based, exclusive)
    end type daqp_node

    type :: daqp_bnb_data
        !! Branch and bound (upstream's `DAQPBnB`).
        integer(ip), allocatable :: bin_ids(:)    !! the binary constraints
        integer(ip) :: nb = 0                     !! number of binary constraints
        integer(ip) :: neq = 0                    !! length of the fixed (immutable) prefix of the working set
        type(daqp_node), allocatable :: tree(:)   !! the stack of nodes
        integer(ip) :: n_nodes = 0
        integer(ip), allocatable :: tree_ws(:)    !! warm starts of the nodes
        integer(ip) :: nws = 0
        integer(ip) :: n_clean = 0
        integer(ip), allocatable :: fixed_ids(:)
        integer(ip) :: nodecount = 0
        integer(ip) :: itercount = 0
        integer(ip), allocatable :: root_ws(:)    !! working set of the latest root relaxation
        integer(ip) :: n_root_ws = 0
    end type daqp_bnb_data

    type :: daqp_avi_data
        !! Affine variational inequality (upstream's `DAQPAVI`). The `n x n`
        !! matrices are in upstream's (row-major) memory order.
        logical :: is_symmetric = .false.
        logical :: retry_rho_needed = .false.
        real(wp), allocatable :: Hsym(:,:), Hs_rho(:,:), H_rho(:,:), LU_H(:,:)
        integer(ip), allocatable :: P_H2(:), P_H(:), P_S(:)
        real(wp), allocatable :: kkt_buffer(:)
        real(wp), allocatable :: xtemp(:), Hx(:), x(:), y(:)
        real(wp) :: rho = 0.0_wp
    end type daqp_avi_data

    type :: daqp_vec
        real(wp), allocatable :: v(:)
    end type daqp_vec

    type :: daqp_eq_data
        !! Elimination of equality constraints (upstream's `DAQPEqElim`).
        integer(ip) :: n = 0, m = 0, ms = 0     !! dimensions of the original problem
        integer(ip) :: neq = 0                  !! number of eliminated equalities
        integer(ip) :: nz = 0                   !! number of variables of the reduced problem
        integer(ip) :: mr = 0                   !! number of constraints of the reduced problem
        integer(ip) :: ndrop = 0                !! number of constraints that the equalities imply
        integer(ip) :: ncand = 0                !! number of equality candidates
        integer(ip) :: path = eq_path_ldp
        logical :: metric = .false.             !! diagonal Hessian (QR in the metric of H)
        logical :: active = .false.             !! a reduction is formed
        integer(ip) :: error = 0                !! exit flag if the latest right-hand side is infeasible
        logical :: installed = .false.          !! the reduced problem is in the workspace
        logical :: allocated_dims = .false.
        integer(ip), allocatable :: eq_ids(:), cand_ids(:), keep(:), drop_ids(:)
        real(wp), allocatable :: V(:,:)         !! Householder vectors (columns 1:neq), then Z
        real(wp), allocatable :: tau(:), s_eq(:)
        real(wp), allocatable :: R(:)           !! triangular factor, packed (row k: `R(arsum(k)+1:...)`)
        real(wp), allocatable :: dsq(:)         !! H^{-1/2} for a diagonal Hessian
        real(wp), allocatable :: W(:,:)         !! null-space basis: `W(:,i)` is row i (`nz x n`)
        real(wp), allocatable :: xp(:)          !! particular solution
        real(wp) :: fp = 0.0_wp                 !! objective at xp
        real(wp), allocatable :: tmp(:)         !! scratch (3n)
        type(daqp_vec), allocatable :: cols(:)  !! responses to the equalities' right-hand sides
        integer(ip) :: ncols = 0
        real(wp), allocatable :: xf(:), gf(:), df(:), sh(:)
        logical :: f_valid = .false.
        real(wp), allocatable :: rho_r(:)       !! soft weights of the reduced constraints (4 mr)
        type(daqp_ldp_data) :: other            !! the LDP (with its problem) that is not in the workspace
    end type daqp_eq_data

    type :: daqp_workspace
        !! The workspace of the solver (upstream's `DAQPWorkspace`).
        !! All arrays are allocated in `daqp_setup`, so a solve allocates nothing
        !! (except the lazily formed caches of the equality elimination).
        type(daqp_settings) :: settings !! the settings
        type(daqp_problem) :: qp        !! the problem (the reduced one while an elimination is installed)
        logical :: has_qp = .true.      !! (false for `daqp_minrep`)
        ! the LDP
        integer(ip) :: n  = 0  !! number of variables
        integer(ip) :: m  = 0  !! number of constraints (including the simple bounds)
        integer(ip) :: ms = 0  !! number of simple bounds
        integer(ip) :: rmode = rinv_none      !! what `R` holds
        real(wp), allocatable :: R(:)         !! packed upper triangular `R^{-1}`
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
        real(wp), allocatable :: xold(:)      !! the previous primal iterate
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
        ! extensions
        type(daqp_bnb_data), allocatable :: bnb  !! branch and bound (allocated if there are binary constraints)
        integer(ip) :: nh = 1                     !! levels (hierarchical QP), or a counter of outer iterations
        logical :: has_bp = .false.               !! break points are given (hierarchical QP)
        integer(ip), allocatable :: break_points(:)
        type(daqp_avi_data), allocatable :: avi   !! AVI (allocated for `daqp_problem_avi`)
        type(daqp_eq_data), allocatable :: eq     !! equality elimination (allocated while used)
        logical :: timer_on = .false.             !! the time limit is checked
        integer(int64) :: timer_start = 0         !! start of the solve (`system_clock` count)
        ! individual weights of the soft constraints (zero selects the settings)
        logical :: has_weights = .false.
        real(wp), allocatable :: rho_ls(:), rho_us(:), w_ls(:), w_us(:)
        integer(ip) :: state = 0              !! state mask
        logical :: is_setup = .false.         !! set up successfully
    end type daqp_workspace

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
!  Whether the workspace holds a hierarchical QP (upstream's `DAQP_IS_HIERARCHICAL`).

    pure logical function is_hierarchical(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    is_hierarchical = work%has_bp .and. work%nh > 1

    end function is_hierarchical
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the equality constraints of the workspace are eliminated
!  (upstream's `DAQP_IS_REDUCED`).

    pure logical function is_reduced(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    is_reduced = .false.
    if (allocated(work%eq)) is_reduced = work%eq%active

    end function is_reduced
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the workspace holds a nonsymmetric AVI.

    pure logical function is_avi_nonsym(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    is_avi_nonsym = .false.
    if (allocated(work%avi)) is_avi_nonsym = .not. work%avi%is_symmetric

    end function is_avi_nonsym
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
    allocate(a(max(0_ip,n)), stat=stat)
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
    allocate(a(max(0_ip,n1),max(0_ip,n2)), stat=stat)
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
    allocate(a(max(0_ip,n)), stat=stat)
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
    allocate(a(max(0_ip,n)), stat=stat)
    if (stat /= 0) istat = stat

    end subroutine resize1l
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two real vectors (no copy).

    subroutine swap_r1(a,b)

    real(wp), allocatable, intent(inout) :: a(:) !! first
    real(wp), allocatable, intent(inout) :: b(:) !! second

    real(wp), allocatable :: t(:)

    call move_alloc(a, t)
    call move_alloc(b, a)
    call move_alloc(t, b)

    end subroutine swap_r1
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two real matrices (no copy).

    subroutine swap_r2(a,b)

    real(wp), allocatable, intent(inout) :: a(:,:) !! first
    real(wp), allocatable, intent(inout) :: b(:,:) !! second

    real(wp), allocatable :: t(:,:)

    call move_alloc(a, t)
    call move_alloc(b, a)
    call move_alloc(t, b)

    end subroutine swap_r2
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two integer vectors (no copy).

    subroutine swap_i1(a,b)

    integer(ip), allocatable, intent(inout) :: a(:) !! first
    integer(ip), allocatable, intent(inout) :: b(:) !! second

    integer(ip), allocatable :: t(:)

    call move_alloc(a, t)
    call move_alloc(b, a)
    call move_alloc(t, b)

    end subroutine swap_i1
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two integers.

    pure subroutine swap_int(a,b)

    integer(ip), intent(inout) :: a !! first
    integer(ip), intent(inout) :: b !! second

    integer(ip) :: t

    t = a; a = b; b = t

    end subroutine swap_int
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two logicals.

    pure subroutine swap_log(a,b)

    logical, intent(inout) :: a !! first
    logical, intent(inout) :: b !! second

    logical :: t

    t = a; a = b; b = t

    end subroutine swap_log
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange two problems (no copy of the arrays).

    subroutine swap_problem(a,b)

    type(daqp_problem), intent(inout) :: a !! first
    type(daqp_problem), intent(inout) :: b !! second

    call swap_int(a%n, b%n)
    call swap_int(a%m, b%m)
    call swap_int(a%ms, b%ms)
    call swap_log(a%has_H, b%has_H)
    call swap_log(a%has_f, b%has_f)
    call swap_log(a%has_sense, b%has_sense)
    call swap_int(a%problem_type, b%problem_type)
    call swap_r2(a%Hc, b%Hc)
    call swap_r1(a%Rf, b%Rf)
    call swap_r1(a%f, b%f)
    call swap_r2(a%At, b%At)
    call swap_r1(a%bupper, b%bupper)
    call swap_r1(a%blower, b%blower)
    call swap_i1(a%sense, b%sense)
    call swap_int(a%nh, b%nh)
    call swap_i1(a%break_points, b%break_points)

    end subroutine swap_problem
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_types
!*****************************************************************************************
