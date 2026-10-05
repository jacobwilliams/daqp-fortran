!*****************************************************************************************
!> author: Jacob Williams
!
!  The Fortran port of the DAQP solver (dual active-set method for dense
!  convex quadratic programs, with branch and bound for binary constraints,
!  hierarchical QPs, affine variational inequalities, and a proximal-point
!  outer loop for semidefinite problems and LPs).
!
!  This module is a modern Fortran translation of the C code of
!  [DAQP](https://github.com/darnstrom/daqp) v0.10.3
!  (Copyright (c) 2022 Daniel Arnström, MIT licence): `daqp.c`,
!  `auxiliary.c`, `factorization.c`, `utils.c`, `api.c`, `daqp_prox.c`,
!  `bnb.c`, `hierarchical.c`, and `avi.c` (the elimination of equalities,
!  `eq_elim.c`, is in [[daqp_eq_elim]]). Changed from the original: translated
!  to Fortran, 1-based indexing, allocatable components instead of pointers
!  (see [[daqp_types]] for the storage).
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

    module daqp_core

    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use daqp_types, ridx_t => ridx, lidx_t => lidx, has_t => has ! (local copies below, which can be inlined)
    use daqp_eq_elim
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    private

    public :: daqp_setup, daqp_solve, daqp_update_ldp, daqp_destroy, daqp_quadprog, daqp_avi
    public :: daqp_ldp, daqp_ldp2qp_solution
    public :: daqp_activate_constraints, daqp_deactivate_constraints
    public :: daqp_set_working_set, daqp_refresh_soft_weights
    public :: daqp_allocate_soft_weights, daqp_set_soft_weights
    public :: daqp_set_primal_start, daqp_primal_init_active, daqp_dual_init_active
    public :: daqp_minrep, daqp_first_violating, daqp_extract_active_duals
    public :: daqp_bnb, daqp_hiqp, daqp_solve_avi, daqp_prox

    contains

!*****************************************************************************************
!>
!  Position of element `(i,j)`, `i<=j`, of an `n x n` upper triangle packed by
!  rows (a copy of [[daqp_types:ridx]] in this module, so that it is inlined
!  in the hot loops).

    pure integer(ip) function ridx(i,j,n)

    integer(ip), intent(in) :: i !! row
    integer(ip), intent(in) :: j !! column
    integer(ip), intent(in) :: n !! dimension

    ridx = ((i-1)*(2*n-i))/2 + j

    end function ridx
!*****************************************************************************************

!*****************************************************************************************
!>
!  Position of element `(i,j)`, `j<=i`, of a lower triangle packed by rows
!  (a copy of [[daqp_types:lidx]], so that it is inlined).

    pure integer(ip) function lidx(i,j)

    integer(ip), intent(in) :: i !! row
    integer(ip), intent(in) :: j !! column

    lidx = ((i-1)*i)/2 + j

    end function lidx
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `flag` is set in `s` (a copy of [[daqp_types:has]], so that it is inlined).

    pure logical function has(s,flag)

    integer(ip), intent(in) :: s    !! flags
    integer(ip), intent(in) :: flag !! flag to test

    has = iand(s,flag) /= 0

    end function has
!*****************************************************************************************
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
        if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) then
          ns_active = ns_active + 1
        end if
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
!  Reciprocal quadratic weight of the active side of soft constraint `id`
!  (a zero individual weight selects `settings%rho_soft`, which is given in
!  the normalized formulation).

    pure real(wp) function soft_rho(work,id)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id            !! constraint

    real(wp) :: rho

    if (work%has_weights) then
        if (has(work%sense(id),daqp_lower)) then
            rho = work%rho_ls(id)
        else
            rho = work%rho_us(id)
        end if
        if (rho /= 0.0_wp) then
            soft_rho = rho*work%scaling(id)*work%scaling(id)
            return
        end if
    end if
    soft_rho = work%settings%rho_soft

    end function soft_rho
!*****************************************************************************************

!*****************************************************************************************
!>
!  Linear weight of the active side of soft constraint `id`, i.e. what the
!  multiplier has to exceed for the slack to become nonzero (a zero individual
!  weight selects `settings%w_soft`).

    pure real(wp) function soft_w(work,id)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id            !! constraint

    real(wp) :: w

    if (work%has_weights) then
        if (has(work%sense(id),daqp_lower)) then
            w = work%w_ls(id)
        else
            w = work%w_us(id)
        end if
        if (w /= 0.0_wp) then
            soft_w = w/work%scaling(id)
            return
        end if
    end if
    soft_w = work%settings%w_soft

    end function soft_w
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the soft constraints may have a linear penalty, or individual weights
!  (else every soft constraint has the same, purely quadratic penalty).

    pure logical function has_l1(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    has_l1 = work%settings%w_soft /= 0.0_wp .or. work%has_weights

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

    w = soft_w(work,id)
    if (w == 0.0_wp) then
        if (lam == 0.0_wp) then
            soft_residual = 0.0_wp
        else
            soft_residual = soft_rho(work,id)*lam
        end if
    else
        if (has(work%sense(id),daqp_lower)) then
            soft_residual = soft_rho(work,id)*(lam + w)
        else
            soft_residual = soft_rho(work,id)*(lam - w)
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

    w = soft_w(work,id)
    if (w == 0.0_wp) then
        soft_penalty = soft_rho(work,id)*lam*lam
    else if (has(work%sense(id),daqp_slack_fixed)) then
        soft_penalty = 0.0_wp
    else
        soft_penalty = soft_rho(work,id)*(lam*lam-w*w)
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
!  refactorization). Returns false if the cycling is to be reported (for the
!  proximal method and hierarchical problems, which handle cycling themselves,
!  for a nonsymmetric AVI, or if the floor is already used).

    logical function set_noise_floor(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (work%n_prox > 0 .or. is_hierarchical(work) .or. is_avi_nonsym(work) .or. &
        has(work%state,state_noise_floor)) then
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
    if (has(work%sense(add_ind),daqp_soft) .and. .not. has(work%sense(add_ind),daqp_slack_fixed)) then
      rho = soft_rho(work,add_ind)
    end if
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
        w = soft_w(work,add_ind)
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
    work%x(1:n) = 0.0_wp
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
            w = soft_w(work,ind)
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
        rho = soft_rho(work,ind)
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
            if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) then
              s = s - soft_residual(work, id, 0.0_wp)
            end if
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
        w = soft_w(work,id)
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
        if (allocated(work%bnb)) then
            if (rm_ind-1 < work%bnb%n_clean) return
        end if
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
    if (has(work%sense(id),daqp_soft)) w = soft_w(work,id)
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
                dependency_residual < -work%settings%primal_tol*dependency_scale) then
              exitflag = daqp_exit_overdetermined_initial
            end if
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
                if (.not. has(work%sense(j),daqp_immutable)) then
                  work%sense(j) = iand(work%sense(j), not(daqp_active))
                end if
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

    if (allocated(work%bnb)) work%bnb%n_root_ws = 0 ! also drop the BnB warm start
    do i = 1, work%n_active
        id = work%WS(i)
        if (has(work%sense(id),daqp_immutable)) cycle
        work%sense(id) = iand(work%sense(id), not(daqp_active))
    end do
    call daqp_reset_workspace(work) ! the next update activates the remaining ones

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
        if (has(work%sense(id),daqp_soft) .and. .not. has(work%sense(id),daqp_slack_fixed)) then
          work%xldl(i) = work%xldl(i) - soft_residual(work, id, work%lam_star(i))
        end if
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
        val = dot_seq(work%n, work%qp%At(:,id-work%ms), work%x)
    end if
    if (has(work%sense(id),daqp_lower)) then
        active_residual = val - work%qp%blower(id)
    else
        active_residual = val - work%qp%bupper(id)
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
    if (na == 0 .or. .not. work%has_qp .or. work%sing_ind /= empty_ind) return

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
!  Time since the start of the solve [s].

    real(wp) function elapsed_time(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(int64) :: count, rate

    call system_clock(count, rate)
    elapsed_time = real(count-work%timer_start,wp)/real(rate,wp)

    end function elapsed_time
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
                    ! (not in BnB, which relies on the order of the working set)
                    if (work%n_active > 2 .and. .not. tried_repair .and. .not. allocated(work%bnb) .and. &
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
                        if (tried_repair .or. allocated(work%bnb)) then
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
        if (work%timer_on .and. mod(iter,32_ip) == 0) then
            if (elapsed_time(work) > work%settings%time_limit) then
                exitflag = daqp_exit_timelimit
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

    pure real(wp) function hessian_scale_of(n,H) result(hs)

    integer(ip), intent(in) :: n     !! dimension
    real(wp), intent(in) :: H(n,n)   !! Hessian

    integer(ip) :: i

    hs = 0.0_wp
    do i = 1, n
        if (abs(H(i,i)) > hs) hs = abs(H(i,i))
    end do

    end function hessian_scale_of
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form the packed Cholesky factor `R` of `H` (with reciprocal diagonal), or
!  copy a given factor `Rf`. Returns 1, 0 if an unfactored Hessian needs a
!  regularization, or `daqp_exit_nonconvex`.

    integer(ip) function form_R(work,regularize_all,eps,min_pivot,max_pivot,H,Rf) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: regularize_all    !! regularize all the directions
    real(wp), intent(in) :: eps              !! the regularization
    real(wp), intent(out) :: min_pivot       !! smallest pivot (of the unregularized directions)
    real(wp), intent(out) :: max_pivot       !! largest pivot
    real(wp), intent(in), optional :: H(work%n,work%n) !! Hessian (`Hc` convention)
    real(wp), intent(in), optional :: Rf(:)  !! Cholesky factor, packed by rows

    integer(ip) :: i, j, k, n, di, kik, kij
    real(wp) :: pivot, inv_diag, s, s0, s1, s2, s3, c

    n = work%n
    min_pivot = daqp_inf
    max_pivot = 0.0_wp
    flag = 1

    if (present(Rf)) then ! the factor is given
        do i = 1, n
            di = ridx(i,i,n)
            if (Rf(di) <= work%settings%zero_tol) then
                flag = daqp_exit_nonconvex
                return
            end if
            work%R(di) = 1.0_wp/Rf(di)
            do j = 1, n-i
                work%R(di+j) = Rf(di+j)
            end do
        end do
        return
    end if

    ! pack (symmetrized) H
    k = 1
    do i = 1, n
        work%R(k) = H(i,i)
        if (regularize_all) then
            work%R(k) = work%R(k) + eps
        else if (work%n_prox > 0) then
            if (work%prox_mask(i)) work%R(k) = work%R(k) + eps
        end if
        k = k + 1
        do j = i+1, n
            work%R(k) = 0.5_wp*(H(j,i) + H(i,j))
            k = k + 1
        end do
    end do

    flag = 0
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
    flag = 1

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

    logical function finish_Rinv(work,is_factored,H)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: is_factored          !! the factor was given (then `H` is not needed)
    real(wp), intent(in), optional :: H(work%n,work%n) !! Hessian (`Hc` convention)

    integer(ip) :: i, j, k, n
    real(wp) :: hinv_max, hmax, hii, s2, cnd
    logical :: installed

    n = work%n
    work%state = iand(work%state, not(state_cholesky_pending))
    call invert_R(work)
    ! cond(H) >= max (H^-1)_ii * max H_ii (H_ii >= R_ii^2 if H is factored)
    hinv_max = 0.0_wp
    hmax = 0.0_wp
    k = 1
    do i = 1, n
        if (is_factored) then
            hii = 1.0_wp/(work%R(k)*work%R(k))
        else
            hii = H(i,i)
        end if
        s2 = 0.0_wp
        do j = i, n
            s2 = s2 + work%R(k)*work%R(k)
            k = k + 1
        end do
        if (s2 > hinv_max) hinv_max = s2
        if (hii > hmax) hmax = hii
    end do
    ! regularize an ill-conditioned Hessian, or mark it for refinement
    installed = .false.
    if (allocated(work%eq)) installed = work%eq%installed
    cnd = hinv_max*hmax
    if (.not. is_factored .and. work%n_prox == 0 .and. .not. allocated(work%avi)) then
        if ((installed .and. cnd > hessian_cond_max) .or. &
            real(n,wp)*epsilon(1.0_wp)*hinv_max*hmax > hessian_cond_eps) then
            finish_Rinv = .false.
            return
        end if
    end if
    if (cnd > refine_cond) work%state = ior(work%state, state_ill_conditioned)
    finish_Rinv = .true.

    end function finish_Rinv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Factor `H` (upstream's `daqp_update_R`): diagonal, dense, or regularized
!  (proximal) when it is singular. With `defer_inverse`, the inverse of a
!  well-conditioned factor is deferred until a constrained solve needs it.
!  `H` absent and `Rf` absent: an LP; `Rf`: the Cholesky factor is given.

    integer(ip) function update_R(work,defer_inverse,H,Rf) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: defer_inverse !! defer the inverse of the factor
    real(wp), intent(in), optional :: H(work%n,work%n) !! Hessian (`Hc` convention)
    real(wp), intent(in), optional :: Rf(:) !! Cholesky factor of H, packed by rows

    integer(ip) :: i, j, n, regularization_tries, di
    real(wp) :: eps, zero_tol, factor_tol, hessian_scale, acceptance_tol, &
                dmin, dmax, hi, min_pivot, max_pivot, d, cond, eps_mach, ptol, abs_diag
    logical :: regularize_all, force_prox, is_diagonal, ok, zero_row, is_factored, installed

    n = work%n
    is_factored = present(Rf)
    eps = work%settings%eps_prox
    zero_tol = work%settings%zero_tol
    factor_tol = zero_tol
    hessian_scale = 0.0_wp
    regularization_tries = 0
    eps_mach = epsilon(1.0_wp)
    flag = 1

    force_prox = work%settings%eps_prox > 0.0_wp .and. .not. is_factored .and. &
                 .not. is_avi_nonsym(work)
    regularize_all = force_prox

    ! reset the semi-proximal mask for this factorization
    work%prox_mask(1:n) = .false.
    work%n_prox = 0
    work%state = iand(work%state, not(state_rinv_normalized + state_ill_conditioned + &
                                      state_cholesky_pending))

    if (.not. present(H) .and. .not. is_factored) then ! LP: all directions need proximal regularization
        if (work%has_qp) then
            if (work%qp%has_f) work%n_prox = n
        end if
        work%scaling(1:work%ms) = 1.0_wp
        return
    end if

    ! check if diagonal
    is_diagonal = .true.
    if (.not. is_factored) then
        do i = 1, n
            abs_diag = H(i,i)
            if (abs_diag < 0.0_wp) abs_diag = -abs_diag
            if (abs_diag > hessian_scale) hessian_scale = abs_diag
            do j = i+1, n
                if (H(j,i) > zero_tol .or. H(j,i) < -zero_tol) then
                    is_diagonal = .false.
                    exit
                end if
            end do
            if (.not. is_diagonal) exit
        end do
    else
        do i = 1, n
            di = ridx(i,i,n)
            do j = 1, n-i
                if (Rf(di+j) > zero_tol .or. Rf(di+j) < -zero_tol) then
                    is_diagonal = .false.
                    exit
                end if
            end do
            if (.not. is_diagonal) exit
        end do
    end if

    if (force_prox) then
        if (.not. is_diagonal) hessian_scale = hessian_scale_of(n, H)
        eps = prox_reg_scaled(work, hessian_scale)
        if (eps <= 0.0_wp) then
            flag = daqp_exit_nonconvex
            return
        end if
        work%n_prox = n
        work%prox_mask(1:n) = .true.
    end if

    ! diagonal case
    if (is_diagonal) then
        if (.not. is_factored) then
            if (hessian_scale > 0.0_wp) factor_tol = zero_tol*hessian_scale
            eps = prox_reg_scaled(work, hessian_scale)
        end if
        ! allow small-scale Hessians without tightening the absolute
        ! acceptance threshold for large-scale Hessians
        acceptance_tol = min(factor_tol, zero_tol)
        work%rmode = rinv_diag
        dmin = daqp_inf
        dmax = 0.0_wp
        do i = 1, n
            if (is_factored) then
                hi = Rf(ridx(i,i,n))
                if (hi <= zero_tol) then
                    flag = daqp_exit_nonconvex
                    return
                end if
            else
                hi = H(i,i)
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
            end if
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
    if (.not. is_factored .and. .not. regularize_all .and. .not. allocated(work%avi)) then
        ! zero rows of H are decoupled => regularize only them (semi-proximal)
        hessian_scale = hessian_scale_of(n, H)
        do i = 1, n
            if (H(i,i) /= 0.0_wp) cycle
            zero_row = .true.
            do j = 1, n
                if (H(j,i) /= 0.0_wp .or. H(i,j) /= 0.0_wp) then
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
        i = form_R(work, regularize_all, eps, min_pivot, max_pivot, H, Rf)
        if (i < 0) then
            flag = i
            return
        end if
        ok = i == 1
        if (ok .and. .not. is_factored) then
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
                installed = .false.
                if (allocated(work%eq)) installed = work%eq%installed
                if (ieee_is_finite(dmin) .and. ieee_is_finite(dmax) .and. dmin > 0.0_wp) then
                    if (ieee_is_finite(cond)) then
                        if (real(n,wp)*eps_mach*cond <= hessian_cond_eps .and. &
                            .not. (installed .and. cond > hessian_cond_max)) then
                            work%state = ior(work%state, state_cholesky_pending)
                            return
                        end if
                    end if
                end if
            end if
            if (finish_Rinv(work, is_factored, H)) return
        end if
        if (is_factored) then ! (not reached: a factor is never regularized)
            flag = daqp_exit_nonconvex
            return
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
            hessian_scale = hessian_scale_of(n, H)
            eps = prox_reg_scaled(work, hessian_scale)
            if (eps <= 0.0_wp) then
                flag = daqp_exit_nonconvex
                return
            end if
            regularize_all = .true.
            work%n_prox = n
            work%prox_mask(1:n) = .true.
        end if
    end do

    end function update_R
!*****************************************************************************************

!*****************************************************************************************
!>
!  Factor the Hessian of the problem in the workspace (`daqp_update_Rinv` on
!  `qp->H`).

    integer(ip) function update_R_qp(work,defer_inverse) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: defer_inverse !! defer the inverse of the factor

    if (work%qp%problem_type == daqp_problem_factored) then
        flag = update_R(work, defer_inverse, Rf=work%qp%Rf)
    else if (work%qp%has_H) then
        flag = update_R(work, defer_inverse, H=work%qp%Hc)
    else
        flag = update_R(work, defer_inverse)
    end if

    end function update_R_qp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Complete the inverse of the factor of the Hessian of the problem in the workspace.

    logical function finish_Rinv_qp(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (work%qp%problem_type == daqp_problem_factored) then
        finish_Rinv_qp = finish_Rinv(work, .true.)
    else
        finish_Rinv_qp = finish_Rinv(work, .false., work%qp%Hc)
    end if

    end function finish_Rinv_qp
!*****************************************************************************************

!*****************************************************************************************
!>
!  The proximal regularization that the factor of `H` was formed with.

    real(wp) function get_proximal_regularization(work) result(eps)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, n
    real(wp) :: recovered, rinv

    eps = 0.0_wp
    n = work%n
    if (work%n_prox == 0 .or. .not. work%has_qp) return
    if (.not. work%qp%has_H) return

    if (work%rmode == rinv_diag .and. work%n_prox < n) then
        i = findloc(work%prox_mask(1:n), .true., dim=1)
        eps = 1.0_wp/(work%R(i)*work%R(i)) - work%qp%Hc(i,i)
        return
    end if
    if (work%rmode == rinv_diag) then
        ! diagonal regularization has no retry loop, so reproduce its
        ! scale-based floor directly
        eps = abs(work%settings%eps_prox)
        if (work%qp%problem_type /= daqp_problem_factored) then
          eps = prox_reg_scaled(work, hessian_scale_of(n, work%qp%Hc))
        end if
        return
    end if

    ! semi-proximal: recover eps from a regularized row of Rinv (e_i/sqrt(eps))
    if (work%n_prox < n) then
        i = findloc(work%prox_mask(1:n), .true., dim=1)
        if (i <= work%ms .and. has(work%state,state_rinv_normalized)) then
            rinv = 1.0_wp/work%scaling(i)
        else
            rinv = work%R(ridx(i,i,n))
        end if
        eps = 1.0_wp/(rinv*rinv)
        return
    end if

    ! handle the eps-shift correctly for simple bounds
    rinv = work%R(1)
    if (work%ms > 0) rinv = rinv/work%scaling(1)
    recovered = 1.0_wp/(rinv*rinv) - work%qp%Hc(1,1)

    eps = prox_reg_scaled(work, hessian_scale_of(n, work%qp%Hc))
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

    if (.not. work%has_v .or. .not. work%qp%has_f) return
    work%v(1:work%n) = work%qp%f(1:work%n)
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
                work%Mr(j,k) = work%qp%At(j,k)/work%scaling(j)
            end do
            do j = ns+1, n
                work%Mr(j,k) = work%qp%At(j,k)
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
                work%Mr(i,k) = work%qp%At(i,k)*work%R(i)
            end do
        end do
    case default ! copy A to M
        do k = 1, mA
            work%Mr(:,k) = work%qp%At(:,k)
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
        work%dupper(i) = work%qp%bupper(i)*work%scaling(i)
        work%dlower(i) = work%qp%blower(i)*work%scaling(i)
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
        diff = work%qp%bupper(i) - work%qp%blower(i)
        if (diff < -work%settings%primal_tol) then ! trivial infeasibility
            do_activate = daqp_exit_infeasible
            return
        else if (diff < work%settings%zero_tol .and. .not. has(work%sense(i),daqp_soft)) then
            ! unmarked equality constraint (blower == bupper)
            if (.not. has(work%sense(i),daqp_auto_equality) .or. .not. has(work%sense(i),daqp_active)) then
              do_activate = 1
            end if
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
            if (work%qp%bupper(i) < -zero_tol .or. work%qp%blower(i) > zero_tol) then
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
    if (work%qp%has_f) then
        vs(1:n) = work%qp%f
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
    ! not a standard QP/AVI
    if (allocated(work%bnb) .or. is_hierarchical(work) .or. work%n_prox > 0) return
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
    else if (is_avi_nonsym(work)) then
        ! AVI: the unconstrained solution is x = -H^{-1} f
        if (work%qp%has_f) then
            call daqp_lu_solve(work%avi%LU_H, work%avi%P_H, work%qp%f, work%x, n)
        else
            work%x(1:n) = 0.0_wp
        end if
        do i = 1, n
            work%x(i) = -work%x(i)
        end do
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
        work%x(1:n) = 0.0_wp ! no linear term: the unconstrained optimum is x = 0
    end if

    ! check the simple bounds
    do i = 1, work%ms
        work%dupper(i) = work%qp%bupper(i) - work%x(i)
        work%dlower(i) = work%qp%blower(i) - work%x(i)
        if (work%dupper(i) < -primal_tol .or. work%dlower(i) > primal_tol) feasible = .false.
    end do
    ! check the general constraints
    do i = work%ms+1, work%m
        s = dot_seq(n, work%qp%At(:,i-work%ms), work%x)
        work%dupper(i) = work%qp%bupper(i) - s
        work%dlower(i) = work%qp%blower(i) - s
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
!  Form the LDP of the problem that is in the workspace (the reduced problem
!  of an equality elimination is passed here as it is), as marked by
!  `mask_in`. Returns 0, or a negative exit flag.

    integer(ip) function update_ldp_core(mask_in,work) result(flag)

    integer(ip), intent(in) :: mask_in          !! the parts to update (`daqp_update_*`)
    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: mask, unconstrained_flag, i, m_tmp
    logical :: do_activate

    do_activate = .false.
    unconstrained_flag = 0

    ! also form what an earlier update left pending (everything stays pending
    ! until this update completes, so an update that fails is redone)
    mask = ior(mask_in, iand(work%state, state_pending))
    work%state = ior(iand(work%state, state_rinv_normalized + state_ill_conditioned + &
                                      state_cholesky_pending), &
                     iand(mask, state_pending))

    ! dimensions of the problem
    work%n = work%qp%n
    work%m = work%qp%m
    work%ms = work%qp%ms

    ! update the constraint flags
    if (has(mask,daqp_update_sense)) then
        if (.not. work%qp%has_sense) then ! all constraints are inequalities
            work%sense(1:work%m) = 0
        else
            work%sense(1:work%m) = work%qp%sense(1:work%m)
            do_activate = .true.
        end if
    end if

    ! check the bounds early
    if (iand(mask, daqp_update_m+daqp_update_v+daqp_update_d+daqp_update_sense) /= 0) then
        flag = check_bounds(work)
        if (flag < 0) return
        if (flag == 1) do_activate = .true.
    end if

    ! form R first; dense QPs defer the inverse until after the candidate check
    if (has(mask,daqp_update_rinv)) then
        if (.not. allocated(work%avi)) then
            flag = update_R_qp(work, .true.)
        else
            call update_avi(work)
            if (work%avi%is_symmetric) then
                flag = update_R(work, .false., H=work%qp%Hc)
            else
                ! early unconstrained check for an AVI: skip the factorization
                ! if x = -H^{-1}f is feasible
                unconstrained_flag = check_unconstrained(work, mask)
                if (unconstrained_flag == unconstrained_optimal) then
                    flag = 0
                    return
                end if
                i = daqp_lu(work%avi%H_rho, work%avi%P_H2, work%n)
                flag = update_R(work, .false., H=work%avi%Hs_rho)
            end if
        end if
        if (flag < 0) return
    end if

    ! update v (if Rinv still holds R, the unconstrained check forms v itself)
    if (.not. has(work%state,state_cholesky_pending) .and. &
        iand(mask, daqp_update_rinv+daqp_update_v) /= 0) call update_v(work)

    if (.not. is_avi_nonsym(work)) unconstrained_flag = check_unconstrained(work, mask)
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
        if (.not. finish_Rinv_qp(work)) then
            ! refactor with regularization => v and d from the check are stale
            work%state = ior(work%state, daqp_update_rinv)
            mask = ior(mask, daqp_update_rinv)
            flag = update_R_qp(work, .false.)
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

    ! update the hierarchy
    if (has(mask,daqp_update_hierarchy)) then
        work%nh = work%qp%nh
        work%has_bp = work%qp%nh > 1
        if (work%has_bp) then
            work%break_points = work%qp%break_points
        else if (allocated(work%break_points)) then
            deallocate(work%break_points)
        end if
    end if

    ! hierarchies are not allowed for prox and nonsymmetric AVIs
    if (is_hierarchical(work) .and. (work%n_prox > 0 .or. is_avi_nonsym(work))) then
        flag = daqp_exit_unsupported
        return
    end if

    flag = 0
    ! an empty working set can be one that a reset left out
    if (do_activate .or. work%n_active == 0) then
        call daqp_reset_workspace(work)
        if (.not. is_hierarchical(work)) then
            flag = daqp_activate_constraints(work)
        else ! activate the first level (since those constraints are hard)
            m_tmp = work%m
            work%m = work%break_points(1)
            flag = daqp_activate_constraints(work)
            work%m = m_tmp
        end if
    end if
    if (flag < 0) return
    work%state = iand(work%state, not(state_pending)) ! everything has been formed
    flag = 0

    end function update_ldp_core
!*****************************************************************************************

!*****************************************************************************************
!>
!  Update the workspace with (changes in) the problem in `work%qp`, as marked
!  by `mask_in` (upstream's `daqp_update_ldp`). If the equality constraints are
!  to be eliminated, the LDP is formed for the reduced problem, and the
!  workspace keeps describing the original problem otherwise.
!  Returns 0, or a negative exit flag.

    integer(ip) function daqp_update_ldp(work,mask_in) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: mask_in !! the parts to update (`daqp_update_*`)

    integer(ip) :: mask
    logical :: was_reduced

    mask = mask_in
    was_reduced = is_reduced(work)
    call daqp_eq_restore(work)
    ! the constraint flags are indexed by the original problem, also while its
    ! equality constraints are eliminated
    if (has(mask,daqp_update_sense)) then
        if (.not. work%qp%has_sense) then
            work%sense(1:work%qp%m) = 0
        else
            work%sense(1:work%qp%m) = work%qp%sense(1:work%qp%m)
        end if
    end if
    if (daqp_eq_wanted(work, mask)) then
        flag = daqp_eq_update(work, mask, update_ldp_core)
        if (flag /= eq_not_reduced) then
            work%n = work%qp%n
            work%m = work%qp%m
            work%ms = work%qp%ms
            return
        end if
    else
        call daqp_eq_deactivate(work)
    end if
    ! the LDP of the problem was not formed while its equalities were eliminated
    if (was_reduced) then
        mask = ior(mask, daqp_update_m + daqp_update_d)
        if (work%qp%has_H .or. work%qp%problem_type == daqp_problem_factored) then
          mask = ior(mask, daqp_update_rinv)
        end if
        if (work%qp%has_f) mask = ior(mask, daqp_update_v)
    end if
    flag = update_ldp_core(mask, work)

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
            val = dot_seq(work%n, work%qp%At(:,i-work%ms), work%x)
        end if
        if (val > work%qp%bupper(i)+tol .or. val < work%qp%blower(i)-tol) then
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
        r = sqrt((work%qp%Hc(i,i)+eps)/(work%qp%Hc(i,i)+eps_new))
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
            r = sqrt((work%qp%Hc(j,j)+eps)/(work%qp%Hc(j,j)+eps_new))
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
            gap = gap - q*work%qp%blower(id)
        else
            gap = gap - q*work%qp%bupper(id)
        end if
        if (lower) then
            wrong_sign = q > 0.0_wp
        else
            wrong_sign = q < 0.0_wp
        end if
        if (.not. has(work%sense(id),daqp_immutable) .and. wrong_sign) then
            if (work%qp%bupper(id) < daqp_inf .and. work%qp%blower(id) > -daqp_inf) then
                gap = gap - abs(q)*(work%qp%bupper(id)-work%qp%blower(id))
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
    if (na >= 9*n/10 .or. work%sing_ind /= empty_ind .or. .not. work%qp%has_H) return
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
    if (work%qp%has_f) then
        do i = 1, n
            gd = gd + work%qp%f(i)*work%xldl(i)
        end do
    end if
    if (work%qp%has_H) then
        do i = 1, n
            hii = abs(work%qp%Hc(i,i))
            s = 0.0_wp
            do j = 1, n
                s = s + work%qp%Hc(j,i)*work%xldl(j)
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
                ax = ax + work%qp%At(j,k)*work%x(j)
                ad = ad + work%qp%At(j,k)*(work%x(j)-work%xold(j))
            end do
        end if
        if (ad > 0.0_wp .and. work%qp%bupper(i) < daqp_inf) then
            sb = (work%qp%bupper(i)-ax)/ad
        else if (ad < 0.0_wp .and. work%qp%blower(i) > -daqp_inf) then
            sb = (work%qp%blower(i)-ax)/ad
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
                    ad = dot_seq(work%n, work%qp%At(:,id-work%ms), work%xldl)
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
            else if (.not. work%qp%has_H .and. .not. moved .and. .not. skipped) then
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
    real(wp), parameter :: relaxation = 1.5_wp

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
    adaptive = .not. is_lp .and. .not. all_pd .and. .not. allocated(work%avi) .and. work%n_prox < nx
    eps_max = eps
    rescaled = .false.
    if (adaptive) then
        hmax = hessian_scale_of(nx, work%qp%Hc)
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
        if (work%settings%dual_tol /= daqp_default_dual_tol .and. 0.1_wp*work%settings%dual_tol < eta) then
          eta = 0.1_wp*work%settings%dual_tol
        end if
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
                work%v(i) = work%qp%f(i)*eps - work%x(i)
            end do
        else
            if (work%n_prox == nx) then ! full shift
                if (work%qp%has_f) then
                    do i = 1, nx
                        work%v(i) = work%qp%f(i) - eps*work%x(i)
                    end do
                else
                    do i = 1, nx
                        work%v(i) = -eps*work%x(i)
                    end do
                end if
            else ! regularize only the singular directions
                do i = 1, nx
                    if (work%prox_mask(i)) then
                        if (work%qp%has_f) then
                            work%v(i) = work%qp%f(i) - eps*work%x(i)
                        else
                            work%v(i) = -eps*work%x(i)
                        end if
                    else
                        if (work%qp%has_f) then
                            work%v(i) = work%qp%f(i) - 0.0_wp*work%x(i)
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
                work%x(1:nx) = work%xold(1:nx) ! the center
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
            if (.not. is_lp .and. .not. allocated(work%avi) .and. work%settings%eps_prox < 0.0_wp .and. &
                work%nh >= prox_face_start .and. work%iterations <= 3) then
                projected = prox_project_step(work, eps)
            end if
        end if
        if ((work%iterations == 1 .or. projected) .and. work%n_active < nx .and. &
            total_iter < work%settings%iter_limit) then
            if (allocated(work%avi)) then
                do i = 1, nx
                    work%x(i) = work%xold(i) + relaxation*(work%x(i) - work%xold(i))
                end do
                center_relaxed = .true.
            else
                step_flag = prox_step(work, s_prev, eps, projected)
                if (step_flag == daqp_exit_unbounded) then
                    exitflag = daqp_exit_unbounded
                    exit
                end if
                center_relaxed = step_flag == 1
            end if
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
!  Deallocate the workspace (the settings are kept).

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
!  Set up the workspace for a problem, and form its LDP (upstream's
!  `setup_daqp_main`, with the problem given in the Fortran layout: `H(n,n)`,
!  `A(m-ms,n)`, or `A` transposed in `At(n,m-ms)`). The settings are taken
!  from `work%settings`, which the caller sets first.
!
!  Returns 1, or a negative exit flag.

    integer(ip) function daqp_setup(work,n,m,ms,bupper,blower,H,f,A,sense,init_mask,At, &
                                    break_points,problem_type,Rf,primal_start,dual_start, &
                                    setup_time) result(flag)

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
    integer(ip), intent(in), optional :: init_mask !! extra update mask (`daqp_update_unconstrained`, `daqp_update_eliminate`)
    real(wp), intent(in), optional :: At(n,m-ms) !! the transpose of `A` (instead of `A`)
    integer(ip), intent(in), optional :: break_points(:) !! the last constraint of each level (hierarchical QP)
    integer(ip), intent(in), optional :: problem_type !! `daqp_problem_qp` (default) or `daqp_problem_avi`
    real(wp), intent(in), optional :: Rf(:)     !! the Cholesky factor of `H`, packed by rows (instead of `H`)
    real(wp), intent(in), optional :: primal_start(n) !! a primal start (its active constraints start the working set)
    real(wp), intent(in), optional :: dual_start(m)   !! a dual start (its nonzero multipliers start the working set)
    real(wp), intent(out), optional :: setup_time    !! time of the setup [s]

    integer(ip) :: istat, nw, ns, nb, mask, i, start
    integer(int64) :: t0, t1, rate

    call system_clock(t0, rate)
    if (present(setup_time)) setup_time = 0.0_wp
    ! (arrays of the right size are kept, for repeated setups)
    work%is_setup = .false.
    if (allocated(work%eq)) deallocate(work%eq)
    if (allocated(work%bnb)) deallocate(work%bnb)
    if (allocated(work%avi)) deallocate(work%avi)
    if (allocated(work%break_points)) deallocate(work%break_points)
    if (allocated(work%qp%break_points)) deallocate(work%qp%break_points)
    if (allocated(work%qp%Rf)) deallocate(work%qp%Rf)
    work%has_weights = .false.

    ! the problem
    work%qp%n = n
    work%qp%m = m
    work%qp%ms = ms
    work%qp%problem_type = daqp_problem_qp
    if (present(problem_type)) work%qp%problem_type = problem_type
    if (present(Rf)) work%qp%problem_type = daqp_problem_factored
    work%qp%has_H = present(H) .or. present(Rf)
    work%qp%has_f = present(f)
    work%qp%has_sense = present(sense) .or. present(primal_start) .or. present(dual_start)
    ! the packed triangles are indexed by default integers
    if ((int(n,int64)+1_int64)*(int(n,int64)+int(m,int64)+2_int64)/2_int64 > int(huge(1_ip),int64) .or. &
        int(n,int64)*int(max(n,m-ms),int64) > int(huge(1_ip),int64)) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if
    istat = 0
    if (present(H)) then
        call resize2(work%qp%Hc, n, n, istat)
    else
        call resize2(work%qp%Hc, 0_ip, 0_ip, istat)
    end if
    if (present(Rf)) call resize1(work%qp%Rf, (n*(n+1))/2, istat)
    call resize1(work%qp%f, n, istat)
    call resize2(work%qp%At, n, m-ms, istat)
    call resize1(work%qp%bupper, m, istat)
    call resize1(work%qp%blower, m, istat)
    call resize1i(work%qp%sense, m, istat)
    if (istat /= 0) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if
    if (present(H)) work%qp%Hc = transpose(H)
    if (present(Rf)) work%qp%Rf = Rf(1:(n*(n+1))/2)
    if (present(f)) then
        work%qp%f = f
    else
        work%qp%f = 0.0_wp
    end if
    if (present(At)) then
        work%qp%At = At
    else if (present(A)) then
        call transpose_into(A, work%qp%At)
    else
        work%qp%At = 0.0_wp
    end if
    work%qp%bupper = bupper
    work%qp%blower = blower
    work%qp%sense = 0
    if (present(sense)) work%qp%sense = sense
    work%qp%nh = 1
    if (present(break_points)) then
        work%qp%nh = int(size(break_points), ip)
        work%qp%break_points = break_points
    end if
    ! the starting working set from a dual or a primal start
    if (present(dual_start)) then
        call daqp_dual_init_active(work%qp, dual_start)
    else if (present(primal_start)) then
        call daqp_primal_init_active(work%qp, primal_start)
    end if

    ! count the soft and binary constraints (to account for them in the allocation)
    ns = count(iand(work%qp%sense, daqp_soft) /= 0)
    nb = count(iand(work%qp%sense, daqp_binary) /= 0)
    if (work%qp%nh > 1) then ! the largest level, if several hierarchies
        ns = 0
        start = 0
        do i = 1, work%qp%nh
            ns = max(ns, work%qp%break_points(i)-start)
            start = work%qp%break_points(i)
        end do
    end if

    ! the iterates
    nw = n + ns
    work%n = n
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
    if (istat /= 0) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if
    ! (the work arrays are written before they are read, as upstream's malloc'ed ones)
    work%x = 0.0_wp  ! an uninitialized iterate is 0
    work%xold = 0.0_wp
    work%D(1) = 0.0_wp
    work%prox_mask = .false.
    work%n_prox = 0
    work%state = 0
    work%soft_slack = 0.0_wp
    work%has_weights = .false.
    work%nh = 1
    work%has_bp = .false.
    work%timer_on = .false.
    work%has_qp = .true.
    call daqp_reset_workspace(work)

    if (work%qp%problem_type == daqp_problem_avi) then
        allocate(work%avi)
        call allocate_avi(work%avi, n)
    end if

    ! branch and bound
    if (nb > n) then
        call daqp_destroy(work)
        flag = daqp_exit_overdetermined_initial
        return
    end if
    if (nb > 0) then
        allocate(work%bnb)
        work%bnb%nb = nb
        allocate(work%bnb%bin_ids(nb), work%bnb%tree(nb+2), work%bnb%tree_ws((nw+1)*(nb+1)), &
                 work%bnb%fixed_ids(nb+1), work%bnb%root_ws(nw+1))
        nb = 0
        do i = 1, m
            if (has(work%qp%sense(i),daqp_binary)) then
                nb = nb + 1
                work%bnb%bin_ids(nb) = i
            end if
        end do
        work%bnb%n_nodes = 0
        work%bnb%nws = 0
        work%bnb%n_root_ws = 0
    end if

    ! the LDP: always update M, d and sense
    mask = daqp_update_m + daqp_update_d + daqp_update_sense
    if (present(init_mask)) mask = ior(mask, init_mask)
    if (work%qp%has_H) mask = ior(mask, daqp_update_rinv)
    if (work%qp%has_f) mask = ior(mask, daqp_update_v)
    ! for an LP, mark all directions as needing proximal regularization
    if (.not. work%qp%has_H .and. work%qp%has_f) work%n_prox = n
    work%m = m
    work%ms = ms
    call resize1(work%scaling, m, istat)
    call resize2(work%Mr, n, m-ms, istat)
    call resize1(work%Mu, m-ms, istat)
    call resize1(work%dupper, m, istat)
    call resize1(work%dlower, m, istat)
    call resize1i(work%sense, m, istat)
    call resize1(work%v, n, istat)
    if (work%qp%has_H) then
        call resize1(work%R, (n*(n+1))/2, istat)
    else
        call resize1(work%R, 0_ip, istat)
    end if
    if (istat /= 0) then
        call daqp_destroy(work)
        flag = daqp_exit_out_of_memory
        return
    end if
    work%scaling = 1.0_wp
    work%v = 0.0_wp
    work%sense = 0
    work%has_v = work%qp%has_f
    if (work%qp%has_H) then
        work%rmode = rinv_dense
    else
        work%rmode = rinv_none
    end if
    if (work%qp%nh > 1) mask = ior(mask, daqp_update_hierarchy)

    istat = daqp_update_ldp(work, mask)
    if (istat < 0) then
        call daqp_destroy(work)
        flag = istat
        return
    end if

    ! a singular quadratic needs v for the proximal linear term even without f
    if (work%n_prox > 0 .and. .not. work%has_v .and. work%qp%has_H) then
        work%has_v = .true.
        work%v = 0.0_wp
    end if

    work%is_setup = .true.
    flag = 1
    call system_clock(t1)
    if (present(setup_time)) setup_time = real(t1-t0,wp)/real(rate,wp)

    if (present(primal_start)) call daqp_set_primal_start(work, primal_start)

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

    integer(ip) :: i, id, c, sm
    logical :: installed

    call daqp_eq_restore(work)
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
    if (is_reduced(work)) then
        ! pass the flags to the constraints of the reduced problem
        sm = daqp_active + daqp_lower
        do c = 1, work%eq%mr
            if (has(work%eq%other%sense(c),daqp_immutable)) cycle
            work%eq%other%sense(c) = ior(iand(work%eq%other%sense(c), not(sm)), &
                                         iand(work%sense(work%eq%keep(c)), sm))
        end do
        installed = daqp_eq_install(work)
        call daqp_reset_workspace(work)
        flag = daqp_activate_constraints(work)
        call daqp_eq_restore(work)
    else
        call daqp_reset_workspace(work)
        flag = daqp_activate_constraints(work)
    end if

    end function daqp_set_working_set
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate the individual weights of the soft constraints (zero, which selects
!  the settings). Returns false if there are no constraints.

    logical function daqp_allocate_soft_weights(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: m, istat

    daqp_allocate_soft_weights = .true.
    if (work%has_weights) return ! already allocated
    ! (the weights are indexed by the original problem, also if its equality
    ! constraints are eliminated)
    m = work%m
    daqp_allocate_soft_weights = .false.
    if (m == 0) return
    istat = 0
    call resize1(work%rho_ls, m, istat)
    call resize1(work%rho_us, m, istat)
    call resize1(work%w_ls, m, istat)
    call resize1(work%w_us, m, istat)
    if (istat /= 0) return
    work%rho_ls = 0.0_wp
    work%rho_us = 0.0_wp
    work%w_ls = 0.0_wp
    work%w_us = 0.0_wp
    work%has_weights = .true.
    daqp_allocate_soft_weights = .true.

    end function daqp_allocate_soft_weights
!*****************************************************************************************

!*****************************************************************************************
!>
!  Refactor the working set after a change of the weights of the soft
!  constraints (upstream's `daqp_refresh_soft_weights`).

    subroutine daqp_refresh_soft_weights(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, flag, m
    logical :: reduced, rebuild

    ! the working set is that of the reduced problem if equalities are eliminated
    reduced = daqp_eq_install(work)
    ! only an active soft constraint makes the factorization stale
    rebuild = .false.
    do i = 1, work%n_active
        if (has(work%sense(work%WS(i)),daqp_soft)) then
            rebuild = .true.
            exit
        end if
    end do
    if (rebuild) then
        call daqp_reset_workspace(work)
        if (is_hierarchical(work)) then
            m = work%m
            work%m = work%break_points(1)
            flag = daqp_activate_constraints(work)
            work%m = m
        else
            flag = daqp_activate_constraints(work)
        end if
    end if
    if (reduced) call daqp_eq_restore(work)

    end subroutine daqp_refresh_soft_weights
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the weights of the soft constraints, one entry per constraint of the
!  original problem (an absent argument leaves that weight untouched; a zero
!  weight selects the settings). A violation `s` of a soft constraint adds
!  `w*s + s^2/(2*rho)` to the objective. Returns false if they could not be
!  allocated.

    logical function daqp_set_soft_weights(work,rho_l,rho_u,w_l,w_u)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in), optional :: rho_l(:) !! reciprocal quadratic weight of the lower side
    real(wp), intent(in), optional :: rho_u(:) !! reciprocal quadratic weight of the upper side
    real(wp), intent(in), optional :: w_l(:)   !! linear weight of the lower side
    real(wp), intent(in), optional :: w_u(:)   !! linear weight of the upper side

    integer(ip) :: m

    daqp_set_soft_weights = daqp_allocate_soft_weights(work)
    if (.not. daqp_set_soft_weights) return
    m = size(work%rho_ls)
    if (present(rho_l)) work%rho_ls(1:m) = rho_l(1:m)
    if (present(rho_u)) work%rho_us(1:m) = rho_u(1:m)
    if (present(w_l)) work%w_ls(1:m) = w_l(1:m)
    if (present(w_u)) work%w_us(1:m) = w_u(1:m)
    call daqp_refresh_soft_weights(work)

    end function daqp_set_soft_weights
!*****************************************************************************************

!*****************************************************************************************
!>
!  Package the result (upstream's `daqp_extract_result`).

    subroutine extract_result(work,x,res,lam)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(inout) :: x(:)             !! the solution
    type(daqp_result), intent(inout) :: res     !! result
    real(wp), intent(inout), optional :: lam(:) !! multipliers

    integer(ip) :: i, id
    real(wp) :: l
    logical :: wrong_sign

    x(1:work%n) = work%x(1:work%n)

    ! multipliers of ordinary QPs (hierarchical QPs form theirs in daqp_hiqp)
    if (present(lam) .and. .not. is_hierarchical(work)) then
        lam(1:work%m) = 0.0_wp
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
    if (work%has_v .and. .not. is_avi_nonsym(work) .and. work%rmode /= rinv_none) then ! QP or symmetric AVI
        res%fval = work%fval
        do i = 1, work%n
            res%fval = res%fval - work%v(i)*work%v(i)
        end do
        res%fval = res%fval*0.5_wp
    else if (work%has_qp .and. work%qp%has_f) then ! LP (or AVI)
        res%fval = 0.0_wp
        do i = 1, work%n
            res%fval = res%fval + work%qp%f(i)*work%x(i)
        end do
    else if (.not. is_hierarchical(work)) then
        ! no linear term: upstream leaves fval unset (which an equality
        ! elimination then forms as 0.5*fval, the objective)
        res%fval = 0.5_wp*work%fval
    end if

    res%soft_slack = work%soft_slack
    res%iter = work%iterations
    if (allocated(work%bnb)) then
        res%nodes = work%bnb%nodecount
    else if (is_hierarchical(work)) then
        res%nodes = 1
    else
        res%nodes = work%nh
    end if

    end subroutine extract_result
!*****************************************************************************************

!*****************************************************************************************
!>
!  The multipliers of the working set, as they are (upstream's
!  `daqp_extract_active_duals`).

    subroutine daqp_extract_active_duals(work,lam)

    type(daqp_workspace), intent(in) :: work !! workspace
    real(wp), intent(out) :: lam(:)          !! multipliers (size `m`)

    integer(ip) :: i

    lam = 0.0_wp
    do i = 1, work%n_active
        lam(work%WS(i)) = work%lam_star(i)
    end do

    end subroutine daqp_extract_active_duals
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the problem from the current working set (upstream's `daqp_solve`).

    subroutine daqp_solve(work,x,res,lam)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(out) :: x(:)               !! the solution (size `n`)
    type(daqp_result), intent(out) :: res       !! exit flag, objective, iterations
    real(wp), intent(out), optional :: lam(:)   !! multipliers (size `m`)

    integer(ip) :: flag
    integer(int64) :: t0, t1, rate
    logical :: reduced

    x = 0.0_wp
    if (present(lam)) lam = 0.0_wp
    if (.not. work%is_setup) then
        res%exitflag = daqp_exit_not_setup
        return
    end if

    call system_clock(t0, rate)
    if (.not. work%has_bp) work%nh = 1
    ! a policy that no longer eliminates takes effect before the solve
    if (is_reduced(work) .and. work%settings%eq_reduction == daqp_eq_reduction_off) then
        flag = daqp_update_ldp(work, 0_ip)
        if (flag < 0) then
            res%exitflag = flag
            return
        end if
    end if
    ! solve the reduced problem if the equalities are eliminated (unless the
    ! latest update found its right-hand side to be infeasible)
    if (is_reduced(work)) then
        if (work%eq%error < 0) then
            res%exitflag = work%eq%error
            return
        end if
    end if
    reduced = daqp_eq_install(work)
    work%timer_on = work%settings%time_limit > 0.0_wp
    work%timer_start = t0

    if (.not. has(work%state,state_unconstrained)) then
        if (work%n_prox == 0) then ! select the algorithm
            if (.not. is_avi_nonsym(work)) then
                if (allocated(work%bnb)) then
                    res%exitflag = daqp_bnb(work)
                else if (is_hierarchical(work)) then
                    res%exitflag = daqp_hiqp(work, lam)
                else
                    res%exitflag = daqp_ldp(work)
                end if
                if (res%exitflag > 0) then
                    call daqp_ldp2qp_solution(work) ! retrieve the QP solution
                    ! refine x (if it might be inaccurate)
                    if (.not. allocated(work%bnb) .and. .not. is_hierarchical(work)) then
                        call refine_primal(work)
                        ! constraints were only added above the rounding level
                        if (has(work%state,state_noise_floor)) then
                            if (violates_hard(work)) res%exitflag = daqp_exit_optimal_inexact
                        end if
                    end if
                end if
            else ! AVI
                res%exitflag = daqp_solve_avi(work)
            end if
        else ! proximal
            if (allocated(work%bnb)) then
                res%exitflag = daqp_exit_nonconvex
            else
                res%exitflag = daqp_prox(work)
            end if
        end if
    else ! unconstrained optimum
        work%iterations = 1
        work%fval = 0.0_wp
        work%soft_slack = 0.0_wp
        res%exitflag = daqp_exit_optimal
    end if
    work%timer_on = .false.

    ! package the result
    call extract_result(work, x, res, lam)
    if (reduced) then
        call daqp_eq_expand(work, x, res%fval, lam)
        call daqp_eq_restore(work)
    end if
    call system_clock(t1)
    res%solve_time = real(t1-t0,wp)/real(rate,wp)

    end subroutine daqp_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up and solve a problem in one call (upstream's `daqp_quadprog`: with the
!  check of the unconstrained optimum, and the automatic elimination of
!  equalities).

    subroutine daqp_quadprog(n,m,ms,bupper,blower,x,res,lam,H,f,A,sense,settings, &
                             break_points,problem_type,Rf,primal_start,dual_start)

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
    integer(ip), intent(in), optional :: break_points(:) !! the last constraint of each level (hierarchical QP)
    integer(ip), intent(in), optional :: problem_type !! `daqp_problem_qp` (default) or `daqp_problem_avi`
    real(wp), intent(in), optional :: Rf(:)     !! the Cholesky factor of `H`, packed by rows (instead of `H`)
    real(wp), intent(in), optional :: primal_start(n) !! a primal start
    real(wp), intent(in), optional :: dual_start(m)   !! a dual start

    type(daqp_workspace) :: work
    integer(ip) :: flag
    real(wp) :: setup_time

    if (present(settings)) work%settings = settings
    flag = daqp_setup(work, n, m, ms, bupper, blower, H, f, A, sense, &
                      init_mask=daqp_update_unconstrained+daqp_update_eliminate, &
                      break_points=break_points, problem_type=problem_type, Rf=Rf, &
                      primal_start=primal_start, dual_start=dual_start, setup_time=setup_time)
    if (flag < 0) then
        res%exitflag = flag
        x = 0.0_wp
        if (present(lam)) lam = 0.0_wp
        return
    end if
    call daqp_solve(work, x, res, lam)
    res%setup_time = setup_time
    call daqp_destroy(work)

    end subroutine daqp_quadprog
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve an affine variational inequality in one call (upstream's `daqp_avi`):
!  find `x` with `blower <= [x(1:ms); A x] <= bupper` and
!  `(H x + f)'(y - x) >= 0` for all feasible `y` (`H` need not be symmetric).

    subroutine daqp_avi(n,m,ms,H,f,bupper,blower,x,res,lam,A,sense,settings)

    integer(ip), intent(in) :: n   !! number of variables
    integer(ip), intent(in) :: m   !! number of constraints (including the simple bounds)
    integer(ip), intent(in) :: ms  !! number of simple bounds
    real(wp), intent(in) :: H(n,n)    !! the matrix of the AVI
    real(wp), intent(in) :: f(n)      !! linear term
    real(wp), intent(in) :: bupper(m) !! upper bounds
    real(wp), intent(in) :: blower(m) !! lower bounds
    real(wp), intent(out) :: x(n)     !! the solution
    type(daqp_result), intent(out) :: res     !! exit flag, iterations
    real(wp), intent(out), optional :: lam(m) !! multipliers
    real(wp), intent(in), optional :: A(m-ms,n) !! constraint matrix
    integer(ip), intent(in), optional :: sense(m) !! constraint flags
    type(daqp_settings), intent(in), optional :: settings !! settings

    call daqp_quadprog(n, m, ms, bupper, blower, x, res, lam, H, f, A, sense, settings, &
                       problem_type=daqp_problem_avi)

    end subroutine daqp_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Mark the constraints that are active at `x` (within 1e-9) as active in the
!  constraint flags of the problem, as a starting working set (upstream's
!  `daqp_primal_init_active`; nothing is done for a problem with binary
!  constraints, for which `x` is only used as an incumbent).

    subroutine daqp_primal_init_active(qp,x)

    type(daqp_problem), intent(inout) :: qp !! the problem
    real(wp), intent(in) :: x(:)            !! primal iterate

    integer(ip) :: i, j
    real(wp) :: ax, slack, tol

    tol = max(1.0e-9_wp, 100.0_wp*epsilon(1.0_wp)) ! (upstream's 1e-9, floored in single precision)
    if (.not. allocated(qp%sense)) then
        allocate(qp%sense(qp%m))
        qp%sense = 0
    end if
    qp%has_sense = .true.
    do i = 1, qp%m
        if (has(qp%sense(i),daqp_binary)) return
    end do
    do i = 1, qp%m
        if (has(qp%sense(i),daqp_immutable)) cycle
        if (i <= qp%ms) then
            ax = x(i)
        else
            ax = 0.0_wp
            do j = 1, qp%n
                ax = ax + x(j)*qp%At(j,i-qp%ms)
            end do
        end if
        slack = ax - qp%bupper(i)
        if (slack < tol .and. slack > -tol) then
            qp%sense(i) = ior(qp%sense(i), daqp_active)
            qp%sense(i) = iand(qp%sense(i), not(daqp_lower))
        else
            slack = ax - qp%blower(i)
            if (slack < tol .and. slack > -tol) then
              qp%sense(i) = ior(qp%sense(i), daqp_active+daqp_lower)
            end if
        end if
    end do

    end subroutine daqp_primal_init_active
!*****************************************************************************************

!*****************************************************************************************
!>
!  Mark the constraints with a nonzero multiplier (beyond 1e-12) as active in
!  the constraint flags of the problem, as a starting working set (upstream's
!  `daqp_dual_init_active`).

    subroutine daqp_dual_init_active(qp,lam)

    type(daqp_problem), intent(inout) :: qp !! the problem
    real(wp), intent(in) :: lam(:)          !! dual iterate

    integer(ip) :: i
    real(wp) :: tol

    tol = max(1.0e-12_wp, 10.0_wp*epsilon(1.0_wp)) ! (upstream's 1e-12, floored in single precision)
    if (.not. allocated(qp%sense)) then
        allocate(qp%sense(qp%m))
        qp%sense = 0
    end if
    qp%has_sense = .true.
    do i = 1, qp%m
        if (has(qp%sense(i),daqp_immutable)) cycle
        if (lam(i) > tol) then
            qp%sense(i) = ior(qp%sense(i), daqp_active)
            qp%sense(i) = iand(qp%sense(i), not(daqp_lower))
        else if (lam(i) < -tol) then
            qp%sense(i) = ior(qp%sense(i), daqp_active+daqp_lower)
        end if
    end do

    end subroutine daqp_dual_init_active
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the starting iterate (upstream's `daqp_set_primal_start`); with binary
!  constraints, `x` is used as an incumbent if it is feasible.

    subroutine daqp_set_primal_start(work,x)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in) :: x(:)                !! iterate of the original problem

    logical :: reduced

    ! x is given for the original problem, also if its equalities are eliminated
    reduced = daqp_eq_install(work)
    if (.not. has(work%state,state_unconstrained)) then
        if (reduced) then
            call daqp_eq_set_primal_start(work, x)
        else
            work%x(1:work%n) = x(1:work%n)
        end if
        if (allocated(work%bnb)) work%state = ior(work%state, state_incumbent)
    end if
    if (reduced) call daqp_eq_restore(work)

    end subroutine daqp_set_primal_start
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determine the redundant constraints of a polyhedron `A x <= b` (upstream's
!  `daqp_minrep`): the first `ms` constraints are `x(1:ms) <= b(1:ms)`, the
!  others `A x <= b(ms+1:)`. `is_redundant(i)` is 1 if constraint `i` is
!  redundant, 0 otherwise.

    subroutine daqp_minrep(A,b,ms,is_redundant,settings)

    real(wp), intent(in) :: A(:,:)              !! general constraints `(m-ms,n)`
    real(wp), intent(in) :: b(:)                !! right-hand sides `(m)`
    integer(ip), intent(in) :: ms               !! number of simple bounds
    integer(ip), intent(out) :: is_redundant(:) !! `(m)`
    type(daqp_settings), intent(in), optional :: settings !! settings

    type(daqp_workspace) :: work
    integer(ip) :: n, m, nw

    n = int(size(A,2), ip)
    m = int(size(b), ip)
    if (present(settings)) work%settings = settings
    work%has_qp = .false.
    work%n = n
    work%m = m
    work%ms = ms
    nw = n
    allocate(work%lam(nw+1), work%lam_star(nw+1), work%WS(nw+1), work%D(nw+1), &
             work%xldl(nw+1), work%zldl(nw+1), work%L(((nw+1)*(nw+2))/2), &
             work%x(n), work%xold(n), work%prox_mask(n))
    work%x = 0.0_wp
    work%xold = 0.0_wp
    work%D(1) = 0.0_wp
    work%prox_mask = .false.
    allocate(work%Mr(n,m-ms), work%Mu(m-ms), work%dupper(m), work%dlower(m), work%sense(m), &
             work%scaling(m), work%v(0), work%R(0))
    call transpose_into(A, work%Mr)
    work%dupper = b
    work%dlower = -daqp_inf
    work%sense = 0
    work%scaling = 1.0_wp  ! (no scaling: exactly as upstream's NULL)
    work%rmode = rinv_none
    work%has_v = .false.
    call daqp_reset_workspace(work)
    call minrep_work(work, is_redundant)

    end subroutine daqp_minrep
!*****************************************************************************************

!*****************************************************************************************
!>
!  Determine the redundant constraints of the LDP in the workspace (upstream's
!  `daqp_minrep_work`).

    subroutine minrep_work(work,is_redundant)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(out) :: is_redundant(:) !! 1 if redundant, 0 otherwise

    integer(ip) :: i, j, exitflag

    is_redundant(1:work%m) = -1
    do i = 1, work%m
        if (is_redundant(i) /= -1 .or. has(work%sense(i),daqp_immutable)) cycle
        call daqp_reset_workspace(work)
        work%sense(i) = daqp_active + daqp_immutable
        call add_constraint(work, i, 1.0_wp)
        exitflag = daqp_ldp(work)
        if (exitflag == daqp_exit_infeasible) then
            is_redundant(i) = 1
            work%sense(i) = iand(work%sense(i), not(daqp_active)) ! (remains immutable -> ignored)
        else
            is_redundant(i) = 0
            work%sense(i) = iand(work%sense(i), not(daqp_immutable))
            if (exitflag == daqp_exit_optimal) then
                do j = 1, work%n_active ! all active constraints are also nonredundant
                    is_redundant(work%WS(j)) = 0
                end do
            end if
        end if
        call daqp_deactivate_constraints(work)
    end do

    end subroutine minrep_work
!*****************************************************************************************

!*****************************************************************************************
!>
!  The first constraint that `x` violates by more than `tol` (upstream's
!  `daqp_first_violating`), or `m+1` if none.

    integer(ip) function daqp_first_violating(x,A,bu,bl,ms,tol) result(ind)

    real(wp), intent(in) :: x(:)    !! point `(n)`
    real(wp), intent(in) :: A(:,:)  !! general constraints `(m-ms,n)`
    real(wp), intent(in) :: bu(:)   !! upper bounds `(m)`
    real(wp), intent(in) :: bl(:)   !! lower bounds `(m)`
    integer(ip), intent(in) :: ms   !! number of simple bounds
    real(wp), intent(in) :: tol     !! tolerance

    integer(ip) :: i, j, m, n
    real(wp) :: ax

    m = int(size(bu), ip)
    n = int(size(x), ip)
    do i = 1, ms
        if (x(i) > bu(i)+tol .or. x(i) < bl(i)-tol) then
            ind = i
            return
        end if
    end do
    do i = ms+1, m
        ax = 0.0_wp
        do j = 1, n
            ax = ax + A(i-ms,j)*x(j)
        end do
        if (ax > bu(i)+tol .or. ax < bl(i)-tol) then
            ind = i
            return
        end if
    end do
    ind = m + 1 ! no constraint is violated

    end function daqp_first_violating
!*****************************************************************************************

!*****************************************************************************************
!>
!  Allocate the data of an AVI (upstream's `allocate_daqp_avi`).

    subroutine allocate_avi(avi,n)

    type(daqp_avi_data), intent(inout) :: avi !! AVI
    integer(ip), intent(in) :: n              !! number of variables

    avi%is_symmetric = .false.
    avi%retry_rho_needed = .false.
    avi%rho = 0.0_wp
    allocate(avi%Hsym(n,n), avi%Hs_rho(n,n), avi%H_rho(n,n), avi%LU_H(n,n), &
             avi%P_H2(n), avi%P_H(n), avi%P_S(n), avi%kkt_buffer(n*n+2*n), &
             avi%Hx(n), avi%x(n), avi%y(n), avi%xtemp(n))

    end subroutine allocate_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  LU factorization with partial pivoting, in place, of the `n x n` row-major
!  matrix `A` (`A(j*n+i+1)` is element `(j,i)`, 0-based). Returns 0, or -1 if
!  a pivot is below 1e-12 (upstream's `daqp_lu`).

    integer(ip) function daqp_lu(A,P,n) result(flag)

    real(wp), intent(inout) :: A(*)      !! matrix; factors on output
    integer(ip), intent(inout) :: P(*)   !! permutation (1-based)
    integer(ip), intent(in) :: n         !! dimension

    integer(ip) :: i, j, k, pivot, tmp_p
    real(wp) :: max_val, pA, tmp

    do i = 1, n
        P(i) = i
    end do
    flag = 0
    do i = 0, n-1
        ! pivot
        max_val = 0.0_wp
        pivot = i
        do j = i, n-1
            pA = A(j*n+i+1)
            if (pA < 0.0_wp) pA = -pA
            if (pA > max_val) then
                max_val = pA
                pivot = j
            end if
        end do
        ! check for singularity
        if (max_val < 1.0e-12_wp) then
            flag = -1
            return
        end if
        ! swap rows
        do k = 0, n-1
            tmp = A(i*n+k+1)
            A(i*n+k+1) = A(pivot*n+k+1)
            A(pivot*n+k+1) = tmp
        end do
        tmp_p = P(i+1)
        P(i+1) = P(pivot+1)
        P(pivot+1) = tmp_p
        ! elimination
        do j = i+1, n-1
            A(j*n+i+1) = A(j*n+i+1)/A(i*n+i+1)
            do k = i+1, n-1
                A(j*n+k+1) = A(j*n+k+1) - A(j*n+i+1)*A(i*n+k+1)
            end do
        end do
    end do

    end function daqp_lu
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve `A x = b` with the factors of [[daqp_lu]] (upstream's `daqp_lu_solve`).

    subroutine daqp_lu_solve(LU,P,b,x,n)

    real(wp), intent(in) :: LU(*)     !! factors
    integer(ip), intent(in) :: P(*)   !! permutation
    real(wp), intent(in) :: b(*)      !! right-hand side
    real(wp), intent(inout) :: x(*)   !! solution
    integer(ip), intent(in) :: n      !! dimension

    integer(ip) :: i, j

    ! solve Ly = Pb
    do i = 0, n-1
        x(i+1) = b(P(i+1))
        do j = 0, i-1
            x(i+1) = x(i+1) - LU(i*n+j+1)*x(j+1)
        end do
    end do
    ! solve Ux = y
    do i = n-1, 0, -1
        do j = i+1, n-1
            x(i+1) = x(i+1) - LU(i*n+j+1)*x(j+1)
        end do
        x(i+1) = x(i+1)/LU(i*n+i+1)
    end do

    end subroutine daqp_lu_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set up the matrices of an AVI (upstream's `daqp_update_avi`): the symmetric
!  part, its shift by `rho` (Douglas-Rachford), and the LU factors of `H`.

    subroutine update_avi(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, n, lu_status
    real(wp) :: val, min_diag, max_row_sum, fro_norm_sq, max_asymmetry, row_sum, &
                asymmetry, hessian_scale, min_lu_pivot, pivot

    n = work%qp%n
    associate (avi => work%avi, H => work%qp%Hc)
    min_diag = daqp_inf
    max_row_sum = 0.0_wp
    fro_norm_sq = 0.0_wp
    max_asymmetry = 0.0_wp
    avi%rho = 0.0_wp
    avi%retry_rho_needed = .false.
    do i = 1, n
        row_sum = 0.0_wp
        do j = 1, n
            if (j > i) then
                asymmetry = abs(H(j,i) - H(i,j))
                if (asymmetry > max_asymmetry) max_asymmetry = asymmetry
            end if
            val = (H(j,i) + H(i,j))*0.5_wp
            avi%Hsym(j,i) = val
            avi%Hs_rho(j,i) = val
            avi%H_rho(j,i) = H(j,i)
            avi%LU_H(j,i) = H(j,i)
            if (val < 0.0_wp) then
                row_sum = row_sum - val
            else
                row_sum = row_sum + val
            end if
            fro_norm_sq = fro_norm_sq + H(j,i)*H(j,i)
            if (i == j .and. val < min_diag) min_diag = val
        end do
        if (row_sum > max_row_sum) max_row_sum = row_sum
    end do
    hessian_scale = sqrt(fro_norm_sq)
    if (hessian_scale < 1.0_wp) hessian_scale = 1.0_wp
    avi%is_symmetric = max_asymmetry <= work%settings%zero_tol*hessian_scale
    if (avi%is_symmetric) return

    ! detect a possibly problematic rho from the LU pivots
    lu_status = daqp_lu(avi%LU_H, avi%P_H, n)
    if (lu_status == 0 .and. min_diag > 0.0_wp) then
        min_lu_pivot = daqp_inf
        do i = 1, n
            pivot = abs(avi%LU_H(i,i))
            if (pivot < min_lu_pivot) min_lu_pivot = pivot
        end do
        if (min_lu_pivot < avi_pivot_trigger*min_diag) avi%retry_rho_needed = .true.
    end if

    ! start with the default step length heuristic
    if (min_diag > 0.0_wp .and. max_row_sum > 0.0_wp) then
        avi%rho = sqrt(min_diag*max_row_sum)
    else
        avi%rho = sqrt(fro_norm_sq)/2.0_wp
    end if
    do i = 1, n
        avi%Hs_rho(i,i) = avi%Hs_rho(i,i) + avi%rho
        avi%H_rho(i,i) = avi%H_rho(i,i) + avi%rho
    end do
    ! (the factorization of H_rho is deferred until needed)
    end associate

    end subroutine update_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Retry an AVI with a reduced `rho` (upstream's
!  `daqp_retry_avi_with_reduced_rho`). Returns 0 if no retry is needed, 1, or
!  a negative exit flag.

    integer(ip) function retry_avi_with_reduced_rho(work) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, n

    flag = 0
    if (.not. allocated(work%avi)) return
    if (.not. work%avi%retry_rho_needed) return
    n = work%n
    associate (avi => work%avi)
    avi%retry_rho_needed = .false. ! at most one retry per setup
    avi%rho = avi%rho/avi_retry_rho_reduction
    avi%Hs_rho = avi%Hsym
    avi%H_rho = work%qp%Hc
    do i = 1, n
        avi%Hs_rho(i,i) = avi%Hs_rho(i,i) + avi%rho
        avi%H_rho(i,i) = avi%H_rho(i,i) + avi%rho
    end do
    i = daqp_lu(avi%H_rho, avi%P_H2, n)
    end associate
    flag = update_R(work, .false., H=work%avi%Hs_rho)
    if (flag < 0) return
    call update_v(work)
    flag = update_M(work)
    if (flag < 0) return
    call normalize_Rinv(work)
    call update_d(work)
    flag = 1

    end function retry_avi_with_reduced_rho
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve a nonsymmetric AVI (upstream's `daqp_solve_avi`): Douglas-Rachford
!  splitting, with each step an LDP of the symmetric part, and Newton steps
!  (a KKT solve on the working set) when the working set settles.

    recursive integer(ip) function daqp_solve_avi(work) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: n, i, j, k, tot_iter, counter, terminate_limit, original_limit, &
                   previous_outer_iterations, retry_flag
    real(wp) :: val, s, s2, minimum_newton_residual
    logical :: retry_requested

    n = work%n
    exitflag = -10
    tot_iter = 0
    counter = 0
    terminate_limit = 5
    retry_requested = .false.
    minimum_newton_residual = daqp_inf

    work%nh = 0 ! counts the outer iterations
    associate (avi => work%avi)
    avi%x(1:n) = work%x(1:n) ! initial iterate

    k = 0
    do while (k < work%settings%iter_limit)
        work%nh = work%nh + 1
        ! xtemp = H*x + f - (Hsym + rho I)x
        do i = 1, n
            s = 0.0_wp
            s2 = 0.0_wp
            do j = 1, n
                s = s + work%qp%Hc(j,i)*avi%x(j)
                s2 = s2 + avi%Hs_rho(j,i)*avi%x(j)
            end do
            avi%Hx(i) = s
            avi%xtemp(i) = s + work%qp%f(i) - s2
        end do

        ! update the linear term
        if (work%has_v) then
            work%v(1:n) = avi%xtemp(1:n)
            call transform_v(work, work%v)
        end if
        call update_d(work)

        exitflag = daqp_ldp(work)

        if (exitflag < 0) exit
        call daqp_ldp2qp_solution(work)
        tot_iter = tot_iter + work%iterations

        if (counter == terminate_limit) then ! check if the Newton step made progress
            s = 0.0_wp
            do i = 1, n
                val = avi%x(i) - work%x(i)
                s = s + val*val
            end do
            ! no decrease since the last Newton iterate -> revert the Newton step
            if (s > minimum_newton_residual) then
                avi%x(1:n) = work%xold(1:n)
                if (terminate_limit == 30 .and. avi%retry_rho_needed) then
                    retry_requested = .true.
                    exit
                end if
                terminate_limit = terminate_limit + 5 ! give DR more time to converge
                if (terminate_limit > 30) terminate_limit = 30
            else
                minimum_newton_residual = s
                avi%y(1:n) = work%x(1:n)
            end if
        else ! update the y iterate
            avi%y(1:n) = work%x(1:n)
        end if

        ! the working set has not changed -> check the KKT conditions
        if (work%iterations == 1) then
            counter = counter + 1
            if (counter == terminate_limit) then
                work%xold(1:n) = avi%x(1:n) ! in case the Newton step fails
                call daqp_solve_avi_kkt(work) ! find a KKT point
                if (daqp_check_optimal_avi(work)) then
                    work%x(1:n) = avi%x(1:n)
                    exitflag = 1
                    exit
                end if
                k = k + 1
                cycle
            end if
        else
            counter = 0
        end if

        do i = 1, n
            avi%xtemp(i) = avi%rho*avi%y(i) + avi%Hx(i)
            avi%y(i) = avi%y(i) - avi%x(i)
        end do
        do i = 1, n
            avi%xtemp(i) = avi%xtemp(i) + 0.5_wp*avi%Hsym(i,i)*avi%y(i) ! diagonal
            do j = i+1, n
                val = 0.5_wp*avi%Hsym(j,i)
                avi%xtemp(i) = avi%xtemp(i) + val*avi%y(j)
                avi%xtemp(j) = avi%xtemp(j) + val*avi%y(i)
            end do
        end do
        call daqp_lu_solve(avi%H_rho, avi%P_H2, avi%xtemp, avi%x, n)
        k = k + 1
    end do
    end associate

    if (retry_requested) then
        original_limit = work%settings%iter_limit
        previous_outer_iterations = work%nh
        retry_flag = retry_avi_with_reduced_rho(work)
        if (retry_flag < 0) then
            exitflag = retry_flag
            return
        end if
        if (retry_flag > 0 .and. k+1 < original_limit) then
            work%settings%iter_limit = original_limit - (k+1)
            retry_flag = daqp_solve_avi(work)
            work%settings%iter_limit = original_limit
            work%iterations = work%iterations + tot_iter
            work%nh = work%nh + previous_outer_iterations
            exitflag = retry_flag
            return
        end if
        work%iterations = tot_iter
        exitflag = daqp_exit_iterlimit
        return
    end if
    if (k == work%settings%iter_limit) exitflag = daqp_exit_iterlimit
    work%iterations = tot_iter

    end function daqp_solve_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  The KKT point of an AVI on the working set (upstream's `daqp_solve_avi_kkt`):
!  `S lam = -A_W H^{-1} f - b_W` with `S = A_W H^{-1} A_W'`, then
!  `H x = -f - A_W' lam`.

    subroutine daqp_solve_avi_kkt(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, k, nAS, n, row_idx, roff, flag
    real(wp) :: s, lj

    nAS = work%n_active
    n = work%n
    roff = nAS*nAS ! the right-hand side follows S in kkt_buffer
    associate (avi => work%avi, qp => work%qp, temp => work%avi%xtemp)

    ! S = A_WS * H^-1 * A_WS^T
    do i = 0, nAS-1
        ! temp = H^-1 * A_row_WS(i)^T
        row_idx = work%WS(i+1)
        if (row_idx <= work%ms) then ! simple bound
            avi%kkt_buffer(roff+1:roff+n) = 0.0_wp
            avi%kkt_buffer(roff+row_idx) = 1.0_wp
            call daqp_lu_solve(avi%LU_H, avi%P_H, avi%kkt_buffer(roff+1:roff+n), temp, n)
        else
            call daqp_lu_solve(avi%LU_H, avi%P_H, qp%At(:,row_idx-work%ms), temp, n)
        end if
        do j = 0, nAS-1
            row_idx = work%WS(j+1)
            if (row_idx <= work%ms) then ! simple bound
                s = temp(row_idx)
            else ! general constraint
                s = 0.0_wp
                do k = 1, n
                    s = s + qp%At(k,row_idx-work%ms)*temp(k)
                end do
            end if
            avi%kkt_buffer(j*nAS+i+1) = s
        end do
    end do

    ! the right-hand side: -A_WS * H^-1 * f - b_WS
    call daqp_lu_solve(avi%LU_H, avi%P_H, qp%f, temp, n)
    do i = 0, nAS-1
        row_idx = work%WS(i+1)
        if (has(work%sense(row_idx),daqp_lower)) then
            s = qp%blower(row_idx)
        else
            s = qp%bupper(row_idx)
        end if
        if (row_idx <= work%ms) then
            s = s + temp(row_idx)
        else
            do k = 1, n
                s = s + qp%At(k,row_idx-work%ms)*temp(k)
            end do
        end if
        avi%kkt_buffer(roff+i+1) = -s
        ! soft constraints -> the diagonal of S is regularized
        if (has(work%sense(row_idx),daqp_soft)) then
          avi%kkt_buffer(i*(nAS+1)+1) = avi%kkt_buffer(i*(nAS+1)+1) + &
                work%settings%rho_soft/(work%scaling(row_idx)*work%scaling(row_idx))
        end if
    end do

    ! lambda: S * lambda = rhs
    flag = daqp_lu(avi%kkt_buffer, avi%P_S, nAS)
    call daqp_lu_solve(avi%kkt_buffer, avi%P_S, avi%kkt_buffer(roff+1:roff+nAS), work%lam_star, nAS)

    ! x: H * x = -f - A_WS^T * lambda
    do i = 1, n
        temp(i) = -qp%f(i)
    end do
    do j = 1, nAS
        lj = work%lam_star(j)
        row_idx = work%WS(j)
        if (row_idx <= work%ms) then
            temp(row_idx) = temp(row_idx) - lj
        else
            do i = 1, n
                temp(i) = temp(i) - qp%At(i,row_idx-work%ms)*lj
            end do
        end if
    end do
    call daqp_lu_solve(avi%LU_H, avi%P_H, temp, avi%x, n)
    end associate

    end subroutine daqp_solve_avi_kkt
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the KKT point of an AVI is optimal (upstream's `daqp_check_optimal_avi`).

    logical function daqp_check_optimal_avi(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, j
    real(wp) :: dual_tol, primal_tol, ax

    daqp_check_optimal_avi = .false.
    dual_tol = work%settings%dual_tol
    primal_tol = work%settings%primal_tol
    ! the dual variables
    do i = 1, work%n_active
        if (has(work%sense(work%WS(i)),daqp_immutable)) cycle
        if (has(work%sense(work%WS(i)),daqp_lower)) then
            if (work%lam_star(i) > dual_tol) return
        else
            if (work%lam_star(i) < -dual_tol) return
        end if
    end do
    ! simple constraints
    do i = 1, work%ms
        if (has(work%sense(i),daqp_active)) cycle
        if (work%avi%x(i) > work%qp%bupper(i) + primal_tol) return
        if (work%avi%x(i) < work%qp%blower(i) - primal_tol) return
    end do
    ! general constraints
    do i = work%ms+1, work%m
        if (has(work%sense(i),daqp_active)) cycle
        ax = 0.0_wp
        do j = 1, work%n
            ax = ax + work%qp%At(j,i-work%ms)*work%avi%x(j)
        end do
        if (ax > work%qp%bupper(i) + primal_tol) return
        if (ax < work%qp%blower(i) - primal_tol) return
    end do
    daqp_check_optimal_avi = .true. ! an optimal KKT point is found

    end function daqp_check_optimal_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Branch-and-bound constraint id helpers: the bit `bnb_lower_bit` marks the
!  lower bound.

    pure integer(ip) function add_lower_flag(x)
    integer(ip), intent(in) :: x !! constraint id
    add_lower_flag = ibset(x, bnb_lower_bit)
    end function add_lower_flag

    pure integer(ip) function remove_lower_flag(x)
    integer(ip), intent(in) :: x !! constraint id
    remove_lower_flag = ibclr(x, bnb_lower_bit)
    end function remove_lower_flag

    pure integer(ip) function toggle_lower_flag(x)
    integer(ip), intent(in) :: x !! constraint id
    toggle_lower_flag = ieor(x, ishft(1_ip, bnb_lower_bit))
    end function toggle_lower_flag

    pure logical function extract_lower_flag(x)
    integer(ip), intent(in) :: x !! constraint id
    extract_lower_flag = btest(x, bnb_lower_bit)
    end function extract_lower_flag
!*****************************************************************************************

!*****************************************************************************************
!>
!  Signed distance of binary constraint `id` from the midpoint of its bounds.

    real(wp) function binary_diff(work,id) result(diff)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id            !! constraint

    integer(ip) :: j, disp

    diff = 0.5_wp*(work%dupper(id)+work%dlower(id))
    if (id <= work%ms) then ! simple bound
        if (work%rmode /= rinv_dense) then ! Hessian is identity (or diagonal)
            diff = diff - work%x(id)
        else
            disp = ridx(id,id,work%n)
            do j = id, work%n
                diff = diff - work%R(disp)*work%x(j)
                disp = disp + 1
            end do
        end if
    else ! general bound (add_infeasible already computed M*u)
        diff = diff - work%Mu(id-work%ms)
    end if

    end function binary_diff
!*****************************************************************************************

!*****************************************************************************************
!>
!  Store the free part of the working set in `ids`. Returns the number stored.

    integer(ip) function bnb_store_ws(work,ids) result(n_ids)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(inout) :: ids(:)     !! storage

    integer(ip) :: i, id

    n_ids = 0
    do i = work%bnb%neq+1, work%n_active
        id = work%WS(i)
        if (iand(work%sense(id), daqp_immutable+daqp_binary) /= daqp_immutable+daqp_binary) then
            n_ids = n_ids + 1
            if (has(work%sense(id),daqp_lower)) then
                ids(n_ids) = add_lower_flag(id)
            else
                ids(n_ids) = id
            end if
        end if
    end do

    end function bnb_store_ws
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add the constraints in `ids` to the working set (aborted if the basis gets singular).

    subroutine bnb_load_ws(work,ids,n_ids)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: ids(:)           !! constraints (with the lower-bound bit)
    integer(ip), intent(in) :: n_ids            !! number of constraints

    integer(ip) :: i, id

    do i = 1, n_ids
        call add_upper_lower(work, ids(i))
        if (work%sing_ind /= empty_ind) then
            id = work%WS(work%n_active)
            work%n_active = work%n_active - 1
            work%sense(id) = iand(work%sense(id), not(daqp_active))
            work%sing_ind = empty_ind
            exit
        end if
    end do

    end subroutine bnb_load_ws
!*****************************************************************************************

!*****************************************************************************************
!>
!  The immutable constraints (e.g., equalities) are kept fixed as a prefix of
!  the working set throughout the tree. Returns the length of that prefix.

    integer(ip) function bnb_setup_root(work) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, nfixed

    nfixed = work%n_active
    do i = 1, work%n_active
        if (.not. has(work%sense(work%WS(i)),daqp_immutable)) then
            nfixed = i - 1
            exit
        end if
    end do
    do j = nfixed+1, work%n_active
        if (has(work%sense(work%WS(j)),daqp_immutable)) exit
    end do
    if (j > work%n_active) then ! the mutable constraints are a warm start
        flag = nfixed
        return
    end if

    ! immutable after mutable => only activate the immutable constraints
    do i = 1, work%n_active
        if (.not. has(work%sense(work%WS(i)),daqp_immutable)) then
          work%sense(work%WS(i)) = iand(work%sense(work%WS(i)), not(daqp_active))
        end if
    end do
    call daqp_reset_workspace(work)
    flag = daqp_activate_constraints(work)
    if (flag >= 0) flag = work%n_active

    end function bnb_setup_root
!*****************************************************************************************

!*****************************************************************************************
!>
!  Use the candidate in `x` as an incumbent. Returns its objective (internal
!  scale), with `u = R*x+v` stored in `xold`, or -1 if it is infeasible.

    real(wp) function bnb_incumbent(work) result(fval)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: i, j, n, disp
    real(wp) :: val, tol

    fval = -1.0_wp
    if (.not. work%has_qp) return
    n = work%n
    tol = work%settings%primal_tol

    ! check feasibility (soft constraints are treated as hard)
    do i = 1, work%m
        if (i <= work%ms) then
            val = work%x(i)
        else
            val = dot_seq(n, work%qp%At(:,i-work%ms), work%x)
        end if
        if (has(work%sense(i),daqp_immutable) .and. .not. has(work%sense(i),daqp_active) .and. &
            .not. has(work%sense(i),daqp_binary)) cycle ! ignored
        if (val > work%qp%bupper(i)+tol .or. val < work%qp%blower(i)-tol) return
        if (has(work%sense(i),daqp_binary) .and. val > work%qp%blower(i)+tol .and. &
            val < work%qp%bupper(i)-tol) return
    end do

    ! invert daqp_ldp2qp_solution: u = R*x + v
    work%xold(1:n) = work%x(1:n)
    if (work%rmode == rinv_dense) then
        do i = 1, work%ms
            work%xold(i) = work%xold(i)*work%scaling(i)
        end do
        do i = n, 1, -1 ! back substitution with the upper triangular Rinv
            disp = ridx(i,i,n)
            do j = i+1, n
                work%xold(i) = work%xold(i) - work%R(disp+j-i)*work%xold(j)
            end do
            work%xold(i) = work%xold(i)/work%R(disp)
        end do
    else if (work%rmode == rinv_diag) then
        do i = 1, n
            work%xold(i) = work%xold(i)/work%R(i)
        end do
    end if
    if (work%has_v) then
        do i = 1, n
            work%xold(i) = work%xold(i) + work%v(i)
        end do
    end if
    fval = 0.0_wp
    do i = 1, n
        fval = fval + work%xold(i)*work%xold(i)
    end do

    end function bnb_incumbent
!*****************************************************************************************

!*****************************************************************************************
!>
!  Branch and bound over the binary constraints (upstream's `daqp_bnb`).
!  Returns the exit flag; `x` (`u`) holds the best solution.

    integer(ip) function daqp_bnb(work) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: branch_id, node
    real(wp) :: fval_bound0, eps_r, fval_inc
    logical :: have_sol

    exitflag = bnb_setup_root(work)
    if (exitflag < 0) return
    work%bnb%neq = exitflag

    ! warm start the root with the root working set of the previous solve
    ! (unless a warm start has been provided)
    if (work%n_active == work%bnb%neq) then
      call bnb_load_ws(work, work%bnb%root_ws, work%bnb%n_root_ws)
    end if

    ! modify the upper bound based on the absolute/relative suboptimality tolerance
    fval_bound0 = work%settings%fval_bound
    eps_r = 1.0_wp/(1.0_wp+work%settings%rel_subopt)
    work%settings%fval_bound = (fval_bound0 - work%settings%abs_subopt)*eps_r
    have_sol = .false.

    ! start from a user-provided integer-feasible solution
    if (has(work%state,state_incumbent)) then
        work%state = iand(work%state, not(state_incumbent))
        fval_inc = 0.5_wp*bnb_incumbent(work)
        if (fval_inc >= 0.0_wp .and. fval_inc < fval_bound0) then
            work%settings%fval_bound = (fval_inc - work%settings%abs_subopt)*eps_r
            have_sol = .true. ! a feasible solution is stored in xold
        end if
    end if

    associate (bnb => work%bnb)
    bnb%itercount = 0
    bnb%nodecount = 0
    ! the root node
    bnb%tree(1) = daqp_node(bin_id=0, depth=-1, ws_start=0, ws_end=0)
    bnb%n_nodes = 1
    bnb%n_clean = bnb%neq
    bnb%nws = 0

    exitflag = daqp_exit_infeasible
    ! tree exploration
    do while (bnb%n_nodes > 0)
        bnb%n_nodes = bnb%n_nodes - 1
        node = bnb%n_nodes + 1
        exitflag = process_node(work, node) ! solve the relaxation
        if (bnb%tree(node)%depth < 0 .and. exitflag > 0) then
            bnb%n_root_ws = bnb_store_ws(work, bnb%root_ws)
        end if
        ! individual relaxations are often too short to reach the timer check in
        ! daqp_ldp, so also enforce the limit across the tree
        if (work%timer_on .and. iand(bnb%nodecount, 31_ip) == 0) then
            if (elapsed_time(work) > work%settings%time_limit) then
                exitflag = daqp_exit_timelimit
                exit
            end if
        end if
        ! cut conditions
        if (exitflag == daqp_exit_infeasible) cycle ! dominance cut
        if (exitflag < 0) exit ! the inner solver failed

        ! find an index to branch over
        branch_id = get_branch_id(work)
        if (branch_id == empty_ind) then ! nothing to branch over => integer feasible
            work%settings%fval_bound = (0.5_wp*work%fval - work%settings%abs_subopt)*eps_r
            call swap_x(work) ! store the feasible solution
            have_sol = .true.
        else
            call spawn_children(work, node, branch_id)
        end if
    end do

    ! exploration completed
    work%iterations = bnb%itercount
    ! restore the root state (unfix the binaries etc.) so that the workspace can
    ! be reused for subsequent solves
    call node_cleanup_workspace(work, bnb%neq)
    bnb%n_clean = bnb%neq
    end associate
    if (.not. have_sol) then
        work%settings%fval_bound = fval_bound0
        if (exitflag >= 0) exitflag = daqp_exit_infeasible
    else
        ! invert fval_bound = (0.5*fval_best - abs_subopt)*eps_r to recover fval_best
        work%fval = 2.0_wp*work%settings%fval_bound/eps_r + 2.0_wp*work%settings%abs_subopt
        work%settings%fval_bound = fval_bound0
        call swap_x(work) ! x (u) holds the best feasible solution
        if (exitflag >= daqp_exit_infeasible) exitflag = daqp_exit_optimal
    end if

    end function daqp_bnb
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve the relaxation of a node (upstream's `daqp_process_node`).

    integer(ip) function process_node(work,node) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: node             !! node (index in the tree)

    integer(ip) :: depth
    logical :: cleanup

    work%bnb%nodecount = work%bnb%nodecount + 1
    depth = work%bnb%tree(node)%depth
    if (depth >= 0) then
        ! fix a binary constraint
        work%bnb%fixed_ids(depth+1) = work%bnb%tree(node)%bin_id
        ! set up the relaxation
        cleanup = work%bnb%n_nodes == 0
        if (.not. cleanup) cleanup = work%bnb%tree(node-1)%depth /= depth
        if (cleanup) then
            ! the sibling has been processed => fix the workspace state
            work%bnb%n_clean = work%bnb%n_clean + (depth - work%bnb%tree(node+1)%depth)
            call node_cleanup_workspace(work, work%bnb%n_clean)
            call warmstart_node(work, node)
        else
            call add_upper_lower(work, work%bnb%tree(node)%bin_id)
            work%sense(remove_lower_flag(work%bnb%tree(node)%bin_id)) = &
                ior(work%sense(remove_lower_flag(work%bnb%tree(node)%bin_id)), daqp_immutable) ! equality
            if (work%sing_ind /= empty_ind) call setup_cold_bnb(work, node) ! cold start, not to miss integer feasible
        end if
    end if
    ! solve the relaxation
    exitflag = daqp_ldp(work)
    work%bnb%itercount = work%bnb%itercount + work%iterations

    if (exitflag == daqp_exit_cycle) then ! try to repair (cold start)
        ! a cycle can be caused by stale cached forward-substitution data
        work%reuse_ind = 0
        call setup_cold_bnb(work, node)
        exitflag = daqp_ldp(work)
        work%bnb%itercount = work%bnb%itercount + work%iterations
    end if

    end function process_node
!*****************************************************************************************

!*****************************************************************************************
!>
!  The binary constraint to branch over (with the lower-bound bit for the
!  lower endpoint), or `empty_ind` if the relaxation is integer feasible.

    integer(ip) function get_branch_id(work) result(branch)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, id
    real(wp) :: diff, dist, tol, ad

    branch = empty_ind
    do i = 1, work%bnb%nb
        id = work%bnb%bin_ids(i)
        if (has(work%sense(id),daqp_active)) cycle ! skip fixed binary constraints
        ! signed distance from the midpoint between the bounds
        diff = binary_diff(work, id)
        ! a zero-dual binary constraint can lie at an endpoint without being
        ! active: it is already integer feasible
        ad = diff
        if (diff < 0.0_wp) ad = -diff
        dist = 0.5_wp*(work%dupper(id)-work%dlower(id)) - ad
        tol = work%settings%primal_tol*work%scaling(id)
        if (dist <= tol) cycle
        ! explore the endpoint nearest to the relaxation first
        if (diff < 0.0_wp) then
            branch = id
        else
            branch = add_lower_flag(id)
        end if
        return
    end do

    end function get_branch_id
!*****************************************************************************************

!*****************************************************************************************
!>
!  Replace a node by its two children.

    subroutine spawn_children(work,node,branch_id)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: node             !! node
    integer(ip), intent(in) :: branch_id        !! constraint to branch over

    call save_warmstart(work, node)
    associate (tree => work%bnb%tree)
    ! child 1 (reuses the current node)
    tree(node)%bin_id = toggle_lower_flag(branch_id)
    tree(node)%depth = tree(node)%depth + 1
    ! child 2
    tree(node+1)%bin_id = branch_id
    tree(node+1)%depth = tree(node)%depth
    tree(node+1)%ws_start = tree(node)%ws_start
    tree(node+1)%ws_end = tree(node)%ws_end
    end associate
    work%bnb%n_nodes = work%bnb%n_nodes + 2

    end subroutine spawn_children
!*****************************************************************************************

!*****************************************************************************************
!>
!  Restore the working set to its first `n_clean` constraints.

    subroutine node_cleanup_workspace(work,n_clean)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: n_clean          !! constraints to keep

    integer(ip) :: i, id

    do i = n_clean+1, work%n_active
        id = work%WS(i)
        if (has(work%sense(id),daqp_binary)) then
            work%sense(id) = iand(work%sense(id), not(daqp_active+daqp_immutable))
        else
            work%sense(id) = iand(work%sense(id), not(daqp_active))
        end if
    end do
    work%sing_ind = empty_ind
    work%n_active = n_clean
    ! only the retained prefix can still have a valid cached substitution
    if (work%reuse_ind > n_clean) work%reuse_ind = n_clean

    end subroutine node_cleanup_workspace
!*****************************************************************************************

!*****************************************************************************************
!>
!  Warm start a node: its fixed constraints, then its stored working set.

    subroutine warmstart_node(work,node)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: node             !! node

    integer(ip) :: i, ws_start, ws_end

    ! add the fixed constraints
    do i = work%bnb%n_clean - work%bnb%neq, work%bnb%tree(node)%depth
        call add_upper_lower(work, work%bnb%fixed_ids(i+1))
        work%sense(remove_lower_flag(work%bnb%fixed_ids(i+1))) = &
            ior(work%sense(remove_lower_flag(work%bnb%fixed_ids(i+1))), daqp_immutable)
    end do
    work%bnb%n_clean = work%bnb%neq + work%bnb%tree(node)%depth
    ! add the free constraints
    ws_start = work%bnb%tree(node)%ws_start
    ws_end = work%bnb%tree(node)%ws_end
    call bnb_load_ws(work, work%bnb%tree_ws(ws_start+1:), ws_end-ws_start)
    work%bnb%nws = ws_start ! always move up the tree after a warm start

    end subroutine warmstart_node
!*****************************************************************************************

!*****************************************************************************************
!>
!  Store the working set of a node, as the warm start of its children.

    subroutine save_warmstart(work,node)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: node             !! node

    integer(ip) :: nstored

    work%bnb%tree(node)%ws_start = work%bnb%nws
    nstored = bnb_store_ws(work, work%bnb%tree_ws(work%bnb%nws+1:))
    work%bnb%nws = work%bnb%nws + nstored
    work%bnb%tree(node)%ws_end = work%bnb%nws

    end subroutine save_warmstart
!*****************************************************************************************

!*****************************************************************************************
!>
!  Add a constraint at its upper bound, or at its lower one if `add_id` has
!  the lower-bound bit.

    subroutine add_upper_lower(work,add_id)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: add_id           !! constraint (with the lower-bound bit)

    integer(ip) :: id

    id = remove_lower_flag(add_id)
    if (extract_lower_flag(add_id)) then
        work%sense(id) = ior(work%sense(id), daqp_lower)
        call add_constraint(work, id, -1.0_wp)
    else
        work%sense(id) = iand(work%sense(id), not(daqp_lower))
        call add_constraint(work, id, 1.0_wp)
    end if

    end subroutine add_upper_lower
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cold start a node: only its fixed constraints.

    subroutine setup_cold_bnb(work,node)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: node             !! node

    integer(ip) :: i

    call node_cleanup_workspace(work, work%bnb%n_clean)
    do i = work%bnb%n_clean - work%bnb%neq, work%bnb%tree(node)%depth
        call add_upper_lower(work, work%bnb%fixed_ids(i+1))
        work%sense(remove_lower_flag(work%bnb%fixed_ids(i+1))) = &
            ior(work%sense(remove_lower_flag(work%bnb%fixed_ids(i+1))), daqp_immutable)
    end do
    work%bnb%n_clean = work%bnb%neq + work%bnb%tree(node)%depth

    end subroutine setup_cold_bnb
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solve a hierarchical (lexicographic) QP (upstream's `daqp_hiqp`): the first
!  level is hard; each following level is soft, and is made hard at the slacks
!  it ends up with. `lam` (if present) is set to the slacks of the soft levels.

    integer(ip) function daqp_hiqp(work,lam) result(exitflag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(inout), optional :: lam(:) !! slacks of the soft levels

    integer(ip) :: i, j, jj, id, start, iend, iterations, nfree, n_active_old
    real(wp) :: w

    iterations = 0
    exitflag = 0
    ! one hierarchy -> just solve a normal LDP
    if (.not. is_hierarchical(work)) then
        exitflag = daqp_ldp(work)
        return
    end if

    ! a previous solve shifted d by the slacks of the soft levels, so it has to
    ! be reformed by an update first
    if (has(work%state,daqp_update_d) .and. work%has_qp) then
        exitflag = daqp_exit_unsupported
        return
    end if

    if (present(lam)) lam(1:work%m) = 0.0_wp

    ! move down the hierarchy ((0-based) constraints start..end-1 form a level)
    start = work%break_points(1)
    ! a previous solve leaves constraints of the soft levels in the working set;
    ! restart from the (hard) first level, whose active set is kept
    do i = 1, work%n_active
        if (work%WS(i) > start) exit
    end do
    if (i <= work%n_active) then
        work%m = start
        call daqp_reset_workspace(work)
        exitflag = daqp_activate_constraints(work)
        if (exitflag < 0) return
        exitflag = 0
    end if
    nfree = work%n
    do i = 2, work%nh
        ! initialize the current level
        iend = work%break_points(i)
        work%m = iend
        ! soften the constraints and activate
        do j = start+1, iend
            work%sense(j) = ior(work%sense(j), daqp_soft)
            if (has(work%sense(j),daqp_active)) then
                if (has(work%sense(j),daqp_lower)) then
                    call add_constraint(work, j, -1.0_wp)
                else
                    call add_constraint(work, j, 1.0_wp)
                end if
                if (work%sing_ind /= empty_ind) then
                    ! dependent constraint (e.g., from a warm start): leave it out,
                    ! and make it mutable so that it is not ignored
                    work%sense(j) = iand(work%sense(j), not(daqp_active))
                    work%sense(j) = iand(work%sense(j), not(daqp_immutable))
                    work%n_active = work%n_active - 1
                    work%sing_ind = empty_ind
                    if (work%reuse_ind > work%n_active) work%reuse_ind = work%n_active
                end if
            end if
        end do

        ! save the best solution in case daqp_ldp fails
        work%xold(1:work%n) = work%x(1:work%n)
        ! solve the LDP
        exitflag = daqp_ldp(work)
        iterations = iterations + work%iterations
        if (exitflag < 0) exit

        if (iterations >= work%settings%iter_limit) then
            exitflag = daqp_exit_iterlimit
            exit
        end if

        ! perturb the right-hand side with the slacks of the level
        do j = 1, work%n_active
            id = work%WS(j)
            if (has(work%sense(id),daqp_soft)) then
                w = soft_slack(work, j)
                if (w < -work%settings%primal_tol) then
                    work%dlower(id) = work%dlower(id) + w
                else if (w > work%settings%primal_tol) then
                    work%dupper(id) = work%dupper(id) + w
                end if
                if (present(lam)) then
                    if (has(work%sense(id),daqp_lower)) then ! for weakly active
                        w = w - 1.0e-14_wp
                    else
                        w = w + 1.0e-14_wp
                    end if
                    lam(id) = w
                end if
            end if
        end do

        ! make the constraints of the current level hard
        do j = start+1, iend
            work%sense(j) = iand(work%sense(j), not(daqp_soft))
        end do

        if (i == work%nh) exit

        ! find the first active constraint of the current level
        do j = 1, work%n_active
            if (work%WS(j) > start) exit
        end do

        ! reactivate the constraints of the current level (to address soft->hard)
        n_active_old = min(work%n_active, work%n)
        do jj = n_active_old+1, work%n_active
            work%sense(work%WS(jj)) = iand(work%sense(work%WS(jj)), not(daqp_active+daqp_immutable))
        end do
        work%n_active = j - 1
        work%reuse_ind = j - 1
        work%sing_ind = empty_ind
        do jj = j, n_active_old
            call add_constraint(work, work%WS(jj), work%lam_star(jj))
            ! skip if the working set becomes overdetermined
            if (work%sing_ind /= empty_ind) then
                call remove_constraint(work, jj)
                work%sing_ind = empty_ind
                work%sense(work%WS(jj)) = iand(work%sense(work%WS(jj)), not(daqp_immutable))
            else
                if (has(work%sense(work%WS(jj)),daqp_immutable)) nfree = nfree - 1
            end if
        end do

        if (nfree <= 0) exit ! no degrees of freedom left
        ! move up the hierarchy
        start = iend
    end do
    ! finalize
    if (exitflag < 0) then ! restore a point that was good before it failed
        work%x(1:work%n) = work%xold(1:work%n)
        exitflag = daqp_exit_no_freedom
    end if
    work%iterations = iterations ! total number of iterations
    work%state = ior(work%state, daqp_update_d) ! the levels have shifted d

    end function daqp_hiqp
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_core
!*****************************************************************************************
