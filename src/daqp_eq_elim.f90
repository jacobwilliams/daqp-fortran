!*****************************************************************************************
!> author: Jacob Williams
!
!  Elimination of equality constraints, before the QP is turned into an LDP:
!  a translation of `eq_elim.c` of [DAQP](https://github.com/darnstrom/daqp)
!  v0.10.3 (Copyright (c) 2022 Daniel Arnström, MIT licence). Changed from the
!  original: translated to Fortran, 1-based indexing.
!
!  The equality constraints `A_E x = b_E` are eliminated through
!  `x = xp + W w`, where the columns of `W` span the null space of `A_E`
!  (from the QR factorization `A_E' = Q [R; 0]` of the normalized rows) and
!  `xp` is a particular solution. With `x = xp + W w` the remaining constraints
!  become constraints on `w` (the simple bounds turning into general
!  constraints), and the reduced problem is posed in one of three ways
!  (`eq%path`):
!
!  * `eq_path_ldp`: if `Z'HZ = L L'` is positive definite, `W = Z L^{-T}` gives
!    `W'HW = I`, and `xp` is the minimizer over the equality constraints, so
!    that the reduced problem is `min 0.5||w||^2`. For a diagonal (positive)
!    Hessian, the QR is formed in the metric of `H` instead.
!  * `eq_path_qp`: otherwise, `W = Z` and the reduced problem keeps the
!    Hessian `Z'HZ` and the linear term `Z'(H xp + f)`.
!  * `eq_path_lp`: for an LP, `W = Z` and the linear term is `Z'f`.
!
!  The reduced problem is an ordinary problem, whose LDP is formed and solved
!  by the usual routines: it is swapped into the workspace (installed) while
!  it is formed or solved, and the workspace describes the original problem
!  otherwise. Matrices are in upstream's (row-major) memory order: `V(:,k)` is
!  the k-th Householder vector (then `Z`), `W(:,i)` is the i-th row of `W`.

    module daqp_eq_elim

    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use daqp_types

    implicit none

    private

    abstract interface
        integer(ip) function update_ldp_proc(mask, work)
            !! Forms the LDP of the problem in the workspace.
            import :: ip, daqp_workspace
            implicit none
            integer(ip), intent(in) :: mask
            type(daqp_workspace), intent(inout) :: work
        end function update_ldp_proc
    end interface

    public :: update_ldp_proc
    public :: daqp_eq_wanted, daqp_eq_update, daqp_eq_deactivate, daqp_eq_install, daqp_eq_restore
    public :: daqp_eq_set_primal_start, daqp_eq_expand, free_daqp_eq

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y <-- (I - tau v v') y` on the indices `>= k` (`v(k) = 1` implicitly).

    pure subroutine eq_reflect(v, tau, k, n, y)

    real(wp), intent(in) :: v(:)       !! Householder vector
    real(wp), intent(in) :: tau        !! Householder scalar
    integer(ip), intent(in) :: k       !! first index
    integer(ip), intent(in) :: n       !! length
    real(wp), intent(inout) :: y(:)    !! vector

    integer(ip) :: i
    real(wp) :: w

    w = y(k)
    do i = k+1, n
        w = w + v(i)*y(i)
    end do
    w = w*tau
    y(k) = y(k) - w
    do i = k+1, n
        y(i) = y(i) - w*v(i)
    end do

    end subroutine eq_reflect
!*****************************************************************************************

!*****************************************************************************************
!>
!  `Y <-- Q'Y` for the `cnt` columns of `Y` (length `n`), with the reflectors
!  `1..nr` (four vectors per pass over a reflector).

    subroutine eq_apply_QT_many(V, tau, nr, n, Y, cnt)

    real(wp), intent(in) :: V(:,:)     !! Householder vectors (columns)
    real(wp), intent(in) :: tau(:)     !! Householder scalars
    integer(ip), intent(in) :: nr      !! number of reflectors
    integer(ip), intent(in) :: n       !! length of the vectors
    real(wp), intent(inout) :: Y(:,:)  !! vectors (columns)
    integer(ip), intent(in) :: cnt     !! number of vectors

    integer(ip) :: i, j, k
    real(wp) :: tk, vi, w0, w1, w2, w3

    j = 1
    do while (j+3 <= cnt)
        do k = 1, nr
            tk = tau(k)
            w0 = Y(k,j); w1 = Y(k,j+1); w2 = Y(k,j+2); w3 = Y(k,j+3)
            do i = k+1, n
                vi = V(i,k)
                w0 = w0 + vi*Y(i,j); w1 = w1 + vi*Y(i,j+1)
                w2 = w2 + vi*Y(i,j+2); w3 = w3 + vi*Y(i,j+3)
            end do
            w0 = w0*tk; w1 = w1*tk; w2 = w2*tk; w3 = w3*tk
            Y(k,j) = Y(k,j) - w0; Y(k,j+1) = Y(k,j+1) - w1
            Y(k,j+2) = Y(k,j+2) - w2; Y(k,j+3) = Y(k,j+3) - w3
            do i = k+1, n
                vi = V(i,k)
                Y(i,j) = Y(i,j) - w0*vi; Y(i,j+1) = Y(i,j+1) - w1*vi
                Y(i,j+2) = Y(i,j+2) - w2*vi; Y(i,j+3) = Y(i,j+3) - w3*vi
            end do
        end do
        j = j + 4
    end do
    do while (j <= cnt)
        do k = 1, nr
            call eq_reflect(V(:,k), tau(k), k, n, Y(:,j))
        end do
        j = j + 1
    end do

    end subroutine eq_apply_QT_many
!*****************************************************************************************

!*****************************************************************************************
!>
!  `Z = Q(:,nr+1:n)`, accumulated in the columns `nr+1..n` of `V` (four at a time).

    subroutine eq_accumulate_Z(V, tau, nr, n)

    real(wp), intent(inout) :: V(:,:)  !! Householder vectors; Z on output
    real(wp), intent(in) :: tau(:)     !! Householder scalars
    integer(ip), intent(in) :: nr      !! number of reflectors
    integer(ip), intent(in) :: n       !! dimension

    integer(ip) :: i, j, k
    real(wp) :: tk, vi, w0, w1, w2, w3

    do j = nr+1, n
        V(1:n,j) = 0.0_wp
        V(j,j) = 1.0_wp
    end do
    do k = nr, 1, -1
        tk = tau(k)
        j = nr + 1
        do while (j+3 <= n)
            w0 = V(k,j); w1 = V(k,j+1); w2 = V(k,j+2); w3 = V(k,j+3)
            do i = k+1, n
                vi = V(i,k)
                w0 = w0 + vi*V(i,j); w1 = w1 + vi*V(i,j+1)
                w2 = w2 + vi*V(i,j+2); w3 = w3 + vi*V(i,j+3)
            end do
            w0 = w0*tk; w1 = w1*tk; w2 = w2*tk; w3 = w3*tk
            V(k,j) = V(k,j) - w0; V(k,j+1) = V(k,j+1) - w1
            V(k,j+2) = V(k,j+2) - w2; V(k,j+3) = V(k,j+3) - w3
            do i = k+1, n
                vi = V(i,k)
                V(i,j) = V(i,j) - w0*vi; V(i,j+1) = V(i,j+1) - w1*vi
                V(i,j+2) = V(i,j+2) - w2*vi; V(i,j+3) = V(i,j+3) - w3*vi
            end do
            j = j + 4
        end do
        do while (j <= n)
            call reflect_col(k, j)
            j = j + 1
        end do
    end do

    contains

        subroutine reflect_col(kk, jj)
            !! column `jj` of `V` <-- reflector `kk` applied to it
            integer(ip), intent(in) :: kk, jj
            integer(ip) :: ii
            real(wp) :: w
            w = V(kk,jj)
            do ii = kk+1, n
                w = w + V(ii,kk)*V(ii,jj)
            end do
            w = w*tau(kk)
            V(kk,jj) = V(kk,jj) - w
            do ii = kk+1, n
                V(ii,jj) = V(ii,jj) - w*V(ii,kk)
            end do
        end subroutine reflect_col

    end subroutine eq_accumulate_Z
!*****************************************************************************************

!*****************************************************************************************
!>
!  `C(j,i) = X(:,i).Y(:,j)` for vectors of length `n`. With `upper`, only
!  `j >= i` is computed, and then mirrored.

    subroutine eq_gemm_tn(n, p, q, X, Y, C, upper)

    integer(ip), intent(in) :: n       !! length of the vectors
    integer(ip), intent(in) :: p       !! number of vectors in X
    integer(ip), intent(in) :: q       !! number of vectors in Y
    real(wp), intent(in) :: X(:,:)     !! vectors (columns)
    real(wp), intent(in) :: Y(:,:)     !! vectors (columns)
    real(wp), intent(inout) :: C(:,:)  !! products
    logical, intent(in) :: upper       !! only the upper triangle (mirrored)

    integer(ip) :: i, j, k, j0
    real(wp) :: s

    do i = 1, p
        j0 = 1
        if (upper) j0 = i
        do j = j0, q
            s = 0.0_wp
            do k = 1, n
                s = s + X(k,i)*Y(k,j)
            end do
            C(j,i) = s
        end do
    end do
    if (upper) then
        do i = 1, p
            do j = 1, i-1
                C(j,i) = C(i,j)
            end do
        end do
    end if

    end subroutine eq_gemm_tn
!*****************************************************************************************

!*****************************************************************************************
!>
!  `B <-- P_k B P_k` for `k = 1..nr`, for a symmetric `B` (lower triangle,
!  row major: `B(j,i)`, `j <= i`). Only the trailing blocks are updated.

    subroutine eq_twoside(B, n, V, tau, nr, p)

    real(wp), intent(inout) :: B(:,:)  !! symmetric matrix (lower triangle)
    integer(ip), intent(in) :: n       !! dimension
    real(wp), intent(in) :: V(:,:)     !! Householder vectors
    real(wp), intent(in) :: tau(:)     !! Householder scalars
    integer(ip), intent(in) :: nr      !! number of reflectors
    real(wp), intent(inout) :: p(:)    !! scratch (n)

    integer(ip) :: i, j, k
    real(wp) :: tk, pv, s, vi, wi

    do k = 1, nr
        tk = tau(k)
        pv = 0.0_wp
        ! p = tk*B*v (v(k) = 1), a symmetric product with the lower triangle
        do i = k, n
            p(i) = 0.0_wp
        end do
        do i = k, n
            vi = vk(i)
            s = B(k,i)
            do j = k+1, i-1
                s = s + B(j,i)*V(j,k)
                p(j) = p(j) + B(j,i)*vi
            end do
            if (i > k) then
                s = s + B(i,i)*vi
                p(k) = p(k) + B(k,i)*vi
            end if
            p(i) = p(i) + s
        end do
        do i = k, n
            p(i) = p(i)*tk
            pv = pv + p(i)*vk(i)
        end do
        ! w = p - (tk/2)(p'v) v, B <-- B - v w' - w v'
        pv = pv*(0.5_wp*tk)
        do i = k, n
            p(i) = p(i) - pv*vk(i)
        end do
        do i = k, n
            vi = vk(i)
            wi = p(i)
            B(k,i) = B(k,i) - (vi*p(k) + wi)
            do j = k+1, i
                B(j,i) = B(j,i) - (vi*p(j) + wi*V(j,k))
            end do
        end do
    end do

    contains

        pure real(wp) function vk(ii)
            !! element `ii` of reflector `k` (1 at `k`)
            integer(ip), intent(in) :: ii
            if (ii == k) then
                vk = 1.0_wp
            else
                vk = V(ii,k)
            end if
        end function vk

    end subroutine eq_twoside
!*****************************************************************************************

!*****************************************************************************************
!>
!  Cholesky factorization `A = L L'` (lower, row major: `A(j,i)`, `j <= i`, in
!  place), with the rows in blocks of four. Returns false if `A` is not
!  positive definite (a pivot below `zero_tol`, or relative to the largest pivot).

    logical function eq_chol(A, n, zero_tol)

    real(wp), intent(inout) :: A(:,:)  !! matrix; factor on output
    integer(ip), intent(in) :: n       !! dimension
    real(wp), intent(in) :: zero_tol   !! tolerance

    integer(ip) :: i, j, k, ib, be
    real(wp) :: min_pivot, max_pivot, dj, s, s0, s1, s2, s3, l

    ! (0-based loop variables, as upstream; A(j+1,i+1) is upstream's A[i*n+j])
    eq_chol = .false.
    min_pivot = daqp_inf
    max_pivot = 0.0_wp
    do ib = 0, n-1, 4
        be = min(ib+4, n)
        do j = 0, ib-1
            dj = 1.0_wp/A(j+1,j+1)
            if (be-ib == 4) then
                s0 = A(j+1,ib+1); s1 = A(j+1,ib+2); s2 = A(j+1,ib+3); s3 = A(j+1,ib+4)
                do k = 0, j-1
                    l = A(k+1,j+1)
                    s0 = s0 - A(k+1,ib+1)*l; s1 = s1 - A(k+1,ib+2)*l
                    s2 = s2 - A(k+1,ib+3)*l; s3 = s3 - A(k+1,ib+4)*l
                end do
                A(j+1,ib+1) = s0*dj; A(j+1,ib+2) = s1*dj
                A(j+1,ib+3) = s2*dj; A(j+1,ib+4) = s3*dj
            else
                do i = ib, be-1
                    s = A(j+1,i+1)
                    do k = 0, j-1
                        s = s - A(k+1,i+1)*A(k+1,j+1)
                    end do
                    A(j+1,i+1) = s*dj
                end do
            end if
        end do
        do i = ib, be-1
            do j = ib, i
                s = A(j+1,i+1)
                do k = 0, j-1
                    s = s - A(k+1,i+1)*A(k+1,j+1)
                end do
                if (j < i) then
                    A(j+1,i+1) = s/A(j+1,j+1)
                else
                    if (s <= zero_tol) return
                    if (s < min_pivot) min_pivot = s
                    if (s > max_pivot) max_pivot = s
                    A(i+1,i+1) = sqrt(s)
                end if
            end do
        end do
    end do
    eq_chol = min_pivot > zero_tol*max_pivot

    end function eq_chol
!*****************************************************************************************

!*****************************************************************************************
!>
!  Rows `X(:,r) <-- X(:,r) L^{-T}` (solves `L x' = x'`), for `cnt` rows
!  (`L` lower, row major: `L(k,j)`, `k <= j`).

    subroutine eq_trsm_rows(L, nz, X, cnt)

    real(wp), intent(in) :: L(:,:)     !! Cholesky factor
    integer(ip), intent(in) :: nz      !! dimension
    real(wp), intent(inout) :: X(:,:)  !! rows (columns of X)
    integer(ip), intent(in) :: cnt     !! number of rows

    integer(ip) :: i, j, k
    real(wp) :: dj, s, s0, s1, s2, s3, lk

    i = 1
    do while (i+3 <= cnt)
        do j = 1, nz
            dj = 1.0_wp/L(j,j)
            s0 = X(j,i); s1 = X(j,i+1); s2 = X(j,i+2); s3 = X(j,i+3)
            do k = 1, j-1
                lk = L(k,j)
                s0 = s0 - lk*X(k,i); s1 = s1 - lk*X(k,i+1)
                s2 = s2 - lk*X(k,i+2); s3 = s3 - lk*X(k,i+3)
            end do
            X(j,i) = s0*dj; X(j,i+1) = s1*dj; X(j,i+2) = s2*dj; X(j,i+3) = s3*dj
        end do
        i = i + 4
    end do
    do while (i <= cnt)
        do j = 1, nz
            s = X(j,i)
            do k = 1, j-1
                s = s - L(k,j)*X(k,i)
            end do
            X(j,i) = s/L(j,j)
        end do
        i = i + 1
    end do

    end subroutine eq_trsm_rows
!*****************************************************************************************

!*****************************************************************************************
!>
!  `C(:,i) = W' a_{ids(i)}` for the general constraints `ids` (the rows of
!  `A` times `W`, with `W(:,k)` the k-th row of `W`).

    subroutine eq_rows_times_W(At, ids, ms, cnt, n, W, nz, C)

    real(wp), intent(in) :: At(:,:)    !! rows of A (columns)
    integer(ip), intent(in) :: ids(:)  !! constraints
    integer(ip), intent(in) :: ms      !! number of simple bounds
    integer(ip), intent(in) :: cnt     !! number of constraints
    integer(ip), intent(in) :: n       !! number of variables
    real(wp), intent(in) :: W(:,:)     !! null-space basis (rows as columns)
    integer(ip), intent(in) :: nz      !! number of reduced variables
    real(wp), intent(inout) :: C(:,:)  !! products (columns)

    integer(ip) :: i, j, k, r0, r1, r2, r3
    real(wp) :: b0, b1, b2, b3, wjk

    i = 1
    do while (i+3 <= cnt)
        r0 = ids(i)-ms; r1 = ids(i+1)-ms; r2 = ids(i+2)-ms; r3 = ids(i+3)-ms
        C(1:nz,i:i+3) = 0.0_wp
        do k = 1, n
            b0 = At(k,r0); b1 = At(k,r1); b2 = At(k,r2); b3 = At(k,r3)
            if (b0 == 0.0_wp .and. b1 == 0.0_wp .and. b2 == 0.0_wp .and. b3 == 0.0_wp) cycle
            do j = 1, nz
                wjk = W(j,k)
                C(j,i) = C(j,i) + b0*wjk; C(j,i+1) = C(j,i+1) + b1*wjk
                C(j,i+2) = C(j,i+2) + b2*wjk; C(j,i+3) = C(j,i+3) + b3*wjk
            end do
        end do
        i = i + 4
    end do
    do while (i <= cnt)
        r0 = ids(i)-ms
        C(1:nz,i) = 0.0_wp
        do k = 1, n
            b0 = At(k,r0)
            if (b0 == 0.0_wp) cycle
            do j = 1, nz
                C(j,i) = C(j,i) + b0*W(j,k)
            end do
        end do
        i = i + 1
    end do

    end subroutine eq_rows_times_W
!*****************************************************************************************

!*****************************************************************************************
!>
!  Dot product, summed in order.

    pure real(wp) function eq_dot(a, b, n)

    real(wp), intent(in) :: a(:)  !! first vector
    real(wp), intent(in) :: b(:)  !! second vector
    integer(ip), intent(in) :: n  !! length

    integer(ip) :: i

    eq_dot = 0.0_wp
    do i = 1, n
        eq_dot = eq_dot + a(i)*b(i)
    end do

    end function eq_dot
!*****************************************************************************************

!*****************************************************************************************
!>
!  `y <-- H x` (dense, or diagonal if `metric`).

    subroutine eq_hess_times(qp, metric, x, y)

    type(daqp_problem), intent(in) :: qp !! the problem
    logical, intent(in) :: metric        !! H is diagonal
    real(wp), intent(in) :: x(:)         !! vector
    real(wp), intent(inout) :: y(:)      !! H x

    integer(ip) :: i, n

    n = qp%n
    if (.not. qp%has_H) then
        y(1:n) = 0.0_wp
        return
    end if
    if (metric) then
        do i = 1, n
            y(i) = qp%Hc(i,i)*x(i)
        end do
        return
    end if
    do i = 1, n
        y(i) = eq_dot(qp%Hc(:,i), x, n)
    end do

    end subroutine eq_hess_times
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the general constraint `i` is an equality that can be eliminated
!  (the same criteria as the detection of equal bounds).

    pure logical function eq_is_candidate(work, i)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: i             !! constraint

    integer(ip) :: s

    s = work%sense(i)
    eq_is_candidate = .false.
    if (iand(s, daqp_soft+daqp_binary) /= 0) return
    if (work%qp%bupper(i) - work%qp%blower(i) < work%settings%zero_tol) then
        eq_is_candidate = .true.
        return
    end if
    if (has(s,daqp_auto_equality)) return ! its bounds are no longer equal
    eq_is_candidate = iand(s, daqp_active+daqp_immutable) == daqp_active+daqp_immutable

    end function eq_is_candidate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Number of equality candidates.

    pure integer(ip) function eq_count_candidates(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i

    eq_count_candidates = 0
    do i = work%qp%ms+1, work%qp%m
        if (eq_is_candidate(work,i)) eq_count_candidates = eq_count_candidates + 1
    end do

    end function eq_count_candidates
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `H` is diagonal.

    pure logical function eq_is_diagonal(qp, zero_tol)

    type(daqp_problem), intent(in) :: qp !! the problem
    real(wp), intent(in) :: zero_tol     !! tolerance

    integer(ip) :: i, j

    eq_is_diagonal = .false.
    if (.not. qp%has_H) return
    do i = 1, qp%n
        do j = 1, qp%n
            if (i /= j .and. (qp%Hc(j,i) > zero_tol .or. qp%Hc(j,i) < -zero_tol)) return
        end do
    end do
    eq_is_diagonal = .true.

    end function eq_is_diagonal
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether eliminating `n_eq` equalities is expected to pay off.

    pure logical function eq_is_worthwhile(work, n_eq)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: n_eq          !! number of equality candidates

    integer(ip) :: n, n_ineq

    n = work%qp%n
    n_ineq = work%qp%m - work%qp%ms - n_eq
    eq_is_worthwhile = .true.
    if (n_eq >= n) return
    eq_is_worthwhile = .false.
    if (n < eq_min_dim .or. n_eq <= eq_min_count .or. eq_min_ratio*n_eq <= n) return
    if (eq_is_diagonal(work%qp, work%settings%zero_tol) .and. &
        (n_ineq == 0 .or. eq_diag_min_ratio*n_eq < n)) return
    eq_is_worthwhile = .true.

    end function eq_is_worthwhile
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the problem in the workspace is to be reduced for an update with
!  `mask` (`settings%eq_reduction`; `auto` only with `daqp_update_eliminate`).

    logical function daqp_eq_wanted(work, mask)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: mask          !! update mask

    integer(ip) :: n_eq, policy

    daqp_eq_wanted = .false.
    policy = work%settings%eq_reduction
    if (policy == daqp_eq_reduction_off) return
    if (policy /= daqp_eq_reduction_on .and. .not. has(mask,daqp_update_eliminate)) return
    if (work%qp%m <= work%qp%ms) return
    ! a hierarchy refers to the constraints by their index, and a factored
    ! Hessian is not available as a Hessian
    if (work%qp%nh > 1 .or. is_hierarchical(work) .or. &
        work%qp%problem_type == daqp_problem_factored) return
    n_eq = eq_count_candidates(work)
    if (n_eq == 0) return
    if (policy == daqp_eq_reduction_on) then
        daqp_eq_wanted = .true.
        return
    end if
    daqp_eq_wanted = eq_is_worthwhile(work, n_eq)

    end function daqp_eq_wanted
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the equality candidates are the ones that the reduction was formed for.

    logical function eq_same_candidates(work)

    type(daqp_workspace), intent(in) :: work !! workspace

    integer(ip) :: i, k

    eq_same_candidates = .false.
    k = 0
    do i = work%qp%ms+1, work%qp%m
        if (.not. eq_is_candidate(work,i)) cycle
        if (k == work%eq%ncand) return
        if (work%eq%cand_ids(k+1) /= i) return
        k = k + 1
    end do
    eq_same_candidates = k == work%eq%ncand

    end function eq_same_candidates
!*****************************************************************************************

!*****************************************************************************************
!>
!  Free the LDP of the reduced problem (not the reduced problem itself).

    subroutine eq_free_reduced_ldp(eq)

    type(daqp_eq_data), intent(inout) :: eq !! elimination

    associate (d => eq%other)
        if (allocated(d%Mr)) deallocate(d%Mr)
        if (allocated(d%dupper)) deallocate(d%dupper)
        if (allocated(d%dlower)) deallocate(d%dlower)
        if (allocated(d%scaling)) deallocate(d%scaling)
        if (allocated(d%Mu)) deallocate(d%Mu)
        if (allocated(d%sense)) deallocate(d%sense)
        if (allocated(d%R)) deallocate(d%R)
        if (allocated(d%v)) deallocate(d%v)
        if (allocated(d%bin_ids)) deallocate(d%bin_ids)
        if (allocated(d%rho_ls)) deallocate(d%rho_ls)
        if (allocated(d%rho_us)) deallocate(d%rho_us)
        if (allocated(d%w_ls)) deallocate(d%w_ls)
        if (allocated(d%w_us)) deallocate(d%w_us)
        d%n = 0; d%m = 0; d%ms = 0
        d%rmode = rinv_none
        d%has_v = .false.
        d%has_weights = .false.
        d%state = 0
        d%n_prox = 0
        d%nb = 0
    end associate
    if (allocated(eq%rho_r)) deallocate(eq%rho_r)

    end subroutine eq_free_reduced_ldp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Storage of the LDP of the reduced problem (which has no simple bounds).

    subroutine eq_allocate_reduced_ldp(eq, nb)

    type(daqp_eq_data), intent(inout) :: eq !! elimination
    integer(ip), intent(in) :: nb           !! number of binary constraints

    integer(ip) :: nz, mr

    call eq_free_reduced_ldp(eq)
    nz = eq%nz
    mr = eq%mr
    associate (d => eq%other)
        d%n = nz; d%m = mr; d%ms = 0
        allocate(d%Mr(nz,mr), d%dupper(mr), d%dlower(mr), d%scaling(mr), d%Mu(mr), d%sense(mr))
        d%sense = 0
        d%scaling = 1.0_wp
        if (d%qp%has_H) then
            allocate(d%R((nz*(nz+1))/2))
            d%rmode = rinv_dense
        else
            allocate(d%R(0))
            d%rmode = rinv_none
        end if
        allocate(d%v(nz))
        d%v = 0.0_wp
        d%has_v = d%qp%has_f
        if (nb > 0) then
            allocate(d%bin_ids(nb))
        else
            allocate(d%bin_ids(0))
        end if
    end associate

    end subroutine eq_allocate_reduced_ldp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Free the cached responses to the right-hand side.

    subroutine eq_free_rhs_cache(eq)

    type(daqp_eq_data), intent(inout) :: eq !! elimination

    if (allocated(eq%cols)) deallocate(eq%cols)
    eq%ncols = 0
    if (allocated(eq%xf)) deallocate(eq%xf)
    if (allocated(eq%gf)) deallocate(eq%gf)
    if (allocated(eq%df)) deallocate(eq%df)
    if (allocated(eq%sh)) deallocate(eq%sh)
    eq%f_valid = .false.

    end subroutine eq_free_rhs_cache
!*****************************************************************************************

!*****************************************************************************************
!>
!  Free the reduction (and the reduced problem).

    subroutine eq_free_reduction(eq)

    type(daqp_eq_data), intent(inout) :: eq !! elimination

    call eq_free_rhs_cache(eq)
    if (allocated(eq%eq_ids)) deallocate(eq%eq_ids)
    if (allocated(eq%cand_ids)) deallocate(eq%cand_ids)
    if (allocated(eq%keep)) deallocate(eq%keep)
    if (allocated(eq%drop_ids)) deallocate(eq%drop_ids)
    if (allocated(eq%V)) deallocate(eq%V)
    if (allocated(eq%tau)) deallocate(eq%tau)
    if (allocated(eq%s_eq)) deallocate(eq%s_eq)
    if (allocated(eq%R)) deallocate(eq%R)
    if (allocated(eq%dsq)) deallocate(eq%dsq)
    if (allocated(eq%W)) deallocate(eq%W)
    if (allocated(eq%xp)) deallocate(eq%xp)
    if (allocated(eq%tmp)) deallocate(eq%tmp)
    eq%other%qp = daqp_problem()
    eq%allocated_dims = .false.

    end subroutine eq_free_reduction
!*****************************************************************************************

!*****************************************************************************************
!>
!  Storage that only depends on the dimensions of the original problem.

    subroutine eq_allocate_reduction(eq, n, m, ms)

    type(daqp_eq_data), intent(inout) :: eq !! elimination
    integer(ip), intent(in) :: n, m, ms     !! dimensions of the original problem

    if (eq%allocated_dims .and. eq%n == n .and. eq%m == m .and. eq%ms == ms) return
    call eq_free_reduction(eq)
    eq%n = n; eq%m = m; eq%ms = ms
    allocate(eq%eq_ids(m), eq%cand_ids(m), eq%keep(m), eq%drop_ids(m))
    allocate(eq%V(n,n), eq%tau(n), eq%s_eq(n), eq%R((n*(n+1))/2), eq%xp(n), eq%tmp(3*n))
    allocate(eq%other%qp%bupper(m), eq%other%qp%blower(m), eq%other%qp%sense(m))
    eq%allocated_dims = .true.

    end subroutine eq_allocate_reduction
!*****************************************************************************************

!*****************************************************************************************
!>
!  Householder QR of `A_E'` (`A_E H^{-1/2}` if `metric`), left-looking in
!  panels of four candidates. Candidates that are (numerically) linearly
!  dependent on the ones before are not eliminated.

    subroutine eq_build_qr(work, zero_tol)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in) :: zero_tol            !! tolerance

    real(wp), allocatable :: panel(:,:)
    real(wp) :: pn(4), tol, nrm, alpha, beta, d
    integer(ip) :: c, i, p, q, neq, cnt, k0, n, ms, r0

    associate (eq => work%eq, qp => work%qp)
    n = eq%n
    ms = eq%ms
    tol = sqrt(zero_tol)
    allocate(panel(n,4))
    neq = 0
    c = 0
    do while (c < eq%ncand .and. neq < n)
        cnt = min(4_ip, eq%ncand-c)
        k0 = neq
        do p = 1, cnt
            r0 = eq%cand_ids(c+p) - ms
            if (eq%metric) then
                do i = 1, n
                    panel(i,p) = qp%At(i,r0)*eq%dsq(i)
                end do
            else
                panel(1:n,p) = qp%At(1:n,r0)
            end if
            nrm = 0.0_wp
            do i = 1, n
                nrm = nrm + panel(i,p)*panel(i,p)
            end do
            if (nrm <= zero_tol) then
                pn(p) = 0.0_wp
            else
                pn(p) = 1.0_wp/sqrt(nrm)
            end if
            do i = 1, n
                panel(i,p) = panel(i,p)*pn(p)
            end do
        end do
        call eq_apply_QT_many(eq%V, eq%tau, k0, n, panel, cnt)
        do p = 1, cnt
            if (neq >= n) exit
            if (pn(p) == 0.0_wp) cycle ! empty constraint
            do q = k0+1, neq
                call eq_reflect(eq%V(:,q), eq%tau(q), q, n, panel(:,p))
            end do
            alpha = 0.0_wp
            do i = neq+1, n
                alpha = alpha + panel(i,p)*panel(i,p)
            end do
            alpha = sqrt(alpha)
            if (alpha <= tol) cycle ! linearly dependent
            if (panel(neq+1,p) > 0.0_wp) then
                beta = -alpha
            else
                beta = alpha
            end if
            d = panel(neq+1,p) - beta
            eq%tau(neq+1) = -d/beta
            do i = neq+2, n
                panel(i,p) = panel(i,p)/d
            end do
            panel(neq+1,p) = beta
            do i = 1, neq+1
                eq%R((neq*(neq+1))/2+i) = panel(i,p)
            end do
            eq%V(1:n,neq+1) = panel(1:n,p)
            eq%s_eq(neq+1) = pn(p)
            neq = neq + 1
            eq%eq_ids(neq) = eq%cand_ids(c+p)
        end do
        c = c + 4
    end do
    eq%neq = neq
    end associate

    end subroutine eq_build_qr
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `H` is symmetric (for an AVI).

    pure logical function eq_is_symmetric(qp, zero_tol)

    type(daqp_problem), intent(in) :: qp !! the problem
    real(wp), intent(in) :: zero_tol     !! tolerance

    integer(ip) :: i, j
    real(wp) :: scale, d, diff

    scale = 0.0_wp
    do i = 1, qp%n
        d = qp%Hc(i,i)
        if (d < 0.0_wp) d = -d
        if (d > scale) scale = d
    end do
    if (scale < 1.0_wp) scale = 1.0_wp
    eq_is_symmetric = .false.
    do i = 1, qp%n
        do j = i+1, qp%n
            diff = qp%Hc(j,i) - qp%Hc(i,j)
            if (diff > zero_tol*scale .or. diff < -zero_tol*scale) return
        end do
    end do
    eq_is_symmetric = .true.

    end function eq_is_symmetric
!*****************************************************************************************

!*****************************************************************************************
!>
!  Split `Z` (`n x nz`) into `[Z1 Z2]`, with `Z2` zero in the rows `cid` (the
!  variables with curvature), so that `H*Z2 = 0` exactly. Forms
!  `Hr = blockdiag(Z1'HZ1, 0)` and returns the dimension of `Z1`.

    integer(ip) function eq_split_flat(qp, Z, n, nz, cid, nc, zero_tol, Hr) result(r)

    type(daqp_problem), intent(in) :: qp   !! the problem
    real(wp), intent(inout) :: Z(:,:)      !! null-space basis (columns)
    integer(ip), intent(in) :: n           !! number of variables
    integer(ip), intent(in) :: nz          !! dimension of the null space
    integer(ip), intent(in) :: cid(:)      !! variables with curvature
    integer(ip), intent(in) :: nc          !! number of them
    real(wp), intent(in) :: zero_tol       !! tolerance
    real(wp), intent(inout) :: Hr(:,:)     !! reduced Hessian

    real(wp), allocatable :: Y(:,:), tau(:), Zr(:,:), G(:,:)
    real(wp) :: tol, alpha, beta, d, sm
    integer(ip) :: i, j, k, c

    tol = sqrt(zero_tol)
    allocate(Y(nz,max(1_ip,nc)), tau(nc+1), Zr(nz,n))
    Y = 0.0_wp
    tau = 0.0_wp
    r = 0

    ! Householder QR of Z_C' (the rows cid of Z), with its rank r
    do c = 1, nc
        if (r >= nz) exit
        do j = 1, nz
            Y(j,r+1) = Z(cid(c),j)
        end do
        do k = 1, r
            call eq_reflect(Y(:,k), tau(k), k, nz, Y(:,r+1))
        end do
        alpha = 0.0_wp
        do j = r+1, nz
            alpha = alpha + Y(j,r+1)*Y(j,r+1)
        end do
        alpha = sqrt(alpha)
        if (alpha <= tol) cycle ! dependent on the rows before
        if (Y(r+1,r+1) > 0.0_wp) then
            beta = -alpha
        else
            beta = alpha
        end if
        d = Y(r+1,r+1) - beta
        tau(r+1) = -d/beta
        do j = r+2, nz
            Y(j,r+1) = Y(j,r+1)/d
        end do
        Y(r+1,r+1) = beta
        r = r + 1
    end do

    ! Z <-- Z Q: rows (Q'z')', the last nz-r columns are Z2
    do i = 1, n
        do j = 1, nz
            Zr(j,i) = Z(i,j)
        end do
    end do
    call eq_apply_QT_many(Y, tau, r, nz, Zr, n)
    do c = 1, nc
        do j = r+1, nz
            Zr(j,cid(c)) = 0.0_wp ! roundoff
        end do
    end do
    do i = 1, n
        do j = 1, nz
            Z(i,j) = Zr(j,i)
        end do
    end do

    ! Hr = blockdiag(Z1_C' H_CC Z1_C, 0), with G = H_CC Z1_C
    Hr(1:nz,1:nz) = 0.0_wp
    allocate(G(max(1_ip,r),max(1_ip,nc)))
    do c = 1, nc
        do j = 1, r
            sm = 0.0_wp
            do k = 1, nc
                sm = sm + qp%Hc(cid(k),cid(c))*Zr(j,cid(k))
            end do
            G(j,c) = sm
        end do
    end do
    do i = 1, r
        do j = i, r
            sm = 0.0_wp
            do c = 1, nc
                sm = sm + Zr(i,cid(c))*G(j,c)
            end do
            Hr(j,i) = sm
            Hr(i,j) = sm
        end do
    end do

    end function eq_split_flat
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form the reduction of the problem in the workspace: the factorizations, the
!  reduced constraints and the storage of the reduced problem. Returns 0 if
!  nothing can be eliminated, or if a soft constraint would be left out.

    integer(ip) function eq_build_reduction(work) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: n, m, ms, i, j, k, c, neq, nz, mI, mtot, mr, nb, nc, id, istat
    logical :: use_twoside, use_refl, symmetric, split, chol_ok
    integer(ip), allocatable :: gen_ids(:), cid(:)
    real(wp), allocatable :: L(:,:), B(:,:), HZ(:,:), Ak(:,:)
    real(wp) :: zero_tol, scale, s, rN, rNZ, rNE, rMI, f_gemm, f_two, f_w, f_refl, &
                eps_mach, wmax, hmax, s2
    logical :: have_L

    n = work%qp%n
    m = work%qp%m
    ms = work%qp%ms
    zero_tol = work%settings%zero_tol
    mI = 0
    nb = 0
    nc = 0
    use_twoside = .false.
    use_refl = .false.
    symmetric = .true.
    split = .false.
    have_L = .false.
    flag = 0
    istat = 0

    call eq_allocate_reduction(work%eq, n, m, ms)
    call eq_free_rhs_cache(work%eq) ! formed for the previous reduction
    work%eq%active = .false.
    if (allocated(work%bnb)) work%bnb%n_root_ws = 0 ! refers to other constraints

    associate (eq => work%eq, qp => work%qp)

    ! equality candidates
    eq%ncand = 0
    do i = ms+1, m
        if (eq_is_candidate(work,i)) then
            eq%ncand = eq%ncand + 1
            eq%cand_ids(eq%ncand) = i
        end if
    end do

    ! a positive diagonal Hessian is used as a metric in the QR
    eq%metric = .false.
    if (qp%has_H) then
        if (eq_is_diagonal(qp, zero_tol)) then
            scale = 0.0_wp
            do i = 1, n
                if (qp%Hc(i,i) > scale) scale = qp%Hc(i,i)
            end do
            eq%metric = scale > 0.0_wp
            do i = 1, n
                if (.not. eq%metric) exit
                if (qp%Hc(i,i) <= zero_tol*scale) eq%metric = .false.
            end do
            if (eq%metric) then
                call resize1(eq%dsq, n, istat)
                do i = 1, n
                    eq%dsq(i) = 1.0_wp/sqrt(qp%Hc(i,i))
                end do
            end if
        end if
        if (.not. eq%metric .and. qp%problem_type == daqp_problem_avi) then
            symmetric = eq_is_symmetric(qp, zero_tol)
        end if
    end if

    call eq_build_qr(work, zero_tol)
    neq = eq%neq
    if (neq == 0) return ! nothing eliminated
    nz = n - neq
    eq%nz = nz

    ! the general constraints that are kept (dependent equalities included)
    allocate(gen_ids(max(1_ip,m-ms)))
    k = 0
    do i = ms+1, m
        if (k < neq) then
            if (eq%eq_ids(k+1) == i) then
                k = k + 1
                cycle
            end if
        end if
        mI = mI + 1
        gen_ids(mI) = i
    end do

    ! the equalities determine x = xp: all other constraints only have to be
    ! consistent with xp
    if (nz == 0) then
        eq%ndrop = 0
        do c = 1, ms+mI
            if (c <= ms) then
                id = c
            else
                id = gen_ids(c-ms)
            end if
            if (has(work%sense(id),daqp_soft)) return
            eq%ndrop = eq%ndrop + 1
            eq%drop_ids(eq%ndrop) = id
        end do
        eq%path = eq_path_ldp
        if (allocated(eq%other%qp%Hc)) deallocate(eq%other%qp%Hc)
        if (allocated(eq%other%qp%f)) deallocate(eq%other%qp%f)
        if (allocated(eq%other%qp%At)) deallocate(eq%other%qp%At)
        if (allocated(eq%W)) deallocate(eq%W)
        allocate(eq%W(0,n), eq%other%qp%At(0,0), eq%other%qp%f(0))
        eq%mr = 0
        call set_reduced_problem(0_ip, 0_ip)
        eq%other%qp%has_H = .false.
        eq%other%qp%has_f = .false.
        call eq_allocate_reduced_ldp(eq, 0_ip)
        flag = 1
        return
    end if

    ! pick the kernels by their flop counts (the two-sided Householder
    ! product is memory bound, so its count is weighted)
    if (qp%has_H .and. .not. eq%metric .and. symmetric) then
        rN = real(n,wp); rNZ = real(nz,wp)
        f_gemm = rN*rN*rNZ + rN*rNZ*rNZ/2.0_wp
        f_two = 2.0_wp/3.0_wp*(rN*rN*rN-rNZ*rNZ*rNZ)
        use_twoside = 1.7_wp*f_two < f_gemm
    end if
    rN = real(n,wp); rNZ = real(nz,wp); rNE = real(neq,wp); rMI = real(mI,wp)
    f_w = rMI*rN*rNZ
    f_refl = rMI*(2.0_wp*rN*rNE + rNZ*rNZ/2.0_wp)
    use_refl = f_refl < f_w

    call eq_accumulate_Z(eq%V, eq%tau, neq, n)
    ! (column j of Z is eq%V(:,neq+j))

    ! variables with curvature (nonzero rows of H), for eq_split_flat
    if (qp%has_H .and. .not. eq%metric .and. symmetric) then
        allocate(cid(n))
        do i = 1, n
            do j = 1, n
                if (qp%Hc(j,i) /= 0.0_wp .or. qp%Hc(i,j) /= 0.0_wp) exit
            end do
            if (j <= n) then
                nc = nc + 1
                cid(nc) = i
            end if
        end do
        split = nc < nz ! then some null directions have no curvature
    end if

    ! the reduced Hessian, and how the reduced problem is posed
    if (allocated(eq%other%qp%Hc)) deallocate(eq%other%qp%Hc)
    if (allocated(eq%other%qp%f)) deallocate(eq%other%qp%f)
    if (.not. qp%has_H) then
        eq%path = eq_path_lp
    else if (eq%metric) then
        eq%path = eq_path_ldp
    else
        allocate(eq%other%qp%Hc(nz,nz))
        if (split) then
            k = eq_split_flat(qp, eq%V(:,neq+1:n), n, nz, cid, nc, zero_tol, eq%other%qp%Hc)
            use_refl = .false. ! the rows (A Q)_2 refer to the unsplit Z
        else if (use_twoside) then
            allocate(B(n,n))
            B = qp%Hc
            call eq_twoside(B, n, eq%V, eq%tau, neq, eq%tmp)
            do i = 1, nz
                do j = 1, i
                    eq%other%qp%Hc(j,i) = B(neq+j,neq+i)
                    eq%other%qp%Hc(i,j) = B(neq+j,neq+i)
                end do
            end do
            deallocate(B)
        else
            allocate(HZ(n,nz)) ! column j = H z_j
            call eq_gemm_tn(n, nz, n, eq%V(:,neq+1:n), qp%Hc, HZ, .false.)
            call eq_gemm_tn(n, nz, nz, eq%V(:,neq+1:n), HZ, eq%other%qp%Hc, symmetric)
            deallocate(HZ)
        end if
        eq%path = eq_path_qp
        if (symmetric) then
            allocate(L(nz,nz))
            L = eq%other%qp%Hc
            chol_ok = eq_chol(L, nz, zero_tol)
            if (chol_ok) then
                eq%path = eq_path_ldp
                have_L = .true.
            else
                deallocate(L)
            end if
        end if
    end if

    ! W: H^{-1/2} Z, Z L^{-T}, or Z
    if (allocated(eq%W)) deallocate(eq%W)
    allocate(eq%W(nz,n))
    do ! form W
        do i = 1, n
            s = 1.0_wp
            if (eq%metric) s = eq%dsq(i)
            do j = 1, nz
                eq%W(j,i) = s*eq%V(i,neq+j)
            end do
        end do
        if (have_L) then
            call eq_trsm_rows(L, nz, eq%W, n)
            ! ill-conditioned Z'HZ => PATH_QP (where it is regularized)
            eps_mach = epsilon(1.0_wp)
            wmax = 0.0_wp
            hmax = 0.0_wp
            do j = 1, nz
                s2 = 0.0_wp
                do i = 1, n
                    s2 = s2 + eq%W(j,i)*eq%W(j,i)
                end do
                if (s2 > wmax) wmax = s2
                if (eq%other%qp%Hc(j,j) > hmax) hmax = eq%other%qp%Hc(j,j)
            end do
            if ((nc < n .and. wmax*hmax > hessian_cond_max) .or. &
                real(nz,wp)*eps_mach*wmax*hmax > hessian_cond_eps) then
                deallocate(L)
                have_L = .false.
                eq%path = eq_path_qp
                cycle
            end if
        end if
        exit
    end do
    if (eq%path == eq_path_ldp) then
        if (.not. allocated(eq%other%qp%Hc)) allocate(eq%other%qp%Hc(nz,nz))
        eq%other%qp%Hc = 0.0_wp
        do i = 1, nz
            eq%other%qp%Hc(i,i) = 1.0_wp
        end do
    else
        allocate(eq%other%qp%f(nz))
    end if

    ! reduced constraints: rows of W for the simple bounds, then A_I W
    mtot = ms + mI
    if (allocated(eq%other%qp%At)) deallocate(eq%other%qp%At)
    allocate(eq%other%qp%At(nz,max(1_ip,mtot)))
    do i = 1, ms
        eq%other%qp%At(1:nz,i) = eq%W(1:nz,i)
    end do
    if (mI > 0) then
        if (.not. use_refl) then
            call eq_rows_times_W(qp%At, gen_ids, ms, mI, n, eq%W, nz, eq%other%qp%At(:,ms+1:mtot))
        else
            allocate(Ak(n,mI))
            do i = 1, mI
                if (eq%metric) then
                    do k = 1, n
                        Ak(k,i) = qp%At(k,gen_ids(i)-ms)*eq%dsq(k)
                    end do
                else
                    Ak(1:n,i) = qp%At(1:n,gen_ids(i)-ms)
                end if
            end do
            call eq_apply_QT_many(eq%V, eq%tau, neq, n, Ak, mI) ! rows <-- (Q'a')'
            do i = 1, mI
                do j = 1, nz
                    eq%other%qp%At(j,ms+i) = Ak(neq+j,i)
                end do
            end do
            deallocate(Ak)
            if (have_L) call eq_trsm_rows(L, nz, eq%other%qp%At(:,ms+1:mtot), mI)
        end if
    end if

    ! keep the constraints that the reduced variables affect; the others are
    ! implied by the equalities: they only have to be consistent with them
    mr = 0
    eq%ndrop = 0
    do c = 1, mtot
        if (c <= ms) then
            id = c
        else
            id = gen_ids(c-ms)
        end if
        if (eq_dot(eq%other%qp%At(:,c), eq%other%qp%At(:,c), nz) <= zero_tol) then
            if (has(work%sense(id),daqp_soft)) return
            eq%ndrop = eq%ndrop + 1
            eq%drop_ids(eq%ndrop) = id
            cycle
        end if
        mr = mr + 1
        if (mr /= c) eq%other%qp%At(1:nz,mr) = eq%other%qp%At(1:nz,c)
        eq%keep(mr) = id
    end do
    eq%mr = mr
    do c = 1, mr
        eq%other%qp%sense(c) = work%sense(eq%keep(c))
        if (has(eq%other%qp%sense(c),daqp_binary)) nb = nb + 1
    end do

    ! the reduced problem
    call set_reduced_problem(nz, mr)
    eq%other%qp%has_H = allocated(eq%other%qp%Hc)
    eq%other%qp%has_f = allocated(eq%other%qp%f)
    if (.not. eq%other%qp%has_f) allocate(eq%other%qp%f(nz))
    if (.not. eq%other%qp%has_H) allocate(eq%other%qp%Hc(0,0))
    if (.not. eq%other%qp%has_f) eq%other%qp%f = 0.0_wp

    call eq_allocate_reduced_ldp(eq, nb)
    flag = 1

    end associate

    contains

        subroutine set_reduced_problem(nr, mrr)
            !! dimensions and flags of the reduced problem
            integer(ip), intent(in) :: nr, mrr
            work%eq%other%qp%n = nr
            work%eq%other%qp%m = mrr
            work%eq%other%qp%ms = 0
            work%eq%other%qp%has_sense = .true.
            work%eq%other%qp%nh = 1
            work%eq%other%qp%problem_type = work%qp%problem_type
            if (allocated(work%eq%other%qp%break_points)) deallocate(work%eq%other%qp%break_points)
            if (allocated(work%eq%other%qp%Rf)) deallocate(work%eq%other%qp%Rf)
        end subroutine set_reduced_problem

    end function eq_build_reduction
!*****************************************************************************************

!*****************************************************************************************
!>
!  Value of equality `id` (the upper bound, unless it is active at its lower one).

    pure real(wp) function eq_value(work, id)

    type(daqp_workspace), intent(in) :: work !! workspace
    integer(ip), intent(in) :: id            !! constraint

    if (has(work%sense(id),daqp_lower)) then
        eq_value = work%qp%blower(id)
    else
        eq_value = work%qp%bupper(id)
    end if

    end function eq_value
!*****************************************************************************************

!*****************************************************************************************
!>
!  Response of `xp` to a unit right-hand side of eliminated equality `k` (in
!  `eq_path_ldp`), stored in `eq%cols(k)%v = [xp_k, H xp_k, shifts]`.

    subroutine eq_rhs_column(work, k)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: k                !! eliminated equality

    integer(ip) :: n, ms, neq, nz, mr, i, j, c, id, rs
    real(wp) :: s, hi

    associate (eq => work%eq, qp => work%qp)
    n = eq%n; ms = eq%ms; neq = eq%neq; nz = eq%nz; mr = eq%mr
    if (.not. allocated(eq%cols)) then
        allocate(eq%cols(neq))
        eq%ncols = neq
    end if
    if (allocated(eq%cols(k)%v)) return
    allocate(eq%cols(k)%v(2*n+mr))
    associate (col => eq%cols(k)%v, t => eq%tmp(2*n+1:3*n))
    ! R'y = s_k e_k, forward substitution from k
    col(1:n) = 0.0_wp
    do i = k, neq
        rs = ((i-1)*i)/2
        if (i == k) then
            s = eq%s_eq(k)
        else
            s = 0.0_wp
        end if
        do j = k, i-1
            s = s - eq%R(rs+j)*col(j)
        end do
        col(i) = s/eq%R(rs+i)
    end do
    do i = neq, 1, -1
        call eq_reflect(eq%V(:,i), eq%tau(i), i, n, col(1:n))
    end do
    if (eq%metric) then
        do i = 1, n
            col(i) = col(i)*eq%dsq(i)
        end do
    end if
    ! xk <-- xk - W W'H xk, hk = H xk
    call eq_hess_times(qp, eq%metric, col(1:n), col(n+1:2*n))
    t(1:nz) = 0.0_wp
    do i = 1, n
        hi = col(n+i)
        if (hi == 0.0_wp) cycle
        do j = 1, nz
            t(j) = t(j) + eq%W(j,i)*hi
        end do
    end do
    do i = 1, n
        col(i) = col(i) - eq_dot(eq%W(:,i), t, nz)
    end do
    call eq_hess_times(qp, eq%metric, col(1:n), col(n+1:2*n))
    do c = 1, mr
        id = eq%keep(c)
        if (id <= ms) then
            col(2*n+c) = -col(id)
        else
            col(2*n+c) = -eq_dot(qp%At(:,id-ms), col(1:n), n)
        end if
    end do
    end associate
    end associate

    end subroutine eq_rhs_column
!*****************************************************************************************

!*****************************************************************************************
!>
!  Response of `xp` to `f` (`eq_path_ldp`): `xf = -W W'f`, `gf = H xf + f`,
!  `df = -A xf`.

    subroutine eq_rhs_f(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: n, nz, mr, i, j, c, istat
    real(wp) :: fi

    istat = 0
    associate (eq => work%eq, qp => work%qp)
    n = eq%n; nz = eq%nz; mr = eq%mr
    call resize1(eq%xf, n, istat)
    call resize1(eq%gf, n, istat)
    call resize1(eq%df, mr, istat)
    call resize1(eq%sh, mr, istat)
    associate (t => eq%tmp(2*n+1:3*n))
    t(1:nz) = 0.0_wp
    if (qp%has_f) then
        do i = 1, n
            fi = qp%f(i)
            if (fi == 0.0_wp) cycle
            do j = 1, nz
                t(j) = t(j) + eq%W(j,i)*fi
            end do
        end do
    end if
    do i = 1, n
        eq%xf(i) = -eq_dot(eq%W(:,i), t, nz)
    end do
    call eq_hess_times(qp, eq%metric, eq%xf, eq%gf)
    if (qp%has_f) then
        do i = 1, n
            eq%gf(i) = eq%gf(i) + qp%f(i)
        end do
    end if
    ! -a_c xf = a_c W W'f = Ar_c W'f (and the rows of W for the simple bounds)
    do c = 1, mr
        eq%df(c) = eq_dot(eq%other%qp%At(:,c), t, nz)
    end do
    end associate
    eq%f_valid = .true.
    end associate

    end subroutine eq_rhs_f
!*****************************************************************************************

!*****************************************************************************************
!>
!  The part that depends on `b` and `f`: the particular solution `xp`, the
!  linear term of the reduced problem (`eq_path_qp`, `eq_path_lp`), and the
!  shifted bounds. With `warm`, the responses to `b_E` and `f` are used if
!  few equalities have a nonzero right-hand side.
!  Returns 1, or a negative exit flag if the implied constraints are violated.

    integer(ip) function eq_reduce_rhs(work, warm) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    logical, intent(in) :: warm                 !! the reduction was not formed anew

    integer(ip) :: n, ms, neq, nz, mr, i, j, k, c, nnz, id, rs
    real(wp) :: primal_tol, b, s, gi, shift, val, fi

    flag = 1
    associate (eq => work%eq, qp => work%qp)
    n = eq%n; ms = eq%ms; neq = eq%neq; nz = eq%nz; mr = eq%mr
    primal_tol = work%settings%primal_tol
    associate (y => eq%tmp(1:n), g => eq%tmp(n+1:2*n), fz => eq%tmp(2*n+1:3*n))
    nnz = 0
    if (warm .and. eq%path == eq_path_ldp) then
        do k = 1, neq
            if (eq_value(work, eq%eq_ids(k)) /= 0.0_wp) nnz = nnz + 1
        end do
    end if

    if (warm .and. eq%path == eq_path_ldp .and. 4*nnz <= neq) then
        if (.not. eq%f_valid) call eq_rhs_f(work)
        do i = 1, n
            eq%xp(i) = eq%xf(i)
            g(i) = eq%gf(i)
        end do
        do c = 1, mr
            eq%sh(c) = eq%df(c)
        end do
        do k = 1, neq
            b = eq_value(work, eq%eq_ids(k))
            if (b == 0.0_wp) cycle
            call eq_rhs_column(work, k)
            associate (col => eq%cols(k)%v)
            do i = 1, n
                eq%xp(i) = eq%xp(i) + b*col(i)
                g(i) = g(i) + b*col(n+i)
            end do
            do c = 1, mr
                eq%sh(c) = eq%sh(c) + b*col(2*n+c)
            end do
            end associate
        end do
        eq%fp = 0.0_wp
        do i = 1, n
            fi = 0.0_wp
            if (qp%has_f) fi = qp%f(i)
            eq%fp = eq%fp + 0.5_wp*eq%xp(i)*(g(i) + fi)
        end do
        do c = 1, mr
            id = eq%keep(c)
            if (qp%bupper(id) >= daqp_inf) then
                eq%other%qp%bupper(c) = daqp_inf
            else
                eq%other%qp%bupper(c) = qp%bupper(id) + eq%sh(c)
            end if
            if (qp%blower(id) <= -daqp_inf) then
                eq%other%qp%blower(c) = -daqp_inf
            else
                eq%other%qp%blower(c) = qp%blower(id) + eq%sh(c)
            end if
        end do
    else
        ! xp = Q1 R^{-T} b_E (scaled by H^{-1/2} in the metric case)
        do k = 1, neq
            rs = ((k-1)*k)/2
            s = eq%s_eq(k)*eq_value(work, eq%eq_ids(k))
            do i = 1, k-1
                s = s - eq%R(rs+i)*y(i)
            end do
            y(k) = s/eq%R(rs+k)
        end do
        do i = 1, n
            if (i <= neq) then
                eq%xp(i) = y(i)
            else
                eq%xp(i) = 0.0_wp
            end if
        end do
        do k = neq, 1, -1
            call eq_reflect(eq%V(:,k), eq%tau(k), k, n, eq%xp)
        end do
        if (eq%metric) then
            do i = 1, n
                eq%xp(i) = eq%xp(i)*eq%dsq(i)
            end do
        end if

        ! g = H xp + f, f(xp) = 0.5 xp'(g + f)
        call eq_hess_times(qp, eq%metric, eq%xp, g)
        if (qp%has_f) then
            do i = 1, n
                g(i) = g(i) + qp%f(i)
            end do
        end if
        eq%fp = 0.0_wp
        do i = 1, n
            fi = 0.0_wp
            if (qp%has_f) fi = qp%f(i)
            eq%fp = eq%fp + 0.5_wp*eq%xp(i)*(g(i) + fi)
        end do

        ! W'g, the linear term of the reduced problem
        fz(1:nz) = 0.0_wp
        do i = 1, n
            gi = g(i)
            if (gi == 0.0_wp) cycle
            do j = 1, nz
                fz(j) = fz(j) + eq%W(j,i)*gi
            end do
        end do
        if (eq%path == eq_path_ldp) then
            ! move xp to the minimizer over the equalities: xp - W W'g
            do i = 1, n
                eq%xp(i) = eq%xp(i) - eq_dot(eq%W(:,i), fz, nz)
            end do
            eq%fp = eq%fp - 0.5_wp*eq_dot(fz, fz, nz)
        else
            do j = 1, nz
                eq%other%qp%f(j) = fz(j)
            end do
        end if

        ! shift the bounds of the kept constraints by their value at xp
        do c = 1, mr
            id = eq%keep(c)
            if (id <= ms) then
                shift = eq%xp(id)
            else
                shift = eq_dot(qp%At(:,id-ms), eq%xp, n)
            end if
            if (qp%bupper(id) >= daqp_inf) then
                eq%other%qp%bupper(c) = daqp_inf
            else
                eq%other%qp%bupper(c) = qp%bupper(id) - shift
            end if
            if (qp%blower(id) <= -daqp_inf) then
                eq%other%qp%blower(c) = -daqp_inf
            else
                eq%other%qp%blower(c) = qp%blower(id) - shift
            end if
        end do
    end if
    ! the constraints that were left out only have to be consistent
    do c = 1, eq%ndrop
        id = eq%drop_ids(c)
        if (id <= ms) then
            val = eq%xp(id)
        else
            val = eq_dot(qp%At(:,id-ms), eq%xp, n)
        end if
        if (qp%bupper(id)-val < -primal_tol .or. qp%blower(id)-val > primal_tol) then
            ! an inconsistent dependent equality is reported the same way as
            ! when it is detected while forming the working set
            if (eq_is_candidate(work,id)) then
                flag = daqp_exit_overdetermined_initial
            else
                flag = daqp_exit_infeasible
            end if
            return
        end if
    end do
    end associate
    end associate

    end function eq_reduce_rhs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Exchange the LDP of the workspace with `work%eq%other` (no copy).

    subroutine eq_swap_ldp(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    associate (d => work%eq%other)
        call swap_problem(work%qp, d%qp)
        call swap_int(work%n, d%n)
        call swap_int(work%m, d%m)
        call swap_int(work%ms, d%ms)
        call swap_r2(work%Mr, d%Mr)
        call swap_r1(work%dupper, d%dupper)
        call swap_r1(work%dlower, d%dlower)
        call swap_r1(work%R, d%R)
        call swap_int(work%rmode, d%rmode)
        call swap_r1(work%v, d%v)
        call swap_log(work%has_v, d%has_v)
        call swap_r1(work%scaling, d%scaling)
        call swap_r1(work%Mu, d%Mu)
        call swap_i1(work%sense, d%sense)
        call swap_log(work%has_weights, d%has_weights)
        call swap_r1(work%rho_ls, d%rho_ls)
        call swap_r1(work%rho_us, d%rho_us)
        call swap_r1(work%w_ls, d%w_ls)
        call swap_r1(work%w_us, d%w_us)
        call swap_int(work%state, d%state)
        call swap_int(work%n_prox, d%n_prox)
        if (allocated(work%bnb)) then
            call swap_i1(work%bnb%bin_ids, d%bin_ids)
            call swap_int(work%bnb%nb, d%nb)
        end if
    end associate

    end subroutine eq_swap_ldp
!*****************************************************************************************

!*****************************************************************************************
!>
!  Swap the reduced problem into the workspace. Returns whether it is installed
!  (by this call).

    logical function daqp_eq_install(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    integer(ip) :: c, id, mr, istat

    daqp_eq_install = .false.
    if (.not. allocated(work%eq)) return
    if (.not. work%eq%active .or. work%eq%installed) return
    ! the soft weights are indexed by the original problem; the reduced
    ! constraints are the original ones (only shifted)
    associate (d => work%eq%other)
        d%has_weights = .false.
        if (work%has_weights) then
            mr = work%eq%mr
            istat = 0
            call resize1(d%rho_ls, mr, istat)
            call resize1(d%rho_us, mr, istat)
            call resize1(d%w_ls, mr, istat)
            call resize1(d%w_us, mr, istat)
            do c = 1, mr
                id = work%eq%keep(c)
                d%rho_ls(c) = work%rho_ls(id)
                d%rho_us(c) = work%rho_us(id)
                d%w_ls(c) = work%w_ls(id)
                d%w_us(c) = work%w_us(id)
            end do
            d%has_weights = .true.
        end if
    end associate
    call eq_swap_ldp(work)
    work%eq%installed = .true.
    daqp_eq_install = .true.

    end function daqp_eq_install
!*****************************************************************************************

!*****************************************************************************************
!>
!  Swap the original problem back into the workspace.

    subroutine daqp_eq_restore(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (.not. allocated(work%eq)) return
    if (.not. work%eq%installed) return
    call eq_swap_ldp(work)
    work%eq%installed = .false.

    end subroutine daqp_eq_restore
!*****************************************************************************************

!*****************************************************************************************
!>
!  Give up the reduction (the LDP of the original problem then has to be formed).

    subroutine daqp_eq_deactivate(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (.not. allocated(work%eq)) return
    call free_daqp_eq(work)
    ! the working set refers to the reduced problem
    call daqp_reset_workspace(work)
    if (allocated(work%bnb)) work%bnb%n_root_ws = 0

    end subroutine daqp_eq_deactivate
!*****************************************************************************************

!*****************************************************************************************
!>
!  Form or update the reduction of the problem in the workspace, and the LDP
!  of the reduced problem (with `update_ldp`). Returns `eq_not_reduced` if
!  nothing can be eliminated, a negative exit flag if the equality constraints
!  cannot be satisfied, and the return value of `update_ldp` otherwise.

    integer(ip) function daqp_eq_update(work, mask, update_ldp) result(flag)

    type(daqp_workspace), intent(inout) :: work !! workspace
    integer(ip), intent(in) :: mask             !! update mask
    procedure(update_ldp_proc) :: update_ldp    !! forms the LDP of the installed problem

    integer(ip) :: mask_r, c, nb, keep_mask
    logical :: rebuild, installed

    keep_mask = iand(mask, daqp_update_unconstrained+daqp_update_eliminate)
    call daqp_eq_restore(work)
    if (.not. allocated(work%eq)) allocate(work%eq)

    rebuild = .not. work%eq%active
    if (.not. rebuild) rebuild = iand(mask, daqp_update_rinv+daqp_update_m) /= 0 .or. &
        work%eq%n /= work%qp%n .or. work%eq%m /= work%qp%m .or. work%eq%ms /= work%qp%ms
    if (.not. rebuild) rebuild = .not. eq_same_candidates(work)

    if (rebuild) then
        flag = eq_build_reduction(work)
        if (flag <= 0) then
            call daqp_eq_deactivate(work)
            if (flag == 0) flag = eq_not_reduced
            return
        end if
        work%eq%active = .true.
        ! (Rinv also for an LP, for which it marks the proximal directions)
        mask_r = daqp_update_rinv + daqp_update_m + daqp_update_d + daqp_update_sense
        if (work%eq%other%qp%has_f) mask_r = ior(mask_r, daqp_update_v)
    else
        mask_r = daqp_update_d
        if (work%eq%other%qp%has_f) mask_r = ior(mask_r, daqp_update_v) ! depends on xp
        if (has(mask,daqp_update_sense)) then
            do c = 1, work%eq%mr
                work%eq%other%qp%sense(c) = work%sense(work%eq%keep(c))
            end do
            mask_r = ior(mask_r, daqp_update_sense)
        end if
    end if
    mask_r = ior(mask_r, keep_mask)

    if (has(mask,daqp_update_v)) work%eq%f_valid = .false. ! f has changed
    flag = eq_reduce_rhs(work, .not. rebuild)
    if (flag < 0) then
        ! formed by the next update; reported by a solve before that
        work%eq%other%state = ior(work%eq%other%state, iand(mask_r, state_pending))
        work%eq%error = flag
        return
    end if
    work%eq%error = 0

    installed = daqp_eq_install(work)
    flag = update_ldp(mask_r, work)
    ! a singular reduced Hessian needs v for the proximal linear term
    if (flag >= 0 .and. work%n_prox > 0 .and. .not. work%has_v .and. work%qp%has_H) then
        work%has_v = .true.
        if (allocated(work%v)) deallocate(work%v)
        allocate(work%v(work%n))
        work%v = 0.0_wp
    end if
    if (allocated(work%bnb)) then
        nb = 0
        do c = 1, work%m
            if (has(work%sense(c),daqp_binary)) then
                nb = nb + 1
                work%bnb%bin_ids(nb) = c
            end if
        end do
        work%bnb%nb = nb
    end if
    call daqp_eq_restore(work)

    end function daqp_eq_update
!*****************************************************************************************

!*****************************************************************************************
!>
!  Multipliers of the eliminated equalities, from the stationarity condition
!  `A_E' lam_E = -(H x + f + A_K' lam_K)`.

    subroutine eq_compute_lam_eq(work, x, lam)

    type(daqp_workspace), intent(inout) :: work !! workspace (with the reduced problem installed)
    real(wp), intent(in) :: x(:)                !! solution of the original problem
    real(wp), intent(inout) :: lam(:)           !! multipliers

    integer(ip) :: n, ms, neq, i, k, c, id, rs
    real(wp) :: l, mu

    associate (eq => work%eq, qp => work%eq%other%qp)
    n = eq%n; ms = eq%ms; neq = eq%neq
    associate (g => eq%tmp(1:n))
    call eq_hess_times(qp, eq%metric, x, g)
    if (qp%has_f) then
        do i = 1, n
            g(i) = g(i) + qp%f(i)
        end do
    end if
    do c = 1, eq%mr
        id = eq%keep(c)
        l = lam(id)
        if (l == 0.0_wp) cycle
        if (id <= ms) then
            g(id) = g(id) + l
        else
            do i = 1, n
                g(i) = g(i) + l*qp%At(i,id-ms)
            end do
        end if
    end do
    if (eq%metric) then
        do i = 1, n
            g(i) = g(i)*eq%dsq(i)
        end do
    end if
    do k = 1, neq
        call eq_reflect(eq%V(:,k), eq%tau(k), k, n, g)
    end do
    do k = neq, 1, -1 ! R mu = -Q1'g
        rs = ((k-1)*k)/2
        mu = -g(k)/eq%R(rs+k)
        g(k) = mu
        do i = 1, k-1
            g(i) = g(i) + eq%R(rs+i)*mu
        end do
    end do
    do k = 1, neq
        lam(eq%eq_ids(k)) = eq%s_eq(k)*g(k)
    end do
    end associate
    end associate

    end subroutine eq_compute_lam_eq
!*****************************************************************************************

!*****************************************************************************************
!>
!  Turn the result of the installed reduced problem into the result of the
!  original problem, and carry its working set over to the original constraints.

    subroutine daqp_eq_expand(work, x, fval, lam)

    type(daqp_workspace), intent(inout) :: work !! workspace (with the reduced problem installed)
    real(wp), intent(inout) :: x(:)             !! solution (of size n of the original problem)
    real(wp), intent(inout) :: fval             !! objective (of the reduced problem on input)
    real(wp), intent(inout), optional :: lam(:) !! multipliers (of the reduced problem on input)

    integer(ip) :: i, c, n, id, state_mask
    real(wp) :: fval_r, l

    if (.not. allocated(work%eq)) return
    if (.not. work%eq%installed) return
    state_mask = daqp_active + daqp_lower + daqp_slack_fixed
    associate (eq => work%eq, qp => work%eq%other%qp)
    n = eq%n
    associate (xx => eq%tmp(2*n+1:3*n))
    ! x = xp + W w
    do i = 1, n
        xx(i) = eq%xp(i) + eq_dot(eq%W(:,i), work%x, eq%nz)
    end do
    x(1:n) = xx(1:n)

    ! objective (the reduced problem of eq_path_ldp has no linear term)
    fval_r = fval
    if (eq%path == eq_path_ldp .and. .not. work%has_v) fval_r = 0.5_wp*work%fval
    fval = fval_r + eq%fp
    if (is_avi_nonsym(work) .and. qp%has_f) fval = eq_dot(qp%f, xx, n)

    ! multipliers of the original constraints (keep is ascending, so they can
    ! be scattered in place from the end)
    if (present(lam) .and. .not. is_hierarchical(work)) then
        do i = eq%m, eq%mr+1, -1
            lam(i) = 0.0_wp
        end do
        do c = eq%mr, 1, -1
            l = lam(c)
            lam(c) = 0.0_wp
            lam(eq%keep(c)) = l
        end do
        call eq_compute_lam_eq(work, xx, lam)
    end if
    end associate

    ! carry the working set over to the original constraints
    do c = 1, eq%mr
        id = eq%keep(c)
        eq%other%sense(id) = ior(iand(eq%other%sense(id), not(state_mask)), &
                                 iand(work%sense(c), state_mask))
    end do
    end associate

    end subroutine daqp_eq_expand
!*****************************************************************************************

!*****************************************************************************************
!>
!  Set the starting iterate of the reduced problem from `x` of the original one.

    subroutine daqp_eq_set_primal_start(work, x)

    type(daqp_workspace), intent(inout) :: work !! workspace
    real(wp), intent(in) :: x(:)                !! iterate of the original problem

    integer(ip) :: n, nz, i, j
    real(wp) :: hi

    associate (eq => work%eq)
    n = eq%n; nz = eq%nz
    associate (t => eq%tmp(1:n), ht => eq%tmp(n+1:2*n))
    ! w = W'H(x-xp) if W'HW = I, and w = W'(x-xp) if W is orthonormal
    do i = 1, n
        t(i) = x(i) - eq%xp(i)
    end do
    if (eq%path == eq_path_ldp) then
        if (eq%installed) then
            call eq_hess_times(eq%other%qp, eq%metric, t, ht)
        else
            call eq_hess_times(work%qp, eq%metric, t, ht)
        end if
    else
        ht(1:n) = t(1:n)
    end if
    work%x(1:nz) = 0.0_wp
    do i = 1, n
        hi = ht(i)
        if (hi == 0.0_wp) cycle
        do j = 1, nz
            work%x(j) = work%x(j) + eq%W(j,i)*hi
        end do
    end do
    end associate
    end associate

    end subroutine daqp_eq_set_primal_start
!*****************************************************************************************

!*****************************************************************************************
!>
!  Free the elimination (restoring the original problem first).

    subroutine free_daqp_eq(work)

    type(daqp_workspace), intent(inout) :: work !! workspace

    if (.not. allocated(work%eq)) return
    call daqp_eq_restore(work)
    deallocate(work%eq)

    end subroutine free_daqp_eq
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_eq_elim
!*****************************************************************************************
