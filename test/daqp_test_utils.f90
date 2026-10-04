!*****************************************************************************************
!>
!  Utilities for the tests and the comparison with upstream: a portable
!  random number generator, generators of random QPs and LPs with a known
!  solution (ports of upstream's `generate_test_QP` and `generate_test_LP`),
!  KKT residuals, and a small test reporter.

    module daqp_test_utils

    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use, intrinsic :: iso_fortran_env, only: int64, output_unit

    implicit none

    private

    integer(int64), save :: rng_state = 12345_int64 !! state of the generator
    integer, save :: n_checks = 0   !! number of checks
    integer, save :: n_failed = 0   !! number of failed checks

    public :: wp, ip
    public :: rng_seed, urand, randn, rand_int, shuffle, orthogonal
    public :: generate_qp, generate_lp, generate_miqp, generate_avi, make_equalities, cholesky_packed
    public :: kkt_residuals, kkt_tolerance
    public :: check, report

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Seed the generator.

    subroutine rng_seed(seed)

    integer, intent(in) :: seed !! seed (any value)

    rng_state = mod(abs(int(seed,int64)), 2147483646_int64) + 1_int64

    end subroutine rng_seed
!*****************************************************************************************

!*****************************************************************************************
!>
!  Uniform random number in (0,1) (Park-Miller minimal standard generator,
!  identical on every compiler).

    real(wp) function urand()

    rng_state = mod(48271_int64*rng_state, 2147483647_int64)
    urand = real(rng_state,wp)/2147483647.0_wp

    end function urand
!*****************************************************************************************

!*****************************************************************************************
!>
!  Standard normal random number (Box-Muller).

    real(wp) function randn()

    real(wp) :: u1, u2

    u1 = urand()
    u2 = urand()
    randn = sqrt(-2.0_wp*log(u1))*cos(2.0_wp*acos(-1.0_wp)*u2)

    end function randn
!*****************************************************************************************

!*****************************************************************************************
!>
!  Random integer in `lo..hi`.

    integer function rand_int(lo,hi)

    integer, intent(in) :: lo !! lower limit
    integer, intent(in) :: hi !! upper limit

    rand_int = lo + min(hi-lo, int(urand()*real(hi-lo+1,wp)))

    end function rand_int
!*****************************************************************************************

!*****************************************************************************************
!>
!  A random permutation of `1..n` (Fisher-Yates).

    subroutine shuffle(p)

    integer, intent(out) :: p(:) !! the permutation

    integer :: i, j, t

    p = [(i, i=1,size(p))]
    do i = size(p), 2, -1
        j = rand_int(1,i)
        t = p(i); p(i) = p(j); p(j) = t
    end do

    end subroutine shuffle
!*****************************************************************************************

!*****************************************************************************************
!>
!  A random orthogonal matrix (modified Gram-Schmidt, twice, of a normal matrix).

    subroutine orthogonal(Q)

    real(wp), intent(out) :: Q(:,:) !! square matrix

    integer :: i, j, k, n, pass
    real(wp) :: r

    n = size(Q,1)
    do j = 1, n
        do i = 1, n
            Q(i,j) = randn()
        end do
    end do
    do j = 1, n
        do pass = 1, 2
            do k = 1, j-1
                r = dot_product(Q(:,k), Q(:,j))
                Q(:,j) = Q(:,j) - r*Q(:,k)
            end do
        end do
        Q(:,j) = Q(:,j)/norm2(Q(:,j))
    end do

    end subroutine orthogonal
!*****************************************************************************************

!*****************************************************************************************
!>
!  A random QP with a known solution and `nact` active constraints
!  (port of upstream's `generate_test_QP`): `cond(H) = kappa`.

    subroutine generate_qp(n, m, ms, nact, kappa, xref, H, f, A, bupper, blower)

    integer, intent(in) :: n      !! number of variables (>= 2)
    integer, intent(in) :: m      !! number of constraints
    integer, intent(in) :: ms     !! number of simple bounds
    integer, intent(in) :: nact   !! number of active constraints at the solution
    real(wp), intent(in) :: kappa !! condition number of `H`
    real(wp), allocatable, intent(out) :: xref(:)   !! the solution
    real(wp), allocatable, intent(out) :: H(:,:)    !! Hessian
    real(wp), allocatable, intent(out) :: f(:)      !! linear term
    real(wp), allocatable, intent(out) :: A(:,:)    !! general constraints `(m-ms,n)`
    real(wp), allocatable, intent(out) :: bupper(:) !! upper bounds
    real(wp), allocatable, intent(out) :: blower(:) !! lower bounds

    real(wp), allocatable :: eigs(:), Q(:,:), T(:,:), Tinv(:,:), Mm(:,:), Ma(:,:), &
                             lam(:), da(:), u(:), v(:), dupper(:), dlower(:)
    integer, allocatable :: perm(:)
    integer :: i, nu, nl, id

    allocate(eigs(n), Q(n,n), T(n,n), Tinv(n,n), Mm(m,n), perm(m), dupper(m), dlower(m), &
             u(n), v(n), lam(nact), Ma(nact,n), da(nact))

    ! H = T'T with cond(H) = kappa
    eigs(1) = 1.0_wp
    eigs(2) = kappa
    do i = 3, n
        eigs(i) = 1.0_wp + (kappa-1.0_wp)*urand()
    end do
    call orthogonal(Q)
    do i = 1, n
        T(i,:) = sqrt(eigs(i))*Q(:,i)
        Tinv(:,i) = Q(:,i)/sqrt(eigs(i))
    end do
    H = matmul(transpose(T), T)
    H = 0.5_wp*(H + transpose(H))

    ! the LDP: min ||u|| s.t. dlower <= M u <= dupper
    Mm(1:ms,:) = Tinv(1:ms,:)
    do i = ms+1, m
        Mm(i,:) = [(randn(), id=1,n)]
    end do
    dupper = 0.0_wp
    dlower = 0.0_wp
    call shuffle(perm)
    nu = rand_int(0,nact)
    nl = nact - nu

    ! active bounds with lam >= 0
    do i = 1, nact
        lam(i) = urand()
        if (i <= nu) then
            Ma(i,:) = Mm(perm(i),:)
        else
            Ma(i,:) = -Mm(perm(i),:)
        end if
    end do
    da = -matmul(Ma, matmul(transpose(Ma), lam))
    u = -matmul(transpose(Ma), lam)
    do i = 1, nact
        id = perm(i)
        if (i <= nu) then
            dupper(id) = da(i)
            dlower(id) = dupper(id) - (0.01_wp + urand())
        else
            dlower(id) = -da(i)
            dupper(id) = dlower(id) + (0.01_wp + urand())
        end if
    end do
    ! inactive constraints are feasible
    do i = nact+1, m
        id = perm(i)
        dupper(id) = dot_product(Mm(id,:),u) + (0.01_wp + urand())
        dlower(id) = dot_product(Mm(id,:),u) - (0.01_wp + urand())
    end do

    ! back to the QP: x = T\(u-v), f = T'v, A = M T, b = d - M v
    v = [(randn(), i=1,n)]
    f = matmul(transpose(T), v)
    xref = matmul(Tinv, u-v)
    A = matmul(Mm(ms+1:m,:), T)
    bupper = dupper - matmul(Mm, v)
    blower = dlower - matmul(Mm, v)

    end subroutine generate_qp
!*****************************************************************************************

!*****************************************************************************************
!>
!  A random LP with a known solution (port of upstream's `generate_test_LP`).

    subroutine generate_lp(n, m, ms, xref, f, A, bupper, blower)

    integer, intent(in) :: n   !! number of variables
    integer, intent(in) :: m   !! number of constraints (>= n)
    integer, intent(in) :: ms  !! number of simple bounds
    real(wp), allocatable, intent(out) :: xref(:)   !! the solution
    real(wp), allocatable, intent(out) :: f(:)      !! linear term
    real(wp), allocatable, intent(out) :: A(:,:)    !! general constraints `(m-ms,n)`
    real(wp), allocatable, intent(out) :: bupper(:) !! upper bounds
    real(wp), allocatable, intent(out) :: blower(:) !! lower bounds

    real(wp), allocatable :: Af(:,:), Aa(:,:), lam(:), ba(:)
    integer, allocatable :: perm(:)
    integer :: i, j, nu, id

    allocate(Af(m,n), Aa(n,n), lam(n), ba(n), perm(m), bupper(m), blower(m))
    Af = 0.0_wp
    do i = 1, ms
        Af(i,i) = 1.0_wp
    end do
    do i = ms+1, m
        Af(i,:) = [(randn(), j=1,n)]
    end do
    bupper = 0.0_wp
    blower = 0.0_wp
    call shuffle(perm)
    nu = rand_int(1,n+1) - 1
    do i = 1, n
        lam(i) = urand()
    end do
    xref = [(randn(), i=1,n)]
    do i = 1, n
        if (i <= nu) then
            Aa(i,:) = Af(perm(i),:)
        else
            Aa(i,:) = -Af(perm(i),:)
        end if
    end do
    f = -matmul(transpose(Aa), lam)
    ba = matmul(Aa, xref)
    do i = 1, n
        id = perm(i)
        if (i <= nu) then
            bupper(id) = ba(i)
            blower(id) = bupper(id) - (0.01_wp + urand())
        else
            blower(id) = -ba(i)
            bupper(id) = blower(id) + (0.01_wp + urand())
        end if
    end do
    do i = n+1, m
        id = perm(i)
        bupper(id) = dot_product(Af(id,:),xref) + (0.01_wp + urand())
        blower(id) = dot_product(Af(id,:),xref) - (0.01_wp + urand())
    end do
    A = Af(ms+1:m,:)

    end subroutine generate_lp
!*****************************************************************************************

!*****************************************************************************************
!>
!  A mixed-integer QP whose first `nb` variables are binary (port of
!  upstream's `generate_test_MIQP`): the origin is feasible, `f` makes it
!  lucrative to leave it, and a cardinality constraint `sum(x(1:nb)) <= nb/2`
!  makes the relaxation fractional. Requires `ms >= nb` and `m > ms`.

    subroutine generate_miqp(n, m, ms, nb, H, f, A, bupper, blower, sense)

    integer, intent(in) :: n, m, ms, nb  !! dimensions, number of binary variables
    real(wp), allocatable, intent(out) :: H(:,:), f(:), A(:,:), bupper(:), blower(:)
    integer(ip), allocatable, intent(out) :: sense(:)

    real(wp), allocatable :: Mt(:,:)
    integer :: i, j

    allocate(Mt(n,n), A(m-ms,n), bupper(m), blower(m), f(n), sense(m))
    do j = 1, n
        do i = 1, n
            Mt(i,j) = randn()
        end do
    end do
    H = matmul(transpose(Mt), Mt)
    do i = 1, n
        H(i,i) = H(i,i) + 1.0_wp
    end do
    do j = 1, n
        do i = 1, m-ms
            A(i,j) = randn()
        end do
    end do
    do i = 1, m
        bupper(i) = 20.0_wp*urand()
    end do
    do i = 1, m
        blower(i) = -20.0_wp*urand()
    end do
    do i = 1, n
        f(i) = 100.0_wp*randn()
    end do
    f(1:nb) = -abs(f(1:nb))
    bupper(1:nb) = 1.0_wp
    blower(1:nb) = 0.0_wp
    sense = 0
    sense(1:nb) = 16 ! binary
    A(1,:) = 0.0_wp
    A(1,1:nb) = 1.0_wp
    bupper(ms+1) = real(nb/2, wp)
    blower(ms+1) = -1.0e30_wp

    end subroutine generate_miqp
!*****************************************************************************************

!*****************************************************************************************
!>
!  An affine variational inequality with a known solution (port of upstream's
!  `generate_test_avi`): `H` positive definite but not symmetric, constraints
!  `A x <= b` (no simple bounds).

    subroutine generate_avi(n, m, xref, H, f, A, b)

    integer, intent(in) :: n, m  !! dimensions
    real(wp), allocatable, intent(out) :: xref(:), H(:,:), f(:), A(:,:), b(:)

    real(wp), allocatable :: Mt(:,:), Nt(:,:), sym(:,:), asym(:,:), lam(:), ax(:)
    integer, allocatable :: perm(:)
    integer :: i, j, nas

    allocate(A(m,n), Mt(n,n), Nt(n,n), lam(m), perm(m), b(m))
    do j = 1, n
        do i = 1, m
            A(i,j) = randn()
        end do
    end do
    call shuffle(perm)
    nas = rand_int(1,n+1) - 1
    lam = 0.0_wp
    do i = 1, nas
        lam(perm(i)) = urand()
    end do
    xref = [(randn(), i=1,n)]
    do j = 1, n
        do i = 1, n
            Mt(i,j) = urand()
        end do
    end do
    do j = 1, n
        do i = 1, n
            Nt(i,j) = randn()
        end do
    end do
    sym = matmul(transpose(Mt), Mt)
    asym = Nt - transpose(Nt)
    H = sym/norm2(sym) + asym/norm2(asym)
    f = -matmul(H, xref) - matmul(transpose(A), lam)
    ax = matmul(A, xref)
    do i = 1, m
        b(i) = ax(i) + 5.0_wp*urand()
    end do
    do i = 1, nas
        b(perm(i)) = ax(perm(i))
    end do

    end subroutine generate_avi
!*****************************************************************************************

!*****************************************************************************************
!>
!  Turn `neq` constraints that are active at `xref` into equalities (as
!  upstream's `generate_test_QP_eq`: general constraints first): their bounds
!  are collapsed onto the active one. Returns the number turned.

    integer function make_equalities(ms, neq, xref, A, bupper, blower) result(k)

    integer, intent(in) :: ms, neq    !! number of simple bounds, of equalities wanted
    real(wp), intent(in) :: xref(:)   !! solution
    real(wp), intent(in) :: A(:,:)    !! general constraints
    real(wp), intent(inout) :: bupper(:), blower(:) !! bounds

    real(wp), allocatable :: r(:)
    integer :: i, pass, m
    real(wp) :: bv

    m = size(bupper)
    r = [xref(1:ms), matmul(A, xref)]
    k = 0
    do pass = 1, 2 ! general constraints, then simple bounds
        do i = 1, m
            if (k >= neq) return
            if ((pass == 1 .and. i <= ms) .or. (pass == 2 .and. i > ms)) cycle
            if (min(abs(r(i)-bupper(i)), abs(r(i)-blower(i))) <= 1.0e-6_wp*(1.0_wp+abs(r(i)))) then
                if (abs(r(i)-bupper(i)) <= abs(r(i)-blower(i))) then
                    bv = bupper(i)
                else
                    bv = blower(i)
                end if
                bupper(i) = bv
                blower(i) = bv
                k = k + 1
            end if
        end do
    end do

    end function make_equalities
!*****************************************************************************************

!*****************************************************************************************
!>
!  The upper Cholesky factor `R` of `H = R'R`, packed by rows.

    function cholesky_packed(H) result(Rp)

    real(wp), intent(in) :: H(:,:)  !! symmetric positive definite matrix
    real(wp), allocatable :: Rp(:)

    real(wp), allocatable :: R(:,:)
    integer :: i, j, k, n

    n = size(H,1)
    allocate(R(n,n), Rp((n*(n+1))/2))
    R = 0.0_wp
    do j = 1, n
        do i = 1, j
            R(i,j) = H(i,j) - dot_product(R(1:i-1,i), R(1:i-1,j))
            if (i < j) then
                R(i,j) = R(i,j)/R(i,i)
            else
                R(j,j) = sqrt(R(j,j))
            end if
        end do
    end do
    k = 0
    do i = 1, n
        do j = i, n
            k = k + 1
            Rp(k) = R(i,j)
        end do
    end do

    end function cholesky_packed
!*****************************************************************************************

!*****************************************************************************************
!>
!  KKT residuals of `x` and `lam` (DAQP's sign convention: `lam > 0` at an
!  upper bound, `lam < 0` at a lower bound), as infinity norms:
!  stationarity `Hx + f + Aall'lam`, primal infeasibility, wrong-signed
!  multipliers, and complementarity.

    subroutine kkt_residuals(ms, bupper, blower, x, lam, stat, prim, dual, comp, H, f, A)

    integer, intent(in) :: ms           !! number of simple bounds
    real(wp), intent(in) :: bupper(:)   !! upper bounds
    real(wp), intent(in) :: blower(:)   !! lower bounds
    real(wp), intent(in) :: x(:)        !! solution
    real(wp), intent(in) :: lam(:)      !! multipliers
    real(wp), intent(out) :: stat       !! stationarity
    real(wp), intent(out) :: prim       !! primal infeasibility
    real(wp), intent(out) :: dual       !! dual infeasibility
    real(wp), intent(out) :: comp       !! complementarity
    real(wp), intent(in), optional :: H(:,:) !! Hessian
    real(wp), intent(in), optional :: f(:)   !! linear term
    real(wp), intent(in), optional :: A(:,:) !! general constraints

    real(wp), allocatable :: g(:), ax(:)
    integer :: i, m, n
    real(wp) :: big

    n = size(x)
    m = size(bupper)
    big = 1.0e29_wp
    allocate(g(n), ax(m))
    g = 0.0_wp
    if (present(H)) g = matmul(H, x)
    if (present(f)) g = g + f
    g(1:ms) = g(1:ms) + lam(1:ms)
    ax(1:ms) = x(1:ms)
    if (present(A)) then
        if (size(A,1) > 0) then
            g = g + matmul(transpose(A), lam(ms+1:m))
            ax(ms+1:m) = matmul(A, x)
        end if
    end if
    stat = maxval(abs(g))
    prim = 0.0_wp
    dual = 0.0_wp
    comp = 0.0_wp
    do i = 1, m
        if (bupper(i) < big) prim = max(prim, ax(i)-bupper(i))
        if (blower(i) > -big) prim = max(prim, blower(i)-ax(i))
        if (lam(i) > 0.0_wp) then
            if (bupper(i) >= big) then
                dual = max(dual, lam(i))
            else
                comp = max(comp, lam(i)*abs(bupper(i)-ax(i)))
            end if
        else if (lam(i) < 0.0_wp) then
            if (blower(i) <= -big) then
                dual = max(dual, -lam(i))
            else
                comp = max(comp, -lam(i)*abs(ax(i)-blower(i)))
            end if
        end if
    end do

    end subroutine kkt_residuals
!*****************************************************************************************

!*****************************************************************************************
!>
!  A tolerance for the KKT residuals: the larger of the solver's primal
!  tolerance and a multiple of `epsilon` scaled by the condition number and
!  the size of the data.

    pure real(wp) function kkt_tolerance(kappa, scale, primal_tol)

    real(wp), intent(in) :: kappa       !! condition number of the problem
    real(wp), intent(in) :: scale       !! size of the data (norm of the solution, multipliers, ...)
    real(wp), intent(in) :: primal_tol  !! the solver's primal tolerance

    kkt_tolerance = max(10.0_wp*primal_tol, 1.0e4_wp*epsilon(1.0_wp)*kappa) * (1.0_wp + scale)

    end function kkt_tolerance
!*****************************************************************************************

!*****************************************************************************************
!>
!  Record a check.

    subroutine check(ok, msg)

    logical, intent(in) :: ok           !! the check passed
    character(len=*), intent(in) :: msg !! description

    n_checks = n_checks + 1
    if (.not. ok) then
        n_failed = n_failed + 1
        write(output_unit,'(A)') 'FAILED: '//msg
    end if

    end subroutine check
!*****************************************************************************************

!*****************************************************************************************
!>
!  Report the checks, and stop with an error if any failed.

    subroutine report(name)

    character(len=*), intent(in) :: name !! name of the test

    write(output_unit,'(A,": ",I0," checks, ",I0," failed")') name, n_checks, n_failed
    if (n_failed > 0) error stop 1

    end subroutine report
!*****************************************************************************************

!*****************************************************************************************
    end module daqp_test_utils
!*****************************************************************************************
