!*****************************************************************************************
!>
!  Compare the Fortran port of DAQP with the upstream C library, for accuracy
!  and speed. Writes a Markdown report to standard output (see
!  `tools/run_compare.sh`, which adds the environment and writes
!  `compare/RESULTS.md`).
!
!  Usage: `compare [quick|accuracy]` (`quick`: a small subset; `accuracy`:
!  no timing).

    program compare

    use iso_c_binding
    use iso_fortran_env, only: int64, real128, output_unit
    use daqp_kinds, only: wp => daqp_wp, ip => daqp_ip
    use daqp_types
    use daqp_core
    use daqp_c_binding
    use daqp_test_utils, only: rng_seed, urand, randn, generate_qp, generate_lp, rand_int, &
                               generate_miqp, generate_avi, make_equalities, cholesky_packed

    implicit none

    type :: qp_problem
        !! A problem in the Fortran layout, plus its upstream (row-major) copy,
        !! and how it is to be set up.
        character(len=:), allocatable :: name
        integer :: n = 0, m = 0, ms = 0
        logical :: has_H = .true.
        logical :: has_f = .true.
        logical :: factored = .false.     !! pass the Cholesky factor instead of H
        logical :: check_kkt = .true.     !! compare the KKT residuals (QPs and LPs)
        integer(ip) :: problem_type = daqp_problem_qp
        integer(ip) :: init_mask = 0      !! 64+128: as upstream's daqp_quadprog
        integer(ip) :: eq_reduction = daqp_eq_reduction_auto
        real(wp) :: kappa = 1.0_wp
        real(wp), allocatable :: H(:,:), f(:), A(:,:), bu(:), bl(:)
        real(wp), allocatable :: At(:,:)   !! transpose of A (the internal layout of the port)
        real(wp), allocatable :: Rf(:)     !! packed Cholesky factor (factored)
        integer(ip), allocatable :: sense(:)
        integer(ip), allocatable :: break_points(:)
        real(wp), allocatable :: primal_start(:), dual_start(:)
        real(wp), allocatable :: rho_l(:), rho_u(:), w_l(:), w_u(:)
        ! row-major copies for C (kept alive while a C workspace points to them)
        real(c_double), allocatable :: Hc(:), Ac(:), fc(:), buc(:), blc(:)
        integer(c_int), allocatable :: sensec(:), bpc(:)
        real(c_double), allocatable :: rho_lc(:), rho_uc(:), w_lc(:), w_uc(:)
    end type qp_problem

    type :: set_stats
        !! Accumulated results of a problem set.
        character(len=:), allocatable :: name
        integer :: n_problems = 0
        integer :: n_flag_equal = 0
        integer :: n_iter_equal = 0
        integer :: n_ws_equal = 0
        integer :: n_x_pass = 0
        integer :: n_lam_pass = 0
        integer :: n_f_pass = 0
        integer :: n_kkt_pass = 0
        integer :: n_solved = 0
        real(wp) :: max_dx = 0.0_wp, max_dlam = 0.0_wp, max_df = 0.0_wp
        character(len=:), allocatable :: mismatches
    end type set_stats

    type(daqp_c_settings), target :: csettings, csettings0
    type(daqp_settings) :: fsettings, fsettings0
    type(set_stats) :: st
    logical :: quick, accuracy_only
    character(len=32) :: arg
    real(wp), parameter :: c_tol = 10.0_wp  !! the constant `c` of the accuracy criteria

    quick = .false.
    accuracy_only = .false.
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        quick = trim(arg) == 'quick'
        accuracy_only = trim(arg) == 'accuracy'
    end if

    ! the same settings for both: the C defaults, copied field by field
    call daqp_default_settings(csettings0)
    fsettings0%primal_tol   = csettings0%primal_tol
    fsettings0%dual_tol     = csettings0%dual_tol
    fsettings0%zero_tol     = csettings0%zero_tol
    fsettings0%pivot_tol    = csettings0%pivot_tol
    fsettings0%progress_tol = csettings0%progress_tol
    fsettings0%cycle_tol    = csettings0%cycle_tol
    fsettings0%iter_limit   = csettings0%iter_limit
    fsettings0%fval_bound   = csettings0%fval_bound
    fsettings0%eps_prox     = csettings0%eps_prox
    fsettings0%eta_prox     = csettings0%eta_prox
    fsettings0%rho_soft     = csettings0%rho_soft
    fsettings0%rel_subopt   = csettings0%rel_subopt
    fsettings0%abs_subopt   = csettings0%abs_subopt
    fsettings0%sing_tol     = csettings0%sing_tol
    fsettings0%refactor_tol = csettings0%refactor_tol
    fsettings0%time_limit   = csettings0%time_limit
    fsettings0%w_soft       = csettings0%w_soft
    fsettings0%eq_reduction = csettings0%eq_reduction
    csettings = csettings0
    fsettings = fsettings0

    write(output_unit,'(A)') '## Accuracy'
    write(output_unit,'(A)') ''
    write(output_unit,'(A)') 'Each problem is set up and solved by both solvers through the split interface'// &
        ' (`setup_daqp_main` + `daqp_solve`, and `daqp_setup` + `daqp_solve`), with the same settings'// &
        ' (upstream''s defaults, field by field) and the same options (`init_mask`, `eq_reduction`,'// &
        ' break points, problem type, primal and dual starts, soft weights).'
    write(output_unit,'(A,ES8.1,A)') 'Pass criteria: `|x_F - x_C|_inf <= c n eps kappa (1 + |x_C|_inf)`'// &
        ' (kappa: the larger of cond(H) and the condition number of the final working set''s Gram'// &
        ' matrix, estimated by 1/min pivot of its LDL'' factors, and 1/rho_soft with soft constraints), the same for the multipliers and'// &
        ' `|f_F - f_C| <= c n eps kappa (1 + |f_C|)`, with c = ', c_tol, &
        '; the KKT residuals of the port (computed in quadruple precision) at most c times those'// &
        ' of the C code (or below the solver''s primal tolerance; QPs and LPs only).'
    write(output_unit,'(A)') ''
    call table_header()

    call set_upstream(); call print_set(st)
    call set_random(); call print_set(st)
    call set_degenerate(); call print_set(st)
    call set_special(); call print_set(st)
    call set_warm(); call print_set(st)
    call set_eq(); call print_set(st)
    call set_bnb(); call print_set(st)
    call set_hierarchical(); call print_set(st)
    call set_avi(); call print_set(st)
    call set_options(); call print_set(st)
    call set_minrep(); call print_set(st)

    if (.not. accuracy_only) call speed()

    contains
!*****************************************************************************************

    subroutine table_header()
        write(output_unit,'(A)') '| Set | Problems | Solved | Same exit flag | Same iterations | Same working set |'// &
            ' x pass | lam pass | f pass | KKT pass | max dx | max dlam | max df |'
        write(output_unit,'(A)') '|---|---|---|---|---|---|---|---|---|---|---|---|---|'
    end subroutine table_header

    subroutine new_set(name)
        character(len=*), intent(in) :: name
        st = set_stats()
        st%name = name
        st%mismatches = ''
    end subroutine new_set

    subroutine print_set(s)
        type(set_stats), intent(in) :: s
        character(len=1024) :: line
        write(line,'("| ",A," | ",I0," | ",I0," | ",I0," | ",I0," | ",I0," | ",I0," | ",I0," | ",I0,'// &
                   '" | ",I0," | ",ES8.1," | ",ES8.1," | ",ES8.1," |")') &
            s%name, s%n_problems, s%n_solved, s%n_flag_equal, s%n_iter_equal, s%n_ws_equal, &
            s%n_x_pass, s%n_lam_pass, s%n_f_pass, s%n_kkt_pass, s%max_dx, s%max_dlam, s%max_df
        write(output_unit,'(A)') trim(line)
        if (len(s%mismatches) > 0) then
            write(output_unit,'(A)') ''
            write(output_unit,'(A)') '<details><summary>'//s%name//': differences</summary>'
            write(output_unit,'(A)') ''
            write(output_unit,'(A)') s%mismatches
            write(output_unit,'(A)') '</details>'
            write(output_unit,'(A)') ''
            call table_header()
        end if
    end subroutine print_set

    subroutine pass_all()
        !! count a problem that neither solver solved (the same exit flag) as passed
        st%n_x_pass = st%n_x_pass + 1
        st%n_lam_pass = st%n_lam_pass + 1
        st%n_f_pass = st%n_f_pass + 1
        st%n_kkt_pass = st%n_kkt_pass + 1
    end subroutine pass_all

    subroutine add_mismatch(msg)
        character(len=*), intent(in) :: msg
        if (len(st%mismatches) < 20000) st%mismatches = st%mismatches//'- '//trim(msg)//new_line('a')
    end subroutine add_mismatch

!*****************************************************************************************
!>
!  Make the row-major copies of a problem for C.

    subroutine finalize_problem(p)
        type(qp_problem), intent(inout) :: p
        integer :: i
        if (.not. allocated(p%A)) allocate(p%A(0,p%n))
        if (.not. allocated(p%sense)) p%sense = [(0_ip, i=1,p%m)]
        if (p%has_H) then
            p%Hc = reshape(transpose(p%H), [p%n*p%n])
        else
            allocate(p%Hc(1))
        end if
        if (p%factored) then
            p%Rf = cholesky_packed(p%H)
            p%Hc = p%Rf
        end if
        p%At = transpose(p%A)
        p%Ac = reshape(p%At, [max(1,size(p%A))], pad=[0.0_wp])
        if (p%has_f) then
            p%fc = p%f
        else
            allocate(p%fc(1))
        end if
        p%buc = p%bu
        p%blc = p%bl
        p%sensec = int(p%sense, c_int)
        if (allocated(p%break_points)) p%bpc = int(p%break_points, c_int)
        if (allocated(p%rho_l)) then
            p%rho_lc = p%rho_l; p%rho_uc = p%rho_u; p%w_lc = p%w_l; p%w_uc = p%w_u
        end if
    end subroutine finalize_problem

    subroutine c_problem(p, cqp)
        type(qp_problem), intent(inout), target :: p
        type(daqp_c_problem), intent(out) :: cqp
        cqp%n = p%n
        cqp%m = p%m
        cqp%ms = p%ms
        cqp%H = c_null_ptr
        if (p%has_H) cqp%H = c_loc(p%Hc)
        cqp%f = c_null_ptr
        if (p%has_f) cqp%f = c_loc(p%fc)
        cqp%A = c_loc(p%Ac)
        cqp%bupper = c_loc(p%buc)
        cqp%blower = c_loc(p%blc)
        cqp%sense = c_loc(p%sensec)
        cqp%nh = 1
        cqp%break_points = c_null_ptr
        if (allocated(p%bpc)) then
            cqp%break_points = c_loc(p%bpc)
            cqp%nh = size(p%bpc)
        end if
        cqp%problem_type = p%problem_type
        if (p%factored) cqp%problem_type = daqp_problem_factored
    end subroutine c_problem

!*****************************************************************************************
!>
!  Set up a problem in both solvers.

    subroutine setup_both(p, work, cqp, cwork, fflag, cflag)
        type(qp_problem), intent(inout), target :: p
        type(daqp_workspace), intent(inout) :: work
        type(daqp_c_problem), intent(inout), target :: cqp
        type(c_ptr), intent(out) :: cwork
        integer(ip), intent(out) :: fflag, cflag
        logical :: ok
        integer(c_int) :: ic
        ! (an unallocated actual argument is an absent optional one)
        work%settings = fsettings
        work%settings%eq_reduction = p%eq_reduction
        if (p%factored) then
            fflag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, f=p%f, sense=p%sense, &
                               init_mask=p%init_mask, At=p%At, break_points=p%break_points, &
                               problem_type=p%problem_type, Rf=p%Rf, &
                               primal_start=p%primal_start, dual_start=p%dual_start)
        else if (p%has_H) then
            fflag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, H=p%H, f=p%f, sense=p%sense, &
                               init_mask=p%init_mask, At=p%At, break_points=p%break_points, &
                               problem_type=p%problem_type, &
                               primal_start=p%primal_start, dual_start=p%dual_start)
        else
            fflag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, f=p%f, sense=p%sense, &
                               init_mask=p%init_mask, At=p%At, break_points=p%break_points, &
                               problem_type=p%problem_type, &
                               primal_start=p%primal_start, dual_start=p%dual_start)
        end if
        if (fflag > 0 .and. allocated(p%rho_l)) &
            ok = daqp_set_soft_weights(work, p%rho_l, p%rho_u, p%w_l, p%w_u)

        call c_problem(p, cqp)
        if (allocated(p%dual_start)) then
            call daqp_c_dual_init_active(cqp, p%dual_start)
        else if (allocated(p%primal_start)) then
            call daqp_c_primal_init_active(cqp, p%primal_start)
        end if
        csettings%eq_reduction = int(p%eq_reduction, c_int)
        cwork = cmp_ws_new(csettings)
        cflag = setup_daqp_main(cqp, cwork, c_null_ptr, int(p%init_mask, c_int))
        if (cflag > 0) then
            if (allocated(p%primal_start)) call daqp_c_set_primal_start(cwork, p%primal_start)
            if (allocated(p%rho_l)) ic = daqp_c_set_soft_weights(cwork, c_loc(p%rho_lc), &
                c_loc(p%rho_uc), c_loc(p%w_lc), c_loc(p%w_uc))
        end if
    end subroutine setup_both

!*****************************************************************************************
!>
!  Solve the problem (already set up) with both solvers, and compare.

    subroutine solve_and_compare(p, work, cwork, label)
        type(qp_problem), intent(inout) :: p
        type(daqp_workspace), intent(inout) :: work
        type(c_ptr), intent(in) :: cwork
        character(len=*), intent(in) :: label

        real(wp), allocatable :: xf(:), lamf(:)
        real(c_double), allocatable, target :: xc(:), lamc(:)
        integer(c_int), allocatable :: wsc(:), lowc(:)
        integer(ip), allocatable :: wsf(:)
        logical, allocatable :: lowf(:)
        type(daqp_result) :: res
        type(daqp_c_result) :: cres
        integer :: na_c, i
        real(wp) :: dx, dlam, df, tolx, kkt_f(4), kkt_c(4), eps, kgram
        logical :: same_ws, kkt_ok
        character(len=512) :: msg

        allocate(xf(p%n), lamf(p%m), xc(p%n), lamc(p%m))
        xc = 0.0_wp
        lamc = 0.0_wp
        call daqp_solve(work, xf, res, lamf)
        cres%x = c_loc(xc)
        cres%lam = c_loc(lamc)
        call daqp_c_solve(cres, cwork)

        st%n_problems = st%n_problems + 1
        if (cres%exitflag > 0) st%n_solved = st%n_solved + 1
        eps = epsilon(1.0_wp)

        ! working sets (as sets of (index, lower))
        na_c = cmp_ws_n_active(cwork)
        allocate(wsc(max(1,na_c)), lowc(max(1,na_c)))
        if (na_c > 0) call cmp_ws_working_set(cwork, wsc, lowc)
        wsf = work%WS(1:work%n_active)
        lowf = [(iand(work%sense(wsf(i)), daqp_lower) /= 0, i=1,work%n_active)]
        same_ws = na_c == work%n_active
        if (same_ws) then
            do i = 1, na_c
                if (.not. any(wsf == wsc(i)+1)) then
                    same_ws = .false.
                else if (any(wsf == wsc(i)+1 .and. (lowf .neqv. (lowc(i) /= 0)))) then
                    same_ws = .false.
                end if
            end do
        end if

        msg = ''
        if (res%exitflag == cres%exitflag) then
            st%n_flag_equal = st%n_flag_equal + 1
        else
            write(msg,'(A,": exit flag ",I0," (Fortran) vs ",I0," (C)")') label, res%exitflag, cres%exitflag
            call add_mismatch(msg)
        end if
        if (res%iter == cres%iter .and. res%nodes == cres%nodes) then
            st%n_iter_equal = st%n_iter_equal + 1
        else
            write(msg,'(A,": iterations (nodes) ",I0," (",I0,") (Fortran) vs ",I0," (",I0,") (C)")') &
                label, res%iter, res%nodes, cres%iter, cres%nodes
            call add_mismatch(msg)
        end if
        if (same_ws) then
            st%n_ws_equal = st%n_ws_equal + 1
        else
            write(msg,'(A,": working sets differ (",I0," vs ",I0," active)")') label, work%n_active, na_c
            call add_mismatch(msg)
        end if

        if (cres%exitflag > 0 .and. res%exitflag > 0) then
            ! the condition number of H, or of the working set's Gram matrix
            kgram = 1.0_wp
            if (work%n_active > 0) kgram = 1.0_wp/max(epsilon(1.0_wp), minval(work%D(1:work%n_active)))
            ! (soft constraints: the penalty s^2/(2 rho) conditions the problem as 1/rho)
            if (any(iand(p%sense, daqp_soft) /= 0)) kgram = max(kgram, 1.0_wp/fsettings%rho_soft)
            tolx = c_tol*real(p%n,wp)*eps*max(1.0_wp, p%kappa, kgram)
            dx = maxval(abs(xf-xc)) / (1.0_wp + maxval(abs(xc)))
            dlam = 0.0_wp
            if (p%m > 0) dlam = maxval(abs(lamf-lamc)) / (1.0_wp + maxval(abs(lamc)))
            df = abs(res%fval-cres%fval) / (1.0_wp + abs(cres%fval))
            st%max_dx = max(st%max_dx, dx)
            st%max_dlam = max(st%max_dlam, dlam)
            st%max_df = max(st%max_df, df)
            if (dx <= tolx) st%n_x_pass = st%n_x_pass + 1
            if (dlam <= tolx) st%n_lam_pass = st%n_lam_pass + 1
            if (df <= tolx) st%n_f_pass = st%n_f_pass + 1
            if (dx > tolx .or. dlam > tolx .or. df > tolx) then
                write(msg,'(A,": dx = ",ES9.2,", dlam = ",ES9.2,", df = ",ES9.2," (tol ",ES9.2,")")') &
                    label, dx, dlam, df, tolx
                call add_mismatch(msg)
            end if
            kkt_ok = .true.
            if (p%check_kkt) then
                call kkt_quad(p, xf, lamf, kkt_f)
                call kkt_quad(p, real(xc,wp), real(lamc,wp), kkt_c)
                do i = 1, 4
                    if (kkt_f(i) > c_tol*kkt_c(i) .and. kkt_f(i) > fsettings%primal_tol) kkt_ok = .false.
                end do
            end if
            if (kkt_ok) then
                st%n_kkt_pass = st%n_kkt_pass + 1
            else
                write(msg,'(A,": KKT (stat, prim, dual, comp) ",4ES9.2," (Fortran) vs ",4ES9.2," (C)")') &
                    label, kkt_f, kkt_c
                call add_mismatch(msg)
            end if
        else if (res%exitflag == cres%exitflag) then ! not solved by either: nothing to compare
            call pass_all()
        end if
    end subroutine solve_and_compare

!*****************************************************************************************
!>
!  KKT residuals in quadruple precision: stationarity, primal infeasibility,
!  wrong-signed multipliers, complementarity (infinity norms, relative to the
!  size of the data). Soft constraints are left out.

    subroutine kkt_quad(p, x, lam, r)
        type(qp_problem), intent(in) :: p
        real(wp), intent(in) :: x(:), lam(:)
        real(wp), intent(out) :: r(4)
        real(real128), allocatable :: g(:), ax(:), xq(:), lq(:)
        real(real128) :: big, scale
        integer :: i, j
        big = 1.0e29_real128
        xq = real(x, real128)
        lq = real(lam, real128)
        allocate(g(p%n), ax(p%m))
        g = 0
        if (p%has_f) g = real(p%f, real128)
        if (p%has_H) then
            do j = 1, p%n
                do i = 1, p%n
                    g(i) = g(i) + real(p%H(i,j),real128)*xq(j)
                end do
            end do
        end if
        do i = 1, p%ms
            g(i) = g(i) + lq(i)
            ax(i) = xq(i)
        end do
        do i = p%ms+1, p%m
            ax(i) = 0
            do j = 1, p%n
                ax(i) = ax(i) + real(p%A(i-p%ms,j),real128)*xq(j)
                g(j) = g(j) + real(p%A(i-p%ms,j),real128)*lq(i)
            end do
        end do
        scale = 1.0_real128 + maxval(abs(xq)) + maxval(abs(g))
        if (p%m > 0) scale = scale + maxval(abs(lq))
        r(1) = real(maxval(abs(g))/scale, wp)
        r(2:4) = 0.0_wp
        do i = 1, p%m
            if (iand(p%sense(i), daqp_soft) /= 0) cycle
            if (p%bu(i) < big) r(2) = max(r(2), real(ax(i)-p%bu(i), wp))
            if (p%bl(i) > -big) r(2) = max(r(2), real(p%bl(i)-ax(i), wp))
            if (iand(p%sense(i), daqp_immutable) /= 0) cycle
            if (lq(i) > 0) then
                if (p%bu(i) >= big) then
                    r(3) = max(r(3), real(lq(i),wp))
                else
                    r(4) = max(r(4), real(lq(i)*abs(p%bu(i)-ax(i))/scale, wp))
                end if
            else if (lq(i) < 0) then
                if (p%bl(i) <= -big) then
                    r(3) = max(r(3), real(-lq(i),wp))
                else
                    r(4) = max(r(4), real(-lq(i)*abs(ax(i)-p%bl(i))/scale, wp))
                end if
            end if
        end do
    end subroutine kkt_quad

!*****************************************************************************************
!>
!  Set up, solve, and compare one problem (`nsolve` solves in a row).

    subroutine run_problem(p, label, nsolve)
        type(qp_problem), intent(inout), target :: p
        character(len=*), intent(in) :: label
        integer, intent(in), optional :: nsolve
        type(daqp_workspace) :: work
        type(daqp_c_problem), target :: cqp
        type(c_ptr) :: cwork
        integer(ip) :: fflag, cflag
        integer :: k, ns
        character(len=256) :: msg
        ns = 1
        if (present(nsolve)) ns = nsolve
        call finalize_problem(p)
        call setup_both(p, work, cqp, cwork, fflag, cflag)
        if (fflag < 0 .or. cflag < 0) then
            st%n_problems = st%n_problems + 1
            if (fflag == cflag) then
                st%n_flag_equal = st%n_flag_equal + 1
                st%n_iter_equal = st%n_iter_equal + 1
                st%n_ws_equal = st%n_ws_equal + 1
                call pass_all()
            else
                write(msg,'(A,": setup flag ",I0," (Fortran) vs ",I0," (C)")') label, fflag, cflag
                call add_mismatch(msg)
            end if
        else
            do k = 1, ns
                if (k == 1) then
                    call solve_and_compare(p, work, cwork, label)
                else
                    write(msg,'(A,", solve ",I0)') label, k
                    call solve_and_compare(p, work, cwork, trim(msg))
                end if
            end do
        end if
        call cmp_ws_free(cwork)
        call daqp_destroy(work)
    end subroutine run_problem

    subroutine random_problem(p, n, m, ms, nact, kappa, name)
        type(qp_problem), intent(out) :: p
        integer, intent(in) :: n, m, ms, nact
        real(wp), intent(in) :: kappa
        character(len=*), intent(in) :: name
        real(wp), allocatable :: xref(:)
        call generate_qp(n, m, ms, nact, kappa, xref, p%H, p%f, p%A, p%bu, p%bl)
        p%n = n; p%m = m; p%ms = ms; p%kappa = kappa; p%name = name
    end subroutine random_problem

!*****************************************************************************************
!>
!  Set 1: upstream's test problems.

    subroutine set_upstream()
        type(qp_problem) :: p
        real(wp), allocatable :: xref(:)
        integer :: k, nqp
        character(len=64) :: label
        call new_set('1. upstream tests')
        call rng_seed(1234)
        nqp = 100
        if (quick) nqp = 10
        do k = 1, nqp
            write(label,'("QP ",I0)') k
            call random_problem(p, 100, 500, 50, 80, 1.0e2_wp, trim(label))
            call run_problem(p, trim(label))
        end do
        do k = 1, nqp/2
            write(label,'("LP ",I0)') k
            p = qp_problem()
            p%n = 100; p%m = 500; p%ms = 50; p%has_H = .false.; p%kappa = 1.0e4_wp
            call generate_lp(p%n, p%m, p%ms, xref, p%f, p%A, p%bu, p%bl)
            call run_problem(p, trim(label))
        end do
        do k = 1, nqp/5 ! (n=20, m=100, ms=0, 16 active), as upstream's Julia tests
            write(label,'("QP small ",I0)') k
            call random_problem(p, 20, 100, 0, 16, 1.0e2_wp, trim(label))
            call run_problem(p, trim(label))
        end do
        do k = 1, nqp/5 ! one-shot (quadprog) semantics: unconstrained check and elimination
            write(label,'("QP quadprog ",I0)') k
            call random_problem(p, 100, 500, 50, 80, 1.0e2_wp, trim(label))
            p%init_mask = daqp_update_unconstrained + daqp_update_eliminate
            call run_problem(p, trim(label))
        end do
    end subroutine set_upstream

!*****************************************************************************************
!>
!  Set 2: random QPs of various sizes and conditioning.

    subroutine set_random()
        type(qp_problem) :: p
        integer :: in, im, ik, ia, n, m, ms, nact
        integer, parameter :: ns(8) = [2, 5, 10, 20, 50, 100, 200, 500]
        real(wp), parameter :: kappas(4) = [1.0_wp, 1.0e2_wp, 1.0e5_wp, 1.0e10_wp]
        character(len=64) :: label
        call new_set('2. random QPs')
        call rng_seed(2026)
        do in = 1, size(ns)
            n = ns(in)
            if (quick .and. n > 50) exit
            do im = 0, 3
                m = (im*n)  ! m = n/2, n, 2n, 3n
                if (im == 0) m = n/2
                do ik = 1, size(kappas)
                    do ia = 1, 2
                        ms = min(m, n)/2
                        nact = min(m, n)
                        if (ia == 1) nact = nact/2
                        write(label,'("n=",I0,", m=",I0,", ms=",I0,", nact=",I0,", kappa=",ES7.0)') &
                            n, m, ms, nact, kappas(ik)
                        call random_problem(p, n, m, ms, nact, kappas(ik), trim(label))
                        call run_problem(p, trim(label))
                    end do
                end do
            end do
        end do
    end subroutine set_random

!*****************************************************************************************
!>
!  Set 3: degenerate QPs: duplicated rows, dependent active rows, more active
!  rows than variables, weakly active constraints.

    subroutine set_degenerate()
        type(qp_problem) :: p
        integer :: k, n, m, ms, nact, nd, i, j, nrep
        real(wp), allocatable :: xs(:), a(:), lam(:)
        integer, allocatable :: dup(:)
        character(len=64) :: label
        call new_set('3. degenerate QPs')
        call rng_seed(77)
        nrep = 20
        if (quick) nrep = 4
        ! duplicated rows (active and inactive)
        do k = 1, nrep
            n = 20; m = 60; ms = 10; nact = 15
            call random_problem(p, n, m, ms, nact, 1.0e3_wp, '')
            nd = 10
            dup = [(rand_int(1, m-ms), i=1,nd)]
            p%A = reshape([transpose(p%A), transpose(p%A(dup,:))], [n, m-ms+nd])
            p%A = transpose(p%A)
            p%bu = [p%bu, p%bu(ms+dup)]
            p%bl = [p%bl, p%bl(ms+dup)]
            p%m = m + nd
            write(label,'("duplicated rows ",I0)') k
            call run_problem(p, trim(label))
        end do
        ! a vertex with more active rows than variables, some weakly active,
        ! some dependent (sums of others): H = I, x* given
        do k = 1, nrep
            n = 10 + mod(k,3)*5
            m = 3*n
            p = qp_problem()
            p%n = n; p%m = m; p%ms = 0; p%kappa = 1.0e2_wp
            allocate(p%A(m,n), p%bu(m), p%bl(m), lam(m))
            xs = [(randn(), i=1,n)]
            do i = 1, m
                if (i > 2*n .and. mod(i,2) == 0) then ! dependent: sum of two earlier rows
                    p%A(i,:) = p%A(i-2*n,:) + p%A(i-2*n+1,:)
                else
                    p%A(i,:) = [(randn(), j=1,n)]
                end if
            end do
            lam = 0.0_wp
            do i = 1, n/2
                lam(i) = urand() ! strictly active
            end do
            p%H = reshape([((merge(1.0_wp, 0.0_wp, i == j), i=1,n), j=1,n)], [n,n])
            p%f = -xs - matmul(transpose(p%A), lam)
            a = matmul(p%A, xs)
            do i = 1, m
                if (i <= 3*n/2 .or. i > 2*n) then ! tight at x* (weakly active if lam = 0)
                    p%bu(i) = a(i)
                else
                    p%bu(i) = a(i) + 0.1_wp + urand()
                end if
                p%bl(i) = -daqp_inf
            end do
            deallocate(lam)
            write(label,'("over-determined vertex ",I0)') k
            call run_problem(p, trim(label))
        end do
    end subroutine set_degenerate

!*****************************************************************************************
!>
!  Set "special": infeasible, semidefinite, LP, soft constraints, diagonal H.

    subroutine set_special()
        type(qp_problem) :: p
        integer :: k, n, m, ms, i, nrep
        real(wp), allocatable :: xref(:)
        character(len=64) :: label
        call new_set('3b. special cases')
        call rng_seed(99)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            ! infeasible: two contradicting rows
            n = 10; m = 30; ms = 5
            call random_problem(p, n, m, ms, 5, 1.0e2_wp, '')
            p%A(2,:) = -p%A(1,:)
            p%bl(ms+2) = -p%bl(ms+1) + 1.0_wp
            p%bu(ms+2) = daqp_inf
            p%bu(ms+1) = daqp_inf
            write(label,'("infeasible ",I0)') k
            call run_problem(p, trim(label))
            ! semidefinite H: zero rows (semi-proximal)
            call random_problem(p, n, m, ms, 5, 1.0e2_wp, '')
            p%H(:,n) = 0.0_wp; p%H(n,:) = 0.0_wp
            p%H(:,n-1) = 0.0_wp; p%H(n-1,:) = 0.0_wp
            p%kappa = 1.0e6_wp
            write(label,'("semidefinite (zero rows) ",I0)') k
            call run_problem(p, trim(label))
            ! semidefinite H: rank deficient, dense (full proximal shift)
            call random_problem(p, n, m, ms, 5, 1.0e2_wp, '')
            p%H = matmul(p%H(:,1:3), transpose(p%H(:,1:3)))
            p%kappa = 1.0e6_wp
            write(label,'("semidefinite (dense) ",I0)') k
            call run_problem(p, trim(label))
            ! diagonal H
            call random_problem(p, n, m, ms, 5, 1.0e2_wp, '')
            p%H = 0.0_wp
            do i = 1, n
                p%H(i,i) = 1.0_wp + 10.0_wp*urand()
            end do
            write(label,'("diagonal H ",I0)') k
            call run_problem(p, trim(label))
            ! soft constraints (on infeasible rows)
            call random_problem(p, n, m, ms, 5, 1.0e2_wp, '')
            p%sense = [(0_ip, i=1,m)]
            p%sense(ms+1:ms+4) = daqp_soft
            call conflicting_soft_rows(p)
            p%check_kkt = .false.
            write(label,'("soft ",I0)') k
            call run_problem(p, trim(label))
            ! small LP
            p = qp_problem()
            p%n = 10; p%m = 40; p%ms = 10; p%has_H = .false.; p%kappa = 1.0e4_wp
            call generate_lp(p%n, p%m, p%ms, xref, p%f, p%A, p%bu, p%bl)
            write(label,'("LP ",I0)') k
            call run_problem(p, trim(label))
        end do
    end subroutine set_special

!*****************************************************************************************
!>
!  Set 4: warm-start sequences (MPC-like): the same QP with a changing linear
!  term and bounds, each solve hot started from the previous working set.

    subroutine set_warm()
        type(qp_problem), target :: p
        type(daqp_workspace) :: work
        type(daqp_c_problem), target :: cqp
        type(c_ptr) :: cwork
        integer(ip) :: fflag, cflag
        integer(c_int) :: mask
        integer :: iseq, k, nseq, nstep, i
        real(wp), allocatable :: f0(:), df(:), db(:)
        character(len=64) :: label
        call new_set('4. warm-start sequences')
        call rng_seed(4)
        nseq = 10; nstep = 20
        if (quick) nseq = 3
        do iseq = 1, nseq
            call random_problem(p, 30, 90, 15, 20, 1.0e2_wp, '')
            call finalize_problem(p)
            f0 = p%f
            df = [(randn(), i=1,p%n)]
            ! the bounds move with a translation of the feasible set (which stays feasible)
            db = [(0.1_wp*randn(), i=1,p%n)]
            db = [db(1:p%ms), matmul(p%A, db)]
            call setup_both(p, work, cqp, cwork, fflag, cflag)
            do k = 0, nstep
                if (k > 0) then
                    p%f = f0 + 0.05_wp*real(k,wp)*df
                    p%bu = p%bu + 0.05_wp*db
                    p%bl = p%bl + 0.05_wp*db
                    p%fc = p%f; p%buc = p%bu; p%blc = p%bl
                    work%qp%f = p%f; work%qp%bupper = p%bu; work%qp%blower = p%bl
                    mask = daqp_update_v + daqp_update_d
                    fflag = daqp_update_ldp(work, int(mask,ip))
                    cflag = daqp_c_update_ldp(mask, cwork, cqp)
                end if
                write(label,'("sequence ",I0,", step ",I0)') iseq, k
                call solve_and_compare(p, work, cwork, trim(label))
            end do
            call cmp_ws_free(cwork)
            call daqp_destroy(work)
        end do
    end subroutine set_warm

!*****************************************************************************************
!>
!  Set 5: elimination of equality constraints (automatic, as `daqp_quadprog`,
!  and forced), with a dense, diagonal, or singular Hessian, an LP, equalities
!  that determine x, and an update of the bounds of a reduced problem.

    subroutine set_eq()
        type(qp_problem), target :: p
        type(daqp_workspace) :: work
        type(daqp_c_problem), target :: cqp
        type(c_ptr) :: cwork
        integer(ip) :: fflag, cflag
        integer(c_int) :: mask
        real(wp), allocatable :: xref(:), shift(:)
        integer :: k, i, nrep, neq, ivar, n, m, ms
        character(len=96) :: label
        call new_set('5. equality elimination')
        call rng_seed(55)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            do ivar = 1, 6
                p = qp_problem()
                n = 30; m = 90; ms = 10
                if (ivar /= 5 .and. ivar /= 6) then
                    call generate_qp(n, m, ms, 25, 1.0e2_wp, xref, p%H, p%f, p%A, p%bu, p%bl)
                    p%n = n; p%m = m; p%ms = ms; p%kappa = 1.0e2_wp
                end if
                p%init_mask = daqp_update_unconstrained + daqp_update_eliminate
                p%eq_reduction = daqp_eq_reduction_on
                select case (ivar)
                case (1) ! dense H, automatic (quadprog)
                    neq = make_equalities(ms, 15, xref, p%A, p%bu, p%bl)
                    p%eq_reduction = daqp_eq_reduction_auto
                    write(label,'("dense H, ",I0," equalities (auto) ",I0)') neq, k
                case (2) ! dense H, forced, split interface (with equality flags)
                    neq = make_equalities(ms, 8, xref, p%A, p%bu, p%bl)
                    p%init_mask = 0
                    p%sense = [(0_ip, i=1,m)]
                    do i = ms+1, m
                        if (p%bu(i) == p%bl(i)) p%sense(i) = daqp_equality
                    end do
                    write(label,'("dense H, ",I0," equalities (on) ",I0)') neq, k
                case (3) ! diagonal H (metric)
                    p%H = 0.0_wp
                    do i = 1, n
                        p%H(i,i) = 1.0_wp + 9.0_wp*urand()
                    end do
                    call generate_qp_bounds_at(p, 15)
                    write(label,'("diagonal H (on) ",I0)') k
                case (4) ! singular H (the reduced Hessian is singular: PATH_QP)
                    p%H(:,n) = 0.0_wp; p%H(n,:) = 0.0_wp
                    p%H(:,n-1) = 0.0_wp; p%H(n-1,:) = 0.0_wp
                    p%kappa = 1.0e6_wp
                    neq = make_equalities(ms, 15, xref, p%A, p%bu, p%bl)
                    write(label,'("singular H (on) ",I0)') k
                case (5) ! LP (PATH_LP)
                    p%n = 20; p%m = 60; p%ms = 5; p%has_H = .false.; p%kappa = 1.0e4_wp
                    call generate_lp(p%n, p%m, p%ms, xref, p%f, p%A, p%bu, p%bl)
                    neq = make_equalities(p%ms, 8, xref, p%A, p%bu, p%bl)
                    write(label,'("LP (on) ",I0)') k
                case (6) ! equalities that determine x (nz = 0)
                    n = 12; m = 30; ms = 0
                    call generate_qp(n, m, ms, n, 1.0e2_wp, xref, p%H, p%f, p%A, p%bu, p%bl)
                    p%n = n; p%m = m; p%ms = ms; p%kappa = 1.0e2_wp
                    neq = make_equalities(ms, n, xref, p%A, p%bu, p%bl)
                    write(label,'("determined x, ",I0," equalities (on) ",I0)') neq, k
                end select
                call run_problem(p, trim(label), nsolve=2)
            end do
            ! an update of the bounds of a reduced problem (warm right-hand side)
            n = 30; m = 90; ms = 10
            p = qp_problem()
            call generate_qp(n, m, ms, 25, 1.0e2_wp, xref, p%H, p%f, p%A, p%bu, p%bl)
            p%n = n; p%m = m; p%ms = ms; p%kappa = 1.0e2_wp
            neq = make_equalities(ms, 15, xref, p%A, p%bu, p%bl)
            p%eq_reduction = daqp_eq_reduction_on
            call finalize_problem(p)
            call setup_both(p, work, cqp, cwork, fflag, cflag)
            write(label,'("update of a reduced problem ",I0)') k
            call solve_and_compare(p, work, cwork, trim(label))
            shift = [(0.1_wp*randn(), i=1,n)]
            shift = [shift(1:ms), matmul(p%A, shift)]
            p%bu = p%bu + shift
            p%bl = p%bl + shift
            p%buc = p%bu; p%blc = p%bl
            work%qp%bupper = p%bu; work%qp%blower = p%bl
            mask = daqp_update_d
            fflag = daqp_update_ldp(work, int(mask,ip))
            cflag = daqp_c_update_ldp(mask, cwork, cqp)
            write(label,'("update of a reduced problem ",I0,", after the update")') k
            call solve_and_compare(p, work, cwork, trim(label))
            call cmp_ws_free(cwork)
            call daqp_destroy(work)
        end do
    end subroutine set_eq

    subroutine conflicting_soft_rows(p)
        !! rows ms+1 and ms+2 (soft): the same row, with a'x >= c+1 and a'x <= c
        type(qp_problem), intent(inout) :: p
        real(wp) :: c
        p%A(2,:) = p%A(1,:)
        c = 0.5_wp*(p%bu(p%ms+1) + p%bl(p%ms+1))
        p%bl(p%ms+1) = c + 1.0_wp
        p%bu(p%ms+1) = daqp_inf
        p%bl(p%ms+2) = -daqp_inf
        p%bu(p%ms+2) = c
    end subroutine conflicting_soft_rows

    subroutine generate_qp_bounds_at(p, neq)
        !! make neq general constraints equalities at the solution of the problem (by the port)
        type(qp_problem), intent(inout) :: p
        integer, intent(in) :: neq
        type(daqp_result) :: r
        real(wp), allocatable :: x(:)
        integer :: k
        allocate(x(p%n))
        call daqp_quadprog(p%n, p%m, p%ms, p%bu, p%bl, x, r, H=p%H, f=p%f, A=p%A, settings=fsettings)
        do k = 1, neq
            p%bu(p%ms+k) = dot_product(p%A(k,:), x)
            p%bl(p%ms+k) = p%bu(p%ms+k)
        end do
    end subroutine generate_qp_bounds_at

!*****************************************************************************************
!>
!  Set 6: branch and bound (binary constraints), as upstream's tests.

    subroutine set_bnb()
        type(qp_problem) :: p
        integer :: k, nrep, i
        character(len=64) :: label
        call new_set('6. branch and bound')
        call rng_seed(66)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            p = qp_problem()
            p%n = 20; p%m = 60; p%ms = 20; p%kappa = 1.0e4_wp
            call generate_miqp(p%n, p%m, p%ms, 10, p%H, p%f, p%A, p%bu, p%bl, p%sense)
            p%check_kkt = .false.
            write(label,'("MIQP ",I0)') k
            call run_problem(p, trim(label), nsolve=2) ! (the second solve is warm started at the root)
            ! with an equality constraint, and a diagonal Hessian
            p = qp_problem()
            p%n = 20; p%m = 60; p%ms = 20; p%kappa = 1.0e4_wp
            call generate_miqp(p%n, p%m, p%ms, 10, p%H, p%f, p%A, p%bu, p%bl, p%sense)
            p%sense(p%m) = daqp_equality
            p%bu(p%m) = 0.0_wp; p%bl(p%m) = 0.0_wp
            if (mod(k,2) == 0) then
                do i = 1, p%n
                    p%H(:,i) = 0.0_wp
                    p%H(i,i) = 1.0_wp + 5.0_wp*urand()
                end do
            end if
            p%check_kkt = .false.
            write(label,'("MIQP with an equality ",I0)') k
            call run_problem(p, trim(label))
            ! with an incumbent (primal start)
            p = qp_problem()
            p%n = 15; p%m = 40; p%ms = 15; p%kappa = 1.0e4_wp
            call generate_miqp(p%n, p%m, p%ms, 8, p%H, p%f, p%A, p%bu, p%bl, p%sense)
            p%primal_start = [(0.0_wp, i=1,p%n)] ! the origin is feasible
            p%check_kkt = .false.
            write(label,'("MIQP with an incumbent ",I0)') k
            call run_problem(p, trim(label))
        end do
        ! upstream's small example
        p = qp_problem()
        p%n = 3; p%m = 5; p%ms = 3
        p%H = reshape([1.0_wp,0.5_wp,0.0_wp, 0.5_wp,1.0_wp,0.5_wp, 0.0_wp,0.5_wp,1.0_wp], [3,3])
        p%f = [1.0_wp, 0.0_wp, 0.0_wp]
        p%A = reshape([1.0_wp,1.0_wp, 2.0_wp,1.0_wp, 3.0_wp,0.0_wp], [2,3])
        p%bu = [1.0_wp, 1.0_wp, 1.0_wp, 1.0e30_wp, 1.0e30_wp]
        p%bl = [0.0_wp, 0.0_wp, 0.0_wp, 4.0_wp, 1.0_wp]
        p%sense = [daqp_binary, daqp_binary, daqp_binary, 0_ip, 0_ip]
        p%check_kkt = .false.
        call run_problem(p, 'small MIQP')
    end subroutine set_bnb

!*****************************************************************************************
!>
!  Set 7: hierarchical QPs, as upstream's tests (also repeated solves after an
!  update).

    subroutine set_hierarchical()
        type(qp_problem), target :: p
        type(daqp_workspace) :: work
        type(daqp_c_problem), target :: cqp
        type(c_ptr) :: cwork
        integer(ip) :: fflag, cflag
        real(wp), allocatable :: L(:,:), x0(:)
        integer :: k, nrep, i, j, n, mh, ml, nl, m
        character(len=64) :: label
        call new_set('7. hierarchical QPs')
        call rng_seed(77)
        nrep = 10
        if (quick) nrep = 3
        ! upstream's examples (no H: minimal norm)
        p = qp_problem()
        p%n = 3; p%m = 6; p%ms = 3; p%has_H = .false.; p%has_f = .false.
        p%A = reshape([1.0_wp,1.0_wp,3.0_wp, 1.0_wp,-1.0_wp,1.0_wp, 1.0_wp,0.0_wp,-1.0_wp], [3,3])
        p%bu = [1.0_wp,1.0_wp,1.0_wp,1.0_wp,0.5_wp,20.0_wp]
        p%bl = [-1.0_wp,-1.0_wp,-1.0_wp,-1.0e30_wp,0.5_wp,10.0_wp]
        p%break_points = [3_ip,4_ip,5_ip,6_ip]
        p%check_kkt = .false.
        call run_problem(p, 'upstream example 1')
        p = qp_problem()
        p%n = 2; p%m = 3; p%ms = 0; p%has_H = .false.; p%has_f = .false.
        p%A = reshape([1.0_wp,1.0_wp,0.0_wp, 0.0_wp,0.0_wp,1.0_wp], [3,2])
        p%bu = [4.0_wp,8.0_wp,1.0_wp]
        p%bl = [4.0_wp,8.0_wp,1.0_wp]
        p%break_points = [0_ip,2_ip,3_ip]
        p%check_kkt = .false.
        call run_problem(p, 'upstream example 2 (degenerate)')
        ! random hierarchies with conflicting soft levels
        do k = 1, nrep
            n = 10; mh = 3; ml = 3; nl = 2
            m = mh + nl*ml
            p = qp_problem()
            p%n = n; p%m = m; p%ms = 0; p%kappa = 1.0e2_wp
            allocate(L(n,n), p%A(m,n))
            do j = 1, n
                do i = 1, n
                    L(i,j) = randn()
                end do
            end do
            p%H = matmul(transpose(L), L)/real(n,wp)
            do i = 1, n
                p%H(i,i) = p%H(i,i) + 1.0_wp
            end do
            p%f = [(randn(), i=1,n)]
            do j = 1, n
                do i = 1, m
                    p%A(i,j) = randn()
                end do
            end do
            x0 = [(randn(), i=1,n)]
            p%bu = matmul(p%A, x0) + [(0.5_wp*urand(), i=1,m)]
            p%bl = matmul(p%A, x0) - [(0.5_wp*urand(), i=1,m)]
            x0 = [(randn(), i=1,n)]
            p%bu(mh+1:) = matmul(p%A(mh+1:,:), x0) + 0.05_wp ! conflicting soft levels
            p%bl(mh+1:) = p%bu(mh+1:) - 0.1_wp
            p%break_points = [(int(mh + i*ml, ip), i=0,nl)]
            p%check_kkt = .false.
            deallocate(L)
            write(label,'("random hierarchy ",I0)') k
            call finalize_problem(p)
            call setup_both(p, work, cqp, cwork, fflag, cflag)
            call solve_and_compare(p, work, cwork, trim(label))
            ! a second solve without an update is unsupported, then after an update
            write(label,'("random hierarchy ",I0,", resolve")') k
            call solve_and_compare(p, work, cwork, trim(label))
            fflag = daqp_update_ldp(work, 0_ip)
            cflag = daqp_c_update_ldp(0_c_int, cwork, cqp)
            write(label,'("random hierarchy ",I0,", after an update")') k
            call solve_and_compare(p, work, cwork, trim(label))
            call cmp_ws_free(cwork)
            call daqp_destroy(work)
        end do
    end subroutine set_hierarchical

!*****************************************************************************************
!>
!  Set 8: affine variational inequalities, as upstream's tests.

    subroutine set_avi()
        type(qp_problem) :: p
        real(wp), allocatable :: xref(:)
        integer :: k, nrep, i
        character(len=64) :: label
        call new_set('8. AVIs')
        call rng_seed(88)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            p = qp_problem()
            p%n = 30; p%m = 100; p%ms = 0; p%kappa = 1.0e4_wp
            call generate_avi(p%n, p%m, xref, p%H, p%f, p%A, p%bu)
            p%bl = [(-1.0e30_wp, i=1,p%m)]
            p%problem_type = daqp_problem_avi
            p%init_mask = daqp_update_unconstrained + daqp_update_eliminate
            p%check_kkt = .false.
            write(label,'("AVI ",I0)') k
            call run_problem(p, trim(label), nsolve=2)
        end do
    end subroutine set_avi

!*****************************************************************************************
!>
!  Set 9: options: a prefactorized Hessian, primal and dual starts,
!  individual soft weights.

    subroutine set_options()
        type(qp_problem) :: p
        real(wp), allocatable :: xref(:), lam(:)
        type(daqp_result) :: r
        integer :: k, nrep, i
        character(len=64) :: label
        call new_set('9. options')
        call rng_seed(99)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            ! prefactorized H
            call random_problem(p, 30, 90, 15, 20, 1.0e2_wp, '')
            p%factored = .true.
            write(label,'("prefactorized H ",I0)') k
            call run_problem(p, trim(label))
            ! primal start
            p = qp_problem()
            call generate_qp(30, 90, 15, 20, 1.0e2_wp, xref, p%H, p%f, p%A, p%bu, p%bl)
            p%n = 30; p%m = 90; p%ms = 15; p%kappa = 1.0e2_wp
            p%primal_start = xref
            write(label,'("primal start ",I0)') k
            call run_problem(p, trim(label))
            ! dual start (the multipliers of the solution)
            deallocate(p%primal_start)
            allocate(lam(p%m))
            call daqp_quadprog(p%n, p%m, p%ms, p%bu, p%bl, xref, r, lam, H=p%H, f=p%f, A=p%A, &
                               settings=fsettings)
            p%dual_start = lam
            deallocate(lam)
            write(label,'("dual start ",I0)') k
            call run_problem(p, trim(label))
            ! individual soft weights (some soft rows violated)
            p = qp_problem()
            call generate_qp(20, 60, 10, 10, 1.0e2_wp, xref, p%H, p%f, p%A, p%bu, p%bl)
            p%n = 20; p%m = 60; p%ms = 10; p%kappa = 1.0e2_wp
            p%sense = [(0_ip, i=1,p%m)]
            p%sense(p%ms+1:p%ms+6) = daqp_soft
            call conflicting_soft_rows(p)
            p%rho_l = [(merge(0.1_wp+urand(), 0.0_wp, mod(i,2) == 0), i=1,p%m)]
            p%rho_u = [(merge(0.1_wp+urand(), 0.0_wp, mod(i,3) == 0), i=1,p%m)]
            p%w_l = [(merge(urand(), 0.0_wp, mod(i,2) == 1), i=1,p%m)]
            p%w_u = [(merge(urand(), 0.0_wp, mod(i,4) == 0), i=1,p%m)]
            p%check_kkt = .false.
            write(label,'("soft weights ",I0)') k
            call run_problem(p, trim(label))
        end do
    end subroutine set_options

!*****************************************************************************************
!>
!  Set 10: redundant constraints of polyhedra (minrep).

    subroutine set_minrep()
        real(wp), allocatable :: A(:,:), b(:), x0(:)
        real(c_double), allocatable :: Ac(:), bc(:)
        integer(ip), allocatable :: rf(:)
        integer(c_int), allocatable :: rc(:)
        integer :: k, nrep, i, j, n, m, ms
        character(len=128) :: msg
        call new_set('10. minrep')
        call rng_seed(10)
        nrep = 10
        if (quick) nrep = 3
        do k = 1, nrep
            n = 5; m = 40; ms = 2
            allocate(A(m-ms,n), b(m), rf(m), rc(m))
            do j = 1, n
                do i = 1, m-ms
                    A(i,j) = randn()
                end do
            end do
            x0 = [(0.1_wp*randn(), i=1,n)]
            b = [x0(1:ms) + 1.0_wp, matmul(A, x0) + [(1.0_wp + urand(), i=1,m-ms)]]
            Ac = reshape(transpose(A), [size(A)])
            bc = b
            call daqp_minrep(A, b, int(ms,ip), rf, fsettings)
            call daqp_c_minrep(rc, Ac, bc, int(n,c_int), int(m,c_int), int(ms,c_int))
            st%n_problems = st%n_problems + 1
            st%n_solved = st%n_solved + 1
            if (all(rf == rc)) then
                st%n_flag_equal = st%n_flag_equal + 1
                st%n_iter_equal = st%n_iter_equal + 1
                st%n_ws_equal = st%n_ws_equal + 1
                call pass_all()
            else
                write(msg,'("polyhedron ",I0,": ",I0," different classifications")') k, count(rf /= rc)
                call add_mismatch(msg)
            end if
            deallocate(A, b, rf, rc)
        end do
    end subroutine set_minrep

!*****************************************************************************************
!>
!  Wall clock time in seconds.

    real(wp) function wtime()
        integer(int64) :: count, rate
        call system_clock(count, rate)
        wtime = real(count,wp)/real(rate,wp)
    end function wtime

    real(wp) function median(t)
        real(wp), intent(in) :: t(:)
        real(wp), allocatable :: s(:)
        integer :: i, j
        real(wp) :: tmp
        s = t
        do i = 2, size(s)
            tmp = s(i)
            j = i - 1
            do while (j >= 1)
                if (s(j) <= tmp) exit
                s(j+1) = s(j)
                j = j - 1
            end do
            s(j+1) = tmp
        end do
        median = s((size(s)+1)/2)
    end function median

!*****************************************************************************************
!>
!  Speed: setup, cold solve, and warm (hot) solve, per QP, as the median of
!  batches of repetitions.

    subroutine speed()
        integer, parameter :: ns(7) = [5, 10, 20, 50, 100, 200, 500]
        integer, parameter :: nbatch = 7
        type(qp_problem), target :: p
        integer :: in, im, n, m, ms, nact, nrep, ib, r
        real(wp) :: tf(3), tc(3), t0, tb(nbatch), tb2(nbatch), target_time
        real(wp), allocatable :: rf(:,:)
        character(len=256) :: line
        integer :: nrow

        target_time = 0.1_wp
        if (quick) target_time = 0.01_wp
        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '## Speed'
        write(output_unit,'(A)') ''
        write(output_unit,'(A)') 'Random QPs (cond(H) = 100, `ms = n/2` simple bounds, `n/2` active'// &
            ' constraints). Time per QP in microseconds: setup (factor H, form the LDP), cold solve'// &
            ' (from the empty working set), and warm solve (`update` of `f` and solve, from the'// &
            ' previous working set). Each is the median of 7 batches of repetitions'// &
            ' (about 0.1 s per batch set). Both workspaces are freed after each setup (a repeated'// &
            ' setup of the same size on one Fortran object reuses its arrays, and is faster).'
        write(output_unit,'(A)') ''
        write(output_unit,'(A)') '| n | m | setup F | setup C | ratio | cold solve F | cold solve C | ratio |'// &
            ' warm solve F | warm solve C | ratio |'
        write(output_unit,'(A)') '|---|---|---|---|---|---|---|---|---|---|---|'
        allocate(rf(3, 2*size(ns)))
        nrow = 0
        call rng_seed(5)
        do in = 1, size(ns)
            n = ns(in)
            if (quick .and. n > 50) exit
            do im = 1, 2
                m = merge(n, 3*n, im == 1)
                ms = n/2
                nact = n/2
                call random_problem(p, n, m, ms, nact, 1.0e2_wp, '')
                call finalize_problem(p)
                ! repetitions per batch: about target_time/nbatch per batch
                t0 = wtime()
                call time_fortran(p, 1, tf)
                nrep = max(1, int(target_time/nbatch/max(1.0e-7_wp, wtime()-t0)))
                do ib = 1, nbatch
                    call time_fortran(p, nrep, tf)
                    tb(ib) = tf(1); tb2(ib) = tf(2)
                end do
                tf(1) = median(tb); tf(2) = median(tb2)
                do ib = 1, nbatch
                    call time_fortran_warm(p, nrep, tb(ib))
                end do
                tf(3) = median(tb)
                do ib = 1, nbatch
                    call time_c(p, nrep, tc)
                    tb(ib) = tc(1); tb2(ib) = tc(2)
                end do
                tc(1) = median(tb); tc(2) = median(tb2)
                do ib = 1, nbatch
                    call time_c_warm(p, nrep, tb(ib))
                end do
                tc(3) = median(tb)
                nrow = nrow + 1
                rf(:,nrow) = tf/tc
                write(line,'("| ",I0," | ",I0,3(" | ",F10.2," | ",F10.2," | ",F5.2),"|")') n, m, &
                    (1.0e6_wp*tf(r), 1.0e6_wp*tc(r), tf(r)/tc(r), r=1,3)
                write(output_unit,'(A)') trim(line)
            end do
        end do
        write(output_unit,'(A)') ''
        write(output_unit,'(A,3F6.2)') 'Median ratio port/C (setup, cold solve, warm solve): ', &
            (median(rf(r,1:nrow)), r=1,3)
    end subroutine speed

    subroutine time_fortran(p, nrep, t)
        !! t(1): setup, t(2): cold solve (per QP)
        type(qp_problem), intent(inout) :: p
        integer, intent(in) :: nrep
        real(wp), intent(out) :: t(2)
        type(daqp_workspace) :: work
        type(daqp_result) :: res
        real(wp), allocatable :: x(:), lam(:)
        integer :: k
        integer(ip) :: flag
        real(wp) :: t0, t1
        allocate(x(p%n), lam(p%m))
        t0 = wtime()
        do k = 1, nrep  ! (freed each time, as the C workspace)
            work%settings = fsettings
            flag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, H=p%H, f=p%f, At=p%At)
            call daqp_destroy(work)
        end do
        t1 = wtime()
        t(1) = (t1-t0)/nrep
        t0 = wtime()
        do k = 1, nrep
            work%settings = fsettings
            flag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, H=p%H, f=p%f, At=p%At)
            call daqp_solve(work, x, res, lam)
            call daqp_destroy(work)
        end do
        t1 = wtime()
        t(2) = (t1-t0)/nrep - t(1)
    end subroutine time_fortran

    subroutine time_fortran_warm(p, nrep, t)
        type(qp_problem), intent(inout) :: p
        integer, intent(in) :: nrep
        real(wp), intent(out) :: t
        type(daqp_workspace) :: work
        type(daqp_result) :: res
        real(wp), allocatable :: x(:), lam(:), f1(:), f2(:)
        integer :: k
        integer(ip) :: flag
        real(wp) :: t0
        allocate(x(p%n), lam(p%m))
        f1 = p%f
        f2 = p%f*1.05_wp
        work%settings = fsettings
        flag = daqp_setup(work, p%n, p%m, p%ms, p%bu, p%bl, H=p%H, f=p%f, At=p%At)
        call daqp_solve(work, x, res, lam)
        t0 = wtime()
        do k = 1, nrep
            if (mod(k,2) == 0) then
                work%qp%f = f1
            else
                work%qp%f = f2
            end if
            flag = daqp_update_ldp(work, daqp_update_v)
            call daqp_solve(work, x, res, lam)
        end do
        t = (wtime()-t0)/nrep
        call daqp_destroy(work)
    end subroutine time_fortran_warm

    subroutine time_c(p, nrep, t)
        type(qp_problem), intent(inout), target :: p
        integer, intent(in) :: nrep
        real(wp), intent(out) :: t(2)
        type(daqp_c_problem), target :: cqp
        type(daqp_c_result) :: cres
        real(c_double), allocatable, target :: x(:), lam(:)
        type(c_ptr) :: cwork
        integer :: k
        integer(c_int) :: flag
        real(wp) :: t0, t1
        allocate(x(p%n), lam(p%m))
        cres%x = c_loc(x)
        cres%lam = c_loc(lam)
        call c_problem(p, cqp)
        csettings = csettings0
        t0 = wtime()
        do k = 1, nrep
            cwork = cmp_ws_new(csettings)
            flag = setup_daqp(cqp, cwork, c_null_ptr)
            call cmp_ws_free(cwork)
        end do
        t1 = wtime()
        t(1) = (t1-t0)/nrep
        t0 = wtime()
        do k = 1, nrep
            cwork = cmp_ws_new(csettings)
            flag = setup_daqp(cqp, cwork, c_null_ptr)
            call daqp_c_solve(cres, cwork)
            call cmp_ws_free(cwork)
        end do
        t1 = wtime()
        t(2) = (t1-t0)/nrep - t(1)
    end subroutine time_c

    subroutine time_c_warm(p, nrep, t)
        type(qp_problem), intent(inout), target :: p
        integer, intent(in) :: nrep
        real(wp), intent(out) :: t
        type(daqp_c_problem), target :: cqp
        type(daqp_c_result) :: cres
        real(c_double), allocatable, target :: x(:), lam(:), f1(:), f2(:)
        type(c_ptr) :: cwork
        integer :: k
        integer(c_int) :: flag
        real(wp) :: t0
        allocate(x(p%n), lam(p%m))
        cres%x = c_loc(x)
        cres%lam = c_loc(lam)
        f1 = p%fc
        f2 = p%fc*1.05_wp
        call c_problem(p, cqp)
        csettings = csettings0
        cwork = cmp_ws_new(csettings)
        flag = setup_daqp(cqp, cwork, c_null_ptr)
        call daqp_c_solve(cres, cwork)
        t0 = wtime()
        do k = 1, nrep
            if (mod(k,2) == 0) then
                p%fc = f1
            else
                p%fc = f2
            end if
            flag = daqp_c_update_ldp(int(daqp_update_v,c_int), cwork, cqp)
            call daqp_c_solve(cres, cwork)
        end do
        t = (wtime()-t0)/nrep
        call cmp_ws_free(cwork)
        p%fc = f1
    end subroutine time_c_warm

!*****************************************************************************************
    end program compare
!*****************************************************************************************
