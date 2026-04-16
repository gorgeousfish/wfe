// ============================================================
// wfe_oneway.mata — One-way FE Estimation Coordinator
//
// Calls: wfe_bridge (I/O), weights_unit/time/fd,
//        vectorize, wwdemean, demean, se_hac, white_test
//
// Workflow:
//   Step 0. Parameter loading + constants + fail-fast validation
//   Step 1. Build exist matrix
//   Step 2. Weight calculation routing (unit/time/fd/unweighted)
//   Step 3. Vectorize weights
//   Step 4. WFE weighted demeaning
//   Step 5. Standard FE demeaning
//   Step 6. Dual OLS regression
//   Step 7. Standard error computation
//   Step 8. Return results
// ============================================================

version 16.0
mata:
mata set matastrict on

// ============================================================
// _wfe_check_balanced() — Panel balance check
//
// Validates that all units have the same number of observations.
// Used for Stock-Watson unbiased SE (panel must be balanced).
//
// Requirement: unit_idx must be sorted ascending.
// ============================================================
void _wfe_check_balanced(real colvector unit_idx,
                         real scalar    NT,
                         real scalar    J_u)
{
    real matrix info
    real colvector counts
    real scalar expected_T

    if (J_u == 0) return

    info = panelsetup(unit_idx, 1)
    counts = info[., 2] - info[., 1] :+ 1

    if (min(counts) != max(counts)) {
        _error("unbiased_se(on) is allowed only when panel is balanced")
    }
}


// ============================================================
// _wfe_oneway_estimate() — One-way WFE Estimation Coordinator
//
// Zero-parameter entry function. Reads data and parameters from
// Stata namespace via _wfe_bridge_construct_views() and
// _wfe_bridge_read_params().
// ============================================================
void _wfe_oneway_estimate()
{
    // ══════════════════════════════════════════════════════════
    // Step 0: Data loading + parameters + constants + validation
    // ══════════════════════════════════════════════════════════

    // 0a. Data views (Bridge)
    real colvector Y, treat, unit_idx, time_idx, cit
    real matrix X
    _wfe_bridge_construct_views(Y, X, treat, unit_idx, time_idx, cit)

    // 0b. Parameters (Bridge)
    string scalar method, qoi, estimator
    real scalar maxdev_did, has_maxdev_did
    real scalar hetero_se, auto_se, df_adj, white, white_alpha
    real scalar unweighted, unbiased_se, verbose, store_wdm, tol
    real scalar N_units, N_times, NT, p
    _wfe_bridge_read_params(method, qoi, estimator,
        maxdev_did, has_maxdev_did,
        hetero_se, auto_se, df_adj, white, white_alpha,
        unweighted, unbiased_se, verbose, store_wdm, tol,
        N_units, N_times, NT, p, X)

    if (estimator == "did" | estimator == "Mdid") {
        errprintf("_wfe_oneway_estimate: estimator '%s' requires weighted two-way FE routing\n",
                  estimator)
        _error(3200)
    }

    // 0c. Runtime constants
    real scalar J_u, n_demean_groups, sigma_df_groups
    real colvector demean_grp_idx, weight_cit

    J_u = N_units

    if (method == "unit") {
        n_demean_groups = N_units
        sigma_df_groups = N_units
        demean_grp_idx  = unit_idx
    }
    else {
        n_demean_groups = N_times
        // R one-way time FE still uses J.u (= number of units) in
        // d.f = fit.final$df - J.u for sigma2/df_r.
        sigma_df_groups = N_units
        demean_grp_idx  = time_idx
    }

    if (rows(Y) != NT) {
        errprintf("_wfe_oneway_estimate: rows(Y)=%g != NT=%g\n", rows(Y), NT)
        _error(3200)
    }

    // 0d. SE combination fail-fast validation
    if (hetero_se == 0 & auto_se == 0) {
        _error("standard errors with independence and homoskedasticity is not supported")
    }
    if (hetero_se == 0 & auto_se == 1) {
        _error("robust standard errors with autocorrelation and homoskedasticity is not supported")
    }

    // 0e. verbose header
    if (verbose) {
        printf("{txt}Weight calculation started\n")
    }


    // Use a detached working copy; integer truncation must not
    // mutate the underlying Stata view before weighted generators consume C_it.
    weight_cit = cit
    if (unweighted == 0 & ((estimator == "" & (method == "unit" | method == "time")) | estimator == "fd")) {
        weight_cit = trunc(cit)
    }

    // ══════════════════════════════════════════════════════════
    // Step 1: Build exist matrix
    // ══════════════════════════════════════════════════════════
    real matrix exist
    exist = wfe_build_exist(unit_idx, time_idx, N_units, N_times, NT)


    // ══════════════════════════════════════════════════════════
    // Step 2: Weight calculation routing
    // ══════════════════════════════════════════════════════════
    real matrix W

    if (unweighted == 1) {
        W = exist
    }
    else if (estimator == "fd") {
        W = wfe_weights_fd(unit_idx, time_idx, treat, weight_cit,
                           N_units, N_times, NT, qoi, exist)
    }
    else if (method == "unit") {
        W = wfe_weights_unit(unit_idx, time_idx, treat, weight_cit,
                             N_units, N_times, NT, qoi, exist)
    }
    else {
        W = wfe_weights_time(unit_idx, time_idx, treat, weight_cit,
                             N_units, N_times, NT, qoi, exist)
    }

    exist = .


    // ══════════════════════════════════════════════════════════
    // Step 3: Vectorize weights
    // ══════════════════════════════════════════════════════════
    real colvector W_vec
    real scalar N_nonzero, N_negative

    W_vec = wfe_vectorize(W, time_idx, unit_idx, NT)

    N_nonzero  = sum(W_vec :!= 0)
    N_negative = sum(W_vec :< 0)

    if (N_nonzero == 0) {
        errprintf("no non-zero one-way weights; estimator not identified\n")
        _error(498)
    }

    if (verbose) {
        printf("{txt}Weight calculation done\n")
        if (estimator == "fd") {
            printf("{txt}Total number observations with non-zero weight: %g\n", N_nonzero)
        }
        if (N_negative > 0 & estimator == "") {
            printf("{err}Warning: negative weights detected in one-way WFE (implementation error suspected)\n")
        }
    }

    if (N_units * N_times > 50000 & verbose) {
        printf("{txt}Warning: large dataset (N*T > 50000), consider using nowhite to reduce memory usage\n")
    }


    // ══════════════════════════════════════════════════════════
    // Step 4-5: Dual demeaning (WFE weighted + standard FE)
    //
    // method="unit": data sorted by unit_idx → use sorted version
    // method="time": data sorted by (unit_idx, time_idx), need time_idx grouping
    //                 → use unsorted version (internal sort/unsort)
    // ══════════════════════════════════════════════════════════
    real colvector Y_tilde, Y_hat
    real matrix    X_tilde, X_hat
    real scalar    k

    X_tilde = J(NT, p, .)
    X_hat   = J(NT, p, .)

    if (method == "unit") {
        Y_tilde = wfe_wwdemean(Y, W_vec, unit_idx, N_units, NT)
        for (k = 1; k <= p; k++) {
            X_tilde[., k] = wfe_wwdemean(X[., k], W_vec, unit_idx,
                                          N_units, NT)
        }

        Y_hat = wfe_demean(Y, unit_idx, N_units, NT)
        for (k = 1; k <= p; k++) {
            X_hat[., k] = wfe_demean(X[., k], unit_idx, N_units, NT)
        }
    }
    else {
        Y_tilde = wfe_wwdemean_unsorted(Y, W_vec, time_idx, N_times, NT)
        for (k = 1; k <= p; k++) {
            X_tilde[., k] = wfe_wwdemean_unsorted(X[., k], W_vec, time_idx,
                                                   N_times, NT)
        }

        Y_hat = wfe_demean_unsorted(Y, time_idx, N_times, NT)
        for (k = 1; k <= p; k++) {
            X_hat[., k] = wfe_demean_unsorted(X[., k], time_idx, N_times, NT)
        }
    }


    // ══════════════════════════════════════════════════════════
    // Step 6: Dual OLS regression
    // ══════════════════════════════════════════════════════════

    struct wfe_ols_result scalar wfe_result, fe_result
    real matrix ginv_XX_tilde, ginv_XX_hat
    real colvector beta_wfe, beta_fe
    real colvector u_tilde, u_hat
    real scalar df_r, sigma2

    wfe_result = _wfe_ols_core(Y_tilde, X_tilde)
    fe_result  = _wfe_ols_core(Y_hat, X_hat)

    ginv_XX_tilde = wfe_result.ginv_XX
    ginv_XX_hat   = fe_result.ginv_XX
    beta_wfe      = wfe_result.beta'
    beta_fe       = fe_result.beta'
    u_tilde       = wfe_result.resid
    u_hat         = fe_result.resid

    _wfe_calc_sigma2(u_tilde, NT, p, sigma_df_groups, sigma2, df_r)

    real colvector diag_ee_tilde, diag_ee_hat
    diag_ee_tilde = u_tilde :^ 2
    diag_ee_hat   = u_hat   :^ 2

    if (verbose) {
        printf("{txt}OLS: p=%g, df_r=%g, sigma2=%g\n", p, df_r, sigma2)
    }


    // ══════════════════════════════════════════════════════════
    // Step 7: Standard error computation
    // ══════════════════════════════════════════════════════════
    real matrix vcov_wfe, vcov_fe, Psi_wfe, Psi_fe
    string scalar se_type

    string scalar se_hetero_str, se_auto_str, se_unbiased_str
    se_hetero_str   = (hetero_se   ? "on" : "off")
    se_auto_str     = (auto_se     ? "on" : "off")
    se_unbiased_str = (unbiased_se ? "on" : "off")

    real scalar is_balanced
    {
        real matrix _info_bal
        real colvector _counts_bal
        _info_bal = panelsetup(unit_idx, 1)
        _counts_bal = _info_bal[., 2] - _info_bal[., 1] :+ 1
        is_balanced = (min(_counts_bal) == max(_counts_bal))
    }

    _wfe_compute_se(
        X_tilde, u_tilde,
        X_hat, u_hat,
        unit_idx,
        J_u,
        p,
        N_times,
        N_nonzero,
        df_adj,
        se_hetero_str,
        se_auto_str,
        se_unbiased_str,
        is_balanced,
        vcov_wfe,
        vcov_fe,
        Psi_wfe,
        Psi_fe,
        ginv_XX_tilde,
        ginv_XX_hat,
        se_type
    )

    if (verbose) {
        printf("{txt}SE type: %s\n", se_type)
    }


    // ══════════════════════════════════════════════════════════
    // Step 8: White specification test
    // ══════════════════════════════════════════════════════════
    real scalar sigma_val
    sigma_val = sqrt(sigma2)

    if (white) {
        real rowvector white_result
        real scalar w_stat, w_pvalue, w_reject
        string scalar w_test_str

        white_result = wfe_white_test_oneway(
            beta_wfe,
            beta_fe,
            Psi_wfe,
            Psi_fe,
            ginv_XX_tilde,
            ginv_XX_hat,
            X_tilde,
            X_hat,
            u_tilde,
            u_hat,
            J_u,
            p,
            N_nonzero,
            white_alpha,
            unit_idx
        )

        w_stat    = white_result[1]
        w_pvalue  = white_result[2]
        w_reject  = white_result[3]
        w_test_str = (w_reject ? "TRUE" : "FALSE")

        if (verbose) {
            printf("{txt}White test: stat=%g, p=%g, reject=%s\n",
                   w_stat, w_pvalue, w_test_str)
        }
    }


    // ══════════════════════════════════════════════════════════
    // Step 9: Return results
    // ══════════════════════════════════════════════════════════

    _wfe_bridge_post_results(
        beta_wfe',
        vcov_wfe,
        W,
        N_nonzero,
        df_r,
        sigma_val,
        se_type,
        beta_fe',
        vcov_fe
    )

    // Store white test after bridge (bridge drops __wfe_white_stat at startup)
    if (white) {
        st_numscalar("__wfe_white_stat", w_stat)
        st_numscalar("__wfe_white_pvalue", w_pvalue)
        st_global("__wfe_white_test", w_test_str)
    }

    st_matrix("__wfe_b_fe", beta_fe')
    st_matrix("__wfe_V_fe", vcov_fe)

    if (store_wdm) {
        st_matrix("__wfe_Y_wdm", Y_tilde)
        st_matrix("__wfe_X_wdm", X_tilde)
    }

    st_matrix("__wfe_Psi_wfe",       Psi_wfe)
    st_matrix("__wfe_Psi_fe",        Psi_fe)
    st_matrix("__wfe_ginv_XX_tilde", ginv_XX_tilde)
    st_matrix("__wfe_ginv_XX_hat",   ginv_XX_hat)
    st_matrix("__wfe_X_tilde",       X_tilde)
    st_matrix("__wfe_X_hat",         X_hat)
    st_matrix("__wfe_u_tilde",       u_tilde)
    st_matrix("__wfe_u_hat",         u_hat)

    st_numscalar("__wfe_sigma2",     sigma2)
    st_numscalar("__wfe_J_u",        J_u)
    st_numscalar("__wfe_N_negative", N_negative)

    if (verbose) {
        printf("\n{txt}=== WFE Estimation Summary ===\n")
        printf("{txt}Method: %s\n", method)
        printf("{txt}QOI: %s\n", qoi)
        printf("{txt}Covariates: %g\n", p)
        printf("{txt}Observations: %g  (non-zero weight: %g)\n", NT, N_nonzero)
        printf("{txt}Units (J_u): %g\n", J_u)
        printf("{txt}Time periods: %g\n", N_times)
        printf("{txt}df_r: %g\n", df_r)
        printf("{txt}sigma: %g\n", sigma_val)
        printf("{txt}SE type: %s\n", se_type)
        printf("{txt}Negative weights: %g\n", N_negative)
    }
}

end
