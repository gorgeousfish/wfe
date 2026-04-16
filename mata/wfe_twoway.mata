// ============================================================
// wfe_twoway.mata -- Two-way FE orchestration for DiD/Mdid estimators
// ============================================================

version 16.0
mata:
mata set matastrict on

void _wfe_tw_gmm_bread_guard(
    real vector Y,
    real matrix X,
    real vector W_it,
    real vector unit_idx,
    real vector time_idx,
    real scalar N_units,
    real scalar N_times
)
{
    real scalar NT, p, g, start, stop, N_units_eff, N_times_eff
    real matrix demeaned, X_dm, U, X_g, X_weighted, X_use
    real matrix panel
    real colvector Y_col, W_col, unit_col, time_col, keep_index
    real colvector unit_support_count, time_support_count
    real colvector Y_sorted, W_sorted, unit_sorted, time_sorted, sort_order

    Y_col = colshape(Y, 1)
    W_col = colshape(W_it, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    X_use = X

    NT = rows(Y_col)
    if (rows(X_use) == 1 & cols(X_use) == NT & NT > 1) {
        X_use = X_use'
    }
    p = cols(X_use)

    keep_index = selectindex(W_col :!= 0)
    if (rows(keep_index) == 0) {
        errprintf("wfe_twoway: at least one non-zero weight is required for two-way covariance\n")
        _error(498)
    }
    if (rows(keep_index) != NT) {
        Y_col = Y_col[keep_index]
        X_use = X_use[keep_index, .]
        W_col = W_col[keep_index]
        unit_col = unit_col[keep_index]
        time_col = time_col[keep_index]
        NT = rows(Y_col)
    }

    unit_support_count = J(N_units, 1, 0)
    time_support_count = J(N_times, 1, 0)
    for (g = 1; g <= NT; g++) {
        unit_support_count[unit_col[g]] = 1
        time_support_count[time_col[g]] = 1
    }
    N_units_eff = sum(unit_support_count)
    N_times_eff = sum(time_support_count)

    unit_col = _wfe_gmm_compact_index(unit_col, N_units)
    time_col = _wfe_gmm_compact_index(time_col, N_times)

    if (NT <= 1) {
        sort_order = 1::NT
    }
    else {
        sort_order = order((unit_col, time_col, (1::NT)), (1, 2, 3))
    }

    Y_sorted = Y_col[sort_order]
    X_use = X_use[sort_order, .]
    W_sorted = W_col[sort_order]
    unit_sorted = unit_col[sort_order]
    time_sorted = time_col[sort_order]

    demeaned = _wfe_twoway_demean(Y_sorted, X_use, unit_sorted, time_sorted,
        N_units_eff, N_times_eff)
    X_dm = demeaned[|1, 2 \ NT, p + 1|]
    panel = panelsetup(unit_sorted, 1)

    U = J(p, p, 0)
    for (g = 1; g <= rows(panel); g++) {
        start = panel[g, 1]
        stop = panel[g, 2]
        X_g = X_dm[|start, 1 \ stop, p|]
        X_weighted = X_g
        X_weighted = X_weighted :* W_sorted[|start \ stop|]
        U = U + cross(X_g, X_weighted)
    }

    if (rank(U / N_units_eff) < p) {
        errprintf("wfe_twoway: weighted GMM moment matrix is singular; two-way covariance is undefined\n")
        _error(498)
    }
}

void _wfe_twoway_estimate()
{
    real colvector Y, treat, unit_idx, time_idx, cit, W_vec, W_kept, beta_fe_col
    real colvector keep_index
    real colvector weight_cit
    real colvector Y_copy, unit_idx_copy, time_idx_copy
    real matrix X, exist, same, W, vcov_wfe, Psi_hat_wfe, Psi_wfe_for_white
    real matrix X_copy
    real scalar maxdev_did, has_maxdev_did, hetero_se, auto_se, df_adj, white, white_alpha
    real scalar unweighted, unbiased_se, verbose, store_wdm, tol
    real scalar N_units, N_times, NT, p, N_nonzero, N_negative, sigma, maxdev_active
    string scalar method, qoi, estimator, vcetype_str, white_test_str
    real rowvector white_result
    real colvector raw_scale, u_w
    struct wfe_twoway_fe_ols_result scalar fe_result
    struct wfe_complex_projection_state scalar cp_state
    struct wfe_complex_ols_result scalar cp_result
    struct wfe_twoway_fe_se_result scalar fe_se_result
    struct wfe_hac_white_result scalar hac_white

    _wfe_bridge_construct_views(Y, X, treat, unit_idx, time_idx, cit)
    _wfe_bridge_read_params(method, qoi, estimator, maxdev_did, has_maxdev_did,
        hetero_se, auto_se, df_adj, white, white_alpha,
        unweighted, unbiased_se, verbose, store_wdm, tol,
        N_units, N_times, NT, p, X)

    // Materialize copies from views for downstream processing
    Y_copy = Y[., .]
    X_copy = X[., .]
    unit_idx_copy = unit_idx[., .]
    time_idx_copy = time_idx[., .]

    if (!hetero_se & auto_se) {
        errprintf("Robust standard errors with autocorrelation and homoskedasticity is not supported\n")
        _error(3200)
    }
    if (!auto_se) {
        errprintf("two-way FE requires hetero_se(on) and auto_se(on)\n")
        _error(3200)
    }
    if (store_wdm) {
        errprintf("wfe_twoway: store_wdm is not supported for two-way estimators because the weighted-demeaned objects are complex-valued\n")
        _error(198)
    }

    exist = wfe_build_exist(unit_idx, time_idx, N_units, N_times, NT)
    same = wfe_build_same(unit_idx, time_idx, treat, N_units, N_times, NT)

    // Truncate C_it to integer before weight generation
    weight_cit = trunc(cit)

    maxdev_active = .
    if (estimator == "Mdid") {
        if (has_maxdev_did) {
            maxdev_active = maxdev_did
            if (verbose & !unweighted) {
                printf(": Matching on Pre-Treatment Outcome Within Maximum Deviation %g\n", maxdev_active)
            }
        }
        else {
            maxdev_active = -1
            if (verbose & !unweighted) {
                printf(": Nearest Neighbor Matching\n")
            }
        }
    }

    if (unweighted) {
        // Full matrix of ones for unweighted case
        W = J(N_times, N_units, 1)
    }
    else if (estimator == "did") {
        W = wfe_weights_did(unit_idx, time_idx, treat, weight_cit, N_units, N_times, NT, qoi, exist, same)
    }
    else if (estimator == "Mdid") {
        W = wfe_weights_mdid(unit_idx, time_idx, treat, weight_cit, Y, maxdev_active,
            N_units, N_times, NT, qoi, exist, same)
    }
    else {
        _error("wfe_twoway: unsupported two-way estimator")
    }

    W_vec = wfe_vectorize(W, time_idx, unit_idx, NT)
    N_nonzero = sum(W_vec :!= 0)
    N_negative = sum(W_vec :< 0)
    if (N_nonzero == 0) {
        errprintf("wfe_twoway: no non-zero DiD/MDiD weights; weighted two-way estimator is unidentified\n")
        _error(498)
    }
    keep_index = selectindex(W_vec :!= 0)

    // Always project on non-zero-weight support for correct coefficient estimation.
    // The real-domain HAC approach (wfe_se_hac_twoway_for_white) computes e_dm and X_w
    // for all NT rows independently, so keep_all=1 is no longer needed for the White test.
    raw_scale = colmax(abs(X_copy[keep_index, .]))'
    cp_state = _wfe_complex_project_core(Y_copy, X_copy, W_vec, unit_idx_copy, time_idx_copy, N_units, N_times, 0, tol)
    W_kept = W_vec[keep_index]
    cp_result = _wfe_complex_ols(cp_state.y_tilde, cp_state.X_tilde, cp_state.w_sqrt, W_kept, 0, raw_scale)
    _wfe_tw_gmm_bread_guard(Y_copy, X_copy, W_vec, unit_idx_copy,
        time_idx_copy, N_units, N_times)

    vcov_wfe = wfe_se_gmm(cp_result.coef_wls, Y_copy, X_copy, W_vec, unit_idx_copy, time_idx_copy,
        N_units, N_times, cols(X), df_adj, NT - sum(W_vec :== 0), Psi_hat_wfe)

    if (cp_result.sigma2 >= 0) {
        sigma = sqrt(cp_result.sigma2)
    }
    else {
        // Complex projection can produce Re(sum(z*z)) < 0 as a numerical
        // artifact when imaginary residual components dominate.  sigma is
        // only a summary statistic (RMSE); the SE uses GMM sandwich and
        // does not depend on sigma2.  Report 0 rather than crashing.
        sigma = 0
        if (verbose) {
            printf("{txt}(note: complex-projection sigma2 = %g < 0; sigma set to 0)\n",
                   cp_result.sigma2)
        }
    }
    vcetype_str = "Heteroscedastic / Autocorrelation Robust Standard Error"

    if (white) {
        // FE-side comparator objects are needed only for White's covariance-difference
        // test and should not constrain nowhite estimation paths.
        fe_result = wfe_twoway_fe_ols(Y_copy, X_copy, unit_idx_copy, time_idx_copy, N_units, N_times)
        fe_se_result = wfe_se_fe_twoway(fe_result.X_hat, fe_result.u_hat, unit_idx_copy,
            N_units, fe_result.ginv_XX_hat, time_idx_copy, N_times)
        // White's covariance-difference test must use the same GMM-side Psi_wfe
        // that defines the posted two-way WFE covariance, while the White helper
        // still needs real-domain X_w/u_w inputs for the cross-HAC Lambda term.
        hac_white = wfe_se_hac_twoway_for_white(cp_result.coef_wls,
            Y_copy, X_copy, W_vec, unit_idx_copy, time_idx_copy, N_units, N_times)
        Psi_wfe_for_white = Psi_hat_wfe
        u_w = hac_white.e_dm :* sqrt(abs(W_vec)) :* sign(W_vec)

        beta_fe_col = fe_result.beta_fe'
        white_result = wfe_white_test_twoway(beta_fe_col, cp_result.coef_wls,
            fe_result.X_hat, fe_result.u_hat, fe_result.ginv_XX_hat,
            hac_white.X_w, u_w, hac_white.ginv_XX_w,
            Psi_wfe_for_white, fe_se_result.Psi_hat_fe,
            unit_idx_copy, N_units, white_alpha,
            W_vec, time_idx_copy, N_times)
        white_test_str = (white_result[3] ? "TRUE" : "FALSE")

        _wfe_bridge_post_results(cp_result.coef_wls', vcov_wfe, W,
            N_nonzero, cp_result.d_f, sigma, vcetype_str,
            fe_result.beta_fe, fe_se_result.var_cov_fe,
            white_result[1], white_result[2], white_test_str)
    }
    else {
        _wfe_bridge_post_results(cp_result.coef_wls', vcov_wfe, W,
            N_nonzero, cp_result.d_f, sigma, vcetype_str)
    }

    st_numscalar("__wfe_sigma2", cp_result.sigma2)
    st_numscalar("__wfe_N_negative", N_negative)
    if (estimator == "Mdid") {
        st_numscalar("__wfe_maxdev_did", maxdev_active)
    }
}

end
