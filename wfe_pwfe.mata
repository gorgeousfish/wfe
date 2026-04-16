// ============================================================
// wfe_pwfe.mata -- PWFE end-to-end orchestration helpers
//
// This layer composes transform, weighting, demeaning, OLS,
// standard-error, and White-test helpers into the pwfe command path.
// ============================================================

version 16.0

mata:
mata set matastrict on

void _wfe_pwfe_len_count(
    real scalar   len_data,
    string scalar caller
)
{
    if (len_data >= . | len_data != floor(len_data) | len_data < 0) {
        errprintf("%s: len_data must be a nonnegative integer\n", caller)
        _error(3200)
    }
}

real matrix _wfe_pwfe_as_design(real matrix X)
{
    if (rows(X) == 1 & cols(X) > 1) {
        return(X')
    }

    return(X)
}

void _wfe_pwfe_validate_dm_idx(
    real vector    unit_idx,
    real scalar    unit_number,
    real scalar    len_data,
    string scalar  caller
)
{
    real scalar i
    real colvector seen

    unit_idx = colshape(unit_idx[., .], 1)

    _wfe_pwfe_len_count(len_data, "_wfe_pwfe_validate_demean_index")

    if (rows(unit_idx) != len_data) {
        errprintf("_wfe_pwfe_validate_demean_index: unit_idx length must equal len_data\n")
        _error(3200)
    }

    if (unit_number >= . | unit_number != floor(unit_number) | unit_number <= 0) {
        errprintf("%s: unit_number must be a positive integer\n", caller)
        _error(3200)
    }

    seen = J(unit_number, 1, 0)
    for (i = 1; i <= len_data; i++) {
        if (unit_idx[i] >= . | unit_idx[i] != floor(unit_idx[i]) | unit_idx[i] < 1 | unit_idx[i] > unit_number) {
            errprintf("%s: unit_idx must contain integer indices in 1..unit_number\n",
                      caller)
            _error(3200)
        }
        seen[unit_idx[i]] = 1
    }

    if (any(seen :== 0)) {
        errprintf("%s: unit_idx must enumerate 1..unit_number without gaps\n",
                  caller)
        _error(3200)
    }
}

void _wfe_pwfe_validate_demean_index(
    real vector    unit_idx,
    real scalar    unit_number,
    real scalar    len_data
)
{
    _wfe_pwfe_validate_dm_idx(unit_idx, unit_number, len_data,
                              "_wfe_pwfe_validate_demean_index")
}

real colvector _wfe_pwfe_demean_bug(
    real vector    var,
    real vector    unit_idx,
    real scalar    unit_number,
    real scalar    len_data
)
{
    return(_wfe_pwfe_demean_bug_core(var, unit_idx, unit_number, len_data,
                                     "_wfe_pwfe_demean_bug"))
}


real colvector _wfe_pwfe_demean_bug_core(
    real vector    var,
    real vector    unit_idx,
    real scalar    unit_number,
    real scalar    len_data,
    string scalar  caller
)
{
    real colvector result, idx
    real scalar i, mean_i

    var = colshape(var[., .], 1)
    unit_idx = colshape(unit_idx[., .], 1)

    _wfe_pwfe_len_count(len_data, caller)

    if (rows(var) != len_data | rows(unit_idx) != len_data) {
        errprintf("%s: var/unit_idx length mismatch\n", caller)
        _error(3200)
    }
    if (len_data == 0) {
        return(J(0, 1, .))
    }

    if (any(var :>= .)) {
        errprintf("%s: var must not contain missing values\n", caller)
        _error(3498)
    }

    if (unit_number >= . | unit_number != floor(unit_number) | unit_number <= 0) {
        errprintf("%s: unit_number must be a positive integer\n", caller)
        _error(3200)
    }

    _wfe_pwfe_validate_dm_idx(unit_idx, unit_number, len_data, caller)

    result = J(len_data, 1, 0)

    for (i = 1; i <= unit_number; i++) {
        idx = selectindex(unit_idx :== i)
        if (rows(idx) > 0) {
            mean_i = sum(var[idx]) / rows(idx)
            result[idx] = var[idx] :- mean_i
        }
    }

    return(result)
}


real matrix _wfe_pwfe_demean_bug_matrix(
    real matrix    data,
    real vector    unit_idx,
    real scalar    unit_number
)
{
    real matrix result
    real scalar j, p

    unit_idx = colshape(unit_idx[., .], 1)
    if (rows(data) == 1 & cols(data) == rows(unit_idx) & rows(unit_idx) > 1) {
        data = data'
    }

    p = cols(data)
    if (rows(data) == 0) {
        return(J(0, p, .))
    }
    if (p == 0) {
        errprintf("_wfe_pwfe_demean_bug_matrix: data must contain at least one column\n")
        _error(3200)
    }

    result = J(rows(data), p, 0)
    for (j = 1; j <= p; j++) {
        result[., j] = _wfe_pwfe_demean_bug_core(data[., j], unit_idx,
                                                 unit_number, rows(data),
                                                 "_wfe_pwfe_demean_bug_matrix")
    }

    return(result)
}


real colvector _wfe_pwfe_compact_tidx(
    real vector    unit_idx,
    real vector    time_idx
)
{
    real colvector compact_time
    real colvector compact_sorted
    real colvector sort_order
    real colvector time_sorted
    real colvector unit_col
    real colvector time_col
    real matrix info
    real scalar i, start_i, stop_i

    unit_col = colshape(unit_idx[., .], 1)
    time_col = colshape(time_idx[., .], 1)

    if (rows(unit_col) != rows(time_col)) {
        errprintf("_wfe_pwfe_compact_tidx: unit_idx/time_idx length mismatch\n")
        _error(3200)
    }
    if (rows(unit_col) == 0) {
        return(J(0, 1, .))
    }

    for (i = 1; i <= rows(unit_col); i++) {
        if (unit_col[i] >= . | unit_col[i] != floor(unit_col[i]) | unit_col[i] < 1) {
            errprintf("_wfe_pwfe_compact_tidx: unit_idx must contain positive integer indices\n")
            _error(3200)
        }
        if (time_col[i] >= . | time_col[i] != floor(time_col[i]) | time_col[i] < 1) {
            errprintf("_wfe_pwfe_compact_tidx: time_idx must contain positive integer indices\n")
            _error(3200)
        }
    }

    sort_order = order((unit_col, (1::rows(unit_col))), (1, 2))
    info = panelsetup(unit_col[sort_order], 1)
    time_sorted = time_col[sort_order]
    compact_time = J(rows(time_col), 1, .)
    compact_sorted = J(rows(time_col), 1, .)

    for (i = 1; i <= rows(info); i++) {
        start_i = info[i, 1]
        stop_i = info[i, 2]
        if (rows(uniqrows(time_sorted[|start_i \ stop_i|])) != (stop_i - start_i + 1)) {
            errprintf("_wfe_pwfe_compact_tidx: time_idx must be unique within unit\n")
            _error(498)
        }
        compact_sorted[|start_i \ stop_i|] = (1::(stop_i - start_i + 1))
    }

    compact_time[sort_order] = compact_sorted
    return(compact_time)
}


void _wfe_pwfe_white_psd_mat(
    real matrix    X,
    string scalar  argname
)
{
    real rowvector evals
    real scalar eig_tol, eig_scale

    symeigensystem(X, ., evals)
    eig_scale = max((1, max(abs(evals))))
    eig_tol = sqrt(epsilon(1)) * eig_scale

    if (min(evals) < -eig_tol) {
        errprintf("_wfe_pwfe_white_command: %s must be positive semidefinite\n",
                  argname)
        _error(3200)
    }
}



real rowvector _wfe_pwfe_white_command(
    real vector    beta_wfe,
    real vector    beta_fe,
    real matrix    Psi_wfe,
    real matrix    Psi_fe,
    real vector    u_hat,
    real vector    u_tilde,
    real matrix    X_hat,
    real matrix    X_tilde,
    real matrix    ginv_XX_hat,
    real matrix    ginv_XX_tilde,
    real scalar    white_alpha
)
{
    real scalar NT, p, white_stat, white_pvalue, test_reject, white_tol
    real scalar psi_wfe_asym, psi_wfe_scale, psi_fe_asym, psi_fe_scale
    real scalar ginv_hat_asym, ginv_hat_scale, ginv_tilde_asym, ginv_tilde_scale
    real scalar ginv_hat_gap, ginv_tilde_gap, ginv_match_scale, ginv_match_tol
    real colvector beta_wfe_col, beta_fe_col, u_hat_col, u_tilde_col
    real colvector diag_ee, beta_diff
    real matrix Lambda1, Lambda2, Phi, Psi_wfe_use, Psi_fe_use
    real matrix ginv_XX_hat_use, ginv_XX_tilde_use
    real matrix ginv_hat_expected, ginv_tilde_expected

    beta_wfe_col = colshape(beta_wfe[., .], 1)
    beta_fe_col = colshape(beta_fe[., .], 1)
    u_hat_col = colshape(u_hat[., .], 1)
    u_tilde_col = colshape(u_tilde[., .], 1)

    // Canonicalize only genuine single-regressor rowmatrices. A true
    // 1 x p design is one observation with p regressors and must keep
    // that geometry instead of being transposed into p observations.
    if (!(rows(X_hat) == 1 & cols(X_hat) > 1 ///
          & rows(u_hat_col) == 1 & rows(u_tilde_col) == 1)) {
        X_hat = _wfe_pwfe_as_design(X_hat)
    }
    if (!(rows(X_tilde) == 1 & cols(X_tilde) > 1 ///
          & rows(u_hat_col) == 1 & rows(u_tilde_col) == 1)) {
        X_tilde = _wfe_pwfe_as_design(X_tilde)
    }

    NT = rows(X_hat)
    p = cols(X_hat)

    if (white_alpha <= 0 | white_alpha >= 1) {
        errprintf("_wfe_pwfe_white_command: white_alpha must lie in (0,1)\n")
        _error(3200)
    }

    if (rows(X_hat) == 0 | rows(X_tilde) == 0) {
        errprintf("_wfe_pwfe_white_command: X_hat and X_tilde must each contain at least one observation\n")
        _error(3200)
    }
    if (rows(X_hat) != rows(X_tilde) | rows(X_hat) != rows(u_hat_col) | rows(X_tilde) != rows(u_tilde_col)) {
        errprintf("_wfe_pwfe_white_command: X_hat, X_tilde, u_hat, and u_tilde must align by row\n")
        _error(3200)
    }
    if (cols(X_hat) == 0 | cols(X_tilde) == 0) {
        errprintf("_wfe_pwfe_white_command: X_hat and X_tilde must each contain at least one regressor\n")
        _error(3200)
    }
    if (cols(X_hat) != cols(X_tilde)) {
        errprintf("_wfe_pwfe_white_command: X_hat and X_tilde must have the same number of regressors\n")
        _error(3200)
    }
    if (rows(beta_wfe_col) != p | cols(beta_wfe_col) != 1 | rows(beta_fe_col) != p | cols(beta_fe_col) != 1) {
        errprintf("_wfe_pwfe_white_command: beta_wfe and beta_fe must be %gx1 vectors\n", p)
        _error(3200)
    }
    if (rows(Psi_wfe) != p | cols(Psi_wfe) != p) {
        errprintf("_wfe_pwfe_white_command: Psi_wfe must be %gx%g to match cols(X_hat)\n", p, p)
        _error(3200)
    }
    if (rows(Psi_fe) != p | cols(Psi_fe) != p) {
        errprintf("_wfe_pwfe_white_command: Psi_fe must be %gx%g to match cols(X_hat)\n", p, p)
        _error(3200)
    }
    if (rows(ginv_XX_hat) != p | cols(ginv_XX_hat) != p) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_hat must be %gx%g to match cols(X_hat)\n", p, p)
        _error(3200)
    }
    if (rows(ginv_XX_tilde) != p | cols(ginv_XX_tilde) != p) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_tilde must be %gx%g to match cols(X_hat)\n", p, p)
        _error(3200)
    }
    if (any(beta_wfe_col :>= .) | any(beta_fe_col :>= .)) {
        errprintf("_wfe_pwfe_white_command: beta vectors must not contain missing values\n")
        _error(3498)
    }
    if (any(Psi_wfe :>= .)) {
        errprintf("_wfe_pwfe_white_command: Psi_wfe must not contain missing values\n")
        _error(3498)
    }
    if (any(Psi_fe :>= .)) {
        errprintf("_wfe_pwfe_white_command: Psi_fe must not contain missing values\n")
        _error(3498)
    }
    psi_wfe_asym = max(abs(vec(Psi_wfe - Psi_wfe')))
    psi_wfe_scale = max((1, max(abs(vec(Psi_wfe)))))
    if (psi_wfe_asym > sqrt(epsilon(1)) * psi_wfe_scale) {
        errprintf("_wfe_pwfe_white_command: Psi_wfe must be symmetric\n")
        _error(3200)
    }
    Psi_wfe_use = 0.5 :* (Psi_wfe + Psi_wfe')
    _wfe_pwfe_white_psd_mat(Psi_wfe_use, "Psi_wfe")

    psi_fe_asym = max(abs(vec(Psi_fe - Psi_fe')))
    psi_fe_scale = max((1, max(abs(vec(Psi_fe)))))
    if (psi_fe_asym > sqrt(epsilon(1)) * psi_fe_scale) {
        errprintf("_wfe_pwfe_white_command: Psi_fe must be symmetric\n")
        _error(3200)
    }
    Psi_fe_use = 0.5 :* (Psi_fe + Psi_fe')
    _wfe_pwfe_white_psd_mat(Psi_fe_use, "Psi_fe")
    if (any(u_hat_col :>= .) | any(u_tilde_col :>= .)) {
        errprintf("_wfe_pwfe_white_command: u_hat and u_tilde must not contain missing values\n")
        _error(3498)
    }
    if (any(X_hat :>= .) | any(X_tilde :>= .)) {
        errprintf("_wfe_pwfe_white_command: X_hat and X_tilde must not contain missing values\n")
        _error(3498)
    }
    if (any(ginv_XX_hat :>= .)) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_hat must not contain missing values\n")
        _error(3498)
    }
    if (any(ginv_XX_tilde :>= .)) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_tilde must not contain missing values\n")
        _error(3498)
    }
    ginv_hat_asym = max(abs(vec(ginv_XX_hat - ginv_XX_hat')))
    ginv_hat_scale = max((1, max(abs(vec(ginv_XX_hat)))))
    if (ginv_hat_asym > sqrt(epsilon(1)) * ginv_hat_scale) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_hat must be symmetric\n")
        _error(3200)
    }
    ginv_XX_hat_use = 0.5 :* (ginv_XX_hat + ginv_XX_hat')
    _wfe_pwfe_white_psd_mat(ginv_XX_hat_use, "ginv_XX_hat")
    ginv_hat_expected = pinv(cross(X_hat, X_hat))
    ginv_hat_gap = max(abs(vec(ginv_XX_hat_use - ginv_hat_expected)))
    ginv_match_scale = max((1, max(abs(vec(ginv_hat_expected))),
        max(abs(vec(ginv_XX_hat_use)))))
    ginv_match_tol = sqrt(epsilon(1)) * rows(X_hat) * ginv_match_scale
    if (ginv_hat_gap > ginv_match_tol) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_hat must match pinv(X_hat'X_hat)\n")
        _error(498)
    }

    ginv_tilde_asym = max(abs(vec(ginv_XX_tilde - ginv_XX_tilde')))
    ginv_tilde_scale = max((1, max(abs(vec(ginv_XX_tilde)))))
    if (ginv_tilde_asym > sqrt(epsilon(1)) * ginv_tilde_scale) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_tilde must be symmetric\n")
        _error(3200)
    }
    ginv_XX_tilde_use = 0.5 :* (ginv_XX_tilde + ginv_XX_tilde')
    _wfe_pwfe_white_psd_mat(ginv_XX_tilde_use, "ginv_XX_tilde")
    ginv_tilde_expected = pinv(cross(X_tilde, X_tilde))
    ginv_tilde_gap = max(abs(vec(ginv_XX_tilde_use - ginv_tilde_expected)))
    ginv_match_scale = max((1, max(abs(vec(ginv_tilde_expected))),
        max(abs(vec(ginv_XX_tilde_use)))))
    ginv_match_tol = sqrt(epsilon(1)) * rows(X_tilde) * ginv_match_scale
    if (ginv_tilde_gap > ginv_match_tol) {
        errprintf("_wfe_pwfe_white_command: ginv_XX_tilde must match pinv(X_tilde'X_tilde)\n")
        _error(498)
    }

    diag_ee = u_hat_col :* u_tilde_col
    Lambda1 = (1 / rows(X_hat)) * cross(X_hat :* diag_ee, X_tilde)
    Lambda2 = (1 / rows(X_tilde)) * cross(X_tilde :* diag_ee, X_hat)
    Phi = Psi_wfe_use + Psi_fe_use ///
        - (rows(X_hat) * ginv_XX_hat_use) * Lambda1 * (rows(X_tilde) * ginv_XX_tilde_use) ///
        - (rows(X_tilde) * ginv_XX_tilde_use) * Lambda2 * (rows(X_hat) * ginv_XX_hat_use)

    beta_diff = beta_fe_col - beta_wfe_col
    white_stat = Re(NT * beta_diff' * pinv(Phi) * beta_diff)
    white_tol = sqrt(epsilon(1)) * max((1, abs(white_stat)))
    if (white_stat < 0) {
        if (abs(white_stat) <= white_tol) {
            white_stat = 0
        }
        else {
            errprintf("Warning: _wfe_pwfe_white_command: White statistic is negative (%g); returning p=1.\n", white_stat)
            return((white_stat, 1, 0))
        }
    }
    white_pvalue = chi2tail(p, white_stat)
    test_reject = (white_pvalue < white_alpha)

    return((white_stat, white_pvalue, test_reject))
}


void _wfe_pwfe_estimate()
{
    real colvector y_star, outcome_raw, treat, unit_idx, time_idx, cit, W_vec
    real colvector weight_cit
    real matrix W, exist, YX_final, Data_wdm, Data_dm
    real matrix X_tilde, X_hat, vcov_ols, vcov_wfe, vcov_fe
    real matrix Psi_wfe, Psi_fe, ginv_XX_tilde, ginv_XX_hat
    real colvector u_tilde, u_hat
    real scalar J_u, J_t, NT, unit_number, fit_final_df, df_r, sigma2, sigma2_fe
    real scalar N_nonzero
    real scalar int_max
    real scalar white_alpha, white_stat, white_pvalue, white_reject
    real rowvector white_result
    string scalar method, qoi, estimator, hetero_se, auto_se, unbiased_se
    string scalar verbose, white, se_type
    struct wfe_ols_result scalar wfe_result, fe_result

    st_view(y_star, ., st_local("pwfe_y_star"), st_local("touse"))
    st_view(outcome_raw, ., st_local("outcome"), st_local("touse"))
    st_view(treat, ., st_local("treat"), st_local("touse"))
    st_view(unit_idx, ., st_local("pwfe_unit_idx"), st_local("touse"))
    st_view(time_idx, ., st_local("pwfe_time_idx"), st_local("touse"))
    st_view(cit, ., st_local("pwfe_cit"), st_local("touse"))

    method = st_local("method")
    qoi = st_local("qoi")
    estimator = st_local("estimator")
    hetero_se = st_local("hetero_se")
    auto_se = st_local("auto_se")
    unbiased_se = st_local("pwfe_unbiased_se")
    verbose = st_local("pwfe_verbose")
    white = st_local("pwfe_white")

    J_u = st_numscalar("__pwfe_J_u")
    J_t = st_numscalar("__pwfe_J_t")
    unit_number = st_numscalar("__pwfe_unit_number")
    white_alpha = st_numscalar("__pwfe_white_alpha")
    NT = rows(y_star)

    // Re-compact the omitted-time clock to the dense within-sample support.
    // Under method(time), downstream demeaning and df accounting must use the
    // same compacted time-group count rather than the caller's stale precompact
    // scalar.
    if (st_local("time") == "") {
        time_idx = _wfe_pwfe_compact_tidx(unit_idx, time_idx)
        J_t = max(time_idx)
        if (method == "time") {
            unit_number = J_t
        }
    }

    // Guard against 32-bit integer overflow for C_it values
    weight_cit = trunc(cit)
    int_max = 2147483647
    if (any(weight_cit :> int_max)) {
        errprintf("_wfe_pwfe_estimate: C_it must lie within the signed 32-bit integer range after truncation\n")
        _error(3498)
    }

    exist = wfe_build_exist(unit_idx, time_idx, J_u, J_t, NT)

    // Keep the public pwfe command contract specific about which one-way path
    // is infeasible, while leaving the low-level weight helpers free to mirror
    // the C reference generators' pure zero-weight output on degenerate
    // support.
    if (estimator == "fd" & sum(rowsum(exist) :> 0) < 2) {
        errprintf("wfe_weights_fd: FD requires at least 2 time periods\n")
        _error(3200)
    }
    if (estimator != "fd" & method == "unit" & sum(rowsum(exist) :> 0) < 2) {
        errprintf("wfe_weights_unit: unit FE requires at least 2 time periods\n")
        _error(3200)
    }
    if (method == "time" & sum(colsum(exist) :> 0) < 2) {
        errprintf("wfe_weights_time: time FE requires at least 2 units\n")
        _error(3200)
    }

    if (estimator == "fd") {
        W = wfe_weights_fd(unit_idx, time_idx, treat, weight_cit, J_u, J_t, NT, qoi, exist)
    }
    else if (method == "unit") {
        W = wfe_weights_unit(unit_idx, time_idx, treat, weight_cit, J_u, J_t, NT, qoi, exist)
    }
    else {
        W = wfe_weights_time(unit_idx, time_idx, treat, weight_cit, J_u, J_t, NT, qoi, exist)
    }

    W_vec = wfe_vectorize(W, time_idx, unit_idx, NT)
    N_nonzero = sum(W_vec :!= 0)

    if (N_nonzero == 0) {
        errprintf("no non-zero pwfe one-way weights; estimator not identified\n")
        _error(498)
    }

    YX_final = y_star, treat

    if (method == "unit") {
        Data_wdm = wfe_wwdemean_matrix(YX_final, W_vec, unit_idx, unit_number, NT)
        Data_dm = wfe_demean_matrix(YX_final, unit_idx, unit_number, NT)
    }
    else {
        Data_wdm = wfe_wwdemean_matrix_unsorted(YX_final, W_vec, time_idx, unit_number, NT)
        /*
           The public FE benchmark must stay on the same fixed-effect axis as
           the weighted PWFE regression. Under method(time), both the reported
           FE comparator and the White cross-covariance therefore use time
           demeaning to maintain consistent comparisons.
        */
        Data_dm = wfe_demean_matrix_unsorted(YX_final, time_idx, unit_number, NT)
    }

    X_tilde = Data_wdm[., 2]
    X_hat = Data_dm[., 2]

    wfe_result = _wfe_ols_core(Data_wdm[., 1], X_tilde)
    fe_result = _wfe_ols_core(Data_dm[., 1], X_hat)

    fit_final_df = NT - cols(X_tilde)
    df_r = fit_final_df - unit_number
    if (df_r <= 0) {
        errprintf("_wfe_pwfe_estimate: fit_final_df - unit_number must be positive\n")
        _error(498)
    }

    sigma2 = quadcross(wfe_result.resid, wfe_result.resid) / df_r
    sigma2_fe = quadcross(fe_result.resid, fe_result.resid) / fit_final_df
    vcov_ols = sigma2_fe * fe_result.ginv_XX

    u_tilde = sqrt(W_vec) :* wfe_result.resid
    u_hat = fe_result.resid

    _wfe_pwfe_compute_se(X_tilde, u_tilde, unit_idx,
                         X_hat, u_hat, unit_idx,
                         J_u, J_t, wfe_result.resid, fit_final_df, unit_number, vcov_ols,
                         hetero_se, auto_se, unbiased_se, verbose,
                         vcov_wfe, vcov_fe, Psi_wfe, Psi_fe,
                         ginv_XX_tilde, ginv_XX_hat, se_type)

    st_matrix("__pwfe_b", wfe_result.beta)
    st_matrix("__pwfe_V", vcov_wfe)
    st_matrix("__pwfe_W", W)
    st_matrix("__pwfe_b_fe", fe_result.beta)
    st_matrix("__pwfe_V_fe", vcov_fe)
    st_numscalar("__pwfe_df_r", df_r)
    st_numscalar("__pwfe_sigma", sqrt(sigma2))
    st_numscalar("__pwfe_sigma2", sigma2)
    st_numscalar("__pwfe_N_units", J_u)
    st_numscalar("__pwfe_N_times", J_t)
    st_global("__pwfe_vcetype", se_type)

    if (white == "on") {
        white_result = _wfe_pwfe_white_command(wfe_result.beta', fe_result.beta',
                                               Psi_wfe, Psi_fe, u_hat, u_tilde,
                                               X_hat, X_tilde, ginv_XX_hat,
                                               ginv_XX_tilde, white_alpha)
        white_stat = white_result[1]
        white_pvalue = white_result[2]
        white_reject = white_result[3]

        st_numscalar("__pwfe_white_stat", white_stat)
        st_numscalar("__pwfe_white_pvalue", white_pvalue)
        st_global("__pwfe_white_test", (white_reject ? "TRUE" : "FALSE"))
    }
}

end
