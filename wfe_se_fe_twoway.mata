// wfe_se_fe_twoway.mata -- FE-side HAC sandwich for two-way FE
// Computes Omega_fe, Psi_hat_fe, and var_cov_fe for two-way fixed effects

version 16.0
mata:
mata set matastrict on

struct wfe_twoway_fe_se_result {
    real matrix Omega_fe
    real matrix Psi_hat_fe
    real matrix var_cov_fe
}

// Struct returned by wfe_se_hac_twoway_for_white — carries all real quantities
// needed for the all-real two-way White test Lambda computation.
struct wfe_hac_white_result {
    real matrix    Psi_wfe      // [p x p]  HAC-based Psi_wfe (sandwich)
    real colvector e_dm         // [NT x 1] FWL-WFE residuals (Y_dm - X_dm*beta)
    real matrix    X_w          // [NT x p] weighted FWL X = X_dm .* sqrt(|W|)
    real matrix    ginv_XX_w    // [p x p]  pinv(X_w' X_w)
}

real colvector _wfe_hac_white_compact_labels(
    real colvector idx,
    real scalar    n_levels,
    string scalar  caller,
    string scalar  argname
)
{
    real colvector levels, compact_idx, mask
    real scalar k

    if (any(idx :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, argname)
        _error(3498)
    }
    if (any(idx :!= floor(idx)) | any(idx :< 1)) {
        errprintf("%s: %s must contain integer indices\n", caller, argname)
        _error(3200)
    }
    levels = uniqrows(sort(idx, 1))
    if (rows(levels) != n_levels) {
        if (argname == "unit_idx") {
            if (caller == "wfe_se_fe_twoway") {
                errprintf("%s: unit_idx must enumerate 1..J_u without gaps\n", caller)
            }
            else {
                errprintf("%s: unit_idx must enumerate 1..N_units without gaps\n", caller)
            }
        }
        else {
            errprintf("%s: time_idx must enumerate 1..N_times without gaps\n", caller)
        }
        _error(3200)
    }

    compact_idx = J(rows(idx), 1, .)
    for (k = 1; k <= rows(levels); k++) {
        mask = selectindex(idx :== levels[k])
        compact_idx[mask] = J(rows(mask), 1, k)
    }

    return(compact_idx)
}

void _wfe_fe_twoway_psd_mat(
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
        errprintf("wfe_se_fe_twoway: %s must be positive semidefinite\n",
                  argname)
        _error(3200)
    }
}

void _wfe_hacw_vunit(
    real colvector unit_idx,
    real scalar    N_units,
    string scalar  caller,
    string scalar  range_arg,
    string scalar  gap_msg
)
{
    real colvector seen_units
    real scalar k

    if (rows(unit_idx) == 0) {
        return
    }
    if (any(unit_idx :>= .)) {
        errprintf("%s: unit_idx must not contain missing values\n", caller)
        _error(3498)
    }
    if (any(unit_idx :!= floor(unit_idx)) | min(unit_idx) < 1 | max(unit_idx) > N_units) {
        errprintf("%s: unit_idx must be integer indices within 1..%s\n",
            caller, range_arg)
        _error(3200)
    }

    seen_units = J(N_units, 1, 0)
    for (k = 1; k <= rows(unit_idx); k++) {
        seen_units[unit_idx[k]] = 1
    }
    if (any(seen_units :== 0)) {
        errprintf("%s: %s\n", caller, gap_msg)
        _error(3200)
    }
}

void _wfe_hac_white_validate_panel(
    real colvector unit_idx,
    real colvector time_idx,
    real scalar    N_units,
    real scalar    N_times,
    string scalar  caller
)
{
    real colvector seen_units, seen_times
    real matrix cell_count
    real scalar k, ui, ti

    if (rows(unit_idx) != rows(time_idx)) {
        errprintf("%s: unit_idx/time_idx length mismatch\n", caller)
        _error(3200)
    }
    if (missing(N_units) | missing(N_times) | ///
        N_units != floor(N_units) | N_times != floor(N_times) | ///
        N_units <= 0 | N_times <= 0) {
        errprintf("%s: N_units and N_times must be positive integers\n", caller)
        _error(3200)
    }
    if (rows(unit_idx) == 0) {
        return
    }
    if (any(unit_idx :>= .) | any(time_idx :>= .)) {
        errprintf("%s: unit_idx/time_idx must not contain missing values\n", caller)
        _error(3498)
    }
    if (any(unit_idx :!= floor(unit_idx)) | any(time_idx :!= floor(time_idx)) | ///
        min(unit_idx) < 1 | max(unit_idx) > N_units | ///
        min(time_idx) < 1 | max(time_idx) > N_times) {
        errprintf("%s: unit_idx/time_idx must be integer indices within 1..N_units and 1..N_times\n",
            caller)
        _error(3200)
    }

    seen_units = J(N_units, 1, 0)
    seen_times = J(N_times, 1, 0)
    cell_count = J(N_units, N_times, 0)

    for (k = 1; k <= rows(unit_idx); k++) {
        ui = unit_idx[k]
        ti = time_idx[k]
        seen_units[ui] = 1
        seen_times[ti] = 1
        cell_count[ui, ti] = cell_count[ui, ti] + 1
    }

    if (max(cell_count) > 1) {
        errprintf("%s: unit-time pair is not unique\n", caller)
        _error(498)
    }

    if (any(seen_units :== 0)) {
        errprintf("%s: unit_idx must enumerate 1..N_units without gaps\n", caller)
        _error(3200)
    }
    if (any(seen_times :== 0)) {
        errprintf("%s: time_idx must enumerate 1..N_times without gaps\n", caller)
        _error(3200)
    }
}


struct wfe_twoway_fe_se_result scalar wfe_se_fe_twoway(
    real matrix    X_hat,
    real vector    u_hat,
    real vector    unit_idx,
    real scalar    J_u,
    real matrix    ginv_XX_hat,
    | real vector  time_idx,
      real scalar  N_times
)
{
    struct wfe_twoway_fe_se_result scalar result
    real matrix bread, Omega_raw, info, X_hat_use, X_sorted, ginv_XX_hat_use
    real matrix cell_count
    real matrix gram_hat, ginv_expected
    real colvector u_col, unit_col, unit_sorted, u_sorted, sort_order, time_col
    real colvector seen_units, seen_times
    real scalar n, p, ginv_asym, ginv_scale, df_fe, k, ui, ti
    real scalar ginv_gap, ginv_match_scale, ginv_match_tol

    // Reshape residuals and unit indices to column vectors
    u_col = colshape(u_hat[., .], 1)
    unit_col = colshape(unit_idx[., .], 1)
    X_hat_use = X_hat
    if (rows(X_hat_use) == 1 & cols(X_hat_use) == rows(u_col) & rows(u_col) > 1) {
        X_hat_use = X_hat_use'
    }
    p = cols(X_hat_use)

    if (args() != 5 & args() != 7) {
        errprintf("wfe_se_fe_twoway: supply time_idx and N_times together\n")
        _error(3200)
    }
    if (rows(X_hat_use) != rows(u_col) | rows(X_hat_use) != rows(unit_col)) {
        errprintf("wfe_se_fe_twoway: input length mismatch\n")
        _error(3200)
    }
    if (p < 1) {
        errprintf("wfe_se_fe_twoway: X_hat must contain at least one regressor\n")
        _error(3200)
    }
    if (cols(ginv_XX_hat) != p | rows(ginv_XX_hat) != p) {
        errprintf("wfe_se_fe_twoway: ginv_XX_hat dimension mismatch\n")
        _error(3200)
    }
    if (any(X_hat_use :>= .) | any(u_col :>= .) | any(ginv_XX_hat :>= .)) {
        errprintf("wfe_se_fe_twoway: inputs must not contain missing values\n")
        _error(3498)
    }
    if (rank(X_hat_use) < p) {
        errprintf("wfe_se_fe_twoway: X_hat must have full column rank\n")
        _error(498)
    }
    // Validate ginv_XX_hat is symmetric
    ginv_asym = max(abs(vec(ginv_XX_hat - ginv_XX_hat')))
    ginv_scale = max((1, max(abs(vec(ginv_XX_hat)))))
    if (ginv_asym > sqrt(epsilon(1)) * ginv_scale) {
        errprintf("wfe_se_fe_twoway: ginv_XX_hat must be symmetric\n")
        _error(3200)
    }
    // Symmetrize ginv_XX_hat to handle numerical roundoff
    ginv_XX_hat_use = 0.5 :* (ginv_XX_hat + ginv_XX_hat')
    _wfe_fe_twoway_psd_mat(ginv_XX_hat_use, "ginv_XX_hat")
    gram_hat = cross(X_hat_use, X_hat_use)
    ginv_expected = pinv(gram_hat)
    ginv_gap = max(abs(vec(ginv_XX_hat_use - ginv_expected)))
    ginv_match_scale = max((1, max(abs(vec(ginv_expected))),
        max(abs(vec(ginv_XX_hat_use)))))
    ginv_match_tol = sqrt(epsilon(1)) * rows(X_hat_use) * ginv_match_scale
    if (ginv_gap > ginv_match_tol) {
        errprintf("wfe_se_fe_twoway: ginv_XX_hat must match pinv(X_hat'X_hat)\n")
        _error(498)
    }
    if (J_u != floor(J_u) | J_u <= 0) {
        errprintf("wfe_se_fe_twoway: J_u must be a positive integer\n")
        _error(3200)
    }
    if (J_u < 2) {
        errprintf("wfe_se_fe_twoway: need at least 2 units for FE-side HAC standard errors\n")
        _error(3200)
    }
    n = rows(unit_col)
    if (n == 0) {
        errprintf("wfe_se_fe_twoway: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(unit_col :>= .)) {
        errprintf("wfe_se_fe_twoway: unit_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | min(unit_col) < 1) {
        errprintf("wfe_se_fe_twoway: unit_idx must contain positive integer labels\n")
        _error(3200)
    }
    if (rows(uniqrows(sort(unit_col, 1))) != J_u) {
        errprintf("wfe_se_fe_twoway: unit_idx must enumerate 1..J_u without gaps\n")
        _error(3200)
    }
    unit_col = _wfe_hac_white_compact_labels(unit_col, J_u,
        "wfe_se_fe_twoway", "unit_idx")
    df_fe = 1
    if (args() >= 7) {
        time_col = colshape(time_idx[., .], 1)
        if (missing(N_times) | N_times != floor(N_times) | N_times < 2) {
            errprintf("wfe_se_fe_twoway: N_times must be an integer >= 2 when time_idx is provided\n")
            _error(3200)
        }
        if (rows(time_col) != n) {
            errprintf("wfe_se_fe_twoway: unit_idx/time_idx length mismatch\n")
            _error(3200)
        }
        if (any(time_col :>= .)) {
            errprintf("wfe_se_fe_twoway: unit_idx/time_idx must not contain missing values\n")
            _error(3498)
        }
        if (any(time_col :!= floor(time_col)) | min(time_col) < 1) {
            errprintf("wfe_se_fe_twoway: time_idx must contain positive integer labels\n")
            _error(3200)
        }
        if (rows(uniqrows(sort(time_col, 1))) != N_times) {
            errprintf("wfe_se_fe_twoway: time_idx must enumerate 1..N_times without gaps\n")
            _error(3200)
        }
        time_col = _wfe_hac_white_compact_labels(time_col, N_times,
            "wfe_se_fe_twoway", "time_idx")
        cell_count = J(J_u, N_times, 0)
        for (k = 1; k <= n; k++) {
            ui = unit_col[k]
            ti = time_col[k]
            cell_count[ui, ti] = cell_count[ui, ti] + 1
        }
        if (max(cell_count) > 1) {
            errprintf("wfe_se_fe_twoway: unit-time pair is not unique\n")
            _error(498)
        }
        if (n - J_u - N_times - p <= 0) {
            errprintf("wfe_se_fe_twoway: FE-side HAC requires NT - J_u - N_times - p > 0\n")
            _error(3352)
        }
        df_fe = (n / (n - 1)) * ((n - p) / (n - J_u - N_times - p))
    }

    // Sort by unit for panelsetup
    if (n <= 1) {
        sort_order = 1::n
    }
    else {
        sort_order = order((unit_col, (1::n)), (1, 2))
    }
    unit_sorted = unit_col[sort_order]
    X_sorted = X_hat_use[sort_order, .]
    u_sorted = u_col[sort_order]

    info = panelsetup(unit_sorted, 1)
    if (rows(info) != J_u) {
        errprintf("wfe_se_fe_twoway: unit_idx must enumerate 1..J_u without gaps\n")
        _error(3200)
    }

    // Compute FE-side HAC components
    Omega_raw = _wfe_omega_hac(X_sorted, u_sorted, unit_sorted, J_u)
    result.Omega_fe = Omega_raw / J_u

    bread = J_u * ginv_XX_hat_use
    result.Psi_hat_fe = df_fe * (bread * result.Omega_fe * bread)
    result.var_cov_fe = result.Psi_hat_fe / J_u

    return(result)
}

// ============================================================
// wfe_se_hac_twoway_for_white() — Real-domain FWL+HAC Psi_wfe for two-way White test
//
// Computes a HAC-based Psi_wfe using FWL-demeaned and weight-scaled variables,
// ensuring scale consistency with wfe_se_fe_twoway's Psi_hat_fe.
//
// Algorithm (plan section 9.1):
//   a) FWL two-way demean Y and X → Y_dm, X_dm
//   b) WFE residuals: e_dm = Y_dm - X_dm * beta_wfe
//   c) Weighted variables: X_w = X_dm .* sqrt(|W|),  u_w = e_dm .* sqrt(|W|) .* sign(W)
//   d) bread = N_units * pinv(X_w' X_w)
//   e) Omega_raw = _wfe_omega_hac(X_w, u_w, unit_idx, N_units)
//   f) Psi = bread * (Omega_raw / N_units) * bread
//
// This is used ONLY for the White misspecification test's Phi matrix.
// The SE path continues to use wfe_se_gmm (GMM sandwich).
//
// @param beta_wfe  real vector   [p x 1]   WFE coefficient vector
// @param Y         real vector   [NT x 1]  Outcome
// @param X         real matrix   [NT x p]  Covariates (same order as Y)
// @param W_vec     real vector   [NT x 1]  DiD/Mdid weights
// @param unit_idx  real vector   [NT x 1]  Unit index (1..N_units)
// @param time_idx  real vector   [NT x 1]  Time index (1..N_times)
// @param N_units   real scalar             Number of units
// @param N_times   real scalar             Number of time periods
// @return          real matrix   [p x p]   HAC-based Psi_wfe (same scale as Psi_hat_fe)
// ============================================================
struct wfe_hac_white_result scalar wfe_se_hac_twoway_for_white(
    real vector    beta_wfe,
    real vector    Y,
    real matrix    X,
    real vector    W_vec,
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times
)
{
    struct wfe_hac_white_result scalar res
    real scalar NT, p, k, ui, ti
    real matrix demeaned, X_dm, Omega_raw, bread, X_use, cell_count
    real colvector Y_col, beta_col, W_col, unit_col, time_col
    real colvector unit_support_count, time_support_count
    real colvector seen_units, seen_times
    real colvector Y_dm, u_w

    Y_col    = colshape(Y, 1)
    beta_col = colshape(beta_wfe, 1)
    W_col    = colshape(W_vec, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    NT = rows(Y_col)
    X_use = X
    if (rows(X_use) == 1 & cols(X_use) == NT & NT > 1) X_use = X_use'
    p = cols(X_use)
    if (p < 1) {
        errprintf("wfe_se_hac_twoway_for_white: X must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(beta_col) != p) {
        errprintf("wfe_se_hac_twoway_for_white: beta_wfe length must match cols(X)\n")
        _error(3200)
    }
    if (any(beta_col :>= .)) {
        errprintf("wfe_se_hac_twoway_for_white: beta_wfe must not contain missing values\n")
        _error(3498)
    }
    if (rows(W_col) != NT) {
        errprintf("wfe_se_hac_twoway_for_white: W_vec length must match Y/X/unit_idx/time_idx\n")
        _error(3200)
    }
    if (rows(unit_col) != NT | rows(time_col) != NT) {
        errprintf("wfe_se_hac_twoway_for_white: unit_idx/time_idx length must match Y and X\n")
        _error(3200)
    }
    if (any(Y_col :>= .) | any(X_use :>= .)) {
        errprintf("wfe_se_hac_twoway_for_white: Y and X must not contain missing values\n")
        _error(3498)
    }
    if (any(W_col :>= .)) {
        errprintf("wfe_se_hac_twoway_for_white: W_vec must not contain missing values\n")
        _error(3498)
    }
    if (missing(N_units) | missing(N_times) | ///
        N_units != floor(N_units) | N_times != floor(N_times) | ///
        N_units <= 0 | N_times <= 0) {
        errprintf("wfe_se_hac_twoway_for_white: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("wfe_se_hac_twoway_for_white: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col)) | ///
        min(unit_col) < 1 | min(time_col) < 1) {
        errprintf("wfe_se_hac_twoway_for_white: unit_idx/time_idx must contain positive integer labels\n")
        _error(3200)
    }
    if (rows(uniqrows(sort(unit_col, 1))) != N_units) {
        errprintf("wfe_se_hac_twoway_for_white: unit_idx must enumerate 1..N_units without gaps\n")
        _error(3200)
    }
    if (rows(uniqrows(sort(time_col, 1))) != N_times) {
        errprintf("wfe_se_hac_twoway_for_white: time_idx must enumerate 1..N_times without gaps\n")
        _error(3200)
    }
    unit_col = _wfe_hac_white_compact_labels(unit_col, N_units,
        "wfe_se_hac_twoway_for_white", "unit_idx")
    time_col = _wfe_hac_white_compact_labels(time_col, N_times,
        "wfe_se_hac_twoway_for_white", "time_idx")
    cell_count = J(N_units, N_times, 0)
    for (k = 1; k <= NT; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ui, ti] = cell_count[ui, ti] + 1
    }
    if (max(cell_count) > 1) {
        errprintf("wfe_se_hac_twoway_for_white: unit-time pair is not unique\n")
        _error(498)
    }
    if (sum(W_col :!= 0) == 0) {
        errprintf("wfe_se_hac_twoway_for_white: at least one non-zero weight is required\n")
        _error(498)
    }

    // FWL two-way demeaning (handles unsorted input internally)
    demeaned = _wfe_twoway_demean(Y_col, X_use, unit_col, time_col, N_units, N_times)
    Y_dm = demeaned[, 1]
    X_dm = demeaned[|1, 2 \ NT, p + 1|]

    unit_support_count = J(N_units, 1, 0)
    time_support_count = J(N_times, 1, 0)
    for (k = 1; k <= NT; k++) {
        if (W_col[k] == 0) continue
        unit_support_count[unit_col[k]] = 1
        time_support_count[time_col[k]] = 1
    }
    if (sum(unit_support_count) < 2) {
        errprintf("wfe_se_hac_twoway_for_white: non-zero weights must span at least 2 units\n")
        _error(498)
    }
    if (sum(time_support_count) < 2) {
        errprintf("wfe_se_hac_twoway_for_white: non-zero weights must span at least 2 time periods\n")
        _error(498)
    }

    // WFE residuals under beta_wfe — stored in struct for White test Lambda
    res.e_dm = Y_dm - X_dm * beta_col

    // Weighted FWL X — stored in struct for White test Lambda
    // X_w[k,j] = X_dm[k,j]*sqrt(|W_k|),  u_w[k] = e_dm[k]*sqrt(|W_k|)*sign(W_k)
    // Zero-weight obs contribute zero (sqrt(0)=0), naturally excluded from clustering
    res.X_w = X_dm :* sqrt(abs(W_col))
    if (rank(res.X_w) < p) {
        errprintf("wfe_se_hac_twoway_for_white: X_w must have full column rank\n")
        _error(498)
    }
    u_w = res.e_dm :* sqrt(abs(W_col)) :* sign(W_col)

    // ginv_XX_w = pinv(X_w' X_w) — stored in struct for White test bread
    res.ginv_XX_w = pinv(res.X_w' * res.X_w)

    // HAC Omega via existing helper (_wfe_omega_hac sorts by unit internally)
    Omega_raw = _wfe_omega_hac(res.X_w, u_w, unit_col, N_units)

    // Bread: N_units * ginv_XX_w — same structure as wfe_se_fe_twoway
    bread = N_units * res.ginv_XX_w

    // Psi_wfe = bread * (Omega_raw / N_units) * bread
    res.Psi_wfe = bread * (Omega_raw / N_units) * bread

    return(res)
}

end
