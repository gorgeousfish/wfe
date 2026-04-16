// wfe_se_gmm.mata -- GMM cluster-robust standard errors for two-way FE

version 16.0
mata:
mata set matastrict on

real colvector _wfe_gmm_compact_index(
    real colvector idx,
    real scalar    max_idx
)
{
    real colvector compact_idx, seen
    real scalar k, level, next_level

    compact_idx = J(rows(idx), 1, .)
    seen = J(max_idx, 1, 0)
    next_level = 0

    for (k = 1; k <= rows(idx); k++) {
        level = idx[k]
        if (seen[level] == 0) {
            next_level = next_level + 1
            seen[level] = next_level
        }
        compact_idx[k] = seen[level]
    }

    return(compact_idx)
}

real matrix _wfe_twoway_demean(
    real vector    Y,
    real matrix    X,
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times
)
{
    real scalar NT, p1, k, i, t
    real matrix Z, demeaned, cell_count, U, Tdummy, counts, P1, Q, X_use
    real colvector Y_col, unit_col, time_col

    Y_col = colshape(Y, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    NT = rows(Y_col)
    X_use = X
    // Handle rowvector input: convert to column vector for consistent processing
    if (rows(X_use) == 1 & cols(X_use) == NT & NT > 1) {
        X_use = X_use'
    }
    if (rows(X_use) != NT | rows(unit_col) != NT | rows(time_col) != NT) {
        errprintf("_wfe_twoway_demean: input length mismatch\n")
        _error(3200)
    }
    if (NT == 0) {
        errprintf("_wfe_twoway_demean: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (N_units != floor(N_units) | N_times != floor(N_times) | N_units <= 0 | N_times <= 0) {
        errprintf("_wfe_twoway_demean: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (N_units < 2) {
        errprintf("_wfe_twoway_demean: two-way demean requires at least 2 units\n")
        _error(3200)
    }
    if (N_times < 2) {
        errprintf("_wfe_twoway_demean: two-way demean requires at least 2 time periods\n")
        _error(3200)
    }
    if (any(Y_col :>= .) | any(X_use :>= .)) {
        errprintf("_wfe_twoway_demean: Y and X must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("_wfe_twoway_demean: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col))) {
        errprintf("_wfe_twoway_demean: unit_idx/time_idx must contain integer indices\n")
        _error(3200)
    }
    if (min(unit_col) < 1 | max(unit_col) > N_units | ///
        min(time_col) < 1 | max(time_col) > N_times) {
        errprintf("_wfe_twoway_demean: unit_idx/time_idx out of range\n")
        _error(3200)
    }
    Z = Y_col, X_use
    p1 = cols(Z)
    cell_count = J(N_units, N_times, 0)
    U = J(NT, N_units, 0)
    Tdummy = J(NT, N_times, 0)

    for (k = 1; k <= NT; k++) {
        i = unit_col[k]
        t = time_col[k]
        cell_count[i, t] = cell_count[i, t] + 1
        if (cell_count[i, t] > 1) {
            errprintf("_wfe_twoway_demean: unit-time pair is not unique\n")
            _error(498)
        }
        U[k, i] = 1
        Tdummy[k, t] = 1
    }

    if (any(colsum(cell_count')' :== 0) | any(colsum(cell_count)' :== 0)) {
        errprintf("_wfe_twoway_demean: unit_idx/time_idx must cover 1..N_units and 1..N_times\n")
        _error(3200)
    }

    // Two-way FE demeaning via FWL projection onto unit and time dummies.
    // The simple mean identity only holds on balanced panels.
    counts = colsum(U)'
    P1 = U * diag(1 :/ counts) * U'
    Q = Tdummy - P1 * Tdummy
    demeaned = Z - P1 * Z - Q * pinv(cross(Q, Q)) * cross(Q, Z)

    return(demeaned)
}


real matrix wfe_se_gmm(
    real vector    beta,
    real vector    Y,
    real matrix    X,
    real vector    W_it,
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times,
    real scalar    nK,
    real scalar    df_adjustment,
    real scalar    Nstar,
    real matrix    Psi_out
)
{
    real scalar NT, p, k, j, g, start, stop, df_correction_denom, ui, ti, actual_nstar
    real scalar support_units, support_times, df_correction, N_units_eff, N_times_eff
    real matrix demeaned, X_dm, U, V, X_g, X_weighted, U_scaled, V_scaled, inv_U, Psi
    real colvector beta_col, Y_col, W_col, unit_col, time_col, keep_index
    real colvector Y_dm, Y_g, W_g, e_g, a_g, sort_order
    real colvector unit_support_count, time_support_count
    real colvector Y_sorted, W_sorted, unit_sorted, time_sorted
    real matrix X_sorted, X_use
    real matrix panel, cell_count

    beta_col = colshape(beta, 1)
    Y_col = colshape(Y, 1)
    W_col = colshape(W_it, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)

    NT = rows(Y_col)
    X_use = X
    // Normalize rowvector to column vector for consistent processing
    if (rows(X_use) == 1 & cols(X_use) == NT & NT > 1) {
        X_use = X_use'
    }
    p = cols(X_use)

    if (p < 1) {
        errprintf("wfe_se_gmm: X must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(beta_col) != p) {
        errprintf("wfe_se_gmm: beta length mismatch\n")
        _error(3200)
    }
    if (rows(X_use) != NT | rows(W_col) != NT) {
        errprintf("wfe_se_gmm: dimension mismatch\n")
        _error(3200)
    }
    if (rows(unit_col) != NT | rows(time_col) != NT) {
        errprintf("wfe_se_gmm: index length mismatch\n")
        _error(3200)
    }
    if (NT == 0) {
        errprintf("wfe_se_gmm: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(beta_col :>= .)) {
        errprintf("wfe_se_gmm: beta must not contain missing values\n")
        _error(3498)
    }
    if (any(Y_col :>= .) | any(X_use :>= .)) {
        errprintf("wfe_se_gmm: Y and X must not contain missing values\n")
        _error(3498)
    }
    if (any(W_col :>= .)) {
        errprintf("wfe_se_gmm: W_it must not contain missing values\n")
        _error(3498)
    }
    actual_nstar = sum(W_col :!= 0)
    if (actual_nstar == 0) {
        errprintf("wfe_se_gmm: at least one non-zero weight is required\n")
        _error(498)
    }
    if (Nstar != floor(Nstar) | Nstar <= 0) {
        errprintf("wfe_se_gmm: Nstar must be a positive integer\n")
        _error(3200)
    }
    if (Nstar != actual_nstar) {
        errprintf("wfe_se_gmm: Nstar must equal sum(W_it != 0)\n")
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("wfe_se_gmm: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (nK != p) {
        errprintf("wfe_se_gmm: nK mismatch\n")
        _error(3200)
    }
    if (N_units != floor(N_units) | N_times != floor(N_times) | N_units <= 0 | N_times <= 0) {
        errprintf("wfe_se_gmm: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (N_times < 2) {
        errprintf("wfe_se_gmm: two-way GMM requires at least 2 time periods\n")
        _error(3200)
    }
    if (N_units < 2) {
        errprintf("wfe_se_gmm: need at least 2 units for cluster standard errors\n")
        _error(3200)
    }
    if (df_adjustment != 0 & df_adjustment != 1) {
        errprintf("wfe_se_gmm: df_adjustment must be 0 or 1\n")
        _error(3200)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col))) {
        errprintf("wfe_se_gmm: unit_idx/time_idx must be integer indices\n")
        _error(3200)
    }
    if (min(unit_col) < 1 | max(unit_col) > N_units) {
        errprintf("wfe_se_gmm: unit_idx out of range\n")
        _error(3200)
    }
    if (min(time_col) < 1 | max(time_col) > N_times) {
        errprintf("wfe_se_gmm: time_idx out of range\n")
        _error(3200)
    }
    unit_support_count = J(N_units, 1, 0)
    time_support_count = J(N_times, 1, 0)
    for (k = 1; k <= NT; k++) {
        if (W_col[k] == 0) continue
        unit_support_count[unit_col[k]] = 1
        time_support_count[time_col[k]] = 1
    }
    // Weighted two-way score requires support spanning both dimensions
    support_units = sum(unit_support_count)
    support_times = sum(time_support_count)
    if (support_units < 2) {
        errprintf("wfe_se_gmm: non-zero weights must span at least 2 units\n")
        _error(498)
    }
    if (support_times < 2) {
        errprintf("wfe_se_gmm: non-zero weights must span at least 2 time periods\n")
        _error(498)
    }

    // Zero-weight observations are outside the weighted two-way score and must
    // not influence the FWL projection or cluster count used by the GMM meat.
    // Re-index the retained support to a dense 1..J_u / 1..J_t panel so the
    // internal projection and cluster aggregation match the effective support.
    keep_index = selectindex(W_col :!= 0)
    if (rows(keep_index) != NT) {
        Y_col = Y_col[keep_index]
        X_use = X_use[keep_index, .]
        W_col = W_col[keep_index]
        unit_col = unit_col[keep_index]
        time_col = time_col[keep_index]
        NT = rows(Y_col)
    }
    unit_col = _wfe_gmm_compact_index(unit_col, N_units)
    time_col = _wfe_gmm_compact_index(time_col, N_times)
    N_units_eff = support_units
    N_times_eff = support_times

    if (NT <= 1) {
        sort_order = 1::NT
    }
    else {
        // Sort by unit and time for cluster aggregation
        sort_order = order((unit_col, time_col, (1::NT)), (1, 2, 3))
    }
    Y_sorted = Y_col[sort_order]
    X_sorted = X_use[sort_order, .]
    W_sorted = W_col[sort_order]
    unit_sorted = unit_col[sort_order]
    time_sorted = time_col[sort_order]

    cell_count = J(N_units_eff, N_times_eff, 0)
    for (k = 1; k <= NT; k++) {
        ui = unit_sorted[k]
        ti = time_sorted[k]
        cell_count[ui, ti] = cell_count[ui, ti] + 1
    }
    if (max(cell_count) > 1) {
        errprintf("wfe_se_gmm: unit-time pair is not unique\n")
        _error(498)
    }
    // Verify time indices are consecutive
    if (any(colsum(cell_count)' :== 0)) {
        errprintf("wfe_se_gmm: time_idx must enumerate 1..N_times without gaps\n")
        _error(3200)
    }

    panel = panelsetup(unit_sorted, 1)
    if (rows(panel) != N_units_eff) {
        errprintf("wfe_se_gmm: unit_idx must enumerate 1..N_units without gaps\n")
        _error(3200)
    }

    demeaned = _wfe_twoway_demean(Y_sorted, X_sorted, unit_sorted, time_sorted,
        N_units_eff, N_times_eff)
    Y_dm = demeaned[, 1]
    X_dm = demeaned[|1, 2 \ NT, p + 1|]

    U = J(p, p, 0)
    V = J(p, p, 0)

    for (g = 1; g <= rows(panel); g++) {
        start = panel[g, 1]
        stop = panel[g, 2]
        X_g = X_dm[|start, 1 \ stop, p|]
        Y_g = Y_dm[|start \ stop|]
        W_g = W_sorted[|start \ stop|]
        e_g = Y_g - X_g * beta_col

        X_weighted = X_g
        for (j = 1; j <= p; j++) {
            X_weighted[, j] = X_g[, j] :* W_g
        }

        U = U + cross(X_g, X_weighted)
        a_g = cross(X_g, W_g :* e_g)
        V = V + a_g * a_g'
    }

    U_scaled = U / N_units_eff
    if (rank(U_scaled) < p) {
        errprintf("wfe_se_gmm: weighted GMM moment matrix is singular; two-way covariance is undefined\n")
        _error(498)
    }

    inv_U = luinv(U_scaled)
    V_scaled = V / N_units_eff
    Psi = inv_U * V_scaled * inv_U

    if (df_adjustment == 1) {
        df_correction_denom = Nstar - support_units - support_times - nK
        if (df_correction_denom <= 0) {
            errprintf("wfe_se_gmm: two-way GMM df_adjustment(on) requires Nstar - N_nonzero_units - N_nonzero_times - p > 0; try df_adjustment(off)\n")
            _error(3351)
        }
        df_correction = (Nstar / (Nstar - 1)) * ((Nstar - nK) / df_correction_denom)
        Psi = df_correction * Psi
    }

    Psi_out = Psi
    return(Psi / N_units_eff)
}

end
