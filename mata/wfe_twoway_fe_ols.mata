// ============================================================
// wfe_twoway_fe_ols.mata -- Two-way fixed effects OLS estimation
//
// Implements FWL projection for two-way fixed effects.
// Demeans Y and X with respect to unit and time effects,
// then applies OLS on the projected data.
// ============================================================

version 16.0
mata:
mata set matastrict on

struct wfe_twoway_fe_ols_result {
    real rowvector beta_fe
    real colvector u_hat
    real matrix    X_hat
    real matrix    ginv_XX_hat
}

real colvector _wfe_twoway_compact_labels(
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
        errprintf("%s: %s must contain positive integer labels\n", caller, argname)
        _error(3200)
    }
    levels = uniqrows(sort(idx, 1))
    if (rows(levels) != n_levels) {
        errprintf("%s: %s must contain exactly %g distinct labels\n",
            caller, argname, n_levels)
        _error(3200)
    }

    compact_idx = J(rows(idx), 1, .)
    for (k = 1; k <= rows(levels); k++) {
        mask = selectindex(idx :== levels[k])
        compact_idx[mask] = J(rows(mask), 1, k)
    }

    return(compact_idx)
}


struct wfe_twoway_fe_ols_result scalar wfe_twoway_fe_ols(
    real vector    Y,
    real matrix    X,
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times
)
{
    struct wfe_twoway_fe_ols_result scalar result
    struct wfe_ols_result scalar ols_result
    real scalar n, p, k, ui, ti, j, scale_j, tol_j, projected_ss
    real colvector Y_col, unit_col, time_col
    real colvector unit_count, time_count, counts
    real matrix X_mat, U, cell_count, Tdummy, P1, Q, YX, P1_YX, YX_dm

    Y_col = colshape(Y[., .], 1)
    X_mat = X[., .]
    unit_col = colshape(unit_idx[., .], 1)
    time_col = colshape(time_idx[., .], 1)

    n = rows(Y_col)
    if (rows(X_mat) == 1 & cols(X_mat) == n & n > 1) {
        X_mat = X_mat'
    }
    p = cols(X_mat)

    if (p < 1) {
        errprintf("wfe_twoway_fe_ols: X must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(X_mat) != n | rows(unit_col) != n | rows(time_col) != n) {
        errprintf("wfe_twoway_fe_ols: input length mismatch\n")
        _error(3200)
    }
    if (n == 0) {
        errprintf("wfe_twoway_fe_ols: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(Y_col :>= .) | any(X_mat :>= .)) {
        errprintf("wfe_twoway_fe_ols: Y and X must not contain missing values\n")
        _error(3498)
    }
    if (N_units != floor(N_units) | N_times != floor(N_times) | ///
        N_units <= 0 | N_times <= 0) {
        errprintf("wfe_twoway_fe_ols: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (N_times < 2) {
        errprintf("wfe_twoway_fe_ols: two-way FE requires at least 2 time periods\n")
        _error(3200)
    }
    if (N_units < 2) {
        errprintf("wfe_twoway_fe_ols: need at least 2 units for two-way FE\n")
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("wfe_twoway_fe_ols: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col))) {
        errprintf("wfe_twoway_fe_ols: unit_idx/time_idx must contain integer indices\n")
        _error(3200)
    }

    unit_col = _wfe_twoway_compact_labels(unit_col, N_units,
        "wfe_twoway_fe_ols", "unit_idx")
    time_col = _wfe_twoway_compact_labels(time_col, N_times,
        "wfe_twoway_fe_ols", "time_idx")

    U = J(n, N_units, 0)
    Tdummy = J(n, N_times, 0)
    cell_count = J(N_units, N_times, 0)
    unit_count = J(N_units, 1, 0)
    time_count = J(N_times, 1, 0)

    for (k = 1; k <= n; k++) {
        ui = unit_col[k]
        ti = time_col[k]

        if (ui < 1 | ui > N_units | ti < 1 | ti > N_times) {
            errprintf("wfe_twoway_fe_ols: index out of range at row %g\n", k)
            _error(3200)
        }

        cell_count[ui, ti] = cell_count[ui, ti] + 1
        unit_count[ui] = unit_count[ui] + 1
        time_count[ti] = time_count[ti] + 1
        U[k, ui] = 1
        Tdummy[k, ti] = 1
    }

    if (max(cell_count) > 1) {
        errprintf("wfe_twoway_fe_ols: unit-time pair is not unique\n")
        _error(498)
    }
    if (any(unit_count :== 0)) {
        errprintf("wfe_twoway_fe_ols: unit_idx must enumerate 1..N_units without gaps\n")
        _error(3200)
    }
    if (any(time_count :== 0)) {
        errprintf("wfe_twoway_fe_ols: time_idx must enumerate 1..N_times without gaps\n")
        _error(3200)
    }

    YX = Y_col, X_mat
    counts = colsum(U)'
    P1 = U * diag(1 :/ counts) * U'
    P1_YX = P1 * YX
    Q = Tdummy - P1 * Tdummy

    YX_dm = YX - P1_YX - Q * pinv(cross(Q, Q)) * cross(Q, YX)

    result.X_hat = YX_dm[., 2..cols(YX_dm)]

    // Check for regressors fully absorbed by FE
    for (j = 1; j <= p; j++) {
        scale_j = max(abs(X_mat[., j]))
        tol_j = sqrt(epsilon(1)) * n * (scale_j^2)
        projected_ss = cross(result.X_hat[., j], result.X_hat[., j])
        if (projected_ss <= tol_j) {
            errprintf("wfe_twoway_fe_ols: regressor fully absorbed by two-way demeaning\n")
            _error(498)
        }
    }
    if (rank(result.X_hat) < p) {
        errprintf("wfe_twoway_fe_ols: regressors are collinear after two-way demeaning\n")
        _error(498)
    }

    ols_result = _wfe_ols_core(YX_dm[., 1], result.X_hat)

    result.beta_fe = ols_result.beta
    result.u_hat = ols_result.resid
    result.ginv_XX_hat = ols_result.ginv_XX

    return(result)
}

end
