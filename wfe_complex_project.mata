// ============================================================
// wfe_complex_project.mata -- Complex DiD projection core
// ============================================================

version 16.0
mata:
mata set matastrict on

struct wfe_twoway_result {
    real colvector beta
    real matrix    vcov
    real matrix    W
    real colvector beta_fe
    real matrix    vcov_fe
    real scalar    sigma
    real scalar    df_r
    real scalar    N_nonzero
    real scalar    white_stat
    real scalar    white_pval
    string scalar  vcetype
    string scalar  white_test
}

struct wfe_complex_projection_state {
    complex colvector y_tilde
    complex matrix    X_tilde
    complex colvector w_sqrt
    real scalar       nz_obs
    real scalar       tn_row
    real scalar       n_Udummy
    real scalar       n_Tdummy
    real matrix       Udummy
    real matrix       Tdummy
}

struct wfe_complex_ols_result {
    real    colvector coef_wls
    complex colvector betaT
    complex colvector e_tilde
    complex colvector resid_vec
    real    scalar    sigma2
    real    scalar    d_f
    complex matrix    ginv_XX_tilde
    complex colvector diag_ee_tilde
}

complex colvector _wfe_complex_im_sqrt(real vector weights)
{
    real scalar n, k
    real colvector weight_col
    complex colvector out

    weight_col = colshape(weights, 1)
    if (any(weight_col :>= .)) {
        errprintf("_wfe_complex_im_sqrt: weights must not contain missing values\n")
        _error(3498)
    }
    n = rows(weight_col)
    out = J(n, 1, C(0, 0))

    for (k = 1; k <= n; k++) {
        if (weight_col[k] > 0) {
            out[k] = C(sqrt(weight_col[k]), 0)
        }
        else if (weight_col[k] < 0) {
            out[k] = C(0, sqrt(-weight_col[k]))
        }
    }

    return(out)
}

real matrix _wfe_complex_unit_dummy(real vector unit_idx, real scalar N_units)
{
    real scalar n, k
    real colvector unit_col
    real matrix U

    if (missing(N_units) | N_units != floor(N_units) | N_units < 1) {
        errprintf("_wfe_complex_unit_dummy: N_units must be a positive integer\n")
        _error(3200)
    }
    unit_col = colshape(unit_idx, 1)
    if (any(unit_col :>= .)) {
        errprintf("_wfe_complex_unit_dummy: unit_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(unit_col :< 1) | any(unit_col :> N_units)) {
        errprintf("_wfe_complex_unit_dummy: unit_idx must be integer indices within 1..N_units\n")
        _error(3200)
    }
    n = rows(unit_col)
    U = J(n, N_units, 0)

    for (k = 1; k <= n; k++) {
        U[k, unit_col[k]] = 1
    }

    return(U)
}

real matrix _wfe_complex_time_dummy(real vector time_idx, real scalar N_times)
{
    real scalar n, k
    real colvector time_col
    real matrix Tdummy

    if (missing(N_times) | N_times != floor(N_times) | N_times < 1) {
        errprintf("_wfe_complex_time_dummy: N_times must be a positive integer\n")
        _error(3200)
    }
    time_col = colshape(time_idx, 1)
    if (any(time_col :>= .)) {
        errprintf("_wfe_complex_time_dummy: time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(time_col :!= floor(time_col)) | any(time_col :< 1) | any(time_col :> N_times)) {
        errprintf("_wfe_complex_time_dummy: time_idx must be integer indices within 1..N_times\n")
        _error(3200)
    }
    n = rows(time_col)
    Tdummy = J(n, N_times, 0)

    for (k = 1; k <= n; k++) {
        Tdummy[k, time_col[k]] = 1
    }

    return(Tdummy)
}

complex matrix _wfe_complex_weight_real_matrix(real matrix X, complex vector w_sqrt)
{
    real scalar n, p, j
    real matrix real_part, imag_part, X_use
    complex colvector w_col

    w_col = colshape(w_sqrt, 1)
    X_use = X
    if (rows(X_use) == 1 & cols(X_use) == rows(w_col) & rows(w_col) > 1) {
        X_use = transposeonly(X_use)
    }

    n = rows(X_use)
    p = cols(X_use)
    if (any(X_use :>= .)) {
        errprintf("_wfe_complex_weight_real_matrix: X must not contain missing values\n")
        _error(3498)
    }
    if (any(Re(w_col) :>= .) | any(Im(w_col) :>= .)) {
        errprintf("_wfe_complex_weight_real_matrix: w_sqrt must not contain missing values\n")
        _error(3498)
    }
    if (rows(w_col) != n) {
        errprintf("_wfe_complex_weight_real_matrix: weight length mismatch\n")
        _error(3200)
    }
    real_part = J(n, p, 0)
    imag_part = J(n, p, 0)

    for (j = 1; j <= p; j++) {
        real_part[., j] = X_use[., j] :* Re(w_col)
        imag_part[., j] = X_use[., j] :* Im(w_col)
    }

    return(C(real_part, imag_part))
}

complex colvector _wfe_complex_general_inv(complex vector weight, real scalar tol)
{
    real scalar k, n
    complex colvector out, weight_col

    if (missing(tol) | tol <= 0) {
        errprintf("_wfe_complex_general_inv: tol must be positive\n")
        _error(3200)
    }
    weight_col = colshape(weight, 1)
    if (any(Re(weight_col) :>= .) | any(Im(weight_col) :>= .)) {
        errprintf("_wfe_complex_general_inv: weight must not contain missing values\n")
        _error(3498)
    }
    n = rows(weight_col)
    out = J(n, 1, C(0, 0))

    for (k = 1; k <= n; k++) {
        if (abs(weight_col[k]) >= tol) {
            out[k] = 1 / weight_col[k]
        }
    }

    return(out)
}

void _wfe_cproj_validate_panel(
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times,
    string scalar  caller
)
{
    real scalar n, k, ui, ti
    real matrix cell_count
    real colvector unit_col, time_col

    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    n = rows(unit_col)
    if (rows(time_col) != n) {
        errprintf("%s: input length mismatch\n", caller)
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("%s: unit_idx/time_idx must not contain missing values\n", caller)
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col)) | ///
        any(unit_col :< 1) | any(unit_col :> N_units) | ///
        any(time_col :< 1) | any(time_col :> N_times)) {
        errprintf("%s: unit_idx/time_idx must be integer indices within 1..N_units and 1..N_times\n", caller)
        _error(3200)
    }

    cell_count = J(N_units, N_times, 0)
    for (k = 1; k <= n; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ui, ti] = cell_count[ui, ti] + 1
    }
    if (max(cell_count) > 1) {
        errprintf("%s: unit-time pair is not unique\n", caller)
        _error(498)
    }
}

struct wfe_complex_projection_state scalar _wfe_complex_project_core(
    real vector    Y,
    real matrix    X,
    real vector    W_vec,
    real vector    unit_idx,
    real vector    time_idx,
    real scalar    N_units,
    real scalar    N_times,
    real scalar    white_flag,
    real scalar    tol
)
{
    struct wfe_complex_projection_state scalar state
    real scalar n, p, k
    real colvector nz_index, keep_mask, col_keep_u, col_keep_t
    real colvector unit_support_count, time_support_count
    real colvector Y_col, W_col, unit_col, time_col
    real colvector Y_nz, W_nz, unit_nz, time_nz
    real matrix X_nz, D, X_use
    complex colvector unit_weight_sum, inv_weight
    complex matrix Y_star, X_star, D_star, D1_star, D2_star
    complex matrix Dginv, P1, Q, QQ_inv, YX_star, transformed

    Y_col = colshape(Y, 1)
    W_col = colshape(W_vec, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    n = rows(Y_col)
    X_use = X
    if (rows(X_use) == 1 & cols(X_use) == n & n > 1) {
        X_use = X_use'
    }
    p = cols(X_use)

    if (rows(X_use) != n | rows(W_col) != n | rows(unit_col) != n | rows(time_col) != n) {
        errprintf("_wfe_complex_project_core: input length mismatch\n")
        _error(3200)
    }
    if (n == 0) {
        errprintf("_wfe_complex_project_core: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (p < 1) {
        errprintf("_wfe_complex_project_core: X must contain at least one regressor\n")
        _error(3200)
    }
    if (N_units != floor(N_units) | N_times != floor(N_times) | N_units <= 0 | N_times <= 0) {
        errprintf("_wfe_complex_project_core: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (N_times < 2) {
        errprintf("_wfe_complex_project_core: two-way projection requires at least 2 time periods\n")
        _error(3200)
    }
    if (N_units < 2) {
        errprintf("_wfe_complex_project_core: two-way projection requires at least 2 units\n")
        _error(3200)
    }
    if (missing(tol) | tol <= 0) {
        errprintf("_wfe_complex_project_core: tol must be positive\n")
        _error(3200)
    }
    if (missing(white_flag) | (white_flag != 0 & white_flag != 1)) {
        errprintf("_wfe_complex_project_core: white_flag must be 0 or 1\n")
        _error(3200)
    }
    if (any(Y_col :>= .)) {
        errprintf("_wfe_complex_project_core: Y must not contain missing values\n")
        _error(3498)
    }
    if (any(X_use :>= .)) {
        errprintf("_wfe_complex_project_core: X must not contain missing values\n")
        _error(3498)
    }
    if (any(W_col :>= .)) {
        errprintf("_wfe_complex_project_core: W_vec must not contain missing values\n")
        _error(3498)
    }
    _wfe_cproj_validate_panel(unit_col, time_col, N_units, N_times,
        "_wfe_complex_project_core")

    state.nz_obs = sum(W_col :!= 0)
    if (state.nz_obs == 0) {
        errprintf("_wfe_complex_project_core: at least one non-zero weight is required\n")
        _error(498)
    }
    // White option only affects postestimation statistics, not identification
    unit_support_count = J(N_units, 1, 0)
    time_support_count = J(N_times, 1, 0)
    for (k = 1; k <= n; k++) {
        if (W_col[k] == 0) continue
        unit_support_count[unit_col[k]] = 1
        time_support_count[time_col[k]] = 1
    }
    if (sum(unit_support_count) < 2) {
        errprintf("_wfe_complex_project_core: non-zero weights must span at least 2 units\n")
        _error(498)
    }
    if (sum(time_support_count) < 2) {
        errprintf("_wfe_complex_project_core: non-zero weights must span at least 2 time periods\n")
        _error(498)
    }

    if (white_flag != 0) {
        keep_mask = J(n, 1, 1)
    }
    else {
        keep_mask = (W_col :!= 0)
    }

    nz_index = selectindex(keep_mask)
    state.tn_row = rows(nz_index)

    Y_nz = Y_col[nz_index]
    X_nz = X_use[nz_index, .]
    W_nz = W_col[nz_index]
    unit_nz = unit_col[nz_index]
    time_nz = time_col[nz_index]

    state.Udummy = _wfe_complex_unit_dummy(unit_nz, N_units)
    col_keep_u = selectindex(colsum(state.Udummy)' :!= 0)
    if (rows(col_keep_u) == 0) {
        state.Udummy = J(rows(state.Udummy), 0, 0)
    }
    else {
        state.Udummy = state.Udummy[., col_keep_u]
    }
    state.n_Udummy = cols(state.Udummy)

    state.Tdummy = _wfe_complex_time_dummy(time_nz, N_times)
    col_keep_t = selectindex(colsum(state.Tdummy)' :!= 0)
    if (rows(col_keep_t) < cols(state.Tdummy)) {
        if (rows(col_keep_t) == 0) {
            state.Tdummy = J(rows(state.Tdummy), 0, 0)
        }
        else {
            state.Tdummy = state.Tdummy[., col_keep_t]
        }
        state.Tdummy = state.Tdummy, J(rows(state.Tdummy), 1, 1000)
    }
    else {
        if (cols(state.Tdummy) > 0) {
            state.Tdummy = state.Tdummy[., 1..(cols(state.Tdummy) - 1)]
        }
        state.Tdummy = state.Tdummy, J(rows(state.Tdummy), 1, 1000)
    }
    state.n_Tdummy = cols(state.Tdummy)

    state.w_sqrt = _wfe_complex_im_sqrt(W_nz)

    Y_star = _wfe_complex_weight_real_matrix(Y_nz, state.w_sqrt)
    X_star = _wfe_complex_weight_real_matrix(X_nz, state.w_sqrt)
    D = state.Udummy, state.Tdummy
    D_star = _wfe_complex_weight_real_matrix(D, state.w_sqrt)

    if (state.n_Udummy > 0) {
        D1_star = D_star[., 1..state.n_Udummy]
        // The unit block projector must invert D1_star'D1_star, which equals
        // the signed within-unit weight totals under transposeonly() rather
        // than the sum of positive square roots alone.
        unit_weight_sum = diagonal(transposeonly(D1_star) * D1_star)
        inv_weight = _wfe_complex_general_inv(unit_weight_sum, tol)
        Dginv = D1_star :* inv_weight'
        P1 = Dginv * transposeonly(D1_star)
    }
    else {
        P1 = J(state.tn_row, state.tn_row, C(0, 0))
    }

    D2_star = D_star[., (state.n_Udummy + 1)..cols(D_star)]
    Q = D2_star - P1 * D2_star
    QQ_inv = pinv(transposeonly(Q) * Q)

    YX_star = Y_star, X_star
    transformed = YX_star - P1 * YX_star - Q * QQ_inv * transposeonly(Q) * YX_star

    state.y_tilde = transformed[., 1]
    if (p > 0) {
        state.X_tilde = transformed[., 2..cols(transformed)]
    }
    else {
        state.X_tilde = J(state.tn_row, 0, C(0, 0))
    }

    return(state)
}

struct wfe_complex_ols_result scalar _wfe_complex_ols(
    complex vector    y_tilde,
    complex matrix    X_tilde,
    complex vector    w_sqrt,
    real    vector    W_kept,
    real    scalar    white_flag,
    | real  vector    raw_scale,
      string scalar  caller
)
{
    struct wfe_complex_ols_result scalar result
    real scalar n, p, j, scale_j, tol_j, projected_ss, root_tol, root_scale
    real scalar effective_n
    real colvector nz_index, W_use, raw_scale_use
    real colvector root_gap
    complex colvector x_j
    complex colvector expected_root, y_use, w_use
    complex matrix Xty, X_use
    string scalar caller_label

    y_use = colshape(y_tilde, 1)
    w_use = colshape(w_sqrt, 1)
    W_use = colshape(W_kept, 1)
    n = rows(y_use)
    X_use = X_tilde
    if (rows(X_use) == 1 & cols(X_use) == n & n > 1) {
        X_use = transposeonly(X_use)
    }
    p = cols(X_use)
    caller_label = "_wfe_complex_ols"

    if (args() == 7) {
        caller_label = caller
    }

    if (rows(X_use) != n | rows(w_use) != n | rows(W_use) != n) {
        errprintf("%s: input length mismatch\n", caller_label)
        _error(3200)
    }
    if (n == 0) {
        errprintf("%s: inputs must contain at least one observation\n", caller_label)
        _error(3200)
    }
    if (p < 1) {
        errprintf("%s: X_tilde must contain at least one regressor\n", caller_label)
        _error(3200)
    }
    if (missing(white_flag) | (white_flag != 0 & white_flag != 1)) {
        errprintf("%s: white_flag must be 0 or 1\n", caller_label)
        _error(3200)
    }
    if (any(Re(y_use) :>= .) | any(Im(y_use) :>= .) | ///
        any(Re(vec(X_use)) :>= .) | any(Im(vec(X_use)) :>= .) | ///
        any(Re(w_use) :>= .) | any(Im(w_use) :>= .) | ///
        any(W_use :>= .)) {
        errprintf("%s: inputs must not contain missing values\n", caller_label)
        _error(3498)
    }
    if (white_flag == 0 & any(W_use :== 0)) {
        errprintf("%s: W_kept must not contain zero weights when white_flag == 0\n",
                  caller_label)
        _error(3498)
    }
    expected_root = _wfe_complex_im_sqrt(W_use)
    root_gap = sqrt((Re(w_use - expected_root) :^ 2) + (Im(w_use - expected_root) :^ 2))
    root_scale = max((1, max(sqrt(abs(W_use)))))
    root_tol = sqrt(epsilon(1)) * root_scale
    if (max(root_gap) > root_tol) {
        errprintf("%s: w_sqrt must match W_kept via _wfe_complex_im_sqrt()\n",
                  caller_label)
        _error(3498)
    }
    if (args() >= 6) {
        raw_scale_use = colshape(raw_scale, 1)
        if (rows(raw_scale_use) != p | cols(raw_scale_use) != 1) {
            errprintf("%s: raw_scale length mismatch\n", caller_label)
            _error(3200)
        }
        if (any(raw_scale_use :>= .) | any(raw_scale_use :< 0)) {
            errprintf("%s: raw_scale must be nonmissing and non-negative\n",
                      caller_label)
            _error(3498)
        }
    }
    // Check for rank deficiency after projection
    for (j = 1; j <= p; j++) {
        x_j = X_use[., j]
        if (args() >= 6) {
            scale_j = raw_scale_use[j]
        }
        else {
            scale_j = max(sqrt((Re(x_j) :^ 2) + (Im(x_j) :^ 2)))
        }
        tol_j = sqrt(epsilon(1)) * n * (scale_j^2)
        projected_ss = sum((Re(x_j) :^ 2) + (Im(x_j) :^ 2))
        if (scale_j == 0 | projected_ss <= tol_j) {
            errprintf("%s: at least one regressor is fully absorbed by weighted two-way projection\n",
                      caller_label)
            _error(498)
        }
    }
    if (rank(X_use) < p) {
        errprintf("%s: regressors are collinear after weighted two-way projection\n",
                  caller_label)
        _error(498)
    }

    if (p > 0) {
        // The complex-path OLS uses plain (non-conjugate) transpose X^T X,
        // not the Hermitian X^H X. This preserves the WLS interpretation
        // where Re(beta) recovers the signed-weight estimator (Imai-Kim 2021).
        result.ginv_XX_tilde = pinv(transposeonly(X_use) * X_use)
        Xty = transposeonly(X_use) * y_use
        result.betaT = result.ginv_XX_tilde * Xty
    }
    else {
        result.ginv_XX_tilde = J(0, 0, C(0, 0))
        result.betaT = J(0, 1, C(0, 0))
    }

    result.coef_wls = Re(result.betaT)
    result.e_tilde = y_use - X_use * result.betaT

    if (white_flag != 0) {
        result.resid_vec = J(n, 1, C(0, 0))
        nz_index = selectindex(W_use :!= 0)
        if (rows(nz_index) > 0) {
            result.resid_vec[nz_index] = result.e_tilde[nz_index] :/ w_use[nz_index]
        }
    }
    else {
        result.resid_vec = result.e_tilde :/ w_use
    }

    result.diag_ee_tilde = result.e_tilde :* result.e_tilde
    effective_n = sum(W_use :!= 0)
    if (effective_n <= 0) {
        errprintf("%s: at least one non-zero weight is required\n", caller_label)
        _error(498)
    }
    // Zero-weight rows are retained in the White path so the reported
    // residual-row count and sigma2 denominator follow the retained
    // projected sample length; the non-White path uses only non-zero
    // weight observations.
    if (white_flag != 0) {
        result.d_f = n
    }
    else {
        result.d_f = effective_n
    }
    result.sigma2 = Re(sum(result.resid_vec :* result.resid_vec) / result.d_f)

    if (result.sigma2 < 0) {
        printf("warning: sigma2 is negative after complex OLS projection\n")
    }

    return(result)
}

struct wfe_twoway_result scalar wfe_complex_project(
    real vector Y,
    real matrix X,
    real vector W_vec,
    real vector unit_idx,
    real vector time_idx,
    real scalar N,
    real scalar T,
    real scalar len_data,
    real scalar tol
)
{
    struct wfe_twoway_result scalar result
    struct wfe_complex_projection_state scalar state
    struct wfe_complex_ols_result scalar ols_result
    real scalar p, N_nonzero, k
    real colvector Y_col, W_col, unit_col, time_col
    real colvector keep_index, W_kept, raw_scale
    real colvector unit_support_count, time_support_count
    real matrix X_use

    Y_col = colshape(Y, 1)
    W_col = colshape(W_vec, 1)
    unit_col = colshape(unit_idx, 1)
    time_col = colshape(time_idx, 1)
    X_use = X
    if (rows(X_use) == 1 & cols(X_use) == rows(Y_col) & rows(Y_col) > 1) {
        X_use = X_use'
    }
    p = cols(X_use)

    if (missing(len_data) | len_data != floor(len_data) | len_data < 0) {
        errprintf("wfe_complex_project: len_data must be a positive integer\n")
        _error(3200)
    }
    if (rows(Y_col) != len_data | rows(X_use) != len_data | rows(W_col) != len_data | ///
        rows(unit_col) != len_data | rows(time_col) != len_data) {
        errprintf("wfe_complex_project: input length mismatch\n")
        _error(3200)
    }
    if (len_data == 0) {
        errprintf("wfe_complex_project: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (p < 1) {
        errprintf("wfe_complex_project: X must contain at least one regressor\n")
        _error(3200)
    }
    if (N != floor(N) | T != floor(T) | N <= 0 | T <= 0) {
        errprintf("wfe_complex_project: N and T must be positive integers\n")
        _error(3200)
    }
    if (T < 2) {
        errprintf("wfe_complex_project: two-way projection requires at least 2 time periods\n")
        _error(3200)
    }
    if (N < 2) {
        errprintf("wfe_complex_project: two-way projection requires at least 2 units\n")
        _error(3200)
    }
    if (missing(tol) | tol <= 0) {
        errprintf("wfe_complex_project: tol must be positive\n")
        _error(3200)
    }
    if (any(Y_col :>= .)) {
        errprintf("wfe_complex_project: Y must not contain missing values\n")
        _error(3498)
    }
    if (any(X_use :>= .)) {
        errprintf("wfe_complex_project: X must not contain missing values\n")
        _error(3498)
    }
    if (any(W_col :>= .)) {
        errprintf("wfe_complex_project: W_vec must not contain missing values\n")
        _error(3498)
    }
    _wfe_cproj_validate_panel(unit_col, time_col, N, T, "wfe_complex_project")
    N_nonzero = sum(W_col :!= 0)
    if (N_nonzero == 0) {
        errprintf("wfe_complex_project: at least one non-zero weight is required\n")
        _error(498)
    }
    unit_support_count = J(N, 1, 0)
    time_support_count = J(T, 1, 0)
    for (k = 1; k <= len_data; k++) {
        if (W_col[k] == 0) continue
        unit_support_count[unit_col[k]] = 1
        time_support_count[time_col[k]] = 1
    }
    if (sum(unit_support_count) < 2) {
        errprintf("wfe_complex_project: non-zero weights must span at least 2 units\n")
        _error(498)
    }
    if (sum(time_support_count) < 2) {
        errprintf("wfe_complex_project: non-zero weights must span at least 2 time periods\n")
        _error(498)
    }

    state = _wfe_complex_project_core(Y_col, X_use, W_col, unit_col, time_col, N, T, 0, tol)
    keep_index = selectindex(W_col :!= 0)
    W_kept = W_col[keep_index]
    raw_scale = colmax(abs(X_use[keep_index, .]))'
    ols_result = _wfe_complex_ols(state.y_tilde, state.X_tilde, state.w_sqrt,
        W_kept, 0, raw_scale, "wfe_complex_project")

    result.beta = ols_result.coef_wls
    result.vcov = J(p, p, 0)
    result.W = J(T, N, 0)
    for (k = 1; k <= len_data; k++) {
        result.W[time_col[k], unit_col[k]] = W_col[k]
    }
    result.beta_fe = J(p, 1, 0)
    result.vcov_fe = J(p, p, 0)
    if (ols_result.sigma2 >= 0) {
        result.sigma = sqrt(ols_result.sigma2)
    }
    else {
        result.sigma = .
    }
    result.df_r = ols_result.d_f
    result.N_nonzero = state.nz_obs
    // White diagnostics are not computed by this public helper; leave them
    // unavailable rather than publishing a computed-looking zero/false pair.
    result.white_stat = .
    result.white_pval = .
    result.vcetype = "Complex projection OLS pending SE stages"
    result.white_test = ""

    return(result)
}

end
