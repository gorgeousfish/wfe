// ============================================================
// wfe_se_pwfe.mata -- PWFE standard-error backend
// ============================================================

version 16.0
mata:

void _wfe_pwfe_validate_bal_idx(
    real matrix    unit_idx,
    string scalar  caller
)
{
    real scalar i, J_u

    if (rows(unit_idx) != 1 & cols(unit_idx) != 1) {
        errprintf("%s: unit_idx must be a vector\n", caller)
        _error(3200)
    }

    // Ensure unit_idx is a column vector for uniform processing
    unit_idx = colshape(unit_idx[., .], 1)

    if (rows(unit_idx) == 0) {
        errprintf("%s: unit_idx must contain at least one observation\n", caller)
        _error(3200)
    }

    if (any(unit_idx :>= .)) {
        errprintf("%s: unit_idx must not contain missing values\n", caller)
        _error(3498)
    }

    J_u = max(unit_idx)
    if (J_u != floor(J_u) | J_u <= 0) {
        errprintf("%s: unit_idx must contain positive integer indices\n", caller)
        _error(3200)
    }

    for (i = 1; i <= rows(unit_idx); i++) {
        if (unit_idx[i] != floor(unit_idx[i]) | unit_idx[i] < 1) {
            errprintf("%s: unit_idx must contain positive integer indices\n", caller)
            _error(3200)
        }
    }
}

real colvector _wfe_pwfe_unit_sort_order(real matrix unit_idx)
{
    real scalar n

    if (rows(unit_idx) != 1 & cols(unit_idx) != 1) {
        errprintf("_wfe_pwfe_unit_sort_order: unit_idx must be a vector\n")
        _error(3200)
    }

    // Ensure unit_idx is a column vector for uniform processing
    unit_idx = colshape(unit_idx[., .], 1)

    n = rows(unit_idx)
    if (n == 0) {
        return(J(0, 1, .))
    }

    _wfe_pwfe_validate_bal_idx(unit_idx, "_wfe_pwfe_unit_sort_order")

    if (n == 1) {
        return(1::n)
    }

    // Sort by unit index, then by original position to ensure stable ordering
    return(order((unit_idx, (1::n)), (1, 2)))
}

real scalar _wfe_pwfe_is_balanced(real matrix unit_idx)
{
    real matrix info
    real colvector counts, sort_order, unit_sorted

    if (rows(unit_idx) != 1 & cols(unit_idx) != 1) {
        errprintf("_wfe_pwfe_is_balanced: unit_idx must be a vector\n")
        _error(3200)
    }
    unit_idx = colshape(unit_idx[., .], 1)
    _wfe_pwfe_validate_bal_idx(unit_idx, "_wfe_pwfe_is_balanced")

    sort_order = _wfe_pwfe_unit_sort_order(unit_idx)
    unit_sorted = unit_idx[sort_order]
    info = panelsetup(unit_sorted, 1)
    if (rows(info) == 0) {
        return(1)
    }

    counts = info[,2] - info[,1] :+ 1
    return(all(counts :== counts[1]))
}


real scalar _wfe_pwfe_balanced_block_length(real matrix unit_idx)
{
    real matrix info
    real colvector counts, sort_order, unit_sorted

    if (rows(unit_idx) != 1 & cols(unit_idx) != 1) {
        errprintf("_wfe_pwfe_balanced_block_length: unit_idx must be a vector\n")
        _error(3200)
    }
    unit_idx = colshape(unit_idx[., .], 1)
    _wfe_pwfe_validate_bal_idx(unit_idx, "_wfe_pwfe_balanced_block_length")

    sort_order = _wfe_pwfe_unit_sort_order(unit_idx)
    unit_sorted = unit_idx[sort_order]
    info = panelsetup(unit_sorted, 1)
    if (rows(info) == 0) {
        return(0)
    }

    counts = info[,2] - info[,1] :+ 1
    if (!all(counts :== counts[1])) {
        return(.)
    }

    return(counts[1])
}


real colvector _wfe_pwfe_compact_unit_idx_core(
    real vector    unit_idx,
    real scalar    J_u,
    string scalar  argname,
    string scalar  caller
)
{
    real scalar i, j, next_idx, found_slot
    real scalar label_j
    real colvector compact_idx, uniq_idx

    unit_idx = colshape(unit_idx[., .], 1)

    if (J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("%s: J_u must be a positive integer\n", caller)
        _error(3200)
    }
    if (any(unit_idx :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, argname)
        _error(3498)
    }
    if (rows(unit_idx) == 0) {
        errprintf("%s: %s must contain at least one observation\n", caller, argname)
        _error(3200)
    }

    for (i = 1; i <= rows(unit_idx); i++) {
        if (unit_idx[i] != floor(unit_idx[i]) | unit_idx[i] < 1) {
            errprintf("%s: %s must contain positive integer unit labels\n",
                      caller, argname)
            _error(3200)
        }
    }

    uniq_idx = uniqrows(sort(unit_idx, 1))
    if (rows(uniq_idx) != J_u) {
        errprintf("%s: %s must identify exactly J_u unique units\n",
                  caller, argname)
        _error(3200)
    }

    compact_idx = J(rows(unit_idx), 1, .)
    uniq_idx = J(J_u, 1, .)
    next_idx = 0

    for (j = 1; j <= rows(unit_idx); j++) {
        label_j = unit_idx[j]
        found_slot = 0
        for (i = 1; i <= next_idx; i++) {
            if (uniq_idx[i] == label_j) {
                found_slot = i
                break
            }
        }
        if (found_slot == 0) {
            next_idx = next_idx + 1
            uniq_idx[next_idx] = label_j
            found_slot = next_idx
        }
        compact_idx[j] = found_slot
    }

    return(compact_idx)
}

void _wfe_pwfe_validate_unit_idx_core(
    real vector    unit_idx,
    real scalar    J_u,
    string scalar  argname,
    string scalar  caller
)
{
    _wfe_pwfe_compact_unit_idx_core(unit_idx, J_u, argname, caller)
}

void _wfe_pwfe_validate_unit_idx(
    real vector    unit_idx,
    real scalar    J_u,
    string scalar  argname
)
{
    _wfe_pwfe_validate_unit_idx_core(unit_idx, J_u, argname,
                                     "_wfe_pwfe_validate_unit_idx")
}

void _wfe_pwfe_nm_mat_core(
    real matrix    X,
    string scalar  argname,
    string scalar  caller
)
{
    if (any(X :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, argname)
        _error(3498)
    }
}

void _wfe_pwfe_nm_mat(
    real matrix    X,
    string scalar  argname
)
{
    _wfe_pwfe_nm_mat_core(X, argname, "_wfe_pwfe_nm_mat")
}


real matrix _wfe_pwfe_sym_mat_core(
    real matrix    X,
    string scalar  argname,
    string scalar  caller
)
{
    real scalar asym, scale

    asym = max(abs(vec(X - X')))
    scale = max((1, max(abs(vec(X)))))
    if (asym > sqrt(epsilon(1)) * scale) {
        errprintf("%s: %s must be symmetric\n", caller, argname)
        _error(3200)
    }

    return(0.5 :* (X + X'))
}

real matrix _wfe_pwfe_sym_mat(
    real matrix    X,
    string scalar  argname
)
{
    return(_wfe_pwfe_sym_mat_core(X, argname, "_wfe_pwfe_sym_mat"))
}


void _wfe_pwfe_psd_mat_core(
    real matrix    X,
    string scalar  argname,
    string scalar  caller
)
{
    real rowvector evals
    real scalar eig_tol, eig_scale

    symeigensystem(X, ., evals)
    eig_scale = max((1, max(abs(evals))))
    eig_tol = sqrt(epsilon(1)) * eig_scale

    if (min(evals) < -eig_tol) {
        errprintf("%s: %s must be positive semidefinite\n", caller, argname)
        _error(3200)
    }
}

void _wfe_pwfe_psd_mat(
    real matrix    X,
    string scalar  argname
)
{
    _wfe_pwfe_psd_mat_core(X, argname, "_wfe_pwfe_psd_mat")
}


void _wfe_pwfe_nm_vec_core(
    real matrix    x,
    string scalar  argname,
    string scalar  caller
)
{
    if (rows(x) != 1 & cols(x) != 1) {
        errprintf("%s: %s must be a vector\n", caller, argname)
        _error(3200)
    }

    x = colshape(x[., .], 1)

    if (any(x :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, argname)
        _error(3498)
    }
}

void _wfe_pwfe_nm_vec(
    real matrix    x,
    string scalar  argname
)
{
    _wfe_pwfe_nm_vec_core(x, argname, "_wfe_pwfe_nm_vec")
}


real matrix _wfe_pwfe_omega_hc(
    real matrix    X,
    real vector    u,
    real scalar    J_u
)
{
    real scalar NT, p, denom
    real colvector u_col

    // Ensure one-dimensional inputs are treated as column vectors
    u_col = colshape(u[., .], 1)
    if (rows(X) == 1 & cols(X) == rows(u_col) & rows(u_col) > 1) {
        // Convert single-regressor row matrix to column vector
        X = colshape(X[., .], 1)
    }

    if (rows(X) != rows(u_col)) {
        errprintf("_wfe_pwfe_omega_hc: X (%g rows) and u (%g rows) mismatch\n",
                  rows(X), rows(u_col))
        _error(3200)
    }
    if (rows(X) == 0) {
        errprintf("_wfe_pwfe_omega_hc: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("_wfe_pwfe_omega_hc: X must not contain missing values\n")
        _error(3498)
    }
    if (any(u_col :>= .)) {
        errprintf("_wfe_pwfe_omega_hc: u must not contain missing values\n")
        _error(3498)
    }
    if (J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("_wfe_pwfe_omega_hc: J_u must be a positive integer\n")
        _error(3200)
    }

    NT = rows(X)
    p = cols(X)
    if (p <= 0) {
        errprintf("_wfe_pwfe_omega_hc: X must contain at least one regressor\n")
        _error(3200)
    }
    if (J_u > NT) {
        errprintf("_wfe_pwfe_omega_hc: J_u must not exceed the number of observations\n")
        _error(3200)
    }
    denom = NT - J_u - p

    if (denom <= 0) {
        errprintf("_wfe_pwfe_omega_hc: NT - J_u - p must be positive; got %g\n",
                  denom)
        _error(498)
    }

    return((1 / denom) * cross(X :* (u_col :^ 2), X))
}


void _wfe_pwfe_sw_refchk(
    real matrix    X,
    real colvector u_col,
    real scalar    J_u,
    real matrix    Omega_HC,
    real matrix    ginv_XX
)
{
    real matrix omega_expected, ginv_expected
    real scalar rel_omega, rel_ginv
    real scalar omega_scale, ginv_scale

    omega_expected = _wfe_pwfe_omega_hc(X, u_col, J_u)
    ginv_expected = pinv(cross(X, X))

    omega_scale = max((1, max(abs(vec(omega_expected)))))
    rel_omega = max(abs(vec(Omega_HC - omega_expected))) / omega_scale
    if (rel_omega > sqrt(epsilon(1))) {
        errprintf("_wfe_pwfe_stockwatson_bias: Omega_HC must match _wfe_pwfe_omega_hc(X, u, J_u)\n")
        _error(498)
    }

    ginv_scale = max((1, max(abs(vec(ginv_expected)))))
    rel_ginv = max(abs(vec(ginv_XX - ginv_expected))) / ginv_scale
    if (rel_ginv > sqrt(epsilon(1))) {
        errprintf("_wfe_pwfe_stockwatson_bias: ginv_XX must match pinv(X'X)\n")
        _error(498)
    }
}


void _wfe_pwfe_stockwatson_bias(
    real matrix    X,
    real vector    u,
    real vector    unit_idx,
    real scalar    J_u,
    real scalar    J_t,
    real matrix    Omega_HC,
    real matrix    ginv_XX,
    real matrix    Psi_sw
)
{
    real scalar p, i, omega_asym, omega_scale, ginv_asym, ginv_scale
    real scalar omega_eig_tol, omega_eig_scale, ginv_eig_tol, ginv_eig_scale
    real matrix B_hat, info, Xi, XX_i, Sigma_HRFE, Bread, X_sorted, Omega_HC_use, ginv_XX_use
    real colvector ui, block_lengths, sort_order, unit_sorted, u_sorted
    real colvector u_col, unit_col
    real rowvector omega_evals, ginv_evals

    // Ensure one-dimensional inputs are treated as column vectors
    u_col = colshape(u[., .], 1)
    unit_col = colshape(unit_idx[., .], 1)
    if (rows(X) == 1 ///
        & cols(X) == rows(u_col) ///
        & cols(X) == rows(unit_col) ///
        & rows(u_col) > 1) {
        // Convert single-regressor row matrix to column format
        X = _wfe_precompute_as_design(X, rows(u_col))
    }

    p = cols(X)

    if (p <= 0) {
        errprintf("_wfe_pwfe_stockwatson_bias: X must contain at least one regressor\n")
        _error(3200)
    }

    if (rows(X) != rows(u_col) | rows(X) != rows(unit_col)) {
        errprintf("_wfe_pwfe_stockwatson_bias: X, u, and unit_idx must align\n")
        _error(3200)
    }
    if (rows(X) == 0) {
        errprintf("_wfe_pwfe_stockwatson_bias: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("_wfe_pwfe_stockwatson_bias: X must not contain missing values\n")
        _error(3498)
    }
    if (any(u_col :>= .)) {
        errprintf("_wfe_pwfe_stockwatson_bias: u must not contain missing values\n")
        _error(3498)
    }
    if (J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("_wfe_pwfe_stockwatson_bias: J_u must be a positive integer\n")
        _error(3200)
    }
    if (J_u < 2) {
        errprintf("_wfe_pwfe_stockwatson_bias: need at least 2 units for Stock-Watson standard errors\n")
        _error(498)
    }
    if (J_t >= . | J_t != floor(J_t) | J_t <= 0) {
        errprintf("_wfe_pwfe_stockwatson_bias: J_t must be a positive integer\n")
        _error(3200)
    }
    if (rows(Omega_HC) != p | cols(Omega_HC) != p) {
        errprintf("_wfe_pwfe_stockwatson_bias: Omega_HC must be %gx%g to match cols(X)\n",
                  p, p)
        _error(3200)
    }
    if (rows(ginv_XX) != p | cols(ginv_XX) != p) {
        errprintf("_wfe_pwfe_stockwatson_bias: ginv_XX must be %gx%g to match cols(X)\n",
                  p, p)
        _error(3200)
    }
    if (any(Omega_HC :>= .)) {
        errprintf("_wfe_pwfe_stockwatson_bias: Omega_HC must not contain missing values\n")
        _error(3498)
    }
    if (any(ginv_XX :>= .)) {
        errprintf("_wfe_pwfe_stockwatson_bias: ginv_XX must not contain missing values\n")
        _error(3498)
    }
    omega_asym = max(abs(vec(Omega_HC - Omega_HC')))
    omega_scale = max((1, max(abs(vec(Omega_HC)))))
    if (omega_asym > sqrt(epsilon(1)) * omega_scale) {
        errprintf("_wfe_pwfe_stockwatson_bias: Omega_HC must be symmetric\n")
        _error(3200)
    }
    // Symmetrize the covariance matrix to remove numerical roundoff
    Omega_HC_use = 0.5 :* (Omega_HC + Omega_HC')
    symeigensystem(Omega_HC_use, ., omega_evals)
    omega_eig_scale = max((1, max(abs(omega_evals))))
    omega_eig_tol = sqrt(epsilon(1)) * omega_eig_scale
    if (min(omega_evals) < -omega_eig_tol) {
        errprintf("_wfe_pwfe_stockwatson_bias: Omega_HC must be positive semidefinite\n")
        _error(3200)
    }
    ginv_asym = max(abs(vec(ginv_XX - ginv_XX')))
    ginv_scale = max((1, max(abs(vec(ginv_XX)))))
    if (ginv_asym > sqrt(epsilon(1)) * ginv_scale) {
        errprintf("_wfe_pwfe_stockwatson_bias: ginv_XX must be symmetric\n")
        _error(3200)
    }
    // Symmetrize the matrix to remove numerical roundoff
    ginv_XX_use = 0.5 :* (ginv_XX + ginv_XX')
    symeigensystem(ginv_XX_use, ., ginv_evals)
    ginv_eig_scale = max((1, max(abs(ginv_evals))))
    ginv_eig_tol = sqrt(epsilon(1)) * ginv_eig_scale
    if (min(ginv_evals) < -ginv_eig_tol) {
        errprintf("_wfe_pwfe_stockwatson_bias: ginv_XX must be positive semidefinite\n")
        _error(3200)
    }
    _wfe_pwfe_sw_refchk(X, u_col, J_u, Omega_HC_use, ginv_XX_use)
    unit_col = _wfe_pwfe_compact_unit_idx_core(unit_col, J_u, "unit_idx",
                                             "_wfe_pwfe_stockwatson_bias")

    if (J_t < 3) {
        errprintf("_wfe_pwfe_stockwatson_bias: need at least 3 time periods for Stock-Watson standard errors\n")
        _error(498)
    }

    sort_order = _wfe_pwfe_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    X_sorted = X[sort_order, .]
    u_sorted = u_col[sort_order]

    B_hat = J(p, p, 0)
    info = panelsetup(unit_sorted, 1)
    block_lengths = info[, 2] - info[, 1] :+ 1

    if (any(block_lengths :!= J_t)) {
        errprintf("_wfe_pwfe_stockwatson_bias: unit_idx must define a balanced panel with exactly J_t observations per unit\n")
        _error(498)
    }

    for (i = 1; i <= rows(info); i++) {
        Xi = panelsubmatrix(X_sorted, i, info)
        if (rows(Xi) > 1) {
            ui = panelsubmatrix(u_sorted, i, info)
            XX_i = cross(Xi, Xi)
            B_hat = B_hat + (1 / J_t) * XX_i * (1 / (J_t - 1)) * sum(ui :^ 2)
        }
    }

    B_hat = B_hat * (1 / J_u)
    Sigma_HRFE = ((J_t - 1) / (J_t - 2)) * (Omega_HC_use - (1 / (J_t - 1)) * B_hat)
    Bread = rows(X) * ginv_XX_use
    Psi_sw = Bread * Sigma_HRFE * Bread
}


void _wfe_pwfe_compute_se(
    real matrix    X_tilde,
    real vector    u_tilde,
    real vector    unit_idx_tilde,
    real matrix    X_hat,
    real vector    u_hat,
    real vector    unit_idx_hat,
    real scalar    J_u,
    real scalar    J_t,
    real vector    resid_final,
    real scalar    fit_final_df,
    real scalar    unit_number,
    real matrix    vcov_ols,
    string scalar  hetero_se,
    string scalar  auto_se,
    string scalar  unbiased_se,
    string scalar  verbose,
    real matrix    vcov_wfe,
    real matrix    vcov_fe,
    real matrix    Psi_hat_wfe,
    real matrix    Psi_hat_fe,
    real matrix    ginv_XX_tilde,
    real matrix    ginv_XX_hat,
    string scalar  se_type
)
{
    real scalar NT, NT_fe, sigma2, actual_J_t_tilde, actual_J_t_hat
    real matrix Omega_raw_wfe, Omega_raw_fe, Omega_wfe, Omega_fe
    real matrix Bread_wfe, Bread_fe, Psi_sw_wfe, Psi_sw_fe

    // Ensure one-dimensional inputs are treated as column vectors
    u_tilde = colshape(u_tilde[., .], 1)
    unit_idx_tilde = colshape(unit_idx_tilde[., .], 1)
    u_hat = colshape(u_hat[., .], 1)
    unit_idx_hat = colshape(unit_idx_hat[., .], 1)
    resid_final = colshape(resid_final[., .], 1)

    // Canonicalize only when a 1 x NT row matrix is provably a single-regressor
    // observation stream. A genuine 1 x p single-observation design must stay
    // 1 x p to preserve the crossprod() contract.
    X_tilde = _wfe_precompute_as_design(X_tilde, rows(u_tilde))
    X_hat   = _wfe_precompute_as_design(X_hat,   rows(u_hat))

    NT = rows(X_tilde)
    NT_fe = rows(X_hat)

    if (rows(X_tilde) != rows(u_tilde) | rows(X_tilde) != rows(unit_idx_tilde)) {
        errprintf("_wfe_pwfe_compute_se: X_tilde, u_tilde, and unit_idx_tilde must align\n")
        _error(3200)
    }
    if (rows(X_hat) != rows(u_hat) | rows(X_hat) != rows(unit_idx_hat)) {
        errprintf("_wfe_pwfe_compute_se: X_hat, u_hat, and unit_idx_hat must align\n")
        _error(3200)
    }
    if (rows(X_tilde) != rows(X_hat)) {
        errprintf("_wfe_pwfe_compute_se: X_tilde and X_hat must have the same number of rows\n")
        _error(3200)
    }
    if (cols(X_tilde) != cols(X_hat)) {
        errprintf("_wfe_pwfe_compute_se: X_tilde (%g cols) and X_hat (%g cols) mismatch\n",
                  cols(X_tilde), cols(X_hat))
        _error(3200)
    }
    if (NT == 0) {
        errprintf("_wfe_pwfe_compute_se: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (cols(X_tilde) == 0) {
        errprintf("_wfe_pwfe_compute_se: X_tilde and X_hat must each contain at least one regressor\n")
        _error(3200)
    }
    _wfe_pwfe_nm_mat_core(X_tilde, "X_tilde", "_wfe_pwfe_compute_se")
    _wfe_pwfe_nm_vec_core(u_tilde, "u_tilde", "_wfe_pwfe_compute_se")
    _wfe_pwfe_nm_mat_core(X_hat, "X_hat", "_wfe_pwfe_compute_se")
    _wfe_pwfe_nm_vec_core(u_hat, "u_hat", "_wfe_pwfe_compute_se")
    if (hetero_se != "on" & hetero_se != "off") {
        errprintf("_wfe_pwfe_compute_se: hetero_se must be on or off\n")
        _error(3200)
    }
    if (auto_se != "on" & auto_se != "off") {
        errprintf("_wfe_pwfe_compute_se: auto_se must be on or off\n")
        _error(3200)
    }
    if (unbiased_se != "on" & unbiased_se != "off") {
        errprintf("_wfe_pwfe_compute_se: unbiased_se must be on or off\n")
        _error(3200)
    }
    if (verbose != "on" & verbose != "off") {
        errprintf("_wfe_pwfe_compute_se: verbose must be on or off\n")
        _error(3200)
    }
    if (unbiased_se == "on" & (hetero_se != "on" | auto_se != "off")) {
        errprintf("_wfe_pwfe_compute_se: unbiased_se(on) requires hetero_se(on) and auto_se(off)\n")
        _error(198)
    }
    if (J_t >= . | J_t != floor(J_t) | J_t <= 0) {
        errprintf("_wfe_pwfe_compute_se: J_t must be a positive integer\n")
        _error(3200)
    }
    unit_idx_tilde = _wfe_pwfe_compact_unit_idx_core(unit_idx_tilde, J_u,
                                                  "unit_idx_tilde",
                                                  "_wfe_pwfe_compute_se")
    unit_idx_hat = _wfe_pwfe_compact_unit_idx_core(unit_idx_hat, J_u,
                                                "unit_idx_hat",
                                                "_wfe_pwfe_compute_se")
    if (any(unit_idx_tilde :!= unit_idx_hat)) {
        errprintf("_wfe_pwfe_compute_se: unit_idx_tilde and unit_idx_hat must align rowwise\n")
        _error(3200)
    }
    if (hetero_se == "off" & auto_se == "on") {
        errprintf("_wfe_pwfe_compute_se: robust standard errors with autocorrelation and homoskedasticity is not supported\n")
        _error(198)
    }
    if (rank(X_tilde) < cols(X_tilde)) {
        errprintf("_wfe_pwfe_compute_se: X_tilde regressors are collinear\n")
        _error(498)
    }
    if (rank(X_hat) < cols(X_hat)) {
        errprintf("_wfe_pwfe_compute_se: X_hat regressors are collinear\n")
        _error(498)
    }

    _wfe_precompute_bread(X_tilde, X_hat, ginv_XX_tilde, ginv_XX_hat, NT)

    if (hetero_se == "on" & auto_se == "on") {
        se_type = "Heteroscedastic / Autocorrelation Robust Standard Error"

        if (J_u < 2) {
            errprintf("_wfe_pwfe_compute_se: need at least 2 units for HAC standard errors\n")
            _error(498)
        }

        Omega_raw_wfe = _wfe_omega_hac(X_tilde, u_tilde, unit_idx_tilde, J_u)
        Omega_raw_fe = _wfe_omega_hac(X_hat, u_hat, unit_idx_hat, J_u)

        Omega_wfe = (1 / NT) * Omega_raw_wfe
        Omega_fe = (1 / NT_fe) * Omega_raw_fe

        Bread_wfe = NT * ginv_XX_tilde
        Bread_fe = NT_fe * ginv_XX_hat
        Psi_hat_wfe = Bread_wfe * Omega_wfe * Bread_wfe
        Psi_hat_fe = Bread_fe * Omega_fe * Bread_fe
    }
    else if (hetero_se == "on" & auto_se == "off") {
        se_type = "Heteroscedastic Robust Standard Error"

        if (unbiased_se == "on") {
            Omega_wfe = _wfe_pwfe_omega_hc(X_tilde, u_tilde, J_u)
            Omega_fe = _wfe_pwfe_omega_hc(X_hat, u_hat, J_u)

            if (!_wfe_pwfe_is_balanced(unit_idx_tilde)) {
                errprintf("unbiased_se requires a balanced WFE panel\n")
                _error(498)
            }
            if (!_wfe_pwfe_is_balanced(unit_idx_hat)) {
                errprintf("unbiased_se requires a balanced FE panel\n")
                _error(498)
            }
            actual_J_t_tilde = _wfe_pwfe_balanced_block_length(unit_idx_tilde)
            actual_J_t_hat = _wfe_pwfe_balanced_block_length(unit_idx_hat)
            if (actual_J_t_tilde != J_t | actual_J_t_hat != J_t) {
                errprintf("unbiased_se: J_t must match the balanced panel block length\n")
                _error(498)
            }
            if (J_t < 3) {
                errprintf("unbiased_se requires at least 3 time periods\n")
                _error(498)
            }

            se_type = "Heteroskedastic Standard Error (Stock-Watson Bias-Corrected)"
            _wfe_pwfe_stockwatson_bias(X_tilde, u_tilde, unit_idx_tilde, J_u, J_t,
                                       Omega_wfe, ginv_XX_tilde, Psi_sw_wfe)
            _wfe_pwfe_stockwatson_bias(X_hat, u_hat, unit_idx_hat, J_u, J_t,
                                       Omega_fe, ginv_XX_hat, Psi_sw_fe)

            if (verbose == "on") {
                printf("time %g\n", J_t)
            }
        }
        else {
            // HC variance remains clustered on unit panels even with time-demeaned data
            Omega_wfe = _wfe_pwfe_omega_hc(X_tilde, u_tilde, J_u)
            Omega_fe = _wfe_pwfe_omega_hc(X_hat, u_hat, J_u)
        }

        if (unbiased_se == "on") {
            Psi_hat_wfe = Psi_sw_wfe
            Psi_hat_fe = Psi_sw_fe
        }
        else {
            Bread_wfe = NT * ginv_XX_tilde
            Bread_fe = NT_fe * ginv_XX_hat
            Psi_hat_wfe = Bread_wfe * Omega_wfe * Bread_wfe
            Psi_hat_fe = Bread_fe * Omega_fe * Bread_fe
        }
    }
    else if (hetero_se == "off" & auto_se == "off") {
        if (rows(resid_final) == 0) {
            errprintf("_wfe_pwfe_compute_se: resid_final cannot be empty\n")
            _error(3200)
        }
        if (rows(vcov_ols) != cols(X_hat) | cols(vcov_ols) != cols(X_hat)) {
            errprintf("_wfe_pwfe_compute_se: vcov_ols must be %gx%g; got %gx%g\n",
                      cols(X_hat), cols(X_hat), rows(vcov_ols), cols(vcov_ols))
            _error(3200)
        }
        _wfe_pwfe_nm_mat_core(vcov_ols, "vcov_ols", "_wfe_pwfe_compute_se")
        vcov_ols = _wfe_pwfe_sym_mat_core(vcov_ols, "vcov_ols",
                                          "_wfe_pwfe_compute_se")
        _wfe_pwfe_psd_mat_core(vcov_ols, "vcov_ols", "_wfe_pwfe_compute_se")
        _wfe_pwfe_nm_vec_core(resid_final, "resid_final",
                              "_wfe_pwfe_compute_se")
        if (rows(resid_final) != NT) {
            errprintf("_wfe_pwfe_compute_se: resid_final must align with X_tilde in homoskedastic branch\n")
            _error(3200)
        }
        if (fit_final_df >= . | fit_final_df != floor(fit_final_df) | fit_final_df <= 0) {
            errprintf("_wfe_pwfe_compute_se: fit_final_df must be a positive integer\n")
            _error(3200)
        }
        if (fit_final_df != (NT - cols(X_tilde))) {
            errprintf("_wfe_pwfe_compute_se: fit_final_df must equal rows(X_tilde) - cols(X_tilde)\n")
            _error(3200)
        }
        if (unit_number >= . | unit_number != floor(unit_number) | unit_number <= 0) {
            errprintf("_wfe_pwfe_compute_se: unit_number must be a positive integer\n")
            _error(3200)
        }
        if (unit_number != J_u & unit_number != J_t) {
            errprintf("_wfe_pwfe_compute_se: unit_number must equal J_u or J_t\n")
            _error(3200)
        }
        if (fit_final_df <= unit_number) {
            errprintf("_wfe_pwfe_compute_se: fit_final_df - unit_number must be positive\n")
            _error(498)
        }

        se_type = "Homoskedastic Standard Error"
        sigma2 = sum(resid_final :^ 2) / (fit_final_df - unit_number)

        Psi_hat_wfe = NT * (sigma2 * ginv_XX_tilde)
        Psi_hat_fe = NT_fe * vcov_ols
    }
    else {
        errprintf("_wfe_pwfe_compute_se: robust standard errors with autocorrelation and homoskedasticity is not supported\n")
        _error(198)
    }

    vcov_wfe = Psi_hat_wfe * (1 / NT)
    vcov_fe = Psi_hat_fe * (1 / NT_fe)
}

end
