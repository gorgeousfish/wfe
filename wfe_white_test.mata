// ============================================================
// wfe_white_test.mata - White (1980) misspecification tests
//
// White test nests:
//   H0: No misspecification (WFE = Standard FE)
//   H1: Misspecification (WFE != Standard FE)
// ============================================================

version 16.0
mata:
mata set matastrict on

void _wfe_white_psd_mat(
    real matrix    X,
    string scalar  caller,
    string scalar  argname
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

real colvector _wfe_white_compact_labels(
    real colvector  idx,
    real scalar     n_levels,
    string scalar   caller,
    string scalar   argname
)
{
    real colvector levels, compact_idx, mask
    real scalar k

    if (any(idx :!= floor(idx)) | any(idx :< 1)) {
        errprintf("%s: %s must contain integer indices in 1..N_%s\n",
            caller, argname, substr(argname, 1, 1) == "u" ? "units" : "times")
        _error(3200)
    }

    levels = uniqrows(sort(idx, 1))
    if (rows(levels) != n_levels) {
        errprintf("%s: %s must enumerate 1..N_%s without gaps\n",
            caller, argname, substr(argname, 1, 1) == "u" ? "units" : "times")
        _error(3200)
    }

    compact_idx = J(rows(idx), 1, .)
    for (k = 1; k <= rows(levels); k++) {
        mask = selectindex(idx :== levels[k])
        compact_idx[mask] = J(rows(mask), 1, k)
    }

    return(compact_idx)
}

void _wfe_white_require_square_dim(
    real matrix    X,
    real scalar    p,
    string scalar  caller,
    string scalar  argname
)
{
    if (rows(X) != p | cols(X) != p) {
        errprintf("%s: %s must be %gx%g to match cols(X_hat)\n", caller, argname, p, p)
        _error(3200)
    }
}

void _wfe_white_require_ginv_match(
    real matrix    X,
    real matrix    ginv_X_use,
    string scalar  caller,
    string scalar  argname,
    string scalar  design_name
)
{
    real matrix gram_X, ginv_expected
    real scalar ginv_gap, ginv_scale, ginv_tol

    gram_X = cross(X, X)
    ginv_expected = pinv(gram_X)
    ginv_gap = max(abs(vec(ginv_X_use - ginv_expected)))
    ginv_scale = max((1, max(abs(vec(ginv_expected))),
        max(abs(vec(ginv_X_use)))))
    ginv_tol = sqrt(epsilon(1)) * rows(X) * ginv_scale

    if (ginv_gap > ginv_tol) {
        errprintf("%s: %s must match pinv(%s'%s)\n",
            caller, argname, design_name, design_name)
        _error(498)
    }
}

real matrix _wfe_white_phi_guard(
    real matrix    Phi,
    real scalar    p,
    string scalar  caller
)
{
    real matrix Phi_use
    real rowvector phi_evals
    real scalar phi_scale, phi_tol

    Phi_use = 0.5 :* (Phi + Phi')
    phi_evals = symeigenvalues(Phi_use)
    phi_scale = max((1, max(abs(phi_evals))))
    phi_tol = sqrt(epsilon(1)) * phi_scale

    if (min(phi_evals) < -phi_tol) {
        printf("{err}%s: Phi is not positive semidefinite; White test skipped\n", caller)
        return(J(0, 0, .))
    }
    if (rank(Phi_use) < p) {
        printf("{err}%s: Phi is rank-deficient; White test skipped\n", caller)
        return(J(0, 0, .))
    }

    return(Phi_use)
}

// ============================================================
// wfe_white_test_oneway - One-way FE White misspecification test
//
// Parameters:
//   beta_wfe       p x 1 WFE coefficient vector
//   beta_fe        p x 1 standard FE coefficient vector
//   Psi_wfe        p x p WFE sandwich matrix before dividing by J_u
//   Psi_fe         p x p FE sandwich matrix before dividing by J_u
//   ginv_XX_tilde  p x p = pinv(X_tilde'X_tilde)
//   ginv_XX_hat    p x p = pinv(X_hat'X_hat)
//   X_tilde        NT x p WFE-demeaned X
//   X_hat          NT x p standard FE-demeaned X
//   u_tilde        NT x 1 WFE residuals
//   u_hat          NT x 1 standard FE residuals
//   J_u            scalar: number of units (always N_units)
//   p              scalar: number of covariates
//   N_nonzero      scalar: retained for caller compatibility; one-way White
//                  scaling uses NT (the full estimation-sample row count)
//   white_alpha    scalar: significance level (default 0.05)
//
// One-way White uses the estimation-sample row count NT in the bread-Lambda
// products and final quadratic form.
//
// Returns: (white_stat, white_pvalue, test_reject) as a 1 x 3 row vector
// ============================================================
real rowvector wfe_white_test_oneway(
    real vector    beta_wfe,
    real vector    beta_fe,
    real matrix    Psi_wfe,
    real matrix    Psi_fe,
    real matrix    ginv_XX_tilde,
    real matrix    ginv_XX_hat,
    real matrix    X_tilde,
    real matrix    X_hat,
    real vector    u_tilde,
    real vector    u_hat,
    real scalar    J_u,
    real scalar    p,
    real scalar    N_nonzero,
    real scalar    white_alpha,
    | real vector  unit_idx)
{
    real scalar NT, hat_asym, hat_scale, tilde_asym, tilde_scale
    real scalar psi_wfe_asym, psi_wfe_scale, psi_fe_asym, psi_fe_scale
    real colvector beta_diff, diag_ee, unit_col_use
    real colvector beta_wfe_col, beta_fe_col, u_hat_col, u_tilde_col
    real matrix Lambda1, Lambda2, Phi
    real matrix X_hat_use, X_tilde_use
    real matrix ginv_XX_hat_use, ginv_XX_tilde_use
    real matrix Psi_wfe_use, Psi_fe_use
    real scalar white_stat, white_pvalue, white_tol
    real scalar test_reject, denom

    beta_wfe_col = colshape(beta_wfe[., .], 1)
    beta_fe_col = colshape(beta_fe[., .], 1)
    u_hat_col = colshape(u_hat[., .], 1)
    u_tilde_col = colshape(u_tilde[., .], 1)
    X_hat_use = X_hat
    X_tilde_use = X_tilde

    if (rows(X_hat_use) == 1 & cols(X_hat_use) == rows(u_hat_col) & rows(u_hat_col) > 1) {
        X_hat_use = X_hat_use'
    }
    if (rows(X_tilde_use) == 1 & cols(X_tilde_use) == rows(u_tilde_col) & rows(u_tilde_col) > 1) {
        X_tilde_use = X_tilde_use'
    }

    NT = rows(X_hat_use)

    if (rows(beta_wfe_col) != rows(beta_fe_col) | cols(beta_wfe_col) != 1 | cols(beta_fe_col) != 1) {
        errprintf("wfe_white_test_oneway: coefficient dimension mismatch\n")
        _error(3200)
    }
    if (J_u != floor(J_u) | J_u <= 1) {
        errprintf("wfe_white_test_oneway: J_u must be an integer >= 2\n")
        _error(3200)
    }
    if (rows(X_hat_use) != rows(X_tilde_use) | NT != rows(u_hat_col) | NT != rows(u_tilde_col)) {
        errprintf("wfe_white_test_oneway: input length mismatch\n")
        _error(3200)
    }
    if (cols(X_hat_use) != p | cols(X_tilde_use) != p) {
        errprintf("wfe_white_test_oneway: X_hat and X_tilde must have p columns\n")
        _error(3200)
    }
    if (white_alpha <= 0 | white_alpha >= 1) {
        errprintf("wfe_white_test_oneway: white_alpha must lie in (0,1)\n")
        _error(3200)
    }
    _wfe_white_require_square_dim(ginv_XX_hat, p, "wfe_white_test_oneway", "ginv_XX_hat")
    _wfe_white_require_square_dim(ginv_XX_tilde, p, "wfe_white_test_oneway", "ginv_XX_tilde")
    _wfe_white_require_square_dim(Psi_wfe, p, "wfe_white_test_oneway", "Psi_wfe")
    _wfe_white_require_square_dim(Psi_fe, p, "wfe_white_test_oneway", "Psi_fe")

    hat_asym  = max(abs(vec(ginv_XX_hat - ginv_XX_hat')))
    hat_scale = max((1, max(abs(vec(ginv_XX_hat)))))
    if (hat_asym > sqrt(epsilon(1)) * hat_scale) {
        errprintf("wfe_white_test_oneway: ginv_XX_hat must be symmetric\n")
        _error(3200)
    }
    ginv_XX_hat_use = 0.5 :* (ginv_XX_hat + ginv_XX_hat')
    _wfe_white_require_ginv_match(X_hat_use, ginv_XX_hat_use,
        "wfe_white_test_oneway", "ginv_XX_hat", "X_hat")

    tilde_asym  = max(abs(vec(ginv_XX_tilde - ginv_XX_tilde')))
    tilde_scale = max((1, max(abs(vec(ginv_XX_tilde)))))
    if (tilde_asym > sqrt(epsilon(1)) * tilde_scale) {
        errprintf("wfe_white_test_oneway: ginv_XX_tilde must be symmetric\n")
        _error(3200)
    }
    ginv_XX_tilde_use = 0.5 :* (ginv_XX_tilde + ginv_XX_tilde')
    _wfe_white_require_ginv_match(X_tilde_use, ginv_XX_tilde_use,
        "wfe_white_test_oneway", "ginv_XX_tilde", "X_tilde")

    psi_wfe_asym  = max(abs(vec(Psi_wfe - Psi_wfe')))
    psi_wfe_scale = max((1, max(abs(vec(Psi_wfe)))))
    if (psi_wfe_asym > sqrt(epsilon(1)) * psi_wfe_scale) {
        errprintf("wfe_white_test_oneway: Psi_wfe must be symmetric\n")
        _error(3200)
    }
    Psi_wfe_use = 0.5 :* (Psi_wfe + Psi_wfe')

    psi_fe_asym  = max(abs(vec(Psi_fe - Psi_fe')))
    psi_fe_scale = max((1, max(abs(vec(Psi_fe)))))
    if (psi_fe_asym > sqrt(epsilon(1)) * psi_fe_scale) {
        errprintf("wfe_white_test_oneway: Psi_fe must be symmetric\n")
        _error(3200)
    }
    Psi_fe_use = 0.5 :* (Psi_fe + Psi_fe')

    // unit_idx is optional, but a malformed vector must fail fast rather than
    // silently change the White statistic formula.
    if (args() >= 15) {
        real matrix Omega_cross, Lambda12
        unit_col_use = colshape(unit_idx[., .], 1)
        if (rows(unit_col_use) != NT) {
            errprintf("wfe_white_test_oneway: unit_idx length must match rows(X_hat)\n")
            _error(3200)
        }
        unit_col_use = _wfe_validate_se_unit_idx(unit_col_use, J_u, "wfe_white_test_oneway")
        Omega_cross = _wfe_omega_hac_cross(X_tilde_use, u_tilde_col, X_hat_use, u_hat_col,
                                           unit_col_use, J_u)
        Lambda12 = J_u * ginv_XX_tilde_use * Omega_cross * ginv_XX_hat_use
        Phi = Psi_wfe_use + Psi_fe_use - Lambda12 - Lambda12'
    }
    else {
        // Fallback: element-wise residual cross-product (no external deps)
        diag_ee = u_hat_col :* u_tilde_col
        denom   = NT - J_u - p
        if (denom <= 0) {
            errprintf("wfe_white_test_oneway: NT - J_u - p must be positive\n")
            _error(498)
        }
        Lambda1 = (1 / denom) * cross(X_hat_use :* diag_ee, X_tilde_use)
        Lambda2 = (1 / denom) * cross(X_tilde_use :* diag_ee, X_hat_use)
        Phi = Psi_wfe_use + Psi_fe_use ///
            - (NT * ginv_XX_hat_use) * Lambda1 * (NT * ginv_XX_tilde_use) ///
            - (NT * ginv_XX_tilde_use) * Lambda2 * (NT * ginv_XX_hat_use)
    }

    Phi = _wfe_white_phi_guard(Phi, p, "wfe_white_test_oneway")
    if (rows(Phi) == 0) return((., 1, 0))

    beta_diff  = beta_fe_col - beta_wfe_col
    white_stat = NT * (beta_diff' * pinv(Phi) * beta_diff)
    if (rows(white_stat) > 1 | cols(white_stat) > 1) white_stat = white_stat[1,1]

    white_tol = sqrt(epsilon(1)) * max((1, abs(white_stat)))
    if (white_stat < 0) {
        if (abs(white_stat) <= white_tol) {
            white_stat = 0
        }
        else {
            errprintf("Warning: wfe_white_test_oneway: White statistic is negative (%g); returning p=1.\n", white_stat)
            return((white_stat, 1, 0))
        }
    }

    white_pvalue = chi2tail(p, white_stat)
    test_reject  = (white_pvalue < white_alpha)

    return((white_stat, white_pvalue, test_reject))
}


// ============================================================
// wfe_white_test_twoway - Two-way FE White misspecification test
//
// All-real implementation using unit-clustered cross-HAC Lambda, which
// guarantees a PSD Phi matrix and consistent normalization with Psi_wfe/Psi_fe.
//
// Algorithm:
//   Phi = Psi_wfe + Psi_fe - Lambda12 - Lambda21
//   where:
//     Lambda12 = df_white * ginv_XX_w * Omega_cross * ginv_hat
//     Lambda21 = Lambda12'
//     Omega_cross = _wfe_omega_hac_cross(X_w, u_w, X_hat, u_hat, unit, N_u)
//                 = sum_i (sum_{t:i} X_w[t] u_w[t]) * (sum_{t:i} X_hat[t] u_hat[t])'
//
//   df_white = (NT - p) / (Mstar - n_nonzero_units - n_nonzero_times - p)
//   white_stat = (beta_fe - beta_wfe)' * pinv(Phi) * (beta_fe - beta_wfe)
//
// The Phi matrix is guaranteed PSD because the 2p x 2p joint HAC meat
//   [Omega_wfe, Omega_cross; Omega_cross', Omega_fe]
// is PSD by construction (outer product of unit-level scores).
//
// Parameters:
//   beta_fe      p x 1  Standard FE coefficient vector
//   beta_wfe     p x 1  WFE coefficient vector
//   X_hat        NT x p FE-demeaned X (standard FE)
//   u_hat        NT x 1 Standard FE residuals
//   ginv_XX_hat  p x p  pinv(X_hat' X_hat)
//   X_w          NT x p Weighted FWL X = X_dm .* sqrt(|W|)
//   u_w          NT x 1 Weighted WFE residuals = e_dm .* sqrt(|W|) .* sign(W)
//   ginv_XX_w    p x p  pinv(X_w' X_w)
//   Psi_hat_wfe  p x p  HAC Psi_wfe (from wfe_se_hac_twoway_for_white)
//   Psi_hat_fe   p x p  HAC Psi_fe  (from wfe_se_fe_twoway)
//   unit_idx     NT x 1 Unit indices (1..N_units)
//   N_units      scalar Number of units
//   white_alpha  scalar Significance level
//   W_vec        optional NT x 1 observation-level weight vector
//   time_idx     optional NT x 1 time indices (1..N_times)
//   N_times      optional scalar number of time periods
// ============================================================

real rowvector wfe_white_test_twoway(
    real vector   beta_fe,
    real vector   beta_wfe,
    real matrix   X_hat,
    real vector   u_hat,
    real matrix   ginv_XX_hat,
    real matrix   X_w,
    real vector   u_w,
    real matrix   ginv_XX_w,
    real matrix   Psi_hat_wfe,
    real matrix   Psi_hat_fe,
    real vector   unit_idx,
    real scalar   N_units,
    real scalar   white_alpha,
    | real vector W_vec,
      real vector time_idx,
      real scalar N_times)
{
    real scalar NT, p, hat_asym, hat_scale, w_asym, w_scale, k
    real scalar psi_wfe_asym, psi_wfe_scale, psi_fe_asym, psi_fe_scale
    real scalar df_white, df_white_denom, Mstar, support_units, support_times
    real scalar zero_scale, zero_tol
    real colvector beta_fe_col, beta_wfe_col, u_hat_col, u_w_col, beta_diff
    real colvector unit_col, W_col, time_col, unit_support_count, time_support_count
    real matrix X_hat_use, X_w_use
    real matrix ginv_XX_hat_use, ginv_XX_w_use
    real matrix Psi_hat_wfe_use, Psi_hat_fe_use
    real matrix Omega_cross, Lambda12, Phi_hat, cell_count
    real scalar white_stat, white_pvalue, test_reject, white_tol

    beta_fe_col  = colshape(beta_fe[., .], 1)
    beta_wfe_col = colshape(beta_wfe[., .], 1)
    u_hat_col    = colshape(u_hat[., .], 1)
    u_w_col      = colshape(u_w[., .], 1)
    unit_col     = colshape(unit_idx[., .], 1)
    X_hat_use    = X_hat
    X_w_use      = X_w

    if (rows(X_hat_use) == 1 & cols(X_hat_use) == rows(u_hat_col) & rows(u_hat_col) > 1) {
        X_hat_use = X_hat_use'
    }
    if (rows(X_w_use) == 1 & cols(X_w_use) == rows(u_w_col) & rows(u_w_col) > 1) {
        X_w_use = X_w_use'
    }

    if (args() != 13 & args() != 16) {
        errprintf("wfe_white_test_twoway: supply W_vec, time_idx, and N_times together\n")
        _error(3200)
    }

    NT = rows(X_hat_use)
    p  = cols(X_hat_use)

    if (rows(beta_fe_col) != rows(beta_wfe_col) | cols(beta_fe_col) != 1 | cols(beta_wfe_col) != 1) {
        errprintf("wfe_white_test_twoway: coefficient dimension mismatch\n")
        _error(3200)
    }
    if (NT <= 0) {
        errprintf("wfe_white_test_twoway: X_hat must be nonempty\n")
        _error(3200)
    }
    if (p < 1) {
        errprintf("wfe_white_test_twoway: X_hat must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(beta_fe_col) != p) {
        errprintf("wfe_white_test_twoway: coefficient length must match cols(X_hat)\n")
        _error(3200)
    }
    if (NT != rows(u_hat_col) | NT != rows(u_w_col) | NT != rows(X_w_use) | ///
        NT != rows(unit_col)) {
        errprintf("wfe_white_test_twoway: input length mismatch\n")
        _error(3200)
    }
    if (cols(X_w_use) != p) {
        errprintf("wfe_white_test_twoway: X_w must have the same column count as X_hat\n")
        _error(3200)
    }
    if (N_units != floor(N_units) | N_units <= 1) {
        errprintf("wfe_white_test_twoway: N_units must be an integer >= 2\n")
        _error(3200)
    }
    if (white_alpha <= 0 | white_alpha >= 1) {
        errprintf("wfe_white_test_twoway: white_alpha must lie in (0,1)\n")
        _error(3200)
    }
    if (any(beta_fe_col :>= .) | any(beta_wfe_col :>= .) | ///
        any(u_hat_col :>= .) | any(u_w_col :>= .) | any(unit_col :>= .) | ///
        any(vec(X_hat_use) :>= .) | any(vec(X_w_use) :>= .) | ///
        any(vec(ginv_XX_hat) :>= .) | any(vec(ginv_XX_w) :>= .) | ///
        any(vec(Psi_hat_wfe) :>= .) | any(vec(Psi_hat_fe) :>= .)) {
        errprintf("wfe_white_test_twoway: inputs must not contain missing values\n")
        _error(3498)
    }
    unit_col = _wfe_white_compact_labels(unit_col, N_units,
        "wfe_white_test_twoway", "unit_idx")
    if (rank(X_hat_use) < p | rank(X_w_use) < p) {
        errprintf("wfe_white_test_twoway: X_hat and X_w must each have full column rank\n")
        _error(498)
    }

    if (rows(ginv_XX_hat) != p | cols(ginv_XX_hat) != p) {
        errprintf("wfe_white_test_twoway: ginv_XX_hat dimension mismatch\n")
        _error(3200)
    }
    hat_asym  = max(abs(vec(ginv_XX_hat - ginv_XX_hat')))
    hat_scale = max((1, max(abs(vec(ginv_XX_hat)))))
    if (hat_asym > sqrt(epsilon(1)) * hat_scale) {
        errprintf("wfe_white_test_twoway: ginv_XX_hat must be symmetric\n")
        _error(3200)
    }
    ginv_XX_hat_use = 0.5 :* (ginv_XX_hat + ginv_XX_hat')
    _wfe_white_psd_mat(ginv_XX_hat_use, "wfe_white_test_twoway", "ginv_XX_hat")
    _wfe_white_require_ginv_match(X_hat_use, ginv_XX_hat_use,
        "wfe_white_test_twoway", "ginv_XX_hat", "X_hat")

    if (rows(ginv_XX_w) != p | cols(ginv_XX_w) != p) {
        errprintf("wfe_white_test_twoway: ginv_XX_w dimension mismatch\n")
        _error(3200)
    }
    w_asym  = max(abs(vec(ginv_XX_w - ginv_XX_w')))
    w_scale = max((1, max(abs(vec(ginv_XX_w)))))
    if (w_asym > sqrt(epsilon(1)) * w_scale) {
        errprintf("wfe_white_test_twoway: ginv_XX_w must be symmetric\n")
        _error(3200)
    }
    ginv_XX_w_use = 0.5 :* (ginv_XX_w + ginv_XX_w')
    _wfe_white_psd_mat(ginv_XX_w_use, "wfe_white_test_twoway", "ginv_XX_w")
    _wfe_white_require_ginv_match(X_w_use, ginv_XX_w_use,
        "wfe_white_test_twoway", "ginv_XX_w", "X_w")

    if (rows(Psi_hat_wfe) != p | cols(Psi_hat_wfe) != p) {
        errprintf("wfe_white_test_twoway: Psi_hat_wfe dimension mismatch\n")
        _error(3200)
    }
    psi_wfe_asym  = max(abs(vec(Psi_hat_wfe - Psi_hat_wfe')))
    psi_wfe_scale = max((1, max(abs(vec(Psi_hat_wfe)))))
    if (psi_wfe_asym > sqrt(epsilon(1)) * psi_wfe_scale) {
        errprintf("wfe_white_test_twoway: Psi_hat_wfe must be symmetric\n")
        _error(3200)
    }
    Psi_hat_wfe_use = 0.5 :* (Psi_hat_wfe + Psi_hat_wfe')
    _wfe_white_psd_mat(Psi_hat_wfe_use, "wfe_white_test_twoway", "Psi_hat_wfe")

    if (rows(Psi_hat_fe) != p | cols(Psi_hat_fe) != p) {
        errprintf("wfe_white_test_twoway: Psi_hat_fe dimension mismatch\n")
        _error(3200)
    }
    psi_fe_asym  = max(abs(vec(Psi_hat_fe - Psi_hat_fe')))
    psi_fe_scale = max((1, max(abs(vec(Psi_hat_fe)))))
    if (psi_fe_asym > sqrt(epsilon(1)) * psi_fe_scale) {
        errprintf("wfe_white_test_twoway: Psi_hat_fe must be symmetric\n")
        _error(3200)
    }
    Psi_hat_fe_use = 0.5 :* (Psi_hat_fe + Psi_hat_fe')
    _wfe_white_psd_mat(Psi_hat_fe_use, "wfe_white_test_twoway", "Psi_hat_fe")

    if (args() == 16) {
        W_col = colshape(W_vec[., .], 1)
        time_col = colshape(time_idx[., .], 1)
        if (rows(W_col) != NT | rows(time_col) != NT) {
            errprintf("wfe_white_test_twoway: W_vec/time_idx length mismatch\n")
            _error(3200)
        }
        if (N_times != floor(N_times) | N_times <= 1) {
            errprintf("wfe_white_test_twoway: N_times must be an integer >= 2 when W_vec/time_idx are provided\n")
            _error(3200)
        }
        if (any(W_col :>= .) | any(time_col :>= .)) {
            errprintf("wfe_white_test_twoway: W_vec/time_idx must not contain missing values\n")
            _error(3498)
        }
        time_col = _wfe_white_compact_labels(time_col, N_times,
            "wfe_white_test_twoway", "time_idx")
        cell_count = J(N_units, N_times, 0)
        for (k = 1; k <= NT; k++) {
            cell_count[unit_col[k], time_col[k]] = cell_count[unit_col[k], time_col[k]] + 1
        }
        if (max(cell_count) > 1) {
            errprintf("wfe_white_test_twoway: unit-time pair is not unique\n")
            _error(498)
        }
        zero_scale = max((1, max(abs(vec(X_w_use))), max(abs(u_w_col))))
        zero_tol = sqrt(epsilon(1)) * zero_scale
        for (k = 1; k <= NT; k++) {
            if (W_col[k] != 0) continue
            if (max(abs(X_w_use[k, .])) > zero_tol | abs(u_w_col[k]) > zero_tol) {
                errprintf("wfe_white_test_twoway: zero-weight rows in W_vec must map to zero rows in X_w and u_w\n")
                _error(3498)
            }
        }
        Mstar = sum(W_col :!= 0)
        if (Mstar <= 0) {
            errprintf("wfe_white_test_twoway: at least one non-zero weight is required\n")
            _error(498)
        }

        unit_support_count = J(N_units, 1, 0)
        time_support_count = J(N_times, 1, 0)
        support_units = 0
        support_times = 0
        for (k = 1; k <= rows(W_col); k++) {
            if (W_col[k] == 0) continue
            unit_support_count[unit_col[k]] = 1
            time_support_count[time_col[k]] = 1
        }
        support_units = sum(unit_support_count)
        support_times = sum(time_support_count)
        df_white_denom = Mstar - support_units - support_times - p
        if (df_white_denom <= 0) {
            errprintf("wfe_white_test_twoway: df.white requires Mstar - N_nonzero_units - N_nonzero_times - p > 0\n")
            _error(3351)
        }
        df_white = (rows(X_hat_use) - p) / df_white_denom
    }

    // Cross-HAC meat: Omega_cross = sum_i s_wfe_i * s_fe_i'
    // where s_*_i = sum_{t:i} X_*[t] u_*[t]  (unit-level scores)
    Omega_cross = _wfe_omega_hac_cross(X_w_use, u_w_col, X_hat_use, u_hat_col,
                                       unit_col, N_units)

    if (args() == 16) {
        Lambda12 = df_white * ginv_XX_w_use * Omega_cross * ginv_XX_hat_use
    }
    else {
        // Backward-compatible fallback when retained-support metadata is unavailable.
        Lambda12 = N_units * ginv_XX_w_use * Omega_cross * ginv_XX_hat_use
    }

    // Phi = Psi_wfe + Psi_fe - Lambda12 - Lambda12'
    Phi_hat = _wfe_white_phi_guard(
        Psi_hat_wfe_use + Psi_hat_fe_use - Lambda12 - Lambda12',
        p, "wfe_white_test_twoway")
    if (rows(Phi_hat) == 0) return((., 1, 0))

    beta_diff  = beta_fe_col - beta_wfe_col
    // R twoway White uses the covariance-difference normalization in Phi_hat
    // directly and does not apply an extra retained-row multiplier here.
    white_stat = (beta_diff' * pinv(Phi_hat) * beta_diff)
    if (rows(white_stat) > 1 | cols(white_stat) > 1) white_stat = white_stat[1,1]

    white_tol = sqrt(epsilon(1)) * max((1, abs(white_stat)))
    if (white_stat < 0) {
        if (abs(white_stat) <= white_tol) {
            white_stat = 0
        }
        else {
            errprintf("Warning: wfe_white_test_twoway: White statistic is negative (%g); returning p=1.\n", white_stat)
            return((white_stat, 1, 0))
        }
    }
    white_pvalue = chi2tail(p, white_stat)
    test_reject  = (white_pvalue < white_alpha)

    return((white_stat, white_pvalue, test_reject))
}

// ============================================================
// wfe_white_test_pwfe - PWFE White misspecification test
//
// Key differences from the one-way helper:
//   1. Lambda denominators use NT with no J_u or df adjustment
//   2. u_tilde is already the weighted residual sqrt(W) * resid(fit.final)
//   3. The path stays in the real domain; Re() is applied for type consistency
// ============================================================
real rowvector wfe_white_test_pwfe(
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
    real scalar    white_alpha)
{
    real scalar NT, p, psi_wfe_asym, psi_wfe_scale, psi_fe_asym, psi_fe_scale
    real scalar hat_asym, hat_scale, tilde_asym, tilde_scale
    real colvector beta_wfe_col, beta_fe_col, u_hat_col, u_tilde_col
    real colvector diag_ee, beta_diff
    real matrix X_hat_use, X_tilde_use
    real matrix Lambda1, Lambda2, Phi, Psi_wfe_use, Psi_fe_use
    real matrix ginv_XX_hat_use, ginv_XX_tilde_use
    real scalar white_stat, white_pvalue, test_reject, white_tol

    beta_wfe_col = colshape(beta_wfe[., .], 1)
    beta_fe_col = colshape(beta_fe[., .], 1)
    u_hat_col = colshape(u_hat[., .], 1)
    u_tilde_col = colshape(u_tilde[., .], 1)
    X_hat_use = X_hat
    X_tilde_use = X_tilde

    if (!(rows(X_hat_use) == 1 & cols(X_hat_use) > 1 ///
          & rows(u_hat_col) == 1 & rows(u_tilde_col) == 1)) {
        if (rows(X_hat_use) == 1 & cols(X_hat_use) == rows(u_hat_col) & rows(u_hat_col) > 1) {
            X_hat_use = X_hat_use'
        }
    }
    if (!(rows(X_tilde_use) == 1 & cols(X_tilde_use) > 1 ///
          & rows(u_hat_col) == 1 & rows(u_tilde_col) == 1)) {
        if (rows(X_tilde_use) == 1 & cols(X_tilde_use) == rows(u_tilde_col) & rows(u_tilde_col) > 1) {
            X_tilde_use = X_tilde_use'
        }
    }

    NT = rows(X_hat_use)
    p = cols(X_hat_use)

    if (rows(beta_wfe_col) != rows(beta_fe_col) | cols(beta_wfe_col) != 1 | cols(beta_fe_col) != 1) {
        errprintf("wfe_white_test_pwfe: coefficient dimension mismatch\n")
        _error(3200)
    }
    if (p != floor(p) | p <= 0) {
        errprintf("wfe_white_test_pwfe: p must be a positive integer\n")
        _error(3200)
    }
    if (any(beta_wfe_col :>= .) | any(beta_fe_col :>= .) | any(u_hat_col :>= .) | any(u_tilde_col :>= .) | ///
        any(vec(Psi_wfe) :>= .) | any(vec(Psi_fe) :>= .) | ///
        any(vec(ginv_XX_hat) :>= .) | any(vec(ginv_XX_tilde) :>= .) | ///
        any(vec(X_hat_use) :>= .) | any(vec(X_tilde_use) :>= .)) {
        errprintf("wfe_white_test_pwfe: inputs must not contain missing values\n")
        _error(3498)
    }
    if (NT <= 0 | rows(X_tilde_use) <= 0) {
        errprintf("wfe_white_test_pwfe: X_hat and X_tilde must be nonempty\n")
        _error(3200)
    }
    if (rows(beta_fe_col) != p) {
        errprintf("wfe_white_test_pwfe: coefficient length must match cols(X_hat)\n")
        _error(3200)
    }
    if (cols(X_tilde_use) != p) {
        errprintf("wfe_white_test_pwfe: X_hat and X_tilde must have the same column count\n")
        _error(3200)
    }
    if (NT != rows(u_hat_col) | rows(X_tilde_use) != rows(u_tilde_col)) {
        errprintf("wfe_white_test_pwfe: residual length mismatch\n")
        _error(3200)
    }
    if (NT != rows(X_tilde_use)) {
        errprintf("wfe_white_test_pwfe: White path requires rows(X_hat) == rows(X_tilde)\n")
        _error(3200)
    }
    if (rows(ginv_XX_hat) != p | cols(ginv_XX_hat) != p) {
        errprintf("wfe_white_test_pwfe: ginv_XX_hat dimension mismatch\n")
        _error(3200)
    }
    hat_asym = max(abs(vec(ginv_XX_hat - ginv_XX_hat')))
    hat_scale = max((1, max(abs(vec(ginv_XX_hat)))))
    if (hat_asym > sqrt(epsilon(1)) * hat_scale) {
        errprintf("wfe_white_test_pwfe: ginv_XX_hat must be symmetric\n")
        _error(3200)
    }
    ginv_XX_hat_use = 0.5 :* (ginv_XX_hat + ginv_XX_hat')
    _wfe_white_psd_mat(ginv_XX_hat_use, "wfe_white_test_pwfe", "ginv_XX_hat")
    _wfe_white_require_ginv_match(X_hat_use, ginv_XX_hat_use,
        "wfe_white_test_pwfe", "ginv_XX_hat", "X_hat")
    if (rows(ginv_XX_tilde) != p | cols(ginv_XX_tilde) != p) {
        errprintf("wfe_white_test_pwfe: ginv_XX_tilde dimension mismatch\n")
        _error(3200)
    }
    tilde_asym = max(abs(vec(ginv_XX_tilde - ginv_XX_tilde')))
    tilde_scale = max((1, max(abs(vec(ginv_XX_tilde)))))
    if (tilde_asym > sqrt(epsilon(1)) * tilde_scale) {
        errprintf("wfe_white_test_pwfe: ginv_XX_tilde must be symmetric\n")
        _error(3200)
    }
    ginv_XX_tilde_use = 0.5 :* (ginv_XX_tilde + ginv_XX_tilde')
    _wfe_white_psd_mat(ginv_XX_tilde_use, "wfe_white_test_pwfe", "ginv_XX_tilde")
    _wfe_white_require_ginv_match(X_tilde_use, ginv_XX_tilde_use,
        "wfe_white_test_pwfe", "ginv_XX_tilde", "X_tilde")
    if (rows(Psi_wfe) != p | cols(Psi_wfe) != p) {
        errprintf("wfe_white_test_pwfe: Psi_wfe dimension mismatch\n")
        _error(3200)
    }
    psi_wfe_asym = max(abs(vec(Psi_wfe - Psi_wfe')))
    psi_wfe_scale = max((1, max(abs(vec(Psi_wfe)))))
    if (psi_wfe_asym > sqrt(epsilon(1)) * psi_wfe_scale) {
        errprintf("wfe_white_test_pwfe: Psi_wfe must be symmetric\n")
        _error(3200)
    }
    Psi_wfe_use = 0.5 :* (Psi_wfe + Psi_wfe')
    _wfe_white_psd_mat(Psi_wfe_use, "wfe_white_test_pwfe", "Psi_wfe")
    if (rows(Psi_fe) != p | cols(Psi_fe) != p) {
        errprintf("wfe_white_test_pwfe: Psi_fe dimension mismatch\n")
        _error(3200)
    }
    psi_fe_asym = max(abs(vec(Psi_fe - Psi_fe')))
    psi_fe_scale = max((1, max(abs(vec(Psi_fe)))))
    if (psi_fe_asym > sqrt(epsilon(1)) * psi_fe_scale) {
        errprintf("wfe_white_test_pwfe: Psi_fe must be symmetric\n")
        _error(3200)
    }
    Psi_fe_use = 0.5 :* (Psi_fe + Psi_fe')
    _wfe_white_psd_mat(Psi_fe_use, "wfe_white_test_pwfe", "Psi_fe")
    if (white_alpha <= 0 | white_alpha >= 1) {
        errprintf("wfe_white_test_pwfe: white_alpha must lie in (0,1)\n")
        _error(3200)
    }

    diag_ee = u_hat_col :* u_tilde_col

    Lambda1 = (1 / NT) * cross(X_hat_use :* diag_ee, X_tilde_use)
    Lambda2 = (1 / NT) * cross(X_tilde_use :* diag_ee, X_hat_use)

    Phi = Psi_wfe_use + Psi_fe_use ///
        - (NT * ginv_XX_hat_use) * Lambda1 * (NT * ginv_XX_tilde_use) ///
        - (NT * ginv_XX_tilde_use) * Lambda2 * (NT * ginv_XX_hat_use)
    Phi = _wfe_white_phi_guard(Phi, p, "wfe_white_test_pwfe")
    if (rows(Phi) == 0) return((., 1, 0))

    beta_diff = beta_fe_col - beta_wfe_col
    white_stat = Re(NT * beta_diff' * pinv(Phi) * beta_diff)
    white_tol = sqrt(epsilon(1)) * max((1, abs(white_stat)))
    if (white_stat < 0) {
        if (abs(white_stat) <= white_tol) {
            white_stat = 0
        }
        else {
            errprintf("Warning: wfe_white_test_pwfe: White statistic is negative (%g); returning p=1.\n", white_stat)
            return((white_stat, 1, 0))
        }
    }
    white_pvalue = chi2tail(p, white_stat)
    test_reject = (white_pvalue < white_alpha)

    return((white_stat, white_pvalue, test_reject))
}

end
