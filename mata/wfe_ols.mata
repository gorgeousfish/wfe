// wfe_ols.mata — OLS estimation core
// Minimal OLS solver with Moore-Penrose pseudoinverse
// Returns: beta (row vector), residuals, and ginv(XX)

version 16.0
mata:
mata set matastrict on

// OLS result structure: beta (1×p rowvector), residuals (NT×1), ginv(XX) (p×p)
struct wfe_ols_result {
    real rowvector  beta
    real colvector  resid
    real matrix     ginv_XX
}


// _wfe_ols_core() — Core OLS solver using Moore-Penrose pseudoinverse
// Computes: beta = ((X'X)+ X'Y)', residuals = Y - X*beta'
// One-dimensional inputs are canonicalized to column vectors.
//
// @param Y_dm   real vector    Demeaned outcome (NT×1 or 1×NT)
// @param X_dm   real matrix    Demeaned regressors (NT×p)
// @return       struct wfe_ols_result scalar
struct wfe_ols_result scalar _wfe_ols_core(
    real vector    Y_dm,    // Demeaned outcome (NT×1 or 1×NT)
    real matrix    X_dm     // Demeaned regressors (NT×p)
)
{
    struct wfe_ols_result scalar result
    real matrix XX, X_use
    real colvector Y_col
    real scalar p

    // Copy inputs to plain arrays and canonicalize 1D orientation
    Y_col = colshape(Y_dm[., .], 1)
    X_use = X_dm[., .]

    if (rows(X_use) == 1 & cols(X_use) == rows(Y_col) & rows(Y_col) > 1) {
        X_use = X_use'
    }

    // Input validation
    if (rows(Y_col) != rows(X_use)) {
        errprintf("_wfe_ols_core: Y (%g rows) and X (%g rows) dimension mismatch\n",
                  rows(Y_col), rows(X_use))
        _error(3200)
    }
    if (rows(Y_col) < 1) {
        errprintf("_wfe_ols_core: Y and X must contain at least one observation\n")
        _error(3200)
    }
    if (any(Y_col :>= .) | any(X_use :>= .)) {
        errprintf("_wfe_ols_core: Y and X must not contain missing values\n")
        _error(3498)
    }
    p = cols(X_use)
    if (p < 1) {
        errprintf("_wfe_ols_core: X must contain at least one regressor\n")
        _error(3200)
    }

    // Gram matrix X'X
    XX = cross(X_use, X_use)

    // Moore-Penrose pseudoinverse of X'X: rank-deficient designs still return
    // the minimum-norm OLS solution instead of failing at the helper boundary.
    result.ginv_XX = pinv(XX)

    // Coefficients: beta = ((X'X)+ X'Y)' — returned as 1×p row vector
    result.beta = (result.ginv_XX * cross(X_use, Y_col))'

    // Residuals: u = Y - X*beta'
    result.resid = Y_col - X_use * result.beta'

    return(result)
}


// _wfe_calc_sigma2() — Residual variance and degrees of freedom
// df_r = NT - p - J_u, sigma2 = (u'u) / df_r
//
// @param u_tilde  real vector    Residuals (NT×1 or 1×NT)
// @param NT       real scalar    Total observations
// @param p        real scalar    Number of regressors
// @param J_u      real scalar    Number of units
// @param sigma2   real scalar [output] Residual variance
// @param df_r     real scalar [output] Residual degrees of freedom
void _wfe_calc_sigma2(
    real vector    u_tilde,
    real scalar    NT,
    real scalar    p,
    real scalar    J_u,
    real scalar    sigma2,
    real scalar    df_r
)
{
    real colvector u_col

    u_col = colshape(u_tilde[., .], 1)

    if (NT >= . | NT != floor(NT) | NT <= 0 | ///
        p >= . | p != floor(p) | p < 0 | ///
        J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("_wfe_calc_sigma2: NT and J_u must be positive integers; p must be a nonnegative integer\n")
        _error(3200)
    }
    if (p == 0) {
        errprintf("_wfe_calc_sigma2: p must be a positive integer\n")
        _error(3200)
    }
    if (rows(u_col) != NT) {
        errprintf("_wfe_calc_sigma2: rows(u_tilde) must equal NT\n")
        _error(3200)
    }
    if (any(u_col :>= .)) {
        errprintf("_wfe_calc_sigma2: u_tilde must not contain missing values\n")
        _error(3498)
    }

    // Degrees of freedom: df_r = NT - p - J_u
    df_r = NT - p - J_u

    // Check for insufficient degrees of freedom
    if (df_r <= 0) {
        errprintf("_wfe_calc_sigma2: insufficient degrees of freedom: " +
                  "NT=%g, p=%g, J_u=%g, df_r=%g\n",
                  NT, p, J_u, df_r)
        _error(3351)
    }

    // Residual variance: sigma2 = (u'u) / df_r
    sigma2 = cross(u_col, u_col) / df_r
}

end
