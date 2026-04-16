version 16.0
mata:
mata set matastrict on

// Numerically stable inverse logit
// Uses separate branches for positive/negative eta to avoid overflow

real colvector _wfe_pwfe_invlogit(real matrix eta)
{
    real scalar i
    real colvector eta_col, prob
    real scalar exp_eta

    if (rows(eta) != 1 & cols(eta) != 1) {
        errprintf("_wfe_pwfe_invlogit: eta must be a vector\n")
        _error(3200)
    }

    eta_col = colshape(eta[., .], 1)
    if (any(eta_col :>= .)) {
        errprintf("_wfe_pwfe_invlogit: eta must not contain missing values\n")
        _error(3498)
    }
    prob = J(rows(eta_col), 1, .)

    for (i = 1; i <= rows(eta_col); i++) {
        if (eta_col[i] >= 0) {
            prob[i] = 1 / (1 + exp(-eta_col[i]))
        }
        else {
            exp_eta = exp(eta_col[i])
            prob[i] = exp_eta / (1 + exp_eta)
        }
    }

    return(prob)
}


real scalar _wfe_pwfe_logit_nll(
    real matrix    X,
    real vector    y,
    real vector    beta
)
{
    real scalar i, obj
    real colvector y_col, beta_col, eta

    y_col = colshape(y[., .], 1)
    beta_col = colshape(beta[., .], 1)
    if (rows(X) == 1 & cols(X) == rows(y_col) & rows(y_col) > 1) {
        // Canonicalize a single-regressor observation stream to column format
        X = colshape(X[., .], 1)
    }

    if (rows(X) == 0) {
        errprintf("_wfe_pwfe_logit_nll: X must contain at least one observation\n")
        _error(3200)
    }
    if (cols(X) == 0) {
        errprintf("_wfe_pwfe_logit_nll: X must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(X) != rows(y_col)) {
        errprintf("_wfe_pwfe_logit_nll: X and y must align by row\n")
        _error(3200)
    }
    if (cols(X) != rows(beta_col)) {
        errprintf("_wfe_pwfe_logit_nll: beta length must equal cols(X)\n")
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("_wfe_pwfe_logit_nll: X must not contain missing values\n")
        _error(3498)
    }
    if (any(y_col :>= .)) {
        errprintf("_wfe_pwfe_logit_nll: y must not contain missing values\n")
        _error(3498)
    }
    if (any((y_col :!= 0) :& (y_col :!= 1))) {
        errprintf("_wfe_pwfe_logit_nll: y must contain only 0/1 values\n")
        _error(3498)
    }
    if (any(beta_col :>= .)) {
        errprintf("_wfe_pwfe_logit_nll: beta must not contain missing values\n")
        _error(3498)
    }

    eta = X * beta_col
    obj = 0
    for (i = 1; i <= rows(eta); i++) {
        if (eta[i] >= 0) {
            obj = obj + (1 - y_col[i]) * eta[i] + ln(1 + exp(-eta[i]))
        }
        else {
            obj = obj - y_col[i] * eta[i] + ln(1 + exp(eta[i]))
        }
    }

    return(obj)
}


real scalar _wfe_pwfe_unique_count(real vector x)
{
    real scalar i, n_unique
    real colvector x_sorted

    x = colshape(x[., .], 1)
    if (rows(x) == 0) {
        return(0)
    }
    if (any(x :>= .)) {
        errprintf("_wfe_pwfe_unique_count: x must not contain missing values\n")
        _error(3498)
    }

    x_sorted = sort(x, 1)
    n_unique = 1
    for (i = 2; i <= rows(x_sorted); i++) {
        if (x_sorted[i] != x_sorted[i - 1]) {
            n_unique = n_unique + 1
        }
    }

    return(n_unique)
}


real colvector _wfe_pwfe_prior_scale_core(
    real matrix    X,
    real scalar    prior_scale,
    real scalar    intercept_only,
    real scalar    rhs_terms,
    string scalar  caller
)
{
    real scalar j, n_cat, x_scale
    real colvector scales, xj

    if (rows(X) == 0) {
        errprintf("%s: X must contain at least one observation\n", caller)
        _error(3200)
    }
    if (cols(X) == 0) {
        errprintf("%s: X must contain at least one regressor\n", caller)
        _error(3200)
    }
    if (prior_scale >= . | prior_scale <= 0) {
        errprintf("%s: prior_scale must be positive\n", caller)
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("%s: X must not contain missing values\n", caller)
        _error(3498)
    }
    if (rhs_terms >= . | rhs_terms != floor(rhs_terms) | rhs_terms < 0) {
        errprintf("%s: rhs_terms must be a nonnegative integer\n", caller)
        _error(3200)
    }

    scales = J(cols(X), 1, prior_scale)

    /*
       This helper must not inspect ambient Stata locals to infer matrix
       geometry. Ambiguous 1 x k inputs require explicit rhs context from the
       caller; otherwise fail fast instead of silently depending on hidden
       state.
    */
    if (rows(X) == 1 & cols(X) > 1) {
        if (rhs_terms == 0) {
            errprintf("%s: ambiguous 1 x k rowvector requires explicit rhs context\n",
                      caller)
            _error(3200)
        }
        if (rhs_terms == 1) {
            X = X'
        }
    }

    scales = J(cols(X), 1, prior_scale)
    if (cols(X) == 1 & rows(X) > 0 & min(X) == 1 & max(X) == 1 ///
        & intercept_only) {
        scales[1] = 4 * prior_scale
        return(scales)
    }

    for (j = 1; j <= cols(X); j++) {
        xj = X[., j]
        n_cat = _wfe_pwfe_unique_count(xj)
        x_scale = 1

        if (n_cat == 2) {
            x_scale = max(xj) - min(xj)
        }
        else if (n_cat > 2) {
            x_scale = 2 * sqrt(variance(xj))
        }

        scales[j] = scales[j] / x_scale
        if (scales[j] < 1e-12) {
            scales[j] = 1e-12
        }
    }

    return(scales)
}


real colvector _wfe_pwfe_prior_scale(
    real matrix    X,
    real scalar    prior_scale
)
{
    return(_wfe_pwfe_prior_scale_core(X, prior_scale, 0, 0,
                                      "_wfe_pwfe_prior_scale"))
}


real colvector _wfe_pwfe_prior_scale_rhs(
    real matrix    X,
    real scalar    prior_scale,
    real scalar    rhs_terms
)
{
    return(_wfe_pwfe_prior_scale_core(X, prior_scale, 0, rhs_terms,
                                      "_wfe_pwfe_prior_scale_rhs"))
}


real colvector _wfe_pwfe_penlogit_prob_core(
    real matrix X,
    real vector y,
    real scalar prior_scale,
    real scalar intercept_only,
    string scalar caller
)
{
    real scalar n, p, iter, max_iter, conv, dev, devold, dev_scale
    real scalar prior_sd_scale
    real colvector y_col, beta, eta, prob, mu_eta, z, prior_scale_vec
    real colvector prior_sd, prior_sd_new, z_star, w_star, beta_new
    real matrix A, V_beta, x_star

    y_col = colshape(y[., .], 1)
    n = rows(y_col)
    if (rows(X) == 1 & cols(X) == n & n > 1) {
        // Canonicalize 1-D rowvector to column vector
        X = colshape(X[., .], 1)
    }
    p = cols(X)

    if (n == 0 | rows(X) == 0) {
        errprintf("%s: X and y must contain at least one observation\n", caller)
        _error(3200)
    }
    if (rows(X) != n) {
        errprintf("%s: X and y must align by row\n", caller)
        _error(3200)
    }
    if (p == 0) {
        errprintf("%s: X must contain at least one regressor\n", caller)
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("%s: X must not contain missing values\n", caller)
        _error(3498)
    }
    if (any(y_col :>= .)) {
        errprintf("%s: y must not contain missing values\n", caller)
        _error(3498)
    }
    if (any((y_col :!= 0) :& (y_col :!= 1))) {
        errprintf("%s: y must contain only 0/1 values\n", caller)
        _error(3498)
    }
    // All-same y (all-0 or all-1) is allowed: Cauchy prior regularization
    // produces finite probabilities near 0 or 1.
    if (prior_scale >= . | prior_scale <= 0) {
        errprintf("%s: prior_scale must be positive\n", caller)
        _error(3200)
    }

    // Penalized IRLS for logistic regression with data-dependent prior
    beta = J(p, 1, 0)
    prior_scale_vec = _wfe_pwfe_prior_scale_core(X, prior_scale,
                                                 intercept_only, cols(X),
                                                 caller)
    prior_sd = prior_scale_vec
    max_iter = 100
    conv = 0
    devold = 2 * _wfe_pwfe_logit_nll(X, y_col, beta)

    for (iter = 1; iter <= max_iter; iter++) {
        eta = X * beta
        prob = _wfe_pwfe_invlogit(eta)
        mu_eta = prob :* (1 :- prob)

        // Clamp mu_eta away from 0 to avoid IRLS singularity
        mu_eta = rowmax((mu_eta, J(rows(mu_eta), 1, 1e-15)))

        z = eta + (y_col - prob) :/ mu_eta
        x_star = X \ I(p)
        z_star = z \ J(p, 1, 0)
        w_star = sqrt(mu_eta) \ (1 :/ prior_sd)
        A = quadcross(x_star, x_star :* (w_star :^ 2))

        V_beta = invsym(A)
        beta_new = V_beta * quadcross(x_star, z_star :* (w_star :^ 2))
        if (any(beta_new :>= .) | any(diagonal(V_beta) :>= .)) {
            errprintf("%s: penalized IRLS produced non-finite coefficients\n",
                      caller)
            _error(498)
        }

        prior_sd_new = sqrt((beta_new :^ 2 :+ diagonal(V_beta) :+ prior_scale_vec :^ 2) :/ 2)
        dev = 2 * _wfe_pwfe_logit_nll(X, y_col, beta_new)
        dev_scale = 0.1 + abs(dev)
        prior_sd_scale = 0.1 + max(abs(prior_sd_new))

        if (iter > 1 ///
            & abs(dev - devold) / dev_scale < 1e-8 ///
            & max(abs(prior_sd_new - prior_sd)) / prior_sd_scale < 1e-8) {
            beta = beta_new
            conv = 1
            break
        }

        beta = beta_new
        prior_sd = prior_sd_new
        devold = dev
    }

    if (!conv) {
        errprintf("%s: penalized IRLS failed to converge\n", caller)
        _error(498)
    }

    prob = _wfe_pwfe_invlogit(X * beta)
    // Clamp boundary probabilities to epsilon range
    prob = rowmax((prob, J(rows(prob), 1, 1e-15)))
    prob = rowmin((prob, J(rows(prob), 1, 1 - 1e-15)))

    return(prob)
}


real colvector _wfe_pwfe_penlogit_prob(
    real matrix X,
    real vector y,
    real scalar prior_scale
)
{
    return(_wfe_pwfe_penlogit_prob_core(X, y, prior_scale, 0,
                                        "_wfe_pwfe_penlogit_prob"))
}


real colvector _wfe_pwfe_penlogit_prob_int(
    real matrix X,
    real vector y,
    real scalar prior_scale
)
{
    return(_wfe_pwfe_penlogit_prob_core(X, y, prior_scale, 1,
                                        "_wfe_pwfe_penlogit_prob_int"))
}


real scalar _wfe_pwfe_bridge_confirm_rc(string scalar cmd)
{
    real scalar rc

    rc = _stata("capture " + cmd)
    rc = st_numscalar("c(rc)")

    return(rc)
}


void _wfe_pwfe_bridge_singlevar_core(
    string scalar varspec,
    string scalar argname,
    string scalar caller
)
{
    if (strtrim(varspec) == "" | cols(tokens(strtrim(varspec))) != 1) {
        errprintf("%s: %s must name exactly one existing variable\n",
                  caller, argname)
        _error(3200)
    }

    if (_wfe_pwfe_bridge_confirm_rc("confirm variable " + varspec) != 0) {
        errprintf("%s: %s must name exactly one existing variable\n",
                  caller, argname)
        _error(3200)
    }
}


void _wfe_pwfe_bridge_singlevar(
    string scalar varspec,
    string scalar argname
)
{
    _wfe_pwfe_bridge_singlevar_core(varspec, argname,
                                    "_wfe_pwfe_bridge_singlevar")
}


void _wfe_pwfe_bridge_numvar_core(
    string scalar varspec,
    string scalar argname,
    string scalar caller
)
{
    _wfe_pwfe_bridge_singlevar_core(varspec, argname, caller)

    if (_wfe_pwfe_bridge_confirm_rc("confirm numeric variable " + varspec) != 0) {
        errprintf("%s: %s must name exactly one existing numeric variable\n",
                  caller, argname)
        _error(3200)
    }
}


void _wfe_pwfe_bridge_numvar(
    string scalar varspec,
    string scalar argname
)
{
    _wfe_pwfe_bridge_numvar_core(varspec, argname,
                                 "_wfe_pwfe_bridge_numvar")
}


void _wfe_pwfe_bridge_fpvar_core(
    string scalar varspec,
    string scalar argname,
    string scalar caller
)
{
    string scalar vartype

    _wfe_pwfe_bridge_numvar_core(varspec, argname, caller)
    vartype = st_vartype(st_varindex(varspec))

    if (vartype != "float" & vartype != "double") {
        errprintf("%s: %s must be float or double to store continuous transformed outcomes\n",
                  caller, argname)
        _error(3200)
    }
}


void _wfe_pwfe_bridge_fpvar(
    string scalar varspec,
    string scalar argname
)
{
    _wfe_pwfe_bridge_fpvar_core(varspec, argname,
                                "_wfe_pwfe_bridge_fpvar")
}


void _wfe_pwfe_bridge_tousemask_core(
    string scalar varspec,
    string scalar caller
)
{
    real colvector mask

    _wfe_pwfe_bridge_numvar_core(varspec, "touse_var", caller)
    mask = st_data(., varspec)

    if (any(mask :>= .)) {
        errprintf("%s: touse_var must not contain missing values\n",
                  caller)
        _error(3498)
    }

    if (any((mask :!= 0) :& (mask :!= 1))) {
        errprintf("%s: touse_var must contain only 0/1 values\n",
                  caller)
        _error(3498)
    }
}


void _wfe_pwfe_bridge_tousemask(string scalar varspec)
{
    _wfe_pwfe_bridge_tousemask_core(varspec,
                                    "_wfe_pwfe_bridge_tousemask")
}


void _wfe_pwfe_apply_transform(
    string scalar outcome_var,
    string scalar treat_var,
    string scalar pscore_var,
    string scalar touse_var,
    string scalar out_var)
{
    real colvector Y, treat, pscore, y_star

    _wfe_pwfe_bridge_numvar_core(outcome_var, "outcome_var",
                                 "_wfe_pwfe_apply_transform")
    _wfe_pwfe_bridge_numvar_core(treat_var, "treat_var",
                                 "_wfe_pwfe_apply_transform")
    _wfe_pwfe_bridge_numvar_core(pscore_var, "pscore_var",
                                 "_wfe_pwfe_apply_transform")
    _wfe_pwfe_bridge_tousemask_core(touse_var,
                                    "_wfe_pwfe_apply_transform")
    _wfe_pwfe_bridge_fpvar_core(out_var, "out_var",
                                "_wfe_pwfe_apply_transform")

    st_view(Y, ., outcome_var, touse_var)
    st_view(treat, ., treat_var, touse_var)
    st_view(pscore, ., pscore_var, touse_var)

    if (rows(Y) == 0) {
        errprintf("pwfe transform bridge requires at least one selected observation\n")
        _error(3200)
    }

    if (any(Y :>= .)) {
        errprintf("_wfe_pwfe_apply_transform: outcome must not contain missing values\n")
        _error(3498)
    }

    if (any(treat :>= .)) {
        errprintf("_wfe_pwfe_apply_transform: treat must not contain missing values\n")
        _error(3498)
    }

    if (any((treat :!= 0) :& (treat :!= 1))) {
        errprintf("_wfe_pwfe_apply_transform: treat must contain only 0/1 values\n")
        _error(3498)
    }

    if (any(pscore :>= .)) {
        errprintf("_wfe_pwfe_apply_transform: pscore must not contain missing values\n")
        _error(3498)
    }

    if (any(pscore :<= 0 :| pscore :>= 1)) {
        errprintf("_wfe_pwfe_apply_transform: pscore must contain only values strictly between 0 and 1\n")
        _error(3498)
    }

    y_star = wfe_transform(Y, treat, pscore)
    st_store(., out_var, touse_var, y_star)
}

end
