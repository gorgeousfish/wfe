// ============================================================
// wfe_bridge.mata — ADO-Mata Data Bridge
//
// Purpose: Transfer data and parameters from Stata to Mata,
//          and post estimation results back to Stata.
//          No preprocessing, no estimation logic.
//
// Compile: set matastrict on; do mata/wfe_bridge.mata
// Dependencies: wfe_utils
// ============================================================

version 16.0
mata:
mata set matastrict on


void _wfe_bridge_chk_panel(
    real scalar N_units,
    real scalar N_times,
    real scalar NT
)
{
    string scalar touse, unitvar, timevar
    real colvector unit_idx, time_idx, seen_units, seen_times
    real matrix cell_count
    real scalar k

    touse   = st_local("touse")
    unitvar = st_local("wfe_unit_idx")
    timevar = st_local("wfe_time_idx")

    if (strtrim(touse) != "") {
        _wfe_b_validate_touse_mask(touse, "_wfe_bridge_read_params")
    }

    st_view(unit_idx, ., unitvar, touse)
    st_view(time_idx, ., timevar, touse)

    if (rows(unit_idx) != NT | rows(time_idx) != NT) {
        errprintf("_wfe_bridge_read_params: wfe_unit_idx/wfe_time_idx rows must equal NT\n")
        _error(3200)
    }

    if (any(unit_idx :>= .) | any(time_idx :>= .)) {
        errprintf("_wfe_bridge_read_params: wfe_unit_idx/wfe_time_idx must not contain missing values\n")
        _error(3498)
    }

    seen_units = J(N_units, 1, 0)
    seen_times = J(N_times, 1, 0)
    cell_count = J(N_times, N_units, 0)

    for (k = 1; k <= NT; k++) {
        if (unit_idx[k] != floor(unit_idx[k]) | time_idx[k] != floor(time_idx[k]) | ///
            unit_idx[k] < 1 | unit_idx[k] > N_units | ///
            time_idx[k] < 1 | time_idx[k] > N_times) {
            errprintf("_wfe_bridge_read_params: N_units/N_times must match observed dense panel support\n")
            _error(3200)
        }
        // Validate panel structure: unique (unit,time) pairs, sorted by (unit,time).
        // The Mata bridge assumes data arrive in (unit_idx, time_idx) order.
        if (k > 1) {
            if (unit_idx[k] < unit_idx[k-1] | ///
                (unit_idx[k] == unit_idx[k-1] & time_idx[k] < time_idx[k-1])) {
                errprintf("_wfe_bridge_read_params: input must be sorted by (unit_idx, time_idx)\n")
                _error(3200)
            }
        }
        cell_count[time_idx[k], unit_idx[k]] = cell_count[time_idx[k], unit_idx[k]] + 1
        if (cell_count[time_idx[k], unit_idx[k]] > 1) {
            errprintf("_wfe_bridge_read_params: unit-time pair is not unique\n")
            _error(498)
        }
        seen_units[unit_idx[k]] = 1
        seen_times[time_idx[k]] = 1
    }

    if (any(seen_units :== 0) | any(seen_times :== 0)) {
        errprintf("_wfe_bridge_read_params: N_units/N_times must match observed dense panel support\n")
        _error(3200)
    }
}


void _wfe_b_singlevar(
    string scalar varspec,
    string scalar localname
)
{
    if (cols(tokens(strtrim(varspec))) != 1) {
        errprintf("_wfe_bridge_construct_views: %s local must name exactly one existing variable\n",
                  localname)
        _error(3200)
    }
}


void _wfe_b_unique_varlist(
    string scalar varspec,
    string scalar localname
)
{
    string rowvector toks
    real scalar i, j

    toks = tokens(strtrim(varspec))
    for (i = 1; i <= cols(toks); i++) {
        for (j = i + 1; j <= cols(toks); j++) {
            if (toks[i] == toks[j]) {
                errprintf("_wfe_bridge_construct_views: %s local must name unique existing regressors\n",
                          localname)
                _error(3200)
            }
        }
    }
}


real scalar _wfe_b_confirm_rc(string scalar cmd)
{
    real scalar rc

    rc = _stata("capture " + cmd)
    rc = st_numscalar("c(rc)")

    return(rc)
}


void _wfe_b_single_numeric_var(
    string scalar varspec,
    string scalar localname
)
{
    _wfe_b_singlevar(varspec, localname)

    if (_wfe_b_confirm_rc("confirm variable " + varspec) != 0) {
        errprintf("_wfe_bridge_construct_views: %s local must name an existing variable\n",
                  localname)
        _error(3200)
    }

    if (_wfe_b_confirm_rc("confirm numeric variable " + varspec) != 0) {
        errprintf("_wfe_bridge_construct_views: %s local must name an existing numeric variable\n",
                  localname)
        _error(3200)
    }
}


void _wfe_b_numeric_varlist(
    string scalar varspec,
    string scalar localname
)
{
    if (_wfe_b_confirm_rc("confirm variable " + varspec) != 0) {
        errprintf("_wfe_bridge_construct_views: %s local must name only existing regressors\n",
                  localname)
        _error(3200)
    }

    if (_wfe_b_confirm_rc("confirm numeric variable " + varspec) != 0) {
        errprintf("_wfe_bridge_construct_views: %s local must name only numeric regressors\n",
                  localname)
        _error(3200)
    }

    _wfe_b_unique_varlist(varspec, localname)
}


void _wfe_b_req_numlocal(
    string scalar varspec,
    string scalar localname,
    string scalar caller
)
{
    if (cols(tokens(strtrim(varspec))) != 1) {
        errprintf("%s: %s must name an existing numeric variable\n",
                  caller, localname)
        _error(3200)
    }

    if (_wfe_b_confirm_rc("confirm variable " + varspec) != 0) {
        errprintf("%s: %s must name an existing numeric variable\n",
                  caller, localname)
        _error(3200)
    }

    if (_wfe_b_confirm_rc("confirm numeric variable " + varspec) != 0) {
        errprintf("%s: %s must name an existing numeric variable\n",
                  caller, localname)
        _error(3200)
    }
}


void _wfe_b_validate_touse_mask(
    string scalar varspec,
    string scalar caller
)
{
    real colvector touse_mask

    _wfe_b_req_numlocal(varspec, "touse", caller)
    st_view(touse_mask, ., varspec)

    if (rows(touse_mask) == 0) {
        return
    }

    if (any(touse_mask :>= .) | any((touse_mask :!= 0) :& (touse_mask :!= 1))) {
        errprintf("%s: touse must contain only 0/1 values and no missing values\n",
                  caller)
        _error(3498)
    }
}


// ============================================================
// _wfe_bridge_construct_views — Construct read-only views from Stata data
//
// Uses st_view() for zero-copy data transfer, filtered by touse.
//
// Output arguments (passed by reference):
//   Y        NT×1 outcome variable view
//   X        NT×p covariate matrix view
//   treat    NT×1 treatment indicator view
//   unit_idx NT×1 unit index view (1..N_units)
//   time_idx NT×1 time index view (1..N_times)
//   cit      NT×1 weight correction (view or unit vector)
//
// Precondition: Data already sorted and recoded (unit_idx, time_idx ready)
// ============================================================
void _wfe_bridge_construct_views(
    real colvector Y,
    real matrix    X,
    real colvector treat,
    real colvector unit_idx,
    real colvector time_idx,
    real colvector cit)
{
    string scalar touse, depvar, indepvars, treatvar, citvar
    string scalar unitvar, timevar
    touse     = st_local("touse")
    depvar    = st_local("depvar")
    indepvars = st_local("indepvars")
    treatvar  = st_local("treat")
    citvar    = st_local("cit")
    unitvar   = st_local("wfe_unit_idx")
    timevar   = st_local("wfe_time_idx")

    // Formula-side locals must resolve to concrete variables before st_view().
    // Otherwise Mata silently materializes empty views and defers a parser
    // contract violation into later numeric checks.
    if (strtrim(touse) != "") {
        _wfe_b_validate_touse_mask(touse, "_wfe_bridge_construct_views")
    }
    if (strtrim(depvar) == "") {
        errprintf("_wfe_bridge_construct_views: depvar local must name an existing variable\n")
        _error(3200)
    }
    _wfe_b_single_numeric_var(depvar, "depvar")
    if (strtrim(indepvars) == "") {
        errprintf("_wfe_bridge_construct_views: indepvars local must name at least one existing regressor\n")
        _error(3200)
    }
    _wfe_b_numeric_varlist(indepvars, "indepvars")

    // Y: outcome (NT×1)
    st_view(Y, ., depvar, touse)
    if (any(Y :>= .)) {
        errprintf("_wfe_bridge_construct_views: Y must not contain missing values\n")
        _error(3498)
    }

    // X: covariates (NT×p)
    st_view(X, ., tokens(indepvars), touse)
    if (cols(X) > 0 & any(X :>= .)) {
        errprintf("_wfe_bridge_construct_views: X must not contain missing values\n")
        _error(3498)
    }

    // The treatment variable is required input; an empty bridge local would
    // otherwise create an empty view and defer the contract violation.
    if (strtrim(treatvar) == "") {
        errprintf("_wfe_bridge_construct_views: treat local must name an existing variable\n")
        _error(3200)
    }
    _wfe_b_single_numeric_var(treatvar, "treat")
    if (strtrim(unitvar) == "") {
        errprintf("_wfe_bridge_construct_views: wfe_unit_idx local must name an existing variable\n")
        _error(3200)
    }
    _wfe_b_single_numeric_var(unitvar, "wfe_unit_idx")
    if (strtrim(timevar) == "") {
        errprintf("_wfe_bridge_construct_views: wfe_time_idx local must name an existing variable\n")
        _error(3200)
    }
    _wfe_b_single_numeric_var(timevar, "wfe_time_idx")

    // treat: treatment indicator (NT×1), must be 0/1
    st_view(treat, ., treatvar, touse)
    if (any(treat :>= .)) {
        errprintf("_wfe_bridge_construct_views: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat :!= 0) :& (treat :!= 1))) {
        errprintf("_wfe_bridge_construct_views: treat must contain only 0/1 values\n")
        _error(3498)
    }

    // unit_idx, time_idx: panel indices
    st_view(unit_idx, ., unitvar, touse)
    st_view(time_idx, ., timevar, touse)
    if (any(unit_idx :>= .) | any(time_idx :>= .)) {
        errprintf("_wfe_bridge_construct_views: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_idx :< 1) | any(time_idx :< 1) | ///
        any(unit_idx :!= floor(unit_idx)) | any(time_idx :!= floor(time_idx))) {
        errprintf("_wfe_bridge_construct_views: unit_idx/time_idx must contain positive integer values\n")
        _error(3200)
    }

    // cit: weight correction C_it (default = 1)
    if (citvar != "") {
        // An explicit cit() input is part of the public bridge contract.
        // If the caller-side variable name is stale, fail here instead of
        // leaking a generic st_view() variable-not-found error.
        _wfe_b_single_numeric_var(citvar, "cit")
        st_view(cit, ., citvar, touse)
        if (any(cit :>= .)) {
            errprintf("_wfe_bridge_construct_views: cit must not contain missing values\n")
            _error(3498)
        }
        if (any(cit :< 0)) {
            errprintf("_wfe_bridge_construct_views: cit must be non-negative\n")
            _error(3498)
        }
    }
    else {
        cit = J(rows(Y), 1, 1)
    }
}


real matrix _wfe_bridge_as_design(
    real matrix  X,
    real scalar  NT
)
{
    real matrix X_use

    X_use = X[., .]
    if (rows(X_use) == 1 & cols(X_use) == NT & NT > 1) {
        X_use = X_use'
    }

    return(X_use)
}


real matrix _wfe_bridge_as_outcome(
    real matrix  Y,
    real scalar  NT
)
{
    real matrix Y_use

    Y_use = Y[., .]
    if (rows(Y_use) == 1 & cols(Y_use) == NT & NT > 1) {
        Y_use = Y_use'
    }

    return(Y_use)
}


// ============================================================
// _wfe_bridge_read_params — Read and convert Stata locals to Mata types
//
// Boolean conversion rules:
//   - hetero_se/auto_se/df_adj: "on"→1, "off"→0
//   - white/verbose: empty=on(1), "nowhite"/"noverbose"=off(0)
//   - unweighted/unbiased_se/store_wdm: empty=off(0), non-empty=on(1)
// ============================================================
void _wfe_bridge_read_params(
    string scalar method,
    string scalar qoi,
    string scalar estimator,
    real scalar   maxdev_did,
    real scalar   has_maxdev_did,
    real scalar   hetero_se,
    real scalar   auto_se,
    real scalar   df_adj,
    real scalar   white,
    real scalar   white_alpha,
    real scalar   unweighted,
    real scalar   unbiased_se,
    real scalar   verbose,
    real scalar   store_wdm,
    real scalar   tol,
    real scalar   N_units,
    real scalar   N_times,
    real scalar   NT,
    real scalar   p,
    real matrix   X)
{
    string scalar hetero_se_local, auto_se_local, df_adj_local
    string scalar white_local, unweighted_local, unbiased_local
    string scalar verbose_local, store_wdm_local

    // String parameters (direct pass)
    method    = st_local("method")
    qoi       = st_local("qoi")
    estimator = st_local("estimator")
    if (method != "unit" & method != "time") {
        errprintf("_wfe_bridge_read_params: method must be 'unit' or 'time'\n")
        _error(3200)
    }
    if (qoi != "ate" & qoi != "att") {
        errprintf("_wfe_bridge_read_params: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    if (estimator != "" & estimator != "fd" & estimator != "did" & estimator != "Mdid") {
        errprintf("_wfe_bridge_read_params: estimator must be empty, 'fd', 'did', or 'Mdid'\n")
        _error(3200)
    }
    if (method == "time" & estimator != "") {
        errprintf("_wfe_bridge_read_params: estimator '%s' is not compatible with method 'time'\n",
                  estimator)
        _error(3200)
    }
    has_maxdev_did = (st_local("maxdev_did") != "")
    if (has_maxdev_did & estimator != "Mdid") {
        errprintf("_wfe_bridge_read_params: maxdev_did is allowed only with estimator 'Mdid'\n")
        _error(3200)
    }
    if (has_maxdev_did) {
        maxdev_did = strtoreal(st_local("maxdev_did"))
        if (missing(maxdev_did)) {
            errprintf("_wfe_bridge_read_params: maxdev_did must be numeric when provided\n")
            _error(3200)
        }
        if (maxdev_did < 0) {
            errprintf("_wfe_bridge_read_params: maxdev_did must be non-negative when provided\n")
            _error(3200)
        }
    }
    else {
        maxdev_did = .
    }

    // Boolean switches (Stata string convention → Mata 0/1)
    hetero_se_local = st_local("hetero_se")
    auto_se_local   = st_local("auto_se")
    df_adj_local    = st_local("df_adjustment")
    white_local     = st_local("white")
    unweighted_local = st_local("unweighted")
    unbiased_local  = st_local("unbiased_se")
    verbose_local   = st_local("verbose")
    store_wdm_local = st_local("store_wdm")
    if (hetero_se_local != "on" & hetero_se_local != "off") {
        errprintf("_wfe_bridge_read_params: hetero_se must be 'on' or 'off'\n")
        _error(3200)
    }
    if (auto_se_local != "on" & auto_se_local != "off") {
        errprintf("_wfe_bridge_read_params: auto_se must be 'on' or 'off'\n")
        _error(3200)
    }
    if (df_adj_local != "on" & df_adj_local != "off") {
        errprintf("_wfe_bridge_read_params: df_adjustment must be 'on' or 'off'\n")
        _error(3200)
    }
    if (white_local != "" & white_local != "nowhite") {
        errprintf("_wfe_bridge_read_params: white must be empty or 'nowhite'\n")
        _error(3200)
    }
    if (unweighted_local != "" & unweighted_local != "unweighted") {
        errprintf("_wfe_bridge_read_params: unweighted must be empty or 'unweighted'\n")
        _error(3200)
    }
    if (unbiased_local != "" & unbiased_local != "unbiased_se") {
        errprintf("_wfe_bridge_read_params: unbiased_se must be empty or 'unbiased_se'\n")
        _error(3200)
    }
    if (verbose_local != "" & verbose_local != "noverbose") {
        errprintf("_wfe_bridge_read_params: verbose must be empty or 'noverbose'\n")
        _error(3200)
    }
    if (store_wdm_local != "" & store_wdm_local != "store_wdm") {
        errprintf("_wfe_bridge_read_params: store_wdm must be empty or 'store_wdm'\n")
        _error(3200)
    }
    hetero_se   = (hetero_se_local == "on")
    auto_se     = (auto_se_local == "on")
    df_adj      = (df_adj_local == "on")
    white       = (white_local == "")                  // empty=on, "nowhite"=off
    unweighted  = (unweighted_local != "")             // empty=off, non-empty=on
    unbiased_se = (unbiased_local != "")
    verbose     = (verbose_local == "")                // empty=on, "noverbose"=off
    store_wdm   = (store_wdm_local != "")

    if (unbiased_se & (!hetero_se | auto_se)) {
        errprintf("_wfe_bridge_read_params: unbiased_se requires hetero_se(on) and auto_se(off)\n")
        _error(3200)
    }

    // Numeric scalars
    white_alpha = st_numscalar("white_alpha")
    tol         = st_numscalar("tol")
    if (white_alpha >= . | white_alpha <= 0 | white_alpha >= 1) {
        errprintf("_wfe_bridge_read_params: white_alpha must lie in (0,1)\n")
        _error(3200)
    }
    if (tol >= . | tol <= 0) {
        errprintf("_wfe_bridge_read_params: tol must be positive\n")
        _error(3200)
    }

    // Panel dimensions (computed by caller)
    N_units = st_numscalar("N_units")
    N_times = st_numscalar("N_times")
    NT      = st_numscalar("NT")
    if (N_units >= . | N_times >= . | NT >= . | ///
        N_units != floor(N_units) | N_times != floor(N_times) | ///
        NT != floor(NT) | N_units < 1 | N_times < 1 | NT < 1) {
        errprintf("_wfe_bridge_read_params: N_units/N_times/NT must be positive integers\n")
        _error(3200)
    }

    X = _wfe_bridge_as_design(X, NT)

    // Number of covariates (inferred from X)
    p = cols(X)
    if (p < 1) {
        errprintf("_wfe_bridge_read_params: X must contain at least one regressor\n")
        _error(3200)
    }
    if (rows(X) != NT) {
        errprintf("_wfe_bridge_read_params: rows(X) must equal NT\n")
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("_wfe_bridge_read_params: X must not contain missing values\n")
        _error(3498)
    }
    _wfe_bridge_chk_panel(N_units, N_times, NT)
}


real scalar _wfe_bridge_covariance_asymmetry(real matrix V)
{
    real scalar skew, scale

    skew = max(abs(V - V'))
    scale = max(abs(V))
    if (scale > 0) {
        return(skew / scale)
    }

    return(skew)
}


void _wfe_bridge_psd_mat(
    real matrix    V,
    string scalar  name
)
{
    real rowvector evals
    real scalar eig_tol, eig_scale

    symeigensystem(V, ., evals)
    eig_scale = max((1, max(abs(evals))))
    eig_tol = sqrt(epsilon(1)) * eig_scale

    if (min(evals) < -eig_tol) {
        errprintf("_wfe_bridge_post_results: %s must be positive semidefinite\n",
                  name)
        _error(3200)
    }
}


void _wfe_bridge_emit_error(string scalar msg)
{
    stata("display as error " + char(34) + msg + char(34))
}


void _wfe_bridge_validate_post_panel(
    real scalar   N_units,
    real scalar   N_times,
    real scalar   NT,
    string scalar caller
)
{
    string scalar touse, unitvar, timevar
    real colvector unit_idx, time_idx, seen_units, seen_times
    real matrix cell_count
    real scalar k

    touse = st_local("touse")
    unitvar = st_local("wfe_unit_idx")
    timevar = st_local("wfe_time_idx")

    _wfe_b_validate_touse_mask(touse, caller)
    _wfe_b_req_numlocal(unitvar, "wfe_unit_idx", caller)
    _wfe_b_req_numlocal(timevar, "wfe_time_idx", caller)

    st_view(unit_idx, ., unitvar, touse)
    st_view(time_idx, ., timevar, touse)

    if (rows(unit_idx) != NT | rows(time_idx) != NT) {
        errprintf("%s: wfe_unit_idx/wfe_time_idx rows must equal NT\n", caller)
        _error(3200)
    }

    if (any(unit_idx :>= .) | any(time_idx :>= .)) {
        errprintf("%s: wfe_unit_idx/wfe_time_idx must not contain missing values\n", caller)
        _error(3498)
    }

    seen_units = J(N_units, 1, 0)
    seen_times = J(N_times, 1, 0)
    cell_count = J(N_times, N_units, 0)

    for (k = 1; k <= NT; k++) {
        if (unit_idx[k] != floor(unit_idx[k]) | time_idx[k] != floor(time_idx[k]) | ///
            unit_idx[k] < 1 | unit_idx[k] > N_units | ///
            time_idx[k] < 1 | time_idx[k] > N_times) {
            errprintf("%s: N_units/N_times must match observed dense panel support\n", caller)
            _error(3200)
        }
        cell_count[time_idx[k], unit_idx[k]] = cell_count[time_idx[k], unit_idx[k]] + 1
        if (cell_count[time_idx[k], unit_idx[k]] > 1) {
            errprintf("%s: unit-time pair is not unique\n", caller)
            _error(498)
        }
        seen_units[unit_idx[k]] = 1
        seen_times[time_idx[k]] = 1
    }

    if (any(seen_units :== 0) | any(seen_times :== 0)) {
        errprintf("%s: N_units/N_times must match observed dense panel support\n", caller)
        _error(3200)
    }
}


void _wfe_bridge_validate_W_support(
    real matrix W,
    real scalar N_units,
    real scalar N_times,
    real scalar NT
)
{
    string scalar touse, unitvar, timevar, estimator_local, unweighted_local
    real colvector unit_idx, time_idx
    real matrix observed_support
    real scalar k

    estimator_local = st_local("estimator")
    unweighted_local = st_local("unweighted")

    // Two-way unweighted publishing posts a dense all-one W matrix on
    // the declared T x N support (W_it = 1 for all observed cells).
    if (unweighted_local != "" & ///
        (estimator_local == "did" | estimator_local == "Mdid")) {
        if (any(W :!= 1)) {
            errprintf("_wfe_bridge_post_results: W must be a dense all-one grid when estimator is unweighted did/Mdid\n")
            _error(498)
        }
        return
    }

    touse = st_local("touse")
    unitvar = st_local("wfe_unit_idx")
    timevar = st_local("wfe_time_idx")

    // Skip W support validation when full Stata context is unavailable
    // (e.g. direct unit-test calls to _wfe_bridge_post_results)
    if (touse == "" | unitvar == "" | timevar == "") {
        return
    }

    // Post-estimation support validation is set-based: it must reject malformed
    // panel keys, duplicates, or dense-support mismatches, but it must not
    // inherit the pre-estimation row-order requirement from read_params().
    _wfe_bridge_validate_post_panel(
        N_units, N_times, NT, "_wfe_bridge_post_results")

    st_view(unit_idx, ., unitvar, touse)
    st_view(time_idx, ., timevar, touse)

    observed_support = J(N_times, N_units, 0)
    for (k = 1; k <= NT; k++) {
        observed_support[time_idx[k], unit_idx[k]] = 1
    }

    if (unweighted_local != "") {
        // One-way unweighted estimation still assigns regression weight 1 to
        // each observed row only. Because the posted Stata object is a full
        // T x N matrix, the corresponding public contract is the observed-
        // support 0/1 mask rather than arbitrary nonzero values on support.
        if (any(W :!= observed_support)) {
            errprintf("_wfe_bridge_post_results: W must equal the observed-support 0/1 mask when one-way unweighted is requested\n")
            _error(498)
        }
        return
    }

    if (any((W :!= 0) :& (observed_support :== 0))) {
        errprintf("_wfe_bridge_post_results: W must be zero outside observed support unless two-way unweighted\n")
        _error(498)
    }
}


void _wfe_bridge_validate_covariance(
    real matrix    V,
    real scalar    p,
    string scalar  name
)
{
    if (rows(V) != p | cols(V) != p) {
        errprintf("_wfe_bridge_post_results: %s must be %gx%g to match beta\n",
                  name, p, p)
        _error(3200)
    }
    if (any(V :>= .)) {
        errprintf("_wfe_bridge_post_results: %s must not contain missing values\n",
                  name)
        _error(3498)
    }
}


// ============================================================
// _wfe_bridge_post_results — Post estimation results back to Stata
//
// Results are stored in temporary matrices/scalars/macros with __wfe_ prefix,
// then moved to e() by the ado layer after ereturn post.
//
// Required arguments:
//   beta     1D coefficient vector (canonicalized to 1×p row vector)
//   vcov     p×p covariance matrix (symmetric)
//   W        T×N weight matrix
//   N_nonzero  count of non-zero regression weights
//   df_r     residual degrees of freedom
//   sigma    residual standard deviation
//   vcetype  standard error type label
//
// Optional arguments (|):
//   beta_fe, vcov_fe    FE benchmark coefficients and covariance
//   white_stat, white_pvalue, white_test   White test results
//   Y_wdm, X_wdm        Weighted-demeaned data (for store_wdm)
// ============================================================
void _wfe_bridge_post_results(
    real vector    beta,
    real matrix    vcov,
    real matrix    W,
    real scalar    N_nonzero,
    real scalar    df_r,
    real scalar    sigma,
    string scalar  vcetype,
    | real vector    beta_fe,
      real matrix    vcov_fe,
      real scalar    white_stat,
      real scalar    white_pvalue,
      string scalar  white_test,
      real matrix    Y_wdm,
      real matrix    X_wdm)
{
    real scalar p, N_units_declared, N_times_declared
    real scalar actual_n_nonzero, NT_cap, expected_n_nonzero

    // Clear any stale results before validation
    stata("capture matrix drop __wfe_b")
    stata("capture matrix drop __wfe_V")
    stata("capture matrix drop __wfe_W")
    stata("capture matrix drop __wfe_b_fe")
    stata("capture matrix drop __wfe_V_fe")
    stata("capture matrix drop __wfe_Y_wdm")
    stata("capture matrix drop __wfe_X_wdm")
    stata("capture scalar drop __wfe_N_nonzero")
    stata("capture scalar drop __wfe_df_r")
    stata("capture scalar drop __wfe_sigma")
    stata("capture scalar drop __wfe_N_units")
    stata("capture scalar drop __wfe_N_times")
    stata("capture scalar drop __wfe_white_stat")
    stata("capture scalar drop __wfe_white_pvalue")
    stata("capture macro drop __wfe_vcetype")
    stata("capture macro drop __wfe_white_test")

    // beta: ensure 1×p row vector
    if (cols(beta) == 1 & rows(beta) > 1) {
        beta = beta'
    }
    p = cols(beta)
    if (rows(beta) != 1 | p < 1) {
        errprintf("_wfe_bridge_post_results: beta must be a nonempty rowvector\n")
        _error(3200)
    }
    if (any(beta :>= .)) {
        errprintf("_wfe_bridge_post_results: beta must not contain missing values\n")
        _error(3498)
    }

    // Validate vcov shape before symmetry so malformed non-square inputs
    // fail with the bridge contract instead of a deep conformability error.
    _wfe_bridge_validate_covariance(vcov, p, "vcov")
    if (_wfe_bridge_covariance_asymmetry(vcov) > 1e-12) {
        _wfe_bridge_emit_error("_wfe_bridge_post_results: vcov must be symmetric")
        _error(3200)
    }
    vcov = 0.5 :* (vcov + vcov')
    _wfe_bridge_psd_mat(vcov, "vcov")

    // W: T×N weight matrix
    N_units_declared = st_numscalar("N_units")
    N_times_declared = st_numscalar("N_times")
    if (N_units_declared >= . | N_times_declared >= . | ///
        N_units_declared != floor(N_units_declared) | ///
        N_times_declared != floor(N_times_declared) | ///
        N_units_declared < 1 | N_times_declared < 1) {
        errprintf("_wfe_bridge_post_results: N_units and N_times must be positive integers\n")
        _error(3200)
    }
    if (rows(W) != N_times_declared | cols(W) != N_units_declared) {
        errprintf("_wfe_bridge_post_results: W must be %gx%g to match N_times x N_units\n",
                  N_times_declared, N_units_declared)
        _error(3200)
    }
    if (any(W :>= .)) {
        errprintf("_wfe_bridge_post_results: W must not contain missing values\n")
        _error(3498)
    }
    stata("capture confirm scalar NT")
    if (st_numscalar("c(rc)") != 0) {
        _wfe_bridge_emit_error("_wfe_bridge_post_results: NT must be a positive integer when validating N_nonzero")
        _error(3200)
    }
    NT_cap = st_numscalar("NT")
    if (NT_cap >= . | NT_cap != floor(NT_cap) | NT_cap <= 0) {
        _wfe_bridge_emit_error("_wfe_bridge_post_results: NT must be a positive integer when validating N_nonzero")
        _error(3200)
    }
    _wfe_bridge_validate_W_support(W, N_units_declared, N_times_declared, NT_cap)

    // Scalar results
    if (N_nonzero >= . | N_nonzero < 0 | N_nonzero != floor(N_nonzero)) {
        errprintf("_wfe_bridge_post_results: N_nonzero must be a nonnegative integer count\n")
        _error(3200)
    }
    actual_n_nonzero = sum(W :!= 0)
    expected_n_nonzero = actual_n_nonzero
    if (NT_cap < expected_n_nonzero) {
        expected_n_nonzero = NT_cap
    }
    // For two-way unweighted estimation, R publishes a dense T×N all-one W grid
    // even on unbalanced panels. The regression runs only on NT observed rows,
    // so stored N_nonzero is min(sum(W != 0), NT).
    if (N_nonzero != expected_n_nonzero) {
        errprintf("_wfe_bridge_post_results: N_nonzero must equal min(sum(W != 0), NT)\n")
        _error(3200)
    }
    if (df_r >= . | df_r != floor(df_r) | df_r <= 0) {
        errprintf("_wfe_bridge_post_results: df_r must be a positive integer\n")
        _error(3200)
    }
    if (sigma >= . | sigma < 0) {
        errprintf("_wfe_bridge_post_results: sigma must be finite and non-negative\n")
        _error(3200)
    }

    // String result
    if (strtrim(vcetype) == "") {
        errprintf("_wfe_bridge_post_results: vcetype must be nonempty\n")
        _error(3200)
    }

    // Optional: FE benchmark results (must provide both or none)
    if (args() == 8) {
        errprintf("_wfe_bridge_post_results: beta_fe and vcov_fe must be provided together\n")
        _error(3200)
    }
    if (args() >= 9) {
        // beta_fe: ensure 1×p row vector
        if (cols(beta_fe) == 1 & rows(beta_fe) > 1) {
            beta_fe = beta_fe'
        }
        if (rows(beta_fe) != 1 | cols(beta_fe) != p) {
            errprintf("_wfe_bridge_post_results: beta_fe must be 1x%g to match beta\n",
                      p)
            _error(3200)
        }
        if (any(beta_fe :>= .)) {
            errprintf("_wfe_bridge_post_results: beta_fe must not contain missing values\n")
            _error(3498)
        }
        _wfe_bridge_validate_covariance(vcov_fe, p, "vcov_fe")
        if (_wfe_bridge_covariance_asymmetry(vcov_fe) > 1e-12) {
            _wfe_bridge_emit_error("_wfe_bridge_post_results: vcov_fe must be symmetric")
            _error(3200)
        }
        vcov_fe = 0.5 :* (vcov_fe + vcov_fe')
        _wfe_bridge_psd_mat(vcov_fe, "vcov_fe")
    }

    // Optional: White test results (must provide all three or none)
    if (args() > 9 & args() < 12) {
        errprintf("_wfe_bridge_post_results: white_stat, white_pvalue, and white_test must be provided together\n")
        _error(3200)
    }
    if (args() >= 12) {
        if (white_stat >= .) {
            errprintf("_wfe_bridge_post_results: white_stat must be finite\n")
            _error(3200)
        }
        if (white_stat < 0) {
            errprintf("_wfe_bridge_post_results: white_stat must be non-negative\n")
            _error(3200)
        }
        if (white_pvalue >= . | white_pvalue < 0 | white_pvalue > 1) {
            errprintf("_wfe_bridge_post_results: white_pvalue must lie in [0,1]\n")
            _error(3200)
        }
        if (white_test != "TRUE" & white_test != "FALSE") {
            errprintf("_wfe_bridge_post_results: white_test must be TRUE or FALSE\n")
            _error(3200)
        }
    }

    // Optional: weighted-demeaned data (must provide both Y_wdm and X_wdm or none)
    if (args() == 13) {
        errprintf("_wfe_bridge_post_results: Y_wdm and X_wdm must be provided together\n")
        _error(3200)
    }
    if (args() >= 14) {
        real scalar has_Y_wdm, has_X_wdm, NT_declared

        has_Y_wdm = (rows(Y_wdm) > 0 & cols(Y_wdm) > 0)
        has_X_wdm = (rows(X_wdm) > 0 & cols(X_wdm) > 0)

        // Empty matrices as placeholders are malformed; caller should omit entirely
        if (!has_Y_wdm & !has_X_wdm) {
            errprintf("_wfe_bridge_post_results: Y_wdm and X_wdm must have NT rows\n")
            _error(3200)
        }

        if (has_Y_wdm != has_X_wdm) {
            errprintf("_wfe_bridge_post_results: Y_wdm and X_wdm must be provided together\n")
            _error(3200)
        }

        if (has_Y_wdm) {
            NT_declared = st_numscalar("NT")
            if (NT_declared >= . | NT_declared != floor(NT_declared) | NT_declared < 1) {
                errprintf("_wfe_bridge_post_results: NT must be a positive integer when store_wdm is posted\n")
                _error(3200)
            }
            Y_wdm = _wfe_bridge_as_outcome(Y_wdm, NT_declared)
            X_wdm = _wfe_bridge_as_design(X_wdm, NT_declared)
            if (cols(Y_wdm) != 1) {
                errprintf("_wfe_bridge_post_results: Y_wdm must be NT x 1\n")
                _error(3200)
            }
            if (rows(X_wdm) != rows(Y_wdm)) {
                errprintf("_wfe_bridge_post_results: Y_wdm and X_wdm must have the same number of rows\n")
                _error(3200)
            }
            if (rows(Y_wdm) != NT_declared) {
                errprintf("_wfe_bridge_post_results: Y_wdm and X_wdm must have NT rows\n")
                _error(3200)
            }
            if (cols(X_wdm) != p) {
                errprintf("_wfe_bridge_post_results: X_wdm must have %g columns to match beta\n",
                          p)
                _error(3200)
            }
            if (any(Y_wdm :>= .) | any(X_wdm :>= .)) {
                errprintf("_wfe_bridge_post_results: Y_wdm/X_wdm must not contain missing values\n")
                _error(3498)
            }
        }
    }

    st_matrix("__wfe_b", beta)
    st_matrix("__wfe_V", vcov)
    st_matrix("__wfe_W", W)

    st_numscalar("__wfe_N_nonzero", N_nonzero)
    st_numscalar("__wfe_df_r", df_r)
    st_numscalar("__wfe_sigma", sigma)

    // Pass N_units/N_times back for ado ereturn
    st_numscalar("__wfe_N_units", st_numscalar("N_units"))
    st_numscalar("__wfe_N_times", st_numscalar("N_times"))
    st_global("__wfe_vcetype", vcetype)

    if (args() >= 9) {
        st_matrix("__wfe_b_fe", beta_fe)
        st_matrix("__wfe_V_fe", vcov_fe)
    }

    if (args() >= 12) {
        st_numscalar("__wfe_white_stat", white_stat)
        st_numscalar("__wfe_white_pvalue", white_pvalue)
        st_global("__wfe_white_test", white_test)
    }

    if (args() >= 14 & rows(Y_wdm) > 0 & cols(Y_wdm) > 0) {
        st_matrix("__wfe_Y_wdm", Y_wdm)
        st_matrix("__wfe_X_wdm", X_wdm)
    }
}


end
