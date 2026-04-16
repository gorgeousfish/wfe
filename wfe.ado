*! wfe.ado — Weighted Fixed Effects estimator for Stata
*! =================================================================

program define wfe, eclass sortpreserve
    version 16

    /*
       Parser-level failure guard: clear stale estimates before syntax
       so Stata-level syntax errors cannot leave an old e() replayable.
    */
    ereturn clear

    /*
       Stata version check
    */
    if c(stata_version) < 16 {
        display as error "wfe requires Stata 16 or later (Mata complex type support)"
        exit 198
    }

    /*
       Syntax parsing
    */
    syntax varlist(min=1 numeric fv) [if] [in],   ///
        TReat(varname)                             ///
        Unit(varname)                              ///
        [                                          ///
        Time(varname)                              ///
        Method(string)                             ///
        QOI(string)                                ///
        ESTIMator(string)                          ///
        CIT(varname)                               ///
        HETERO_se(string)                          ///
        AUTO_se(string)                            ///
        DF_adjustment(string)                      ///
        noWHITE                                    ///
        WHITE_alpha(string)                        ///
        UNWeighted                                 ///
        UNBIASED_se                                ///
        noVerbose                                  ///
        STORE_wdm                                  ///
        MAXDEV_did(string)                         ///
        TOL(string)                                ///
        DIAGnose                                   ///
        ]

    gettoken depvar indepvars_spec : varlist
    capture confirm numeric variable `depvar'
    if _rc {
        display as error "dependent variable must be a numeric variable name; factor-variable notation is not allowed"
        exit 198
    }

    capture confirm numeric variable `treat'
    if _rc {
        display as error "treat() must specify a numeric variable"
        exit 198
    }

    /*
       Formula-sample boundary for diagnostics
    */
    marksample touse_ifin, novarlist
    marksample touse
    markout `touse' `treat'

    quietly count if `touse_ifin'
    local nt_ifin = r(N)
    quietly count if `touse'
    local nt_formula = r(N)

    if "`verbose'" != "noverbose" & `nt_formula' < `nt_ifin' {
        display as text " " _newline "Missing values are removed"
    }

    if `nt_formula' == 0 {
        display as error "no observations"
        exit 2000
    }

    /*
       Panel identifiers define panel geometry even when y/X missings
       later remove an observation from the estimation sample.
    */
    quietly count if `touse_ifin' & missing(`unit')
    if r(N) > 0 {
        display as error "missing values in unit() are not allowed"
        exit 198
    }
    if "`time'" != "" {
        quietly count if `touse_ifin' & missing(`time')
        if r(N) > 0 {
            display as error "missing values in time() are not allowed"
            exit 198
        }

        /*
           Duplicate unit-time dyads are structural errors only on the
           formula-valid estimation sample (missings are dropped before
           panel indexing).
        */
        tempvar wfe_dup
        quietly duplicates tag `unit' `time' if `touse', gen(`wfe_dup')
        quietly count if `touse' & `wfe_dup' > 0
        if r(N) > 0 {
            display as error "duplicate observations found for unit-time combination"
            exit 198
        }
        drop `wfe_dup'
    }

    /*
       Stata's public interface separates treat() from additional
       regressors. Keep the causal treatment regressor in the scored
       design even when users only spell the extra covariates.
    */
    local score_indepvars_spec `"`indepvars_spec'"'
    if `"`score_indepvars_spec'"' == "" {
        local score_indepvars_spec "`treat'"
    }
    else {
        local __wfe_has_treat_main 0
        local __wfe_treat_suffix ".`treat'"
        fvexpand `score_indepvars_spec' if `touse'
        local __wfe_score_terms `"`r(varlist)'"'
        local __wfe_n_score_terms : word count `__wfe_score_terms'
        forvalues __wfe_i = 1/`__wfe_n_score_terms' {
            local __wfe_term : word `__wfe_i' of `__wfe_score_terms'
            if strpos("`__wfe_term'", "#") == 0 {
                if "`__wfe_term'" == "`treat'" | "`__wfe_term'" == "c.`treat'" {
                    local __wfe_has_treat_main 1
                }
                else if strlen("`__wfe_term'") >= strlen("`__wfe_treat_suffix'") {
                    local __wfe_term_tail = substr("`__wfe_term'", ///
                        strlen("`__wfe_term'") - strlen("`__wfe_treat_suffix'") + 1, .)
                    if "`__wfe_term_tail'" == "`__wfe_treat_suffix'" {
                        local __wfe_has_treat_main 1
                    }
                }
            }
        }
        if !`__wfe_has_treat_main' {
            local score_indepvars_spec `"`treat' `score_indepvars_spec'"'
        }
    }

    /*
       Default value population
    */
    if "`method'" == "" local method "unit"
    if "`qoi'" == "" local qoi "ate"
    if "`hetero_se'" == "" local hetero_se "on"
    if "`auto_se'" == "" local auto_se "on"
    if "`df_adjustment'" == "" local df_adjustment "on"
    if "`white_alpha'" == "" {
        local white_alpha 0.05
    }
    else {
        local __white_alpha_real = real("`white_alpha'")
        if `__white_alpha_real' == . {
            display as error "white_alpha() must lie in (0,1)"
            exit 198
        }
        local white_alpha `__white_alpha_real'
    }

    if "`tol'" == "" {
        local tol = -1
    }
    else {
        local __tol_real = real("`tol'")
        if `__tol_real' == . {
            display as error "tol() must be positive"
            exit 198
        }
        local tol `__tol_real'
    }

    * tol: sentinel -1 means "not specified" → use sqrt(c(epsdouble))
    if `tol' == -1 {
        local tol = sqrt(c(epsdouble))
    }
    else if missing(`tol') | `tol' <= 0 {
        display as error "tol() must be positive"
        exit 198
    }

    /*
       Parameter compatibility guards
    */

    /* Check: time missing + method="time" */
    if "`time'" == "" & "`method'" == "time" {
        display as error "'time.index' should be provided"
        exit 198
    }

    /* Check: method="time" + estimator="fd" */
    if "`method'" == "time" & "`estimator'" == "fd" {
        display as error `"First Difference is not compatible with 'time' method: set method == 'unit'"'
        exit 198
    }

    /* Check: method="time" + estimator="did" */
    if "`method'" == "time" & "`estimator'" == "did" {
        display as error `"Difference-in-Differences is not compatible with 'time' method: set method == 'unit'"'
        exit 198
    }

    /* Check: method="time" + estimator="Mdid" */
    if "`method'" == "time" & "`estimator'" == "Mdid" {
        display as error `"Match-Difference-in-Differences is not compatible with 'time' method: set method == 'unit'"'
        exit 198
    }

    /* Check: time missing + estimator="fd" */
    if "`time'" == "" & "`estimator'" == "fd" {
        display as error `"First Difference cannot calculate when 'time.index' is missing"'
        exit 198
    }

    /* Check: time missing + estimator="did"/"Mdid" */
    if "`time'" == "" & inlist("`estimator'", "did", "Mdid") {
        display as error "'time.index' should be provided"
        exit 198
    }

    /* Check: cit input validity */
    if "`cit'" != "" {
        capture confirm numeric variable `cit'
        if _rc {
            display as error "'C.it' must be a numeric vector with length equal to number of observations"
            exit 198
        }
        quietly count if `touse' & missing(`cit')
        if r(N) > 0 {
            display as error "missing values in cit() are not allowed"
            exit 198
        }
        quietly count if `touse' & `cit' < 0
        if r(N) > 0 {
            display as error "'C.it' must be a non-negative numeric vector"
            exit 198
        }
    }

    /* Check: maxdev_did negative value; zero is exact matching */
    if "`maxdev_did'" != "" {
        local __maxdev_real = real("`maxdev_did'")
        if `__maxdev_real' == . {
            display as error "maxdev_did() must be a numeric value"
            exit 198
        }
        if `__maxdev_real' < 0 {
            display as error "maxdev_did() must be non-negative"
            exit 198
        }
    }

    /* Check: method invalid value */
    if "`method'" != "unit" & "`method'" != "time" {
        display as error "method should be either unit or time"
        exit 198
    }

    /* Check: qoi invalid value */
    if "`qoi'" != "ate" & "`qoi'" != "att" {
        display as error "qoi() must be 'ate' or 'att'"
        exit 198
    }

    /* Check: estimator invalid value */
    if "`estimator'" != "" & "`estimator'" != "fd" & "`estimator'" != "did" & "`estimator'" != "Mdid" {
        display as error "estimator() must be 'fd', 'did', or 'Mdid'"
        exit 198
    }

    /* Check: maxdev_did context guard */
    if "`maxdev_did'" != "" & "`estimator'" != "Mdid" {
        display as error "maxdev_did() is allowed only with estimator(Mdid)"
        exit 198
    }

    /* Check: hetero_se invalid value */
    if "`hetero_se'" != "on" & "`hetero_se'" != "off" {
        display as error "hetero_se() must be 'on' or 'off'"
        exit 198
    }

    /* Check: auto_se invalid value */
    if "`auto_se'" != "on" & "`auto_se'" != "off" {
        display as error "auto_se() must be 'on' or 'off'"
        exit 198
    }

    /* Check: df_adjustment invalid value */
    if "`df_adjustment'" != "on" & "`df_adjustment'" != "off" {
        display as error "df_adjustment() must be 'on' or 'off'"
        exit 198
    }

    /* Check: white_alpha invalid probability */
    if `white_alpha' <= 0 | `white_alpha' >= 1 {
        display as error "white_alpha() must lie in (0,1)"
        exit 198
    }

    /*
       Core treatment-support guard on the formula-valid sample.
       This must run before preprocess/bridge work so invalid public
       treat() input fails fast without leaking internal tempvar chatter.
    */
    quietly levelsof `treat' if `touse', local(__wfe_treat_levels)
    local __wfe_treat_count : word count `__wfe_treat_levels'
    if `__wfe_treat_count' != 2 {
        if `__wfe_treat_count' > 2 {
            display as error "'treat' must be a binary vector: there are more than two values of treatment"
        }
        else {
            display as error "'treat' must contain both 0 and 1 values"
        }
        exit 198
    }
    quietly summarize `treat' if `touse', meanonly
    if r(min) != 0 | r(max) != 1 {
        display as error "'treat' must be a either zero or one where one indicates treatment"
        exit 198
    }

    /*
       SE combination guards (path-aware)

       One-way FE path: estimator="" or estimator="fd"
       Two-way FE path: estimator="did" or estimator="Mdid"
    */

    * Path determination
    local __is_twoway = 0
    if "`estimator'" == "did" | "`estimator'" == "Mdid" {
        local __is_twoway = 1
    }

    /* Check: unbiased_se path validity */
    if "`unbiased_se'" != "" {
        if `__is_twoway' == 1 | ("`hetero_se'" == "on" & "`auto_se'" == "on") {
            display as error "unbiased_se is allowed only for one-way hetero_se(on) auto_se(off)"
            exit 198
        }
    }

    if `__is_twoway' == 0 {
        * One-way FE path

        /* Check: hetero_se=off + auto_se=off */
        if "`hetero_se'" == "off" & "`auto_se'" == "off" {
            display as error "standard errors with independence and homoskedasticity is not supported"
            exit 198
        }

        /* Check: hetero_se=off + auto_se=on */
        if "`hetero_se'" == "off" & "`auto_se'" == "on" {
            display as error "robust standard errors with autocorrelation and homoskedasticity is not supported"
            exit 198
        }
    }
    else {
        * Two-way FE path

        /* Check: hetero_se=on + auto_se=off */
        if "`hetero_se'" == "on" & "`auto_se'" == "off" {
            display as error "two-way FE requires hetero_se(on) and auto_se(on)"
            exit 198
        }

        /* Check: hetero_se=off + auto_se=off */
        if "`hetero_se'" == "off" & "`auto_se'" == "off" {
            display as error "two-way FE requires hetero_se(on) and auto_se(on)"
            exit 198
        }

        /* Check: hetero_se=off + auto_se=on */
        if "`hetero_se'" == "off" & "`auto_se'" == "on" {
            display as error "Robust standard errors with autocorrelation and homoskedasticity is not supported"
            exit 198
        }
    }

    /*
       qoi → causal mapping + unweighted override
    */
    if "`unweighted'" == "unweighted" {
        local __wfe_causal "Unweighted (Standard) Fixed Effect"
    }
    else if "`qoi'" == "ate" {
        local __wfe_causal "ATE (Average Treatment Effect)"
    }
    else if "`qoi'" == "att" {
        local __wfe_causal "ATT (Average Treatment Effect for the Treated)"
    }

    /*
       Diagnostic output
    */
    if "`diagnose'" == "diagnose" {
        * Clear stale estimation results
        ereturn clear
        local __white_state "on"
        if "`white'" == "nowhite" {
            local __white_state "off"
        }
        local __estimator_state "`estimator'"
        if "`__estimator_state'" == "" {
            local __estimator_state "NULL"
        }
        local __verbose_state "on"
        if "`verbose'" == "noverbose" {
            local __verbose_state "off"
        }
        local __unweighted_state "off"
        if "`unweighted'" == "unweighted" {
            local __unweighted_state "on"
        }
        local __unbiased_state "off"
        if "`unbiased_se'" != "" {
            local __unbiased_state "on"
        }
        local __store_wdm_state "off"
        if "`store_wdm'" == "store_wdm" {
            local __store_wdm_state "on"
        }
        display as text "method     = `method'"
        display as text "qoi        = `qoi'"
        display as text "estimator  = `__estimator_state'"
        display as text "hetero_se  = `hetero_se'"
        display as text "auto_se    = `auto_se'"
        display as text "df_adj     = `df_adjustment'"
        display as text "white      = `__white_state'"
        display as text "white_alpha = `white_alpha'"
        display as text "tol        = `tol'"
        display as text "unweighted = `__unweighted_state'"
        display as text "unbiased   = `__unbiased_state'"
        display as text "verbose    = `__verbose_state'"
        display as text "store_wdm  = `__store_wdm_state'"
        display as text "maxdev_did = `maxdev_did'"
        display as text "causal     = `__wfe_causal'"
        display as text "touse N    = " _continue
        quietly count if `touse'
        display as text r(N)
        exit
    }

    /*
       Preprocessing path
    */
    tempvar wfe_unit_idx wfe_time_idx wfe_cit wfe_row_order
    local __wfe_model_varlist `"`depvar' `score_indepvars_spec'"'
    local __preprocess_cmd `"`__wfe_model_varlist' `if' `in', treat(`treat') unit(`unit') method(`method')"'
    if "`time'" != "" {
        local __preprocess_cmd `"`__preprocess_cmd' time(`time')"'
    }
    if "`cit'" != "" {
        local __preprocess_cmd `"`__preprocess_cmd' cit(`cit')"'
    }
    if "`verbose'" == "noverbose" {
        local __preprocess_cmd `"`__preprocess_cmd' noverbose"'
    }
    local __preprocess_cmd `"`__preprocess_cmd' outunit(`wfe_unit_idx') outtime(`wfe_time_idx') outcit(`wfe_cit') outroworder(`wfe_row_order')"'
    _wfe_preprocess_inspect `__preprocess_cmd'

    tempvar wfe_esample
    quietly gen byte `wfe_esample' = !missing(`wfe_unit_idx')
    local touse `wfe_esample'
    local cit `wfe_cit'

    /*
       Sort by dense (unit_idx, time_idx) before the bridge builds
       zero-copy views. The Mata bridge assumes this internal order.
    */
    sort `wfe_unit_idx' `wfe_time_idx' `wfe_row_order'

    local N = r(N)
    local T = r(T)
    local NT = r(NT)
    local unit_number = r(unit_number)
    local panel_balanced = r(panel_balanced)

    scalar N_units = `N'
    scalar N_times = `T'
    scalar NT = `NT'
    scalar panel_balanced = `panel_balanced'

    if `N' * `T' > 50000 & "`verbose'" != "noverbose" {
        display as text "note: N*T = " `N' * `T' ///
            " is large. Consider using nowhite option to reduce memory usage."
    }

    /*
       Bridge pre-contract checks
    */

    * Verify tempvar existence
    capture confirm variable `wfe_unit_idx'
    if _rc {
        display as error "[wfe bridge] Input view construction failed: panel unit index variable does not exist"
        exit 498
    }

    capture confirm variable `wfe_time_idx'
    if _rc {
        display as error "[wfe bridge] Input view construction failed: panel time index variable does not exist"
        exit 498
    }

    * Verify count scalars
    foreach s in N_units N_times NT {
        if scalar(`s') == . | scalar(`s') <= 0 {
            display as error "[wfe bridge] Input view construction failed: scalar `s' invalid (value=`=scalar(`s')')"
            exit 498
        }
    }

    * Verify critical locals
    foreach loc in method qoi {
        if `"``loc''"' == "" {
            display as error "[wfe bridge] Input view construction failed: parse result `loc' is empty"
            exit 498
        }
    }

    /*
       Scalar bridge for Mata access
    */
    local __had_bridge_white_alpha = 0
    local __had_bridge_tol = 0
    local __had_bridge_p = 0
    tempname __wfe_saved_bridge_white_alpha __wfe_saved_bridge_tol __wfe_saved_bridge_p

    capture confirm scalar white_alpha
    if !_rc {
        scalar `__wfe_saved_bridge_white_alpha' = scalar(white_alpha)
        local __had_bridge_white_alpha = 1
    }
    capture confirm scalar tol
    if !_rc {
        scalar `__wfe_saved_bridge_tol' = scalar(tol)
        local __had_bridge_tol = 1
    }
    capture confirm scalar p
    if !_rc {
        scalar `__wfe_saved_bridge_p' = scalar(p)
        local __had_bridge_p = 1
    }

    scalar white_alpha = `white_alpha'
    scalar tol = `tol'

    /*
       Extract depvar and expand factor-variable regressors
    */
    local indepvars_spec `score_indepvars_spec'

    fvexpand `indepvars_spec' if `touse'
    local __wfe_fvexpanded `"`r(varlist)'"'

    fvrevar `indepvars_spec' if `touse'
    local __wfe_fvvars `"`r(varlist)'"'

    local indepvars_names
    local indepvars
    local __wfe_n_terms : word count `__wfe_fvexpanded'
    forvalues __wfe_i = 1/`__wfe_n_terms' {
        local __wfe_term : word `__wfe_i' of `__wfe_fvexpanded'
        local __wfe_var : word `__wfe_i' of `__wfe_fvvars'
        if strpos("`__wfe_term'", "b.") == 0 & strpos("`__wfe_term'", "o.") == 0 {
            local indepvars_names "`indepvars_names' `__wfe_term'"
            local indepvars "`indepvars' `__wfe_var'"
        }
    }

    * Count of independent variables after factor-variable expansion
    local p : word count `indepvars_names'
    scalar p = `p'

    /*
       One-way unit-FE sigma2/df_r uses NT - p - J_u.
       When that denominator is nonpositive, the public command should
       fail-fast at the entry boundary instead of leaking a deep Mata
       helper prefix plus a generic 3351 trailer.
    */
    if "`method'" == "unit" & !inlist("`estimator'", "did", "Mdid") {
        local __wfe_unit_df_r = `NT' - `p' - `N'
        if `__wfe_unit_df_r' <= 0 {
            __wfe_restore_bridge_scalars, ///
                hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
                hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
                hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
            display as error ///
                "wfe: insufficient degrees of freedom: NT=`NT', p=`p', J_u=`N', df_r=`__wfe_unit_df_r'"
            exit 498
        }
    }

    __wfe_ensure_backend_mata

    /*
       Path dispatch to Mata
    */

    if "`store_wdm'" == "store_wdm" & inlist("`estimator'", "did", "Mdid") {
        __wfe_restore_bridge_scalars, ///
            hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
            hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
            hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
        display as error "store_wdm is not supported for two-way estimators (did/Mdid)"
        exit 198
    }

    if inlist("`estimator'", "did", "Mdid") {
        * Two-way FE path: DiD / Multi-DiD
        capture mata: _wfe_twoway_estimate()
        if _rc {
            if _rc == 3351 {
                __wfe_restore_bridge_scalars, ///
                    hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
                    hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
                    hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
                display as error ///
                    "two-way GMM df_adjustment(on) requires Nstar - N_units - N_times - p + 1 > 0; try df_adjustment(off)"
                exit 498
            }
            else if _rc == 3352 {
                * FE-side HAC infeasible for White test; retry with nowhite
                local white "nowhite"
                display as text "(note: White test skipped — FE-side degrees of freedom insufficient)"
                capture mata: _wfe_twoway_estimate()
                if _rc {
                    capture noisily mata: _wfe_twoway_estimate()
                    local __wfe_runtime_rc = _rc
                    __wfe_restore_bridge_scalars, ///
                        hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
                        hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
                        hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
                    exit `__wfe_runtime_rc'
                }
            }
            else {
                capture noisily mata: _wfe_twoway_estimate()
                local __wfe_runtime_rc = _rc
                __wfe_restore_bridge_scalars, ///
                    hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
                    hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
                    hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
                exit `__wfe_runtime_rc'
            }
        }
    }
    else {
        * One-way FE path: standard FE / First-Difference
        capture noisily mata: _wfe_oneway_estimate()
        if _rc {
            local __wfe_runtime_rc = _rc
            __wfe_restore_bridge_scalars, ///
                hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
                hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
                hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
            exit `__wfe_runtime_rc'
        }
    }

    /*
       Result postback & ereturn assembly
    */

    * Step 0: Result integrity checks
    capture confirm matrix __wfe_b
    if _rc {
        __wfe_restore_bridge_scalars, ///
            hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
            hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
            hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
        display as error "[wfe bridge] Result postback failed: coefficient matrix __wfe_b does not exist"
        exit 498
    }
    capture confirm matrix __wfe_V
    if _rc {
        __wfe_restore_bridge_scalars, ///
            hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
            hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
            hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
        display as error "[wfe bridge] Result postback failed: covariance matrix __wfe_V does not exist"
        exit 498
    }

    * Dimension validation
    local p_b = colsof(__wfe_b)
    local r_b = rowsof(__wfe_b)
    local p_V = colsof(__wfe_V)
    local r_V = rowsof(__wfe_V)

    if `p_b' != `p_V' | `r_V' != `p_V' {
        __wfe_restore_bridge_scalars, ///
            hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
            hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
            hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
        display as error "[wfe bridge] Result postback failed: dimension mismatch b(`p_b') vs V(`r_V'x`p_V')"
        exit 498
    }

    if `r_b' != 1 {
        __wfe_restore_bridge_scalars, ///
            hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
            hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
            hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
        display as error "[wfe bridge] Result postback failed: b should be 1xp row vector, got `r_b'x`p_b'"
        exit 498
    }

    * Step 1: Matrix stripe setting
    matrix colnames __wfe_b = `indepvars_names'
    matrix colnames __wfe_V = `indepvars_names'
    matrix rownames __wfe_V = `indepvars_names'
    matrix rownames __wfe_b = `depvar'

    * Step 2: ereturn post (clears all existing e() storage)
    ereturn post __wfe_b __wfe_V, esample(`touse') obs(`=scalar(NT)')

    * Step 3: ereturn local (command metadata)
    local __posted_method "`method'"
    if inlist("`estimator'", "did", "Mdid") {
        local __posted_method "Weighted Two-way"
    }

    local __posted_estimator "`estimator'"
    if "`__posted_estimator'" == "" {
        local __posted_estimator "NULL"
    }

    ereturn local cmd        "wfe"
    ereturn local depvar     "`depvar'"
    ereturn local method     "`__posted_method'"
    ereturn local estimator  "`__posted_estimator'"
    ereturn local properties "b V"
    ereturn local predict    "wfe_predict"
    ereturn local estat_cmd  "wfe_estat"

    * Sign the original score specification
    quietly signestimationsample `depvar' `score_indepvars_spec'
    quietly __wfe_make_value_signature `depvar' `score_indepvars_spec'
    ereturn local datasignature `"`r(datasignature)'"'
    quietly __wfe_make_value_signature `score_indepvars_spec'
    ereturn hidden local _caller_n "`=_N'"
    ereturn hidden local _xb_varspec "`score_indepvars_spec'"
    ereturn hidden local _xb_datasignature "`r(datasignature)'"

    * vcetype from Mata
    local __bridge_vcetype "${__wfe_vcetype}"
    ereturn local vcetype    "`__bridge_vcetype'"

    * qoi storage
    local __posted_qoi "`qoi'"
    if "`unweighted'" != "" {
        local __posted_qoi "unweighted"
    }
    ereturn local qoi "`__posted_qoi'"

    * qoi description
    if "`unweighted'" != "" {
        ereturn local qoi_desc "Unweighted (Standard) Fixed Effect"
    }
    else if "`qoi'" == "ate" {
        ereturn local qoi_desc "ATE (Average Treatment Effect)"
    }
    else if "`qoi'" == "att" {
        ereturn local qoi_desc "ATT (Average Treatment Effect for the Treated)"
    }

    * estimator description
    if "`estimator'" == "" {
        ereturn local estimator_desc "NULL"
    }
    else if "`estimator'" == "fd" {
        ereturn local estimator_desc "FD (First-Difference)"
    }
    else if "`estimator'" == "did" {
        ereturn local estimator_desc "DID (Difference-in-Differences)"
    }
    else if "`estimator'" == "Mdid" {
        ereturn local estimator_desc "DID (Difference-in-Differences) with Matching on Pre-treatment Outcome"
    }

    * Step 4: ereturn scalar (numeric results)
    ereturn scalar N_units    = scalar(__wfe_N_units)
    ereturn scalar N_times    = scalar(__wfe_N_times)
    ereturn scalar df_r       = scalar(__wfe_df_r)
    ereturn scalar sigma      = scalar(__wfe_sigma)
    if "`estimator'" == "Mdid" {
        capture confirm scalar __wfe_maxdev_did
        if _rc {
            display as error "[wfe bridge] Result postback failed: Mdid scalar __wfe_maxdev_did does not exist"
            exit 498
        }
        ereturn scalar maxdev_did = scalar(__wfe_maxdev_did)
    }
    capture confirm scalar __wfe_sigma2
    if !_rc {
        ereturn scalar sigma2 = scalar(__wfe_sigma2)
    }
    ereturn scalar N_nonzero  = scalar(__wfe_N_nonzero)
    capture confirm scalar __wfe_N_negative
    if !_rc {
        ereturn scalar N_negative = scalar(__wfe_N_negative)
    }
    else {
        ereturn scalar N_negative = 0
    }

    * White test results (stored only when White is enabled)
    if "`white'" == "" {
        capture confirm scalar __wfe_white_stat
        if !_rc {
            ereturn scalar white_stat   = scalar(__wfe_white_stat)
            ereturn scalar white_pvalue = scalar(__wfe_white_pvalue)
            ereturn scalar white_alpha  = scalar(white_alpha)

            local __bridge_white_test "${__wfe_white_test}"
            ereturn local  white_test   "`__bridge_white_test'"
        }
    }

    * Step 5: ereturn matrix
    * W: weight matrix T×N
    capture confirm matrix __wfe_W
    if !_rc {
        ereturn matrix W = __wfe_W
    }

    * b_fe, V_fe: standard FE coefficients and covariance (only when White is enabled)
    if "`white'" == "" {
        capture confirm matrix __wfe_b_fe
        if !_rc {
            ereturn matrix b_fe = __wfe_b_fe
        }
        capture confirm matrix __wfe_V_fe
        if !_rc {
            ereturn matrix V_fe = __wfe_V_fe
        }
    }

    * store_wdm: weighted demeaned data (stored only when user-specified)
    if "`store_wdm'" != "" {
        capture confirm matrix __wfe_Y_wdm
        if !_rc {
            ereturn matrix Y_wdm = __wfe_Y_wdm
        }
        capture confirm matrix __wfe_X_wdm
        if !_rc {
            ereturn matrix X_wdm = __wfe_X_wdm
        }
    }

    * Step 6: Formatted output

    * Two-column header (Stata-native style)
    local __wfe_col2 49
    local __wfe_eq   67
    local __wfe_val  69

    * Short display labels (full descriptions stored in e())
    local __qoi_short "`qoi'"
    if "`unweighted'" != "" {
        local __qoi_short "unweighted"
    }
    local __est_short ""
    if "`estimator'" == "fd" {
        local __est_short "FD"
    }
    else if "`estimator'" == "did" {
        local __est_short "DiD"
    }
    else if "`estimator'" == "Mdid" {
        local __est_short "Matched DiD"
    }

    display ""
    display as text "Weighted Fixed Effects Estimation" ///
        _col(`__wfe_col2') "Number of obs"     ///
        _col(`__wfe_eq')   "="                  ///
        _col(`__wfe_val')  as result %9.0g e(N)

    display as text "  Method:      " as result e(method) ///
        _col(`__wfe_col2') as text "Number of units"  ///
        _col(`__wfe_eq')   "="                         ///
        _col(`__wfe_val')  as result %9.0g e(N_units)

    display as text "  Quantity:    " as result "`__qoi_short'" ///
        _col(`__wfe_col2') as text "Time periods"     ///
        _col(`__wfe_eq')   "="                         ///
        _col(`__wfe_val')  as result %9.0g e(N_times)

    if "`estimator'" != "" {
        display as text "  Estimator:   " as result "`__est_short'" ///
            _col(`__wfe_col2') as text "Non-zero wt"      ///
            _col(`__wfe_eq')   "="                         ///
            _col(`__wfe_val')  as result %9.0g e(N_nonzero)

        display ///
            _col(`__wfe_col2') as text "Residual df"      ///
            _col(`__wfe_eq')   "="                         ///
            _col(`__wfe_val')  as result %9.0g e(df_r)
    }
    else {
        display ///
            _col(`__wfe_col2') as text "Non-zero wt"      ///
            _col(`__wfe_eq')   "="                         ///
            _col(`__wfe_val')  as result %9.0g e(N_nonzero)

        display ///
            _col(`__wfe_col2') as text "Residual df"      ///
            _col(`__wfe_eq')   "="                         ///
            _col(`__wfe_val')  as result %9.0g e(df_r)
    }

    display ///
        _col(`__wfe_col2') as text "Neg. weights"     ///
        _col(`__wfe_eq')   "="                         ///
        _col(`__wfe_val')  as result %9.0g e(N_negative)

    display ///
        _col(`__wfe_col2') as text "Sigma"             ///
        _col(`__wfe_eq')   "="                         ///
        _col(`__wfe_val')  as result %9.0g e(sigma)

    * Coefficient table
    ereturn display

    * White test results
    if "`white'" == "" {
        capture confirm scalar e(white_stat)
        if !_rc {
            local __white_test_result "`=e(white_test)'"
            display as text "White (1980) Misspecification Test"
            display as text "  H0: No misspecification (WFE = Standard FE)"
            local __p_display : word count `indepvars_names'
            display as text "  Chi2(" as result `__p_display' as text ")" ///
                _col(15) as text "=" _col(17) as result %9.4f e(white_stat)
            display as text "  P-value" ///
                _col(15) as text "=" _col(17) as result %9.4f e(white_pvalue)
            if "`__white_test_result'" == "TRUE" {
                display as error "  -> Reject H0 at alpha = " %5.3f e(white_alpha)
            }
            else {
                display as text "  -> Fail to reject H0 at alpha = " %5.3f e(white_alpha)
            }
        }
    }

    * Negative weight warning (verbose only)
    if e(N_negative) > 0 & "`verbose'" != "noverbose" {
        display _newline as error ///
            "  Warning: " e(N_negative) " observations have negative weights."
        if inlist("`estimator'", "did", "Mdid") {
            display as text ///
                "  Negative weights are admissible on the weighted two-way DiD path."
        }
        else {
            display as error ///
                "  One-way WFE weights should theoretically be non-negative."
        }
    }

    * Step 7: Cleanup
    foreach __wfe_cleanup_s in __wfe_N_nonzero __wfe_N_negative __wfe_df_r   ///
        __wfe_sigma __wfe_N_units __wfe_N_times __wfe_sigma2 __wfe_J_u       ///
        __wfe_maxdev_did __wfe_white_stat __wfe_white_pvalue                 ///
        N_units N_times NT panel_balanced {                                  
        capture scalar drop `__wfe_cleanup_s'
    }
    capture macro drop __wfe_vcetype __wfe_white_test
    * Restore any caller-visible scalars shadowed by the bridge.
    __wfe_restore_bridge_scalars, ///
        hadwhite(`__had_bridge_white_alpha') savedwhite(`__wfe_saved_bridge_white_alpha') ///
        hadtol(`__had_bridge_tol') savedtol(`__wfe_saved_bridge_tol') ///
        hadp(`__had_bridge_p') savedp(`__wfe_saved_bridge_p')
    * Cleanup intermediate matrices
    foreach __wfe_cleanup_m in __wfe_b_fe __wfe_V_fe __wfe_Psi_wfe __wfe_Psi_fe ///
        __wfe_ginv_XX_tilde __wfe_ginv_XX_hat __wfe_X_tilde __wfe_X_hat      ///
        __wfe_u_tilde __wfe_u_hat {
        capture matrix drop `__wfe_cleanup_m'
    }

end

program define __wfe_restore_bridge_scalars
    version 16

    syntax, HADWHITE(integer) SAVEDWHITE(name) HADTOL(integer) SAVEDTOL(name) HADP(integer) SAVEDP(name)

    if `hadwhite' {
        capture scalar white_alpha = scalar(`savedwhite')
    }
    else {
        capture scalar drop white_alpha
    }

    if `hadtol' {
        capture scalar tol = scalar(`savedtol')
    }
    else {
        capture scalar drop tol
    }

    if `hadp' {
        capture scalar p = scalar(`savedp')
    }
    else {
        capture scalar drop p
    }

    foreach __wfe_rs_s in N_units N_times NT panel_balanced {
        capture scalar drop `__wfe_rs_s'
    }
end

program define __wfe_ensure_backend_mata
    version 16

    capture mata: mata describe _wfe_oneway_estimate()
    local need_oneway = _rc
    capture mata: mata describe _wfe_twoway_estimate()
    local need_twoway = _rc
    capture mata: mata describe wfe_store_weight_summary()
    local need_post = _rc

    if `need_oneway' == 0 & `need_twoway' == 0 & `need_post' == 0 {
        __wfe_load_vshim
        exit 0
    }

    quietly findfile wfe.ado
    if _rc {
        display as error "wfe backend could not locate its ado source path"
        exit 601
    }

    local ado_dir = subinstr(`"`r(fn)'"', "wfe.ado", "", .)

    capture quietly do `"`ado_dir'mata/wfe_utils.mata"'
    capture quietly do `"`ado_dir'mata/wfe_weights_unit.mata"'
    capture quietly do `"`ado_dir'mata/wfe_weights_time.mata"'
    capture quietly do `"`ado_dir'mata/wfe_weights_fd.mata"'
    capture quietly do `"`ado_dir'mata/wfe_weights_did.mata"'
    capture quietly do `"`ado_dir'mata/wfe_weights_mdid.mata"'
    capture quietly do `"`ado_dir'mata/wfe_transform.mata"'
    capture quietly do `"`ado_dir'mata/wfe_wwdemean.mata"'
    capture quietly do `"`ado_dir'mata/wfe_ols.mata"'
    capture quietly do `"`ado_dir'mata/wfe_se_hac.mata"'
    capture quietly do `"`ado_dir'mata/wfe_se_pwfe.mata"'
    capture quietly do `"`ado_dir'mata/wfe_twoway_fe_ols.mata"'
    capture quietly do `"`ado_dir'mata/wfe_se_fe_twoway.mata"'
    capture quietly do `"`ado_dir'mata/wfe_white_test.mata"'
    capture quietly do `"`ado_dir'mata/wfe_pwfe.mata"'
    capture quietly do `"`ado_dir'mata/wfe_se_gmm.mata"'
    capture quietly do `"`ado_dir'mata/wfe_complex_project.mata"'
    capture quietly do `"`ado_dir'mata/wfe_bridge.mata"'
    capture quietly do `"`ado_dir'mata/wfe_postestimation.mata"'
    capture quietly do `"`ado_dir'mata/wfe_oneway.mata"'
    capture quietly do `"`ado_dir'mata/wfe_twoway.mata"'
    __wfe_load_vshim
end

program define __wfe_make_value_signature, rclass
    version 16

    syntax anything(name=varspec)

    local normalized_vars
    quietly fvrevar `varspec' if e(sample)
    local source_vars `"`r(varlist)'"'
    foreach source_var of local source_vars {
        tempvar normalized_var
        quietly generate double `normalized_var' = `source_var' if e(sample)
        local normalized_vars `normalized_vars' `normalized_var'
    }

    quietly _datasignature `normalized_vars', esample nodefault nonames
    return local datasignature `"`r(datasignature)'"'
end

program define __wfe_load_vshim
    version 16

    capture mata: mata describe _wfe_as_colvector()
    local need_ascol = _rc
    capture mata: mata describe __wfe_vshim_version()
    local need_token = _rc
    if !`need_ascol' & !`need_token' {
        exit 0
    }

    quietly findfile __wfe_vshim.ado
    if _rc {
        display as error "wfe backend could not locate __wfe_vshim.ado"
        exit 601
    }

    capture program drop __wfe_vshim
    quietly run `"`r(fn)'"'
end
