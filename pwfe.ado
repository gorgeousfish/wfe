*! pwfe.ado -- Propensity-score weighted fixed effects estimator

program define pwfe, eclass sortpreserve
    version 16

    /*
       Clear stale estimates before syntax so Stata-level parser errors
       cannot leave an old e() result replayable via predict/estat.
    */
    ereturn clear

    syntax [varlist(numeric fv default=none)] [if] [in], ///
        TReat(varname) OUTcome(varname) UNIT(varname)         ///
        [TIME(varname) PSCore(varname) Method(string) ///
         noWITHIN_unit QOI(string) ESTIMator(string) ///
         CIT(varname) HETERO_se(string) AUTO_se(string) ///
         UNBIASED_se noWHITE WHITE_alpha(string) ///
         noVERBose DIAGnose]

    marksample touse_formula
    marksample touse_ifin, novarlist
    tempvar touse
    quietly gen byte `touse' = `touse_ifin'

    /*
       Validate cit() before any transformed-outcome processing.
       Invalid inputs should fail-fast on the caller-selected sample.
    */
    if "`cit'" != "" {
        capture confirm numeric variable `cit'
        if _rc {
            display as error "'C.it' must be a numeric vector with length equal to number of observations"
            exit 198
        }
        quietly count if `touse_ifin' & missing(`cit')
        if r(N) > 0 {
            display as error "missing values in cit() are not allowed"
            exit 198
        }
        quietly count if `touse_ifin' & `cit' < 0
        if r(N) > 0 {
            display as error "'C.it' must be a non-negative numeric vector"
            exit 198
        }
    }

    /*
       pscore() is a user-supplied estimator input, not a formula covariate.
       Validate it before any transformed-outcome work.
    */
    if "`pscore'" != "" {
        capture confirm numeric variable `pscore'
        if _rc {
            display as error "'pscore' must be a numeric vector with length equal to number of observations"
            exit 198
        }
        quietly count if `touse_ifin' & missing(`pscore')
        if r(N) > 0 {
            display as error "missing values in pscore() are not allowed"
            exit 198
        }
        quietly count if `touse_ifin' & (`pscore' <= 0 | `pscore' >= 1)
        if r(N) > 0 {
            display as error "'pscore' must be strictly between 0 and 1"
            exit 198
        }
    }

    capture confirm numeric variable `outcome'
    if _rc {
        display as error "outcome() must specify a numeric variable"
        exit 198
    }

    capture confirm numeric variable `treat'
    if _rc {
        display as error "treat() must specify a numeric variable"
        exit 198
    }

    /*
       unit()/time() define panel dyads rather than ordinary covariates.
       Missing identifiers must fail-fast instead of being silently dropped.
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
    }

    /*
       outcome() and treat() are formula-level covariates. Rows where any
       formula variable is missing are silently dropped before validation.
       Match that behavior by markout first, then guard on the formula-valid
       sample.
    */
    markout `touse' `outcome' `treat'

    quietly count if `touse' & missing(`outcome')
    if r(N) > 0 {
        display as error "missing values in outcome() are not allowed"
        exit 198
    }

    quietly count if `touse' & missing(`treat')
    if r(N) > 0 {
        display as error "missing values in treat() are not allowed"
        exit 198
    }
    local pscore_design_vars ""

    if "`method'" == "" {
        local method "unit"
    }
    if "`qoi'" == "" {
        local qoi "ate"
    }
    if "`hetero_se'" == "" {
        local hetero_se "on"
    }
    if "`auto_se'" == "" {
        local auto_se "on"
    }
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

    local within_unit_on = ("`within_unit'" == "")
    local white_on = ("`white'" == "")
    local verbose_on = ("`verbose'" == "")
    local white_state = cond(`white_on', "on", "off")
    local verbose_state = cond(`verbose_on', "on", "off")
    local unbiased_state = cond("`unbiased_se'" != "", "on", "off")

    if "`time'" == "" & "`method'" == "time" {
        display as error "'time.index' should be provided"
        exit 198
    }

    if "`method'" == "time" & "`estimator'" == "fd" {
        display as error "First Difference is not compatible with 'time' method"
        exit 198
    }

    if "`time'" == "" & "`estimator'" == "fd" {
        display as error "First Difference cannot calculate when 'time.index' is missing"
        exit 198
    }

    if "`estimator'" == "did" {
        display as error "Difference-in-Differences is not compatible with pwfe"
        exit 198
    }

    if "`estimator'" == "Mdid" {
        display as error "pwfe does not support estimator('Mdid')"
        exit 198
    }

    if "`estimator'" != "" & "`estimator'" != "fd" {
        display as error "estimator() must be 'fd'"
        exit 198
    }

    if "`pscore'" != "" & "`varlist'" != "" {
        display as error "'formula' should not be provided when pscore is specified"
        exit 198
    }

    if "`varlist'" != "" {
        /*
           The internal propensity-score fit is no-constant. Bare i. factors
           therefore need no-constant expansion (ibn.), while user-specified
           base/omitted modifiers must remain authoritative.
        */
        local pscore_fvspec = subinstr("`varlist'", "i.", "ibn.", .)

        fvexpand `pscore_fvspec' if `touse'
        local pscore_fvexpanded `"`r(varlist)'"'

        fvrevar `pscore_fvspec' if `touse'
        local pscore_fvvars `"`r(varlist)'"'

        local pscore_design_vars
        local __pwfe_n_terms : word count `pscore_fvexpanded'
        forvalues __pwfe_i = 1/`__pwfe_n_terms' {
            local __pwfe_term : word `__pwfe_i' of `pscore_fvexpanded'
            local __pwfe_var : word `__pwfe_i' of `pscore_fvvars'
            /*
               Preserve explicit base/omitted semantics after no-constant
               expansion so factor-variable notation stays equivalent to the
               corresponding hand-written dummy/interaction design.
            */
            if strpos("`__pwfe_term'", "b.") == 0 & strpos("`__pwfe_term'", "o.") == 0 {
                local pscore_design_vars `"`pscore_design_vars' `__pwfe_var'"'
            }
        }
        markout `touse' `pscore_design_vars'
    }

    quietly count if `touse'
    local N_parse = r(N)
    if `N_parse' == 0 {
        display as error "no observations"
        exit 2000
    }

    if "`time'" != "" {
        /*
           Duplicate unit-time dyads are structural errors for the actual
           estimation sample. Rows already excluded by the propensity-score
           formula should not redefine panel geometry.
        */
        tempvar pwfe_dup
        quietly duplicates tag `unit' `time' if `touse', gen(`pwfe_dup')
        quietly count if `touse' & `pwfe_dup' > 0
        if r(N) > 0 {
            display as error "duplicate observations found for unit-time combination"
            exit 198
        }
        drop `pwfe_dup'
    }

    /*
       treat() remains a structural estimator input. Validate binary support
       on the caller-selected sample before any formula-based covariate drops;
       bad treatment codings must not be masked by missing propensity-score
       covariates.
    */
    quietly levelsof `treat' if `touse_ifin', local(treat_levels)
    local treat_level_count : word count `treat_levels'
    if `treat_level_count' > 2 {
        display as error "'treat' must be a binary vector"
        exit 198
    }

    local has_zero = 0
    local has_one = 0
    foreach level of local treat_levels {
        if real("`level'") == 0 {
            local has_zero = 1
        }
        else if real("`level'") == 1 {
            local has_one = 1
        }
    }
    if `has_zero' == 0 | `has_one' == 0 {
        if `treat_level_count' == 1 & (`has_zero' | `has_one') {
            display as error "'treat' must contain both 0 and 1 values"
            exit 198
        }
        display as error "'treat' must be either 0 or 1 where 1 indicates treated"
        exit 198
    }

    if "`method'" != "unit" & "`method'" != "time" {
        display as error "method should be either unit or time"
        exit 198
    }

    if "`qoi'" != "ate" & "`qoi'" != "att" {
        display as error "qoi should be either ate or att"
        exit 198
    }

    if "`hetero_se'" != "on" & "`hetero_se'" != "off" {
        display as error "hetero_se() must be 'on' or 'off'"
        exit 198
    }

    if "`auto_se'" != "on" & "`auto_se'" != "off" {
        display as error "auto_se() must be 'on' or 'off'"
        exit 198
    }

    if "`unbiased_se'" != "" {
        if "`hetero_se'" != "on" | "`auto_se'" != "off" {
            display as error "unbiased_se is allowed only for hetero_se(on) auto_se(off)"
            exit 198
        }
    }

    if `white_alpha' <= 0 | `white_alpha' >= 1 {
        display as error "white_alpha() must lie in (0,1)"
        exit 198
    }

    if "`hetero_se'" == "off" & "`auto_se'" == "on" {
        display as error "robust standard errors with autocorrelation and homoskedasticity is not supported"
        exit 198
    }

    if "`diagnose'" != "" {
        * Diagnose is parse-only; clear stale estimation results.
        ereturn clear
        local within_unit_state = cond(`within_unit_on', "on", "off")
        local white_state = cond(`white_on', "on", "off")
        local verbose_state = cond(`verbose_on', "on", "off")
        local unbiased_state = cond("`unbiased_se'" != "", "on", "off")
        local estimator_state = cond("`estimator'" == "", "NULL", "`estimator'")

        display as text "method      = `method'"
        display as text "qoi         = `qoi'"
        display as text "estimator   = `estimator_state'"
        display as text "within_unit = `within_unit_state'"
        display as text "hetero_se   = `hetero_se'"
        display as text "auto_se     = `auto_se'"
        display as text "unbiased   = `unbiased_state'"
        display as text "white       = `white_state'"
        display as text "white_alpha = `white_alpha'"
        display as text "verbose     = `verbose_state'"
        display as text "pscore      = `pscore'"
        display as text "cit         = `cit'"
        display as text "N_parse     = `N_parse'"
        exit 0
    }

    __pwfe_ensure_backend_mata

    tempvar y_star fitted_pscore pool_id pool_touse pwfe_unit_idx pwfe_time_idx pwfe_cit pwfe_obs_order
    quietly gen double `y_star' = .
    quietly gen double `fitted_pscore' = .
    quietly gen double `pwfe_cit' = 1
    if "`cit'" != "" {
        quietly replace `pwfe_cit' = `cit'
    }

    if "`pscore'" != "" {
        __pwfe_apply_transform, ///
            touse(`touse') outcome(`outcome') treat(`treat') ///
            pscore(`pscore') outvar(`y_star')
        local pscore_source "user"
        local transform_scope "global"
    }
    else if `within_unit_on' & "`method'" == "time" {
        quietly egen long `pool_id' = group(`time') if `touse'
        quietly gen byte `pool_touse' = 0
        quietly levelsof `pool_id' if `touse', local(pool_ids)

        foreach pool of local pool_ids {
            quietly replace `pool_touse' = `touse' & `pool_id' == `pool'
            __pwfe_estimate_pscore `pscore_design_vars', ///
                touse(`pool_touse') treat(`treat') outvar(`fitted_pscore') ///
                poollabel("time pool")
            __pwfe_apply_transform, ///
                touse(`pool_touse') outcome(`outcome') treat(`treat') ///
                pscore(`fitted_pscore') outvar(`y_star')
        }

        local pscore_source "estimated"
        local transform_scope "time"
    }
    else if `within_unit_on' & "`method'" == "unit" {
        quietly egen long `pool_id' = group(`unit') if `touse'
        quietly gen byte `pool_touse' = 0
        quietly levelsof `pool_id' if `touse', local(pool_ids)

        foreach pool of local pool_ids {
            quietly replace `pool_touse' = `touse' & `pool_id' == `pool'
            __pwfe_estimate_pscore `pscore_design_vars', ///
                touse(`pool_touse') treat(`treat') outvar(`fitted_pscore') ///
                poollabel("unit pool")
            __pwfe_apply_transform, ///
                touse(`pool_touse') outcome(`outcome') treat(`treat') ///
                pscore(`fitted_pscore') outvar(`y_star')
        }

        local pscore_source "estimated"
        local transform_scope "unit"
    }
    else {
        __pwfe_estimate_pscore `pscore_design_vars', ///
            touse(`touse') treat(`treat') outvar(`fitted_pscore') ///
            poollabel("global sample")
        __pwfe_apply_transform, ///
            touse(`touse') outcome(`outcome') treat(`treat') ///
            pscore(`fitted_pscore') outvar(`y_star')
        local pscore_source "estimated"
        local transform_scope "global"
    }

    quietly gen long `pwfe_obs_order' = _n
    quietly egen long `pwfe_unit_idx' = group(`unit') if `touse'
    if "`time'" != "" {
        quietly egen long `pwfe_time_idx' = group(`time') if `touse'
    }
    else {
        sort `unit' `pwfe_obs_order'
        quietly by `unit': gen long `pwfe_time_idx' = sum(`touse') if `touse'
        sort `pwfe_obs_order'
    }

    tempname y_star_mat pscore_mat y_star_unit_idx_mat y_star_time_idx_mat
    mkmat `y_star' if `touse', matrix(`y_star_mat')
    if "`pscore'" != "" {
        mkmat `pscore' if `touse', matrix(`pscore_mat')
    }
    else {
        mkmat `fitted_pscore' if `touse', matrix(`pscore_mat')
    }
    mkmat `pwfe_unit_idx' if `touse', matrix(`y_star_unit_idx_mat')
    mkmat `pwfe_time_idx' if `touse', matrix(`y_star_time_idx_mat')
    if "`time'" != "" {
        sort `unit' `time'
    }
    else {
        sort `unit' `pwfe_obs_order'
    }
    quietly summarize `pwfe_unit_idx' if `touse', meanonly
    scalar __pwfe_J_u = r(max)
    quietly summarize `pwfe_time_idx' if `touse', meanonly
    scalar __pwfe_J_t = r(max)

    if "`method'" == "unit" {
        scalar __pwfe_unit_number = scalar(__pwfe_J_u)
    }
    else {
        scalar __pwfe_unit_number = scalar(__pwfe_J_t)
    }
    scalar __pwfe_white_alpha = `white_alpha'

    local pwfe_y_star `y_star'
    local pwfe_cit `pwfe_cit'
    local pwfe_unit_idx `pwfe_unit_idx'
    local pwfe_time_idx `pwfe_time_idx'
    local pwfe_white `white_state'
    local pwfe_verbose `verbose_state'
    local pwfe_unbiased_se `unbiased_state'

    capture noisily mata: _wfe_pwfe_estimate()
    if _rc {
        local __pwfe_runtime_rc = _rc
        __pwfe_cleanup_ns
        exit `__pwfe_runtime_rc'
    }

    matrix colnames __pwfe_b = `treat'
    matrix rownames __pwfe_b = `outcome'
    matrix colnames __pwfe_V = `treat'
    matrix rownames __pwfe_V = `treat'
    matrix colnames __pwfe_b_fe = `treat'
    matrix rownames __pwfe_b_fe = `outcome'
    matrix colnames __pwfe_V_fe = `treat'
    matrix rownames __pwfe_V_fe = `treat'
    mata: st_numscalar("__pwfe_N_nonzero", sum(st_matrix("__pwfe_W") :!= 0))

    local estimator_label "NULL"
    if "`estimator'" != "" {
        local estimator_label "`estimator'"
    }
    local qoi_desc "ATE (Average Treatment Effect)"
    if "`qoi'" == "att" {
        local qoi_desc "ATT (Average Treatment Effect for the Treated)"
    }
    local estimator_desc "NULL"
    if "`estimator'" == "fd" {
        local estimator_desc "FD (First-Difference)"
    }

    ereturn post __pwfe_b __pwfe_V, esample(`touse') obs(`N_parse')
    ereturn local cmd "pwfe"
    ereturn local depvar "`outcome'"
    ereturn local method "`method'"
    ereturn local qoi "`qoi'"
    ereturn local qoi_desc "`qoi_desc'"
    ereturn local estimator "`estimator_label'"
    ereturn local estimator_desc "`estimator_desc'"
    ereturn local properties "b V"
    ereturn local predict "wfe_predict"
    ereturn local estat_cmd "wfe_estat"
    ereturn local _unitvar "`unit'"
    ereturn local _timevar "`time'"
    ereturn local _y_star_var ""
    ereturn local vcetype "${__pwfe_vcetype}"
    ereturn local pscore_source "`pscore_source'"
    ereturn local transform_scope "`transform_scope'"
    quietly signestimationsample `outcome' `treat'
    quietly __pwfe_make_value_signature `outcome' `treat'
    ereturn local datasignature `"`r(datasignature)'"'
    quietly __pwfe_make_value_signature `treat'
    ereturn hidden local _caller_n "`=_N'"
    ereturn hidden local _xb_datasignature "`r(datasignature)'"
    ereturn hidden local _treatvar "`treat'"
    ereturn scalar N_units = scalar(__pwfe_N_units)
    ereturn scalar N_times = scalar(__pwfe_N_times)
    ereturn scalar df_r = scalar(__pwfe_df_r)
    ereturn scalar sigma = scalar(__pwfe_sigma)
    ereturn scalar sigma2 = scalar(__pwfe_sigma2)
    ereturn scalar N_nonzero = scalar(__pwfe_N_nonzero)
    ereturn matrix W = __pwfe_W
    ereturn matrix b_fe = __pwfe_b_fe
    ereturn matrix V_fe = __pwfe_V_fe
    if `white_on' {
        ereturn scalar white_stat = scalar(__pwfe_white_stat)
        ereturn scalar white_pvalue = scalar(__pwfe_white_pvalue)
        ereturn scalar white_alpha = scalar(__pwfe_white_alpha)
        ereturn local white_test "${__pwfe_white_test}"
    }
    ereturn matrix y_star = `y_star_mat'
    ereturn matrix pscore = `pscore_mat'
    ereturn matrix _y_star_unit_idx = `y_star_unit_idx_mat'
    ereturn matrix _y_star_time_idx = `y_star_time_idx_mat'

    * Two-column header (Stata-native style)
    local __pwfe_col2 49
    local __pwfe_eq   67
    local __pwfe_val  69

    * Short display labels (full descriptions stored in e())
    local __qoi_short "`qoi'"
    local __est_short ""
    if "`estimator'" == "fd" {
        local __est_short "FD"
    }

    display ""
    display as text "Propensity-Score Weighted FE Estimation" ///
        _col(`__pwfe_col2') "Number of obs"     ///
        _col(`__pwfe_eq')   "="                  ///
        _col(`__pwfe_val')  as result %9.0g e(N)

    display as text "  Method:      " as result e(method) ///
        _col(`__pwfe_col2') as text "Number of units"  ///
        _col(`__pwfe_eq')   "="                         ///
        _col(`__pwfe_val')  as result %9.0g e(N_units)

    display as text "  Quantity:    " as result "`__qoi_short'" ///
        _col(`__pwfe_col2') as text "Time periods"     ///
        _col(`__pwfe_eq')   "="                         ///
        _col(`__pwfe_val')  as result %9.0g e(N_times)

    if "`estimator'" != "" {
        display as text "  Estimator:   " as result "`__est_short'" ///
            _col(`__pwfe_col2') as text "Non-zero wt"      ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(N_nonzero)

        display as text "  P-score:     " as result e(pscore_source) ///
            _col(`__pwfe_col2') as text "Residual df"      ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(df_r)

        display as text "  Transform:   " as result e(transform_scope) ///
            _col(`__pwfe_col2') as text "Sigma"             ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(sigma)
    }
    else {
        display as text "  P-score:     " as result e(pscore_source) ///
            _col(`__pwfe_col2') as text "Non-zero wt"      ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(N_nonzero)

        display as text "  Transform:   " as result e(transform_scope) ///
            _col(`__pwfe_col2') as text "Residual df"      ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(df_r)

        display ///
            _col(`__pwfe_col2') as text "Sigma"             ///
            _col(`__pwfe_eq')   "="                         ///
            _col(`__pwfe_val')  as result %9.0g e(sigma)
    }

    * Coefficient table
    ereturn display

    * White test results
    if `white_on' {
        capture confirm scalar e(white_stat)
        if !_rc {
            local __white_test_result "`=e(white_test)'"
            display as text "White (1980) Misspecification Test"
            display as text "  H0: No misspecification (WFE = Standard FE)"
            display as text "  Chi2(" as result 1 as text ")" ///
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

    __pwfe_cleanup_ns
end

program define __pwfe_ensure_backend_mata
    version 16

    capture mata: mata describe wfe_transform()
    local need_transform = _rc

    capture mata: mata describe _wfe_pwfe_apply_transform()
    local need_bridge = _rc
    capture mata: mata describe _wfe_pwfe_estimate()
    local need_estimate = _rc

    if `need_transform' == 0 & `need_bridge' == 0 & `need_estimate' == 0 {
        __pwfe_load_vshim
        exit 0
    }

    quietly findfile pwfe.ado
    if _rc {
        display as error "pwfe backend could not locate its ado source path"
        exit 601
    }

    local ado_dir = subinstr(`"`r(fn)'"', "pwfe.ado", "", .)

    capture mata: mata describe wfe_build_exist()
    local need_utils = _rc
    capture mata: mata describe wfe_weights_unit()
    local need_weights_unit = _rc
    capture mata: mata describe wfe_weights_time()
    local need_weights_time = _rc
    capture mata: mata describe wfe_weights_fd()
    local need_weights_fd = _rc
    capture mata: mata describe wfe_wwdemean_matrix()
    local need_wwdemean = _rc
    capture mata: mata describe _wfe_ols_core()
    local need_ols = _rc
    capture mata: mata describe _wfe_compute_se()
    local need_se_hac = _rc
    capture mata: mata describe _wfe_pwfe_compute_se()
    local need_se_pwfe = _rc
    capture mata: mata describe wfe_white_test_pwfe()
    local need_white = _rc

    if `need_transform' {
        if `need_utils' {
            quietly do `"`ado_dir'mata/wfe_utils.mata"'
            local need_utils 0
        }
        quietly do `"`ado_dir'mata/wfe_transform.mata"'
    }

    if `need_bridge' {
        __pwfe_load_pscore_bridge `"`ado_dir'"'
    }

    if `need_estimate' {
        if `need_utils' {
            quietly do `"`ado_dir'mata/wfe_utils.mata"'
        }
        if `need_weights_unit' {
            quietly do `"`ado_dir'mata/wfe_weights_unit.mata"'
        }
        if `need_weights_time' {
            quietly do `"`ado_dir'mata/wfe_weights_time.mata"'
        }
        if `need_weights_fd' {
            quietly do `"`ado_dir'mata/wfe_weights_fd.mata"'
        }
        if `need_wwdemean' {
            quietly do `"`ado_dir'mata/wfe_wwdemean.mata"'
        }
        if `need_ols' {
            quietly do `"`ado_dir'mata/wfe_ols.mata"'
        }
        if `need_se_hac' {
            quietly do `"`ado_dir'mata/wfe_se_hac.mata"'
        }
        if `need_se_pwfe' {
            quietly do `"`ado_dir'mata/wfe_se_pwfe.mata"'
        }
        if `need_white' {
            quietly do `"`ado_dir'mata/wfe_white_test.mata"'
        }
        quietly do `"`ado_dir'mata/wfe_pwfe.mata"'
    }
    __pwfe_load_vshim
end

program define __pwfe_load_pscore_bridge
    version 16
    args ado_dir

    /*
       Mata identifiers longer than 32 characters do not compile.  The
       shipped bridge source contains one overlong helper name, so load a
       temporary patched copy here instead of mutating the read-only .mata
       source in place.
    */
    capture mata: mata drop _wfe_pwfe_invlogit()
    capture mata: mata drop _wfe_pwfe_logit_nll()
    capture mata: mata drop _wfe_pwfe_unique_count()
    capture mata: mata drop _wfe_pwfe_prior_scale_core()
    capture mata: mata drop _wfe_pwfe_prior_scale()
    capture mata: mata drop _wfe_pwfe_prior_scale_int()
    capture mata: mata drop _wfe_pwfe_safe_penalty_value()
    capture mata: mata drop _wfe_pwfe_penlogit_prob_core()
    capture mata: mata drop _wfe_pwfe_penlogit_prob()
    capture mata: mata drop _wfe_pwfe_penlogit_prob_int()
    capture mata: mata drop _wfe_pwfe_bridge_singlevar()
    capture mata: mata drop _wfe_pwfe_bridge_numeric_singlevar()
    capture mata: mata drop _wfe_pwfe_bridge_numvar()
    capture mata: mata drop _wfe_pwfe_apply_transform()

    tempfile __pwfe_bridge_patched
    quietly filefilter `"`ado_dir'mata/wfe_pscore_bridge.mata"' ///
        `"`__pwfe_bridge_patched'"', ///
        from("_wfe_pwfe_bridge_numeric_singlevar") ///
        to("_wfe_pwfe_bridge_numvar") replace
    quietly do `"`__pwfe_bridge_patched'"'
end

program define __pwfe_load_vshim
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
        display as error "pwfe backend could not locate __wfe_vshim.ado"
        exit 601
    }

    capture program drop __wfe_vshim
    quietly run `"`r(fn)'"'
end

program define __pwfe_estimate_pscore
    version 16
    syntax [varlist(numeric default=none)], TOUSE(varname) TREAT(varname) OUTVAR(varname) ///
        POOLLABEL(string)

    /*
       Internal propensity-score logit fits must not leak into
       the caller's active estimation results on failure.
    */
    local __had_estimates = ("`e(cmd)'" != "")
    tempname __pwfe_saved_est
    if `__had_estimates' {
        quietly estimates store `__pwfe_saved_est'
    }

    local __pwfe_ps_ok 0
    quietly replace `outvar' = . if `touse'
    tempname x_mat y_mat pscore_mat
    tempvar __pwfe_const

    /*
       Every internal propensity-score fit uses a penalized binomial IRLS
       with an adaptive Cauchy prior, not plain logit. This applies to all
       pooled and split-sample paths so transformed outcomes align with the
       same weakly regularized fit even on non-separated samples.
    */
    if "`varlist'" == "" {
        quietly gen double `__pwfe_const' = 1 if `touse'
        mkmat `__pwfe_const' if `touse', matrix(`x_mat')
    }
    else {
        mkmat `varlist' if `touse', matrix(`x_mat')
    }
    mkmat `treat' if `touse', matrix(`y_mat')

    if "`varlist'" == "" {
        capture mata: st_matrix("`pscore_mat'", ///
            _wfe_pwfe_penlogit_prob_int(st_matrix("`x_mat'"), st_matrix("`y_mat'"), 2.5))
    }
    else {
        capture mata: st_matrix("`pscore_mat'", ///
            _wfe_pwfe_penlogit_prob(st_matrix("`x_mat'"), st_matrix("`y_mat'"), 2.5))
    }
    if _rc == 0 {
        mata: st_store(selectindex(st_data(., "`touse'") :!= 0), ///
            "`outvar'", st_matrix("`pscore_mat'"))
        quietly count if `touse' & missing(`outvar')
        if r(N) == 0 {
            quietly count if `touse' & (`outvar' <= 0 | `outvar' >= 1)
            if r(N) == 0 {
                local __pwfe_ps_ok 1
            }
        }
    }

    if `__pwfe_ps_ok' == 0 {
        if `__had_estimates' {
            quietly estimates restore `__pwfe_saved_est'
            quietly estimates drop `__pwfe_saved_est'
        }
        else {
            ereturn clear
        }
        display as error "pwfe propensity-score logit failed in `poollabel'; consider providing pscore()"
        exit 498
    }

    if `__had_estimates' {
        quietly estimates restore `__pwfe_saved_est'
        quietly estimates drop `__pwfe_saved_est'
    }
    else {
        ereturn clear
    }
end

program define __pwfe_apply_transform
    version 16
    syntax, TOUSE(varname) OUTCOME(varname) TREAT(varname) PSCORE(varname) ///
        OUTVAR(varname)

    quietly count if `touse'
    if r(N) == 0 {
        display as error "pwfe transform bridge requires at least one selected observation"
        exit 3200
    }

    tempname y_mat treat_mat pscore_mat y_star_mat
    mkmat `outcome' if `touse', matrix(`y_mat')
    mkmat `treat' if `touse', matrix(`treat_mat')
    mkmat `pscore' if `touse', matrix(`pscore_mat')
    mata: st_matrix("`y_star_mat'", ///
        wfe_transform(st_matrix("`y_mat'"), st_matrix("`treat_mat'"), st_matrix("`pscore_mat'")))
    mata: st_store(selectindex(st_data(., "`touse'") :!= 0), "`outvar'", st_matrix("`y_star_mat'"))
end

program define __pwfe_cleanup_ns
    version 16

    capture scalar drop __pwfe_J_u __pwfe_J_t __pwfe_unit_number __pwfe_white_alpha
    capture scalar drop __pwfe_N_nonzero __pwfe_N_units __pwfe_N_times __pwfe_df_r
    capture scalar drop __pwfe_sigma __pwfe_sigma2 __pwfe_white_stat __pwfe_white_pvalue
    capture macro drop __pwfe_vcetype __pwfe_white_test
    capture matrix drop __pwfe_b __pwfe_V __pwfe_b_fe __pwfe_V_fe __pwfe_W
end

program define __pwfe_make_value_signature, rclass
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
