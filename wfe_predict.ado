*! wfe_predict.ado - Prediction program for wfe and pwfe commands

program define wfe_predict
    version 16

    syntax newvarname [if] [in] [, XB FITted RESiduals]

    if !inlist("`e(cmd)'", "wfe", "pwfe") {
        error 301
    }

    local want_xb = ("`xb'" != "") + ("`fitted'" != "")
    local want_residuals = ("`residuals'" != "")

    if `want_xb' > 1 {
        display as error "xb and fitted may not be combined"
        exit 198
    }

    if `want_xb' == 0 & `want_residuals' == 0 {
        local want_xb = 1
    }

    if `want_residuals' & `want_xb' {
        display as error "xb/fitted and residuals may not be combined"
        exit 198
    }

    /*
       The hidden caller-side row-count key is retained only as metadata.
       Replay validity is governed by e(sample) plus the statistic-specific
       variable/signature guards below, not by this historical snapshot.
    */

    quietly count if e(sample)
    if r(N) != e(N) {
        display as error "wfe_predict: current data no longer match stored estimation sample"
        exit 498
    }

    marksample touse, novarlist
    quietly replace `touse' = 0 if !e(sample)
    local depvar : word 1 of `e(depvar)'
    local data_signature `"`e(datasignature)'"'

    if `want_residuals' & "`depvar'" == "" {
        display as error "wfe_predict: stored dependent variable name is unavailable"
        exit 498
    }

    if "`e(cmd)'" == "wfe" {
        local xb_varspec "`e(_xb_varspec)'"
        local xb_datasignature "`e(_xb_datasignature)'"

        if `want_residuals' {
            capture __wfe_predict_confirm_exact_var `depvar'
            if _rc {
                display as error "wfe_predict: stored dependent variable is no longer available in current data"
                exit 498
            }
            if "`xb_varspec'" == "" {
                display as error "wfe_predict: stored score-design specification is unavailable"
                exit 498
            }
            if "`xb_datasignature'" == "" {
                display as error "wfe_predict: stored score-design signature is unavailable"
                exit 498
            }
            if `"`data_signature'"' == "" {
                display as error "wfe_predict: stored data-signature is unavailable"
                exit 498
            }
            capture quietly __wfe_pred_valsig `depvar' `xb_varspec'
            if _rc {
                display as error "wfe_predict: current score variables no longer match the stored estimation sample"
                exit 498
            }
            if `"`r(datasignature)'"' != `"`e(datasignature)'"' {
                display as error "wfe_predict: current score variables no longer match the stored estimation sample"
                exit 498
            }
        }
        else {
            if "`xb_varspec'" == "" {
                display as error "wfe_predict: stored score-design specification is unavailable"
                exit 498
            }
            if "`xb_datasignature'" == "" {
                display as error "wfe_predict: stored score-design signature is unavailable"
                exit 498
            }
            capture quietly __wfe_pred_valsig `xb_varspec'
            if _rc == 111 {
                display as error "wfe_predict: regressors from estimation are no longer available in current data"
                exit 498
            }
            if _rc {
                display as error "wfe_predict: current score variables no longer match the stored estimation sample"
                exit 498
            }
            if `"`r(datasignature)'"' != `"`xb_datasignature'"' {
                display as error "wfe_predict: current score variables no longer match the stored estimation sample"
                exit 498
            }
        }
    }
    else if "`e(cmd)'" == "pwfe" {
        local treatvar "`e(_treatvar)'"
        if "`treatvar'" == "" {
            display as error "wfe_predict: stored treatment variable name is unavailable"
            exit 498
        }

        if `want_residuals' {
            capture __wfe_predict_confirm_exact_var `depvar'
            if _rc {
                display as error "wfe_predict: stored dependent variable is no longer available in current data"
                exit 498
            }
            capture __wfe_predict_confirm_exact_var `treatvar'
            if _rc {
                display as error "wfe_predict: stored treatment variable is no longer available in current data"
                exit 498
            }
            if `"`data_signature'"' == "" {
                display as error "wfe_predict: stored data-signature is unavailable"
                exit 498
            }
            capture quietly __wfe_pred_valsig `depvar' `treatvar'
            if _rc {
                display as error "wfe_predict: current outcome/treat variables no longer match the stored estimation sample"
                exit 498
            }
            if `"`r(datasignature)'"' != `"`e(datasignature)'"' {
                display as error "wfe_predict: current outcome/treat variables no longer match the stored estimation sample"
                exit 498
            }
        }
        else {
            capture __wfe_predict_confirm_exact_var `treatvar'
            if _rc {
                display as error "wfe_predict: stored treatment variable is no longer available in current data"
                exit 498
            }
            local xb_datasignature "`e(_xb_datasignature)'"
            if "`xb_datasignature'" == "" {
                display as error "wfe_predict: stored treatment-only signature is unavailable"
                exit 498
            }
            capture quietly __wfe_pred_valsig `treatvar'
            if _rc {
                display as error "wfe_predict: current treatment variable no longer matches the stored estimation sample"
                exit 498
            }
            if `"`r(datasignature)'"' != `"`xb_datasignature'"' {
                display as error "wfe_predict: current treatment variable no longer matches the stored estimation sample"
                exit 498
            }
        }
    }

    tempname b
    tempvar xb_hat
    matrix `b' = e(b)
    capture quietly matrix score double `xb_hat' = `b' if `touse'
    if _rc {
        if _rc == 111 {
            display as error "wfe_predict: regressors from estimation are no longer available in current data"
            exit 498
        }
        exit _rc
    }

    if `want_residuals' {
        capture __wfe_predict_confirm_exact_var `depvar'
        if _rc {
            display as error "wfe_predict: stored dependent variable is no longer available in current data"
            exit 498
        }
        quietly generate `typlist' `varlist' = `depvar' - `xb_hat' if `touse'
        label variable `varlist' "Residuals"
    }
    else {
        quietly generate `typlist' `varlist' = `xb_hat' if `touse'
        label variable `varlist' "Linear prediction"
    }
end

program define __wfe_pred_valsig, rclass
    version 16

    syntax anything(name=varspec)

    local normalized_vars
    capture quietly fvrevar `varspec' if e(sample)
    if _rc {
        exit _rc
    }

    local source_vars `"`r(varlist)'"'
    foreach source_var of local source_vars {
        tempvar normalized_var
        capture confirm numeric variable `source_var'
        if _rc {
            exit 498
        }
        capture quietly generate double `normalized_var' = `source_var' if e(sample)
        if _rc {
            exit 498
        }
        local normalized_vars `normalized_vars' `normalized_var'
    }

    capture quietly _datasignature `normalized_vars', esample nodefault nonames
    if _rc {
        exit _rc
    }
    return local datasignature `"`r(datasignature)'"'
end

* Validate exact variable name match
program define __wfe_predict_confirm_exact_var
    version 16

    syntax name(name=target)

    capture unab __wfe_predict_exact : `target'
    if _rc {
        exit 111
    }

    if `"`__wfe_predict_exact'"' != "`target'" {
        exit 111
    }
end

