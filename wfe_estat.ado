*! wfe_estat.ado - Postestimation program for wfe and pwfe commands

program define wfe_estat, rclass
    version 16

    if !inlist("`e(cmd)'", "wfe", "pwfe") {
        error 301
    }

    // Parse subcommand

    local cmdline `"`0'"'
    gettoken subcmd 0 : 0, parse(" ,")
    local subcmd = lower(`"`subcmd'"')
    if "`subcmd'" != "wfe_weights" {
        estat_default `cmdline'
        exit
    }
    local trailing = trim(`"`0'"')
    if !inlist(`"`trailing'"', "", ",") {
        display as error "estat wfe_weights does not accept arguments or options"
        exit 198
    }

    // Validate weight matrix availability and load Mata support
    capture confirm matrix e(W)
    if _rc {
        display as error "e(W) is not available"
        exit 498
    }

    __wfe_postestimation_ensure_mata

    tempname W
    matrix `W' = e(W)
    mata: st_numscalar("__wfe_estat_W_has_missing", any(st_matrix("`W'") :>= .))
    if scalar(__wfe_estat_W_has_missing) {
        capture scalar drop __wfe_estat_W_has_missing
        display as error "estat wfe_weights: e(W) must not contain missing values"
        exit 498
    }

    // Compute and display weight summary statistics
    mata: wfe_store_weight_summary(st_matrix("`W'"))
    if scalar(__wfe_sum_n_nonzero) > 0 & scalar(__wfe_sum_n_negative) == 0 {
        scalar __wfe_sum_neg_ratio = 0
    }

    // Return results in r()
    return clear
    return scalar T = scalar(__wfe_sum_T)
    return scalar N = scalar(__wfe_sum_N)
    return scalar total = scalar(__wfe_sum_total)
    return scalar n_nonzero = scalar(__wfe_sum_n_nonzero)
    return scalar n_positive = scalar(__wfe_sum_n_positive)
    return scalar n_negative = scalar(__wfe_sum_n_negative)
    return scalar w_min = scalar(__wfe_sum_w_min)
    return scalar w_max = scalar(__wfe_sum_w_max)
    return scalar w_mean = scalar(__wfe_sum_w_mean)
    return scalar w_sd = scalar(__wfe_sum_w_sd)
    return scalar nz_min = scalar(__wfe_sum_nz_min)
    return scalar nz_max = scalar(__wfe_sum_nz_max)
    return scalar nz_mean = scalar(__wfe_sum_nz_mean)
    return scalar nz_sd = scalar(__wfe_sum_nz_sd)
    return scalar pos_min = scalar(__wfe_sum_pos_min)
    return scalar pos_max = scalar(__wfe_sum_pos_max)
    return scalar pos_mean = scalar(__wfe_sum_pos_mean)
    return scalar pos_sd = scalar(__wfe_sum_pos_sd)
    return scalar neg_min = scalar(__wfe_sum_neg_min)
    return scalar neg_max = scalar(__wfe_sum_neg_max)
    return scalar neg_mean = scalar(__wfe_sum_neg_mean)
    return scalar neg_sd = scalar(__wfe_sum_neg_sd)
    return scalar neg_ratio = scalar(__wfe_sum_neg_ratio)

    display _newline as text "Weight summary from e(W)"
    display as text "{hline 60}"
    display as text "Rows (T)"            _col(42) as result %10.0g scalar(__wfe_sum_T)
    display as text "Columns (N)"         _col(42) as result %10.0g scalar(__wfe_sum_N)
    display as text "Total elements"      _col(42) as result %10.0g scalar(__wfe_sum_total)
    display as text "Nonzero weights"     _col(42) as result %10.0g scalar(__wfe_sum_n_nonzero)
    display as text "Positive weights"    _col(42) as result %10.0g scalar(__wfe_sum_n_positive)
    display as text "Negative weights"    _col(42) as result %10.0g scalar(__wfe_sum_n_negative)
    display as text "{hline 60}"
    display as text "All weights:   min"  _col(22) as result %12.6f scalar(__wfe_sum_w_min) ///
        _col(40) as text "max" _col(45) as result %12.6f scalar(__wfe_sum_w_max)
    display as text "All weights:  mean"  _col(22) as result %12.6f scalar(__wfe_sum_w_mean) ///
        _col(40) as text "sd"  _col(45) as result %12.6f scalar(__wfe_sum_w_sd)

    if scalar(__wfe_sum_n_nonzero) > 0 {
        display as text "{hline 60}"
        display as text "Nonzero:      min" _col(22) as result %12.6f scalar(__wfe_sum_nz_min) ///
            _col(40) as text "max" _col(45) as result %12.6f scalar(__wfe_sum_nz_max)
        display as text "Nonzero:     mean" _col(22) as result %12.6f scalar(__wfe_sum_nz_mean) ///
            _col(40) as text "sd"  _col(45) as result %12.6f scalar(__wfe_sum_nz_sd)
    }

    if scalar(__wfe_sum_n_positive) > 0 {
        display as text "{hline 60}"
        display as text "Positive:     min" _col(22) as result %12.6f scalar(__wfe_sum_pos_min) ///
            _col(40) as text "max" _col(45) as result %12.6f scalar(__wfe_sum_pos_max)
        display as text "Positive:    mean" _col(22) as result %12.6f scalar(__wfe_sum_pos_mean) ///
            _col(40) as text "sd"  _col(45) as result %12.6f scalar(__wfe_sum_pos_sd)
    }

    if scalar(__wfe_sum_n_nonzero) > 0 {
        display as text "{hline 60}"
        display as text "Negative weight ratio" _col(42) as result %10.6f scalar(__wfe_sum_neg_ratio)
    }

    if scalar(__wfe_sum_n_negative) > 0 {
        display as text "Negative:     min" _col(22) as result %12.6f scalar(__wfe_sum_neg_min) ///
            _col(40) as text "max" _col(45) as result %12.6f scalar(__wfe_sum_neg_max)
        display as text "Negative:    mean" _col(22) as result %12.6f scalar(__wfe_sum_neg_mean) ///
            _col(40) as text "sd"  _col(45) as result %12.6f scalar(__wfe_sum_neg_sd)
    }

    // Clean up temporary scalars
    capture scalar drop __wfe_estat_W_has_missing
    foreach __wfe_scalar in ///
        __wfe_sum_T __wfe_sum_N __wfe_sum_total __wfe_sum_n_nonzero ///
        __wfe_sum_n_positive __wfe_sum_n_negative __wfe_sum_w_min ///
        __wfe_sum_w_max __wfe_sum_w_mean __wfe_sum_w_sd ///
        __wfe_sum_nz_min __wfe_sum_nz_max __wfe_sum_nz_mean ///
        __wfe_sum_nz_sd __wfe_sum_pos_min __wfe_sum_pos_max ///
        __wfe_sum_pos_mean __wfe_sum_pos_sd __wfe_sum_neg_min ///
        __wfe_sum_neg_max __wfe_sum_neg_mean __wfe_sum_neg_sd ///
        __wfe_sum_neg_ratio {
        capture scalar drop `__wfe_scalar'
    }
end

program define __wfe_postestimation_ensure_mata
    version 16
    // Ensure Mata functions are loaded; compile from package directory if needed

    capture mata: mata describe wfe_summarize_weights()
    if !_rc {
        exit
    }

    // findfile works for both flat net install and development layouts
    capture quietly findfile wfe_postestimation.mata
    if !_rc {
        quietly do `"`r(fn)'"'
        exit
    }

    // Fallback: try mata/ subdirectory relative to wfe_estat.ado (dev mode)
    capture quietly findfile wfe_estat.ado
    if !_rc {
        local ado_path `"`r(fn)'"'
        local ado_dir = substr(`"`ado_path'"', 1, ///
            length(`"`ado_path'"') - length("wfe_estat.ado"))
        capture quietly do `"`ado_dir'mata/wfe_postestimation.mata"'
        if !_rc {
            exit
        }
    }

    display as error "wfe_estat: cannot load Mata postestimation functions"
    exit 601
end
