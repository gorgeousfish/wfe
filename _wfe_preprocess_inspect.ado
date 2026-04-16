*! _wfe_preprocess_inspect.ado -- Internal preprocessing and panel validation

program define _wfe_preprocess_inspect, rclass
    version 16

    syntax varlist(min=2 numeric fv) [if] [in],   ///
        TReat(varname numeric)                     ///
        Unit(varname)                              ///
        [                                          ///
        Time(varname)                              ///
        Method(string)                             ///
        CIT(varname numeric)                       ///
        OUTUNIT(name)                              ///
        OUTTIME(name)                              ///
        OUTCIT(name)                               ///
        OUTROWORDER(name)                          ///
        noVerbose                                  ///
        ]

    if "`method'" == "" local method "unit"
    if "`outunit'" == "" local outunit "__wfe_unit_idx"
    if "`outtime'" == "" local outtime "__wfe_time_idx"
    if "`outcit'" == "" local outcit "__wfe_cit"
    if "`outroworder'" == "" local outroworder "__wfe_row_order"

    tempvar wfe_dup
    marksample touse_ifin, novarlist
    marksample touse_formula

    quietly count if `touse_ifin'
    local nt_ifin = r(N)

    quietly count if `touse_formula'
    local nt_formula = r(N)

    if "`verbose'" != "noverbose" & `nt_formula' < `nt_ifin' {
        display as text " " _newline "Missing values are removed"
    }

    if `nt_formula' == 0 {
        display as error "no observations"
        exit 2000
    }

    capture drop `outunit' `outtime' `outcit' `outroworder'

    quietly gen long `outroworder' = .
    quietly replace `outroworder' = _n if `touse_formula'

    /*
       First restrict to the formula sample and only then validate the
       treatment indicator. A row that drops because y/X are missing must
       not invalidate the estimation sample merely because its off-sample
       treat() value is malformed.
    */
    quietly count if `touse_formula' & missing(`treat')
    if r(N) > 0 {
        display as error "missing values in treat() are not allowed"
        exit 198
    }

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
           Duplicate unit-time cells are structural errors only on the
           formula-valid sample that survives model-frame NA removal.
        */
        quietly duplicates tag `unit' `time' if `touse_formula', gen(`wfe_dup')
        quietly count if `touse_formula' & `wfe_dup' > 0
        if r(N) > 0 {
            display as error "duplicate observations found for unit-time combination"
            exit 198
        }
        drop `wfe_dup'
    }

    if "`cit'" != "" {
        quietly count if `touse_formula' & missing(`cit')
        if r(N) > 0 {
            display as error "missing values in cit() are not allowed"
            exit 198
        }
        quietly count if `touse_formula' & `cit' < 0
        if r(N) > 0 {
            display as error "'C.it' must be a non-negative numeric vector"
            exit 198
        }
    }

    quietly gen double `outcit' = .
    if "`cit'" == "" {
        quietly replace `outcit' = 1 if `touse_formula'
    }
    else {
        quietly replace `outcit' = `cit' if `touse_formula'
    }

    quietly egen long `outunit' = group(`unit') if `touse_formula'

    if "`time'" != "" {
        quietly egen long `outtime' = group(`time') if `touse_formula'
    }
    else {
        quietly gen long `outtime' = .
        sort `outunit' `outroworder'
        quietly by `outunit': replace `outtime' = sum(`touse_formula') if `touse_formula'
    }

    quietly levelsof `treat' if `touse_formula', local(__wfe_treat_levels)
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

    quietly summarize `treat' if `touse_formula', meanonly
    if r(min) != 0 | r(max) != 1 {
        display as error "'treat' must be a either zero or one where one indicates treatment"
        exit 198
    }

    quietly levelsof `outunit' if `touse_formula', local(__wfe_unit_levels)
    local N : word count `__wfe_unit_levels'
    quietly levelsof `outtime' if `touse_formula', local(__wfe_time_levels)
    local T : word count `__wfe_time_levels'
    quietly count if `touse_formula'
    local NT = r(N)

    if `N' == 1 {
        display as error "panel has only one unit; estimation not possible"
        exit 198
    }

    /*
       One-way time FE only needs cross-unit variation within each observed
       period. A single observed period remains feasible as long as there are
       at least two units. The single-period guard therefore applies only to
       paths that rely on within-unit time comparisons.
    */
    if "`method'" != "time" & `T' == 1 {
        display as error "panel has only one time period; estimation not possible"
        exit 198
    }

    local panel_balanced = (`NT' == `N' * `T')
    if "`method'" == "time" {
        local unit_number = `T'
    }
    else {
        local unit_number = `N'
    }


    return scalar N = `N'
    return scalar T = `T'
    return scalar NT = `NT'
    return scalar unit_number = `unit_number'
    return scalar panel_balanced = `panel_balanced'
    return local outunit `outunit'
    return local outtime `outtime'
    return local outcit `outcit'
    return local outroworder `outroworder'
end
