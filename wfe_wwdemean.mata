// wfe_wwdemean.mata — Weighted Within Demeaning and Standard FE Demeaning
//
// Functions (8):
//   wfe_wwdemean()                  — Sorted + single column
//   wfe_demean()                    — Sorted + single column (unweighted)
//   wfe_wwdemean_unsorted()         — Unsorted + single column
//   wfe_demean_unsorted()           — Unsorted + single column (unweighted)
//   wfe_wwdemean_matrix()           — Sorted + multiple columns
//   wfe_demean_matrix()             — Sorted + multiple columns (unweighted)
//   wfe_wwdemean_matrix_unsorted()  — Unsorted + multiple columns
//   wfe_demean_matrix_unsorted()    — Unsorted + multiple columns (unweighted)

version 16.0
mata:
mata set matastrict on

void _wfe_wwd_validate_len_data_count(
    real scalar   len_data,
    string scalar caller
)
{
    if (len_data >= . | len_data < 0 | len_data != floor(len_data)) {
        errprintf("%s: len_data must be a nonnegative integer\n", caller)
        _error(3200)
    }
}

real colvector _wfe_wwd_as_colvector(real vector x)
{
    return(colshape(x[., .], 1))
}

real matrix _wfe_wwd_as_design(
    real matrix  data,
    real scalar  len_data
)
{
    real matrix data_use

    data_use = data[., .]
    if (rows(data_use) == 1 & cols(data_use) == len_data & len_data > 1) {
        data_use = data_use'
    }

    return(data_use)
}

void _wfe_validate_empty_grp0(
    real scalar   n_grp,
    real scalar   len_data,
    string scalar caller
)
{
    if (len_data == 0) {
        if (n_grp != 0) {
            errprintf("%s: n_grp must be zero when len_data == 0\n", caller)
            _error(3200)
        }
    }
}

void _wfe_validate_grp_contract(
    real colvector grp_idx,
    real scalar    n_grp,
    string scalar  caller
)
{
    real colvector grp_sorted
    real scalar i, n_unique

    if (n_grp != floor(n_grp) | n_grp < 1) {
        errprintf("%s: n_grp must be a positive integer\n", caller)
        _error(3200)
    }

    if (any(grp_idx :>= .)) {
        errprintf("%s: grp_idx must not contain missing values\n", caller)
        _error(3498)
    }

    if (min(grp_idx) < 1 | any(grp_idx :!= floor(grp_idx))) {
        errprintf("%s: grp_idx must contain positive integer labels\n", caller)
        _error(3200)
    }

    grp_sorted = sort(grp_idx, 1)
    n_unique = 0
    for (i = 1; i <= rows(grp_sorted); i++) {
        if (i == 1) {
            n_unique++
        }
        else if (grp_sorted[i] != grp_sorted[i - 1]) {
            n_unique++
        }
    }

    // Within-demeaning depends only on the group partition, not on whether
    // caller-side labels are literally 1..n_grp. Gapful positive relabeling
    // should therefore be accepted as long as it still identifies n_grp groups.
    if (n_unique != n_grp) {
        errprintf("%s: grp_idx must identify exactly n_grp unique groups\n", caller)
        _error(3200)
    }
}

void _wfe_validate_no_missing_vector(
    real colvector x,
    string scalar  name,
    string scalar  caller
)
{
    if (any(x :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, name)
        _error(3498)
    }
}

void _wfe_validate_no_missing_matrix(
    real matrix    x,
    string scalar  name,
    string scalar  caller
)
{
    if (any(x :>= .)) {
        errprintf("%s: %s must not contain missing values\n", caller, name)
        _error(3498)
    }
}

// wfe_wwdemean() — Weighted within demeaning (sorted + single column)
//
// result[k] = sqrt(weight[k]) * (var[k] - wmean[grp(k)])
// wmean[i]  = sum(weight * var) / sum(weight) within group i
//
// @param var       real colvector [len_data x 1]  Data to demean
// @param weight    real colvector [len_data x 1]  Weights (non-negative)
// @param grp_idx   real colvector [len_data x 1]  Group index (MUST be sorted)
// @param n_grp     real scalar                    Number of unique groups
// @param len_data  real scalar                    Total observations
// @return          real colvector [len_data x 1]  WWDemeaned result
//
// Zero-weight handling: If sum(weight)==0 for a group, wmean=0 and result=0.
// Pre-condition: grp_idx must be sorted ascending for panelsetup().
real colvector wfe_wwdemean(real vector    var,
                            real vector    weight,
                            real vector    grp_idx,
                            real scalar    n_grp,
                            real scalar    len_data)
{
    real colvector result, var_i, w_i, var_col, weight_col, grp_col
    real matrix    info
    real scalar    i, a, b, sum_w, wmean, k

    // Input validation
    _wfe_wwd_validate_len_data_count(len_data, "wfe_wwdemean")
    var_col = _wfe_wwd_as_colvector(var)
    weight_col = _wfe_wwd_as_colvector(weight)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(var_col) != len_data | rows(weight_col) != len_data |
        rows(grp_col) != len_data) {
        errprintf("wfe_wwdemean: var/weight/grp_idx length mismatch\n")
        _error(3200)
    }

    // Empty input handling
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_wwdemean")
    if (len_data == 0) return(J(0, 1, .))

    _wfe_validate_no_missing_vector(var_col, "var", "wfe_wwdemean")
    _wfe_validate_no_missing_vector(weight_col, "weight", "wfe_wwdemean")

    if (any(weight_col :< 0)) {
        errprintf("wfe_wwdemean: weight must be non-negative\n")
        _error(3498)
    }
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_wwdemean")

    // Initialize result
    result = J(len_data, 1, 0)

    for (k = 2; k <= len_data; k++) {
        if (grp_col[k] < grp_col[k-1]) {
            errprintf("wfe_wwdemean: grp_idx must be sorted in contiguous blocks\n")
            _error(3200)
        }
    }

    // Check sort order
    info = panelsetup(grp_col, 1)
    if (rows(info) != n_grp) {
        errprintf("wfe_wwdemean: grp_idx must identify exactly n_grp unique groups\n")
        _error(3200)
    }

    // Panel setup
    for (i = 1; i <= rows(info); i++) {
        a = info[i, 1]
        b = info[i, 2]

        var_i = var_col[a::b]
        w_i   = weight_col[a::b]

        // Pass 1: weighted mean
        sum_w = quadcolsum(w_i)
        if (sum_w != 0) {
            wmean = quadcolsum(w_i :* var_i) / sum_w
        }
        else {
            wmean = 0
        }

        // Pass 2: sqrt(W) * (var - wmean)
        result[a::b] = sqrt(w_i) :* (var_i :- wmean)
    }

    return(result)
}


// wfe_demean() — Standard FE demeaning (sorted + single column)
//
// result[k] = var[k] - mean[grp(k)]
// mean[i]   = sum(var) / count within group i
//
// @param var       real colvector [len_data x 1]  Data to demean
// @param grp_idx   real colvector [len_data x 1]  Group index (MUST be sorted)
// @param n_grp     real scalar                    Number of unique groups
// @param len_data  real scalar                    Total observations
// @return          real colvector [len_data x 1]  Demeaned result
//
// No weight parameter. No zero-division protection needed
// (panelsetup guarantees T_i >= 1 for each panel).
real colvector wfe_demean(real vector    var,
                          real vector    grp_idx,
                          real scalar    n_grp,
                          real scalar    len_data)
{
    real colvector result, var_i, var_col, grp_col
    real matrix    info
    real scalar    i, a, b, T_i, mean_i, k

    // Input validation
    _wfe_wwd_validate_len_data_count(len_data, "wfe_demean")
    var_col = _wfe_wwd_as_colvector(var)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(var_col) != len_data | rows(grp_col) != len_data) {
        errprintf("wfe_demean: var/grp_idx length mismatch\n")
        _error(3200)
    }

    // Empty input handling
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_demean")
    if (len_data == 0) return(J(0, 1, .))

    _wfe_validate_no_missing_vector(var_col, "var", "wfe_demean")
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_demean")

    // Initialize
    result = J(len_data, 1, 0)
    for (k = 2; k <= len_data; k++) {
        if (grp_col[k] < grp_col[k-1]) {
            errprintf("wfe_demean: grp_idx must be sorted in contiguous blocks\n")
            _error(3200)
        }
    }
    info   = panelsetup(grp_col, 1)
    if (rows(info) != n_grp) {
        errprintf("wfe_demean: grp_idx must identify exactly n_grp unique groups\n")
        _error(3200)
    }

    // Per-panel loop
    for (i = 1; i <= rows(info); i++) {
        a = info[i, 1]
        b = info[i, 2]

        var_i = var_col[a::b]
        T_i   = b - a + 1

        // Simple mean
        mean_i = quadcolsum(var_i) / T_i

        // Demean
        result[a::b] = var_i :- mean_i
    }

    return(result)
}


// wfe_wwdemean_unsorted() — Weighted within demeaning for unsorted data (single column)
//
// Uses permutation-unpermutation for cases where data is sorted by unit_idx
// but needs to be grouped by time_idx.
real colvector wfe_wwdemean_unsorted(real vector    var,
                                     real vector    weight,
                                     real vector    grp_idx,
                                     real scalar    n_grp,
                                     real scalar    len_data)
{
    real colvector perm, sorted_result, result
    real colvector var_col, weight_col, grp_col

    _wfe_wwd_validate_len_data_count(len_data, "wfe_wwdemean_unsorted")
    var_col = _wfe_wwd_as_colvector(var)
    weight_col = _wfe_wwd_as_colvector(weight)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(var_col) != len_data | rows(weight_col) != len_data | ///
        rows(grp_col) != len_data) {
        errprintf("wfe_wwdemean_unsorted: var/weight/grp_idx length mismatch\n")
        _error(3200)
    }
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_wwdemean_unsorted")
    if (len_data == 0) return(J(0, 1, .))
    _wfe_validate_no_missing_vector(var_col, "var", "wfe_wwdemean_unsorted")
    _wfe_validate_no_missing_vector(weight_col, "weight", "wfe_wwdemean_unsorted")
    if (any(weight_col :< 0)) {
        errprintf("wfe_wwdemean_unsorted: weight must be non-negative\n")
        _error(3498)
    }
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_wwdemean_unsorted")

    // Step 1: sort by group index
    perm = order(grp_col, 1)

    // Step 2: call sorted version on permuted data
    sorted_result = wfe_wwdemean(var_col[perm], weight_col[perm],
                                 grp_col[perm], n_grp, len_data)

    // Step 3: unpermute back to original row order
    result = J(len_data, 1, .)
    result[perm] = sorted_result

    return(result)
}


// wfe_demean_unsorted() — Standard FE demeaning for unsorted data (single column)
real colvector wfe_demean_unsorted(real vector    var,
                                   real vector    grp_idx,
                                   real scalar    n_grp,
                                   real scalar    len_data)
{
    real colvector perm, sorted_result, result
    real colvector var_col, grp_col

    _wfe_wwd_validate_len_data_count(len_data, "wfe_demean_unsorted")
    var_col = _wfe_wwd_as_colvector(var)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(var_col) != len_data | rows(grp_col) != len_data) {
        errprintf("wfe_demean_unsorted: var/grp_idx length mismatch\n")
        _error(3200)
    }
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_demean_unsorted")
    if (len_data == 0) return(J(0, 1, .))
    _wfe_validate_no_missing_vector(var_col, "var", "wfe_demean_unsorted")
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_demean_unsorted")

    perm = order(grp_col, 1)
    sorted_result = wfe_demean(var_col[perm], grp_col[perm], n_grp, len_data)

    result = J(len_data, 1, .)
    result[perm] = sorted_result

    return(result)
}


// wfe_wwdemean_matrix() — Weighted within demeaning for sorted multi-column data
//
// Optimized: single panelsetup() call, share sum_w and sqrt_w_i
// across all columns within each panel.
real matrix wfe_wwdemean_matrix(real matrix    data,
                                real vector    weight,
                                real vector    grp_idx,
                                real scalar    n_grp,
                                real scalar    len_data)
{
    real scalar    p, i, a, b, sum_w, k
    real matrix    info, result, data_i, data_use
    real colvector w_i, sqrt_w_i, weight_col, grp_col
    real rowvector wmean_row

    _wfe_wwd_validate_len_data_count(len_data, "wfe_wwdemean_matrix")
    data_use = _wfe_wwd_as_design(data, len_data)
    weight_col = _wfe_wwd_as_colvector(weight)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(data_use) != len_data | rows(weight_col) != len_data | ///
        rows(grp_col) != len_data) {
        errprintf("wfe_wwdemean_matrix: data/weight/grp_idx length mismatch\n")
        _error(3200)
    }

    p = cols(data_use)
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_wwdemean_matrix")
    if (len_data == 0) return(J(0, p, .))
    if (p == 0) {
        errprintf("wfe_wwdemean_matrix: data must contain at least one column\n")
        _error(3200)
    }

    _wfe_validate_no_missing_matrix(data_use, "data", "wfe_wwdemean_matrix")
    _wfe_validate_no_missing_vector(weight_col, "weight", "wfe_wwdemean_matrix")

    if (any(weight_col :< 0)) {
        errprintf("wfe_wwdemean_matrix: weight must be non-negative\n")
        _error(3498)
    }
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_wwdemean_matrix")

    result = J(len_data, p, 0)
    for (k = 2; k <= len_data; k++) {
        if (grp_col[k] < grp_col[k-1]) {
            errprintf("wfe_wwdemean_matrix: grp_idx must be sorted in contiguous blocks\n")
            _error(3200)
        }
    }
    info   = panelsetup(grp_col, 1)
    if (rows(info) != n_grp) {
        errprintf("wfe_wwdemean_matrix: grp_idx must identify exactly n_grp unique groups\n")
        _error(3200)
    }

    for (i = 1; i <= rows(info); i++) {
        a = info[i, 1]
        b = info[i, 2]

        w_i      = weight_col[a::b]
        sqrt_w_i = sqrt(w_i)
        sum_w    = quadcolsum(w_i)

        data_i = data_use[a::b, .]

        if (sum_w != 0) {
            // Vectorized weighted mean: [1 x p]
            wmean_row = quadcolsum(w_i :* data_i) / sum_w
        }
        else {
            wmean_row = J(1, p, 0)
        }

        // sqrt(W) * (data - wmean) with broadcast
        result[a::b, .] = sqrt_w_i :* (data_i :- wmean_row)
    }

    return(result)
}


// wfe_demean_matrix() — Standard FE demeaning for sorted multi-column data
real matrix wfe_demean_matrix(real matrix    data,
                              real vector    grp_idx,
                              real scalar    n_grp,
                              real scalar    len_data)
{
    real scalar    p, i, a, b, T_i, k
    real matrix    info, result, data_i, data_use
    real colvector grp_col
    real rowvector mean_row

    _wfe_wwd_validate_len_data_count(len_data, "wfe_demean_matrix")
    data_use = _wfe_wwd_as_design(data, len_data)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(data_use) != len_data | rows(grp_col) != len_data) {
        errprintf("wfe_demean_matrix: data/grp_idx length mismatch\n")
        _error(3200)
    }

    p = cols(data_use)
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_demean_matrix")
    if (len_data == 0) return(J(0, p, .))
    if (p == 0) {
        errprintf("wfe_demean_matrix: data must contain at least one column\n")
        _error(3200)
    }

    _wfe_validate_no_missing_matrix(data_use, "data", "wfe_demean_matrix")
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_demean_matrix")

    result = J(len_data, p, 0)
    for (k = 2; k <= len_data; k++) {
        if (grp_col[k] < grp_col[k-1]) {
            errprintf("wfe_demean_matrix: grp_idx must be sorted in contiguous blocks\n")
            _error(3200)
        }
    }
    info   = panelsetup(grp_col, 1)
    if (rows(info) != n_grp) {
        errprintf("wfe_demean_matrix: grp_idx must identify exactly n_grp unique groups\n")
        _error(3200)
    }

    for (i = 1; i <= rows(info); i++) {
        a = info[i, 1]
        b = info[i, 2]
        T_i = b - a + 1

        data_i   = data_use[a::b, .]
        mean_row = quadcolsum(data_i) / T_i

        result[a::b, .] = data_i :- mean_row
    }

    return(result)
}


// wfe_wwdemean_matrix_unsorted() — Weighted within demeaning for unsorted multi-column data
//
// Single permute → call matrix sorted → single unpermute.
// More efficient than per-column unsorted calls (one sort, not p).
real matrix wfe_wwdemean_matrix_unsorted(real matrix    data,
                                         real vector    weight,
                                         real vector    grp_idx,
                                         real scalar    n_grp,
                                         real scalar    len_data)
{
    real colvector perm
    real matrix    sorted_result, result, data_use
    real colvector weight_col, grp_col

    _wfe_wwd_validate_len_data_count(len_data, "wfe_wwdemean_matrix_unsorted")
    data_use = _wfe_wwd_as_design(data, len_data)
    weight_col = _wfe_wwd_as_colvector(weight)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(data_use) != len_data | rows(weight_col) != len_data | ///
        rows(grp_col) != len_data) {
        errprintf("wfe_wwdemean_matrix_unsorted: data/weight/grp_idx length mismatch\n")
        _error(3200)
    }
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_wwdemean_matrix_unsorted")
    if (len_data == 0) return(J(0, cols(data_use), .))
    if (cols(data_use) == 0) {
        errprintf("wfe_wwdemean_matrix_unsorted: data must contain at least one column\n")
        _error(3200)
    }
    _wfe_validate_no_missing_matrix(data_use, "data", "wfe_wwdemean_matrix_unsorted")
    _wfe_validate_no_missing_vector(weight_col, "weight", "wfe_wwdemean_matrix_unsorted")
    if (any(weight_col :< 0)) {
        errprintf("wfe_wwdemean_matrix_unsorted: weight must be non-negative\n")
        _error(3498)
    }
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_wwdemean_matrix_unsorted")

    perm = order(grp_col, 1)
    sorted_result = wfe_wwdemean_matrix(data_use[perm, .], weight_col[perm],
                                        grp_col[perm], n_grp, len_data)

    result = J(len_data, cols(data_use), .)
    result[perm, .] = sorted_result

    return(result)
}


// wfe_demean_matrix_unsorted() — Standard FE demeaning for unsorted multi-column data
real matrix wfe_demean_matrix_unsorted(real matrix    data,
                                       real vector    grp_idx,
                                       real scalar    n_grp,
                                       real scalar    len_data)
{
    real colvector perm
    real matrix    sorted_result, result, data_use
    real colvector grp_col

    _wfe_wwd_validate_len_data_count(len_data, "wfe_demean_matrix_unsorted")
    data_use = _wfe_wwd_as_design(data, len_data)
    grp_col = _wfe_wwd_as_colvector(grp_idx)
    if (rows(data_use) != len_data | rows(grp_col) != len_data) {
        errprintf("wfe_demean_matrix_unsorted: data/grp_idx length mismatch\n")
        _error(3200)
    }
    _wfe_validate_empty_grp0(n_grp, len_data, "wfe_demean_matrix_unsorted")
    if (len_data == 0) return(J(0, cols(data_use), .))
    if (cols(data_use) == 0) {
        errprintf("wfe_demean_matrix_unsorted: data must contain at least one column\n")
        _error(3200)
    }
    _wfe_validate_no_missing_matrix(data_use, "data", "wfe_demean_matrix_unsorted")
    _wfe_validate_grp_contract(grp_col, n_grp, "wfe_demean_matrix_unsorted")

    perm = order(grp_col, 1)
    sorted_result = wfe_demean_matrix(data_use[perm, .], grp_col[perm],
                                       n_grp, len_data)

    result = J(len_data, cols(data_use), .)
    result[perm, .] = sorted_result

    return(result)
}

end
