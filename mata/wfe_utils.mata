// Utility routines for matrix precomputation.
// Provides helper functions for building existence matrices,
// treatment-stability matrices, and vectorization operations.

version 16.0
mata:
mata set matastrict on

void _wfe_validate_len_data_count(
    real scalar   len_data,
    string scalar caller
)
{
    if (len_data >= . | len_data < 0 | len_data != floor(len_data)) {
        errprintf("%s: len_data must be a nonnegative integer\n", caller)
        _error(3200)
    }
}

void _wfe_validate_empty_panel(
    real scalar   N,
    real scalar   T,
    string scalar caller
)
{
    if (N != 0 | T != 0) {
        errprintf("%s: N and T must both be zero when len_data == 0\n", caller)
        _error(3200)
    }
}

void _wfe_validate_empty_W(
    real matrix    W,
    string scalar  caller
)
{
    if (rows(W) != 0 | cols(W) != 0) {
        errprintf("%s: W must be 0x0 when len_data == 0\n", caller)
        _error(3200)
    }
}

real colvector _wfe_as_colvector(real vector x)
{
    return(colshape(x[., .], 1))
}

real colvector _wfe_require_vector_input(
    real matrix   x,
    string scalar argname,
    string scalar caller
)
{
    if (rows(x) != 1 & cols(x) != 1) {
        errprintf("%s: %s must be a vector\n", caller, argname)
        _error(3200)
    }

    return(colshape(x[., .], 1))
}

real colvector _wfe_trunc_validate_cit(
    real colvector C_it,
    string scalar  caller
)
{
    real colvector weight_cit
    real scalar int_max

    if (any(C_it :>= .) | any(C_it :< 0)) {
        errprintf("%s: C_it must be nonmissing and non-negative\n", caller)
        _error(3498)
    }

    weight_cit = trunc(C_it)
    int_max = 2147483647

    if (any(weight_cit :> int_max)) {
        errprintf("%s: C_it must lie within the signed 32-bit integer range after truncation\n",
                  caller)
        _error(3498)
    }

    return(weight_cit)
}

void _wfe_validate_panel_index_domain(
    real colvector unit_idx,
    real colvector time_idx,
    real scalar    N,
    real scalar    T,
    string scalar  caller
)
{
    if (N != floor(N) | T != floor(T) | N < 1 | T < 1) {
        errprintf("%s: N and T must be positive integers\n", caller)
        _error(3200)
    }

    if (rows(unit_idx) == 0) {
        return
    }

    if (any(unit_idx :>= .) | any(time_idx :>= .)) {
        errprintf("%s: unit_idx/time_idx must not contain missing values\n", caller)
        _error(3498)
    }

    if (min(unit_idx) < 1 | max(unit_idx) > N | ///
        min(time_idx) < 1 | max(time_idx) > T | ///
        any(unit_idx :!= floor(unit_idx)) | ///
        any(time_idx :!= floor(time_idx))) {
        errprintf("%s: unit_idx/time_idx must be integer indices within 1..N and 1..T\n", caller)
        _error(3200)
    }
}

void _wfe_validate_sorted_panel_order(
    real colvector unit_idx,
    real colvector time_idx,
    string scalar  caller
)
{
    real scalar k

    if (rows(unit_idx) <= 1) {
        return
    }

    for (k = 2; k <= rows(unit_idx); k++) {
        if (unit_idx[k] < unit_idx[k-1] | ///
            (unit_idx[k] == unit_idx[k-1] & time_idx[k] < time_idx[k-1])) {
            errprintf("%s: input must be sorted by (unit_idx, time_idx)\n", caller)
            _error(3200)
        }
    }
}

void _wfe_validate_index_gaps(
    real colvector unit_idx,
    real colvector time_idx,
    real scalar    N,
    real scalar    T,
    string scalar  caller
)
{
    real colvector seen_units, seen_times
    real scalar k

    if (rows(unit_idx) == 0) {
        return
    }

    seen_units = J(N, 1, 0)
    seen_times = J(T, 1, 0)

    for (k = 1; k <= rows(unit_idx); k++) {
        seen_units[unit_idx[k]] = 1
        seen_times[time_idx[k]] = 1
    }

    if (any(seen_units :== 0) | any(seen_times :== 0)) {
        errprintf("%s: unit_idx/time_idx must enumerate 1..N and 1..T without gaps\n", caller)
        _error(3200)
    }
}

void _wfe_validate_time_gaps(
    real colvector time_idx,
    real scalar    T,
    string scalar  caller
)
{
    real colvector seen_times
    real scalar k

    if (rows(time_idx) == 0) {
        return
    }

    seen_times = J(T, 1, 0)

    for (k = 1; k <= rows(time_idx); k++) {
        seen_times[time_idx[k]] = 1
    }

    if (any(seen_times :== 0)) {
        errprintf("%s: time_idx must enumerate 1..T without gaps\n", caller)
        _error(3200)
    }
}

void _wfe_validate_exist_contract(
    real matrix    exist,
    real colvector unit_idx,
    real colvector time_idx,
    real scalar    N,
    real scalar    T,
    real scalar    len_data,
    string scalar  caller
)
{
    real matrix expected
    real scalar k

    if (rows(exist) != T | cols(exist) != N) {
        errprintf("%s: exist must match observed unit-time support\n", caller)
        _error(498)
    }

    if (any(exist :>= .) | any(exist :< 0) | any(exist :> 1) | ///
        any(exist :!= floor(exist))) {
        errprintf("%s: exist must match observed unit-time support\n", caller)
        _error(498)
    }

    expected = J(T, N, 0)
    for (k = 1; k <= len_data; k++) {
        if (expected[time_idx[k], unit_idx[k]] != 0) {
            errprintf("%s: unit-time pair is not unique\n", caller)
            _error(498)
        }
        expected[time_idx[k], unit_idx[k]] = 1
    }

    if (any(exist :!= expected)) {
        errprintf("%s: exist must match observed unit-time support\n", caller)
        _error(498)
    }
}

void _wfe_validate_same_contract(
    real matrix    same,
    real colvector unit_idx,
    real colvector time_idx,
    real colvector treat,
    real scalar    N,
    real scalar    T,
    real scalar    len_data,
    string scalar  caller
)
{
    real matrix expected, cell_count, obs_lookup
    real scalar k, ui, ti

    if (rows(same) != T | cols(same) != N) {
        errprintf("%s: same must match adjacent-period treatment stability\n", caller)
        _error(498)
    }

    if (any(same :>= .) | any(same :< 0) | any(same :> 1) | ///
        any(same :!= floor(same))) {
        errprintf("%s: same must match adjacent-period treatment stability\n", caller)
        _error(498)
    }

    if (any(treat :>= .)) {
        errprintf("%s: treat must not contain missing values\n", caller)
        _error(3498)
    }

    if (any((treat :!= 0) :& (treat :!= 1))) {
        errprintf("%s: treat must contain only 0/1 values\n", caller)
        _error(3498)
    }

    expected = J(T, N, 0)
    cell_count = J(T, N, 0)
    obs_lookup = J(T, N, 0)

    for (k = 1; k <= len_data; k++) {
        ui = unit_idx[k]
        ti = time_idx[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("%s: unit-time pair is not unique\n", caller)
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }

    for (ui = 1; ui <= N; ui++) {
        for (ti = 2; ti <= T; ti++) {
            if (obs_lookup[ti, ui] > 0 & obs_lookup[ti - 1, ui] > 0) {
                if (treat[obs_lookup[ti, ui]] == treat[obs_lookup[ti - 1, ui]]) {
                    expected[ti, ui] = 1
                }
            }
        }
    }

    if (any(same :!= expected)) {
        errprintf("%s: same must match adjacent-period treatment stability\n", caller)
        _error(498)
    }
}

// Construct the T x N observation-existence matrix.
real matrix wfe_build_exist(real matrix unit_idx, real matrix time_idx,
                            real scalar N, real scalar T, real scalar len_data)
{
    real matrix exist
    real colvector unit_col, time_col
    real scalar k

    _wfe_validate_len_data_count(len_data, "wfe_build_exist")
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_build_exist")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_build_exist")
    if (rows(unit_col) != len_data | rows(time_col) != len_data) {
        errprintf("wfe_build_exist: input vector length mismatch\n")
        _error(3200)
    }
    if (len_data == 0) {
        if (N != floor(N) | T != floor(T) | N < 0 | T < 0) {
            errprintf("wfe_build_exist: N and T must be nonnegative integers when len_data == 0\n")
            _error(3200)
        }
        return(J(T, N, 0))
    }

    // Support builders operate on the declared grid and may return
    // zero rows/columns for unobserved unit/time levels.
    _wfe_validate_panel_index_domain(unit_col, time_col, N, T, "wfe_build_exist")

    exist = J(T, N, 0)

    for (k = 1; k <= len_data; k++) {
        // A declared panel cell (i,t) is a unique observation slot per
        // Imai-Kim (2021). Repeated rows duplicate the same dyad rather than
        // add new support and must fail-fast at the helper boundary.
        if (exist[time_col[k], unit_col[k]] != 0) {
            errprintf("wfe_build_exist: unit-time pair is not unique\n")
            _error(498)
        }
        exist[time_col[k], unit_col[k]] = 1
    }

    return(exist)
}


// Construct the T x N treatment-stability matrix.
real matrix wfe_build_same(real matrix unit_idx, real matrix time_idx,
                           real matrix treat, real scalar N, real scalar T,
                           real scalar len_data)
{
    real matrix same, obs_lookup, treat_lookup
    real colvector unit_col, time_col, treat_col
    real scalar ui, ti, k

    _wfe_validate_len_data_count(len_data, "wfe_build_same")
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_build_same")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_build_same")
    treat_col = _wfe_require_vector_input(treat, "treat", "wfe_build_same")
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data) {
        errprintf("wfe_build_same: input vector length mismatch\n")
        _error(3200)
    }
    if (len_data == 0) {
        if (N != floor(N) | T != floor(T) | N < 0 | T < 0) {
            errprintf("wfe_build_same: N and T must be nonnegative integers when len_data == 0\n")
            _error(3200)
        }
        return(J(T, N, 0))
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_build_same: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_build_same: treat must contain only 0/1 values\n")
        _error(3498)
    }

    _wfe_validate_panel_index_domain(unit_col, time_col, N, T, "wfe_build_same")

    same = J(T, N, 0)
    obs_lookup = J(T, N, 0)
    treat_lookup = J(T, N, .)

    for (k = 1; k <= len_data; k++) {
        ti = time_col[k]
        ui = unit_col[k]
        if (obs_lookup[ti, ui] != 0) {
            errprintf("wfe_build_same: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
        treat_lookup[ti, ui] = treat_col[k]
    }

    for (ui = 1; ui <= N; ui++) {
        for (ti = 2; ti <= T; ti++) {
            if (obs_lookup[ti, ui] > 0 & obs_lookup[ti - 1, ui] > 0) {
                if (treat_lookup[ti, ui] == treat_lookup[ti - 1, ui]) {
                    same[ti, ui] = 1
                }
            }
        }
    }

    return(same)
}


// Map a T x N weight matrix back to the sorted observation order.
real colvector wfe_vectorize(real matrix W, real matrix time_idx,
                             real matrix unit_idx, real scalar len_data)
{
    real colvector W_vec
    real colvector time_col, unit_col
    real scalar k

    _wfe_validate_len_data_count(len_data, "wfe_vectorize")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_vectorize")
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_vectorize")
    if (rows(time_col) != len_data | rows(unit_col) != len_data) {
        errprintf("wfe_vectorize: input vector length mismatch\n")
        _error(3200)
    }
    if (len_data == 0) {
        if (rows(W) < 0 | cols(W) < 0) {
            errprintf("wfe_vectorize: W dimensions must be nonnegative when len_data == 0\n")
            _error(3200)
        }
        return(J(0, 1, 0))
    }
    if (any(W :>= .)) {
        errprintf("wfe_vectorize: W must not contain missing values\n")
        _error(3498)
    }

    // This is a pure lookup over observed dyads. As long as the
    // observed indices are in-range, declared-grid gaps are legitimate and
    // should map to the corresponding W[t, i] cells without fail-fast.
    _wfe_validate_panel_index_domain(unit_col, time_col, cols(W), rows(W), "wfe_vectorize")

    W_vec = J(len_data, 1, 0)

    for (k = 1; k <= len_data; k++) {
        // Repeated dyads get the same grid lookup instead of being treated
        // as an error. This helper is a pure W[t, i] -> row remap;
        // panel-uniqueness checks belong to higher-level estimators.
        W_vec[k] = W[time_col[k], unit_col[k]]
    }

    return(W_vec)
}

end
