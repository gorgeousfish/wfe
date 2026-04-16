// ============================================================
// wfe_weights_unit.mata
//
// One-way unit fixed effects weight matrix construction.
// Returns T x N weight matrix for unit-based matching with opposite-treatment.
// ============================================================

version 16.0

mata:
mata set matastrict on

/*
 * Compute unit FE weight matrix W[T, N].
 *
 * For each (unit, time) with valid observation, find within-unit
 * observations with opposite treatment status. Set self-weight to 1,
 * matched observations get equal weight 1/count, accumulated into
 * the unit's column of W scaled by C_it and treatment status.
 */
real matrix wfe_weights_unit(
    real matrix unit_idx,
    real matrix time_idx,
    real matrix treat,
    real matrix C_it,
    real scalar N,
    real scalar T_val,
    real scalar len_data,
    string scalar qoi,
    real matrix exist
)
{
    real matrix W, obs_lookup, cell_count
    real colvector w_it
    real colvector weight_cit
    real colvector unit_col, time_col, treat_col, cit_col
    real scalar i, t, t_prime, t_it, c_it_val, count, v
    real scalar ui, ti, k
    real scalar is_ate

    _wfe_validate_len_data_count(len_data, "wfe_weights_unit")

    if (qoi != "ate" & qoi != "att") {
        errprintf("wfe_weights_unit: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_weights_unit")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_weights_unit")
    treat_col = _wfe_require_vector_input(treat, "treat", "wfe_weights_unit")
    cit_col = _wfe_require_vector_input(C_it, "C_it", "wfe_weights_unit")
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data | rows(cit_col) != len_data) {
        errprintf("wfe_weights_unit: input vector length mismatch\n")
        _error(3200)
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_weights_unit: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_weights_unit: treat must contain only 0/1 values\n")
        _error(3498)
    }
    weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_unit")
    if (len_data == 0) {
        if (N != floor(N) | T_val != floor(T_val) | N < 0 | T_val < 0) {
            errprintf("wfe_weights_unit: N and T must be nonnegative integers when len_data == 0\n")
            _error(3200)
        }
        _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                     len_data, "wfe_weights_unit")
        return(J(T_val, N, 0))
    }
    if (N != floor(N) | T_val != floor(T_val) | N < 1 | T_val < 1) {
        errprintf("wfe_weights_unit: N and T must be positive integers\n")
        _error(3200)
    }
    _wfe_validate_panel_index_domain(unit_col, time_col, N, T_val, "wfe_weights_unit")
    W = J(T_val, N, 0)

    obs_lookup = J(T_val, N, 0)
    cell_count = J(T_val, N, 0)

    is_ate = (qoi == "ate")

    for (k = 1; k <= len_data; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("wfe_weights_unit: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                 len_data, "wfe_weights_unit")

    for (i = 1; i <= N; i++) {
        for (t = 1; t <= T_val; t++) {
            if (exist[t, i] == 0) continue

            k = obs_lookup[t, i]
            t_it = treat_col[k]
            c_it_val = weight_cit[k]

            count = 0
            for (t_prime = 1; t_prime <= T_val; t_prime++) {
                if (t_prime == t) continue
                if (exist[t_prime, i] == 0) continue
                if (treat_col[obs_lookup[t_prime, i]] == t_it) continue
                count++
            }

            if (count == 0) continue

            v = 1 / count

            // Build weight vector: self = 1, opposite-treatment matches = 1/count
            w_it = J(T_val, 1, 0)
            w_it[t] = 1
            for (t_prime = 1; t_prime <= T_val; t_prime++) {
                if (t_prime == t) continue
                if (exist[t_prime, i] == 0) continue
                if (treat_col[obs_lookup[t_prime, i]] == t_it) continue
                w_it[t_prime] = v
            }

            if (is_ate) {
                W[., i] = W[., i] + c_it_val * w_it
            }
            else {
                W[., i] = W[., i] + (c_it_val * t_it) * w_it
            }
        }
    }

    return(W)
}

end
