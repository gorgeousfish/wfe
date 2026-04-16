// ============================================================
// wfe_weights_time.mata — Time fixed effects weight computation
//
// Compute weight matrix W[T, N] for time fixed effects estimator
// W[t, i] = cumulative weight for unit i at time t
//
// Matching rule: For observation (i,t), match with other units at same
// time t that have opposite treatment status (X=1-X)
//
// Dependencies: wfe_utils (wfe_build_exist)
// ============================================================

version 16.0

mata:
mata set matastrict on

/*
 * wfe_weights_time() - Compute time FE weights
 *
 * Parameters:
 *   unit_idx  - unit index vector (1..N)
 *   time_idx  - time index vector (1..T)
 *   treat     - treatment indicator (0/1)
 *   C_it      - weight adjustment factor (non-negative, default 1)
 *   N         - number of unique units
 *   T_val     - number of unique time periods (T is Mata keyword)
 *   len_data  - total observations
 *   qoi       - "ate" or "att"
 *   exist     - existence matrix exist[T, N]
 *
 * Returns:
 *   W[T_val, N] weight matrix
 *
 * Algorithm:
 *   For each time t and each existing observation (i,t):
 *   1. Find other units at same time with opposite treatment
 *   2. With n matches: self gets weight 1, each match gets 1/n
 *   3. ATE: W[t,.] += c_it * w_it
 *      ATT: W[t,.] += c_it * t_it * w_it (zeros out control group)
 */
real matrix wfe_weights_time(
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
    real colvector weight_cit
    real colvector unit_col, time_col, treat_col, cit_col
    real rowvector w_it
    real scalar t, i, i_prime, k, k_prime, t_it, c_it_val
    real scalar ui, ti
    real scalar count, v, m, is_ate
    // match_indices: temporary storage for match positions
    real rowvector match_indices

    _wfe_validate_len_data_count(len_data, "wfe_weights_time")

    if (qoi != "ate" & qoi != "att") {
        errprintf("wfe_weights_time: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_weights_time")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_weights_time")
    treat_col = _wfe_require_vector_input(treat, "treat", "wfe_weights_time")
    cit_col = _wfe_require_vector_input(C_it, "C_it", "wfe_weights_time")
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data | rows(cit_col) != len_data) {
        errprintf("wfe_weights_time: input vector length mismatch\n")
        _error(3200)
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_weights_time: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_weights_time: treat must contain only 0/1 values\n")
        _error(3498)
    }
    weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_time")
    if (len_data == 0) {
        if (N != floor(N) | T_val != floor(T_val) | N < 0 | T_val < 0) {
            errprintf("wfe_weights_time: N and T must be nonnegative integers when len_data == 0\n")
            _error(3200)
        }
        _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                     len_data, "wfe_weights_time")
        return(J(T_val, N, 0))
    }
    if (N != floor(N) | T_val != floor(T_val) | N < 1 | T_val < 1) {
        errprintf("wfe_weights_time: N and T must be positive integers\n")
        _error(3200)
    }
    _wfe_validate_panel_index_domain(unit_col, time_col, N, T_val, "wfe_weights_time")
    W = J(T_val, N, 0)
    // Initialize
    cell_count = J(T_val, N, 0)
    is_ate = (qoi == "ate")

    // Build obs_lookup table
    // obs_lookup[t, i] = row index of observation (i,t) in data, 0 if absent
    obs_lookup = J(T_val, N, 0)
    for (k = 1; k <= len_data; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("wfe_weights_time: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                 len_data, "wfe_weights_time")

    // Main loop: iterate over time periods
    for (t = 1; t <= T_val; t++) {

        // Inner loop: iterate over units at current time
        for (i = 1; i <= N; i++) {

            // Check existence
            if (exist[t, i] == 0) continue

            // Get treatment status and C_it
            k = obs_lookup[t, i]
            t_it = treat_col[k]
            c_it_val = weight_cit[k]

            // Find matching units at same time t with opposite treatment
            count = 0
            w_it = J(1, N, 0)
            w_it[1, i] = 1
            match_indices = J(1, N, 0)

            for (i_prime = 1; i_prime <= N; i_prime++) {
                if (i_prime == i) continue
                if (exist[t, i_prime] == 0) continue
                k_prime = obs_lookup[t, i_prime]
                if (treat_col[k_prime] == t_it) continue
                count++
                match_indices[1, count] = i_prime
            }

            // No matches found
            if (count == 0) continue

            // Assign match weights
            v = 1 / count
            for (m = 1; m <= count; m++) {
                w_it[1, match_indices[1, m]] = v
            }

            // Accumulate to W matrix
            if (is_ate) {
                W[t, .] = W[t, .] + c_it_val * w_it
            }
            else {
                // ATT: control group (t_it=0) contributes zero
                W[t, .] = W[t, .] + (c_it_val * t_it) * w_it
            }
        }
    }

    return(W)
}

end
