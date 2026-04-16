// wfe_weights_fd.mata — First-Difference (FD) Weight Matrix
//
// FD estimator eliminates unit fixed effects by differencing adjacent periods
// where treatment changes. All weights are non-negative.
//
// Mata implementation for Stata 16.0+

version 16.0

mata:
mata set matastrict on

/*
 * wfe_weights_fd() — Compute FD weight matrix
 *
 * FD weights for treatment change periods. All weights are 0/1 based.
 *
 * Parameters:
 *   unit_idx, time_idx, treat, C_it, N, T_val, len_data, qoi, exist
 *
 * Returns:
 *   Wfd[T, N] real matrix — time as rows, units as columns
 */
real matrix wfe_weights_fd(
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
    real matrix Wfd, obs_lookup, cell_count
    real colvector weight_cit
    real colvector unit_col, time_col, treat_col, cit_col
    real scalar i, t, t_it, t_it_prev, c_it, k_prev
    real scalar is_ate
    real scalar ui, ti, k

    _wfe_validate_len_data_count(len_data, "wfe_weights_fd")

    if (qoi != "ate" & qoi != "att") {
        errprintf("wfe_weights_fd: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    unit_col = _wfe_require_vector_input(unit_idx, "unit_idx", "wfe_weights_fd")
    time_col = _wfe_require_vector_input(time_idx, "time_idx", "wfe_weights_fd")
    treat_col = _wfe_require_vector_input(treat, "treat", "wfe_weights_fd")
    cit_col = _wfe_require_vector_input(C_it, "C_it", "wfe_weights_fd")
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data | rows(cit_col) != len_data) {
        errprintf("wfe_weights_fd: input vector length mismatch\n")
        _error(3200)
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_weights_fd: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_weights_fd: treat must contain only 0/1 values\n")
        _error(3498)
    }
    weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_fd")
    if (len_data == 0) {
        if (N != floor(N) | T_val != floor(T_val) | N < 0 | T_val < 0) {
            errprintf("wfe_weights_fd: N and T must be nonnegative integers when len_data == 0\n")
            _error(3200)
        }
        _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                     len_data, "wfe_weights_fd")
        return(J(T_val, N, 0))
    }
    if (N != floor(N) | T_val != floor(T_val) | N < 1 | T_val < 1) {
        errprintf("wfe_weights_fd: N and T must be positive integers\n")
        _error(3200)
    }
    _wfe_validate_panel_index_domain(unit_col, time_col, N, T_val, "wfe_weights_fd")
    // FD weight generation compares only adjacent declared time slots j and j-1.
    // A globally empty row therefore stays a zero row on the declared T x N
    // grid instead of bridging a non-adjacent switch into a ghost first
    // difference.
    Wfd = J(T_val, N, 0)
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                 len_data, "wfe_weights_fd")

    // Build observation lookup matrix
    obs_lookup = J(T_val, N, 0)
    cell_count = J(T_val, N, 0)

    is_ate = (qoi == "ate")

    for (k = 1; k <= len_data; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("wfe_weights_fd: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
                                 len_data, "wfe_weights_fd")

    for (t = 2; t <= T_val; t++) {
        for (i = 1; i <= N; i++) {
            if (exist[t, i] == 0 | exist[t - 1, i] == 0) continue

            k = obs_lookup[t, i]
            k_prev = obs_lookup[t - 1, i]
            t_it = treat_col[k]
            t_it_prev = treat_col[k_prev]

            if (t_it == t_it_prev) {
                continue
            }

            c_it = weight_cit[k]

            if (is_ate) {
                Wfd[t, i]     = Wfd[t, i]     + c_it
                Wfd[t - 1, i] = Wfd[t - 1, i] + c_it
            }
            else {
                Wfd[t, i]     = Wfd[t, i]     + c_it * t_it
                Wfd[t - 1, i] = Wfd[t - 1, i] + c_it * t_it
            }
        }
    }

    return(Wfd)
}

end
