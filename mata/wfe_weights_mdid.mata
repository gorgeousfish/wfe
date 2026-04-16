// wfe_weights_mdid.mata — MDiD weight matrix computation (nearest-neighbor matching)

version 16.0

mata:
mata set matastrict on

real matrix wfe_weights_mdid(
    real vector unit_idx,
    real vector time_idx,
    real vector treat,
    real vector C_it,
    real vector y,
    real scalar maxdev,
    real scalar N,
    real scalar T_val,
    real scalar len_data,
    string scalar qoi,
    real matrix exist,
    real matrix same
)
{
    real matrix Wdid, W, obs_lookup, cell_count
    real colvector weight_cit
    real colvector unit_col, time_col, treat_col, cit_col, y_col
    real scalar t, i, k, i_prime, k_it, k_it1, ui, ti
    real scalar t_it, c_it, y_it1, diff, min_diff, scale
    real scalar is_ate, has_match, match_unit, target_treat
    real scalar match_count, v_it

    if (qoi != "ate" & qoi != "att") {
        errprintf("wfe_weights_mdid: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    if (missing(len_data) | len_data != floor(len_data) | len_data < 0) {
        errprintf("wfe_weights_mdid: len_data must be a nonnegative integer\n")
        _error(3200)
    }
    unit_col = _wfe_as_colvector(unit_idx)
    time_col = _wfe_as_colvector(time_idx)
    treat_col = _wfe_as_colvector(treat)
    cit_col = _wfe_as_colvector(C_it)
    y_col = _wfe_as_colvector(y)
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data | rows(cit_col) != len_data | ///
        rows(y_col) != len_data) {
        errprintf("wfe_weights_mdid: input vector length mismatch\n")
        _error(3200)
    }
    if (missing(maxdev)) {
        errprintf("wfe_weights_mdid: maxdev must be nonmissing; use a negative value for nearest-neighbor matching\n")
        _error(3498)
    }
    // Negative maxdev enables nearest-neighbor matching
    if (any(y_col :>= .)) {
        errprintf("wfe_weights_mdid: y must not contain missing values\n")
        _error(3498)
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_weights_mdid: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_weights_mdid: treat must contain only 0/1 values\n")
        _error(3498)
    }
    if (len_data == 0) {
        if (N != floor(N) | T_val != floor(T_val) | N < 0 | T_val < 0) {
            errprintf("wfe_weights_mdid: N and T must be nonnegative integers\n")
            _error(3200)
        }
        if ((N == 0 & T_val != 0) | (N != 0 & T_val == 0)) {
            errprintf("wfe_weights_mdid: N and T must both be zero or both be positive when len_data == 0\n")
            _error(3200)
        }
        _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
            len_data, "wfe_weights_mdid")
        _wfe_validate_same_contract(same, unit_col, time_col, treat_col, N, T_val,
            len_data, "wfe_weights_mdid")
        return(J(T_val, N, 0))
    }
    if (N != floor(N) | T_val != floor(T_val) | N < 1 | T_val < 1) {
        errprintf("wfe_weights_mdid: N and T must be positive integers\n")
        _error(3200)
    }
    if (rows(exist) != T_val | cols(exist) != N) {
        errprintf("wfe_weights_mdid: exist matrix dimension mismatch\n")
        _error(3200)
    }
    if (rows(same) != T_val | cols(same) != N) {
        errprintf("wfe_weights_mdid: same matrix dimension mismatch\n")
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("wfe_weights_mdid: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col)) | ///
        any(unit_col :< 1) | any(unit_col :> N) | ///
        any(time_col :< 1) | any(time_col :> T_val)) {
        errprintf("wfe_weights_mdid: unit_idx/time_idx must be integer indices within 1..N and 1..T\n")
        _error(3200)
    }

    Wdid = J(T_val, N, 0)
    obs_lookup = J(T_val, N, 0)
    cell_count = J(T_val, N, 0)
    weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_mdid")
    is_ate = (qoi == "ate")

    for (k = 1; k <= len_data; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("wfe_weights_mdid: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }
    // Validate support metadata matches the observed panel
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
        len_data, "wfe_weights_mdid")
    _wfe_validate_same_contract(same, unit_col, time_col, treat_col, N, T_val,
        len_data, "wfe_weights_mdid")
    // The reference C generator walks the declared T by N grid and leaves
    // structurally missing rows/columns at zero. This helper should do the
    // same once the observed support itself is validated by exist/same.
    // If the declared grid has only one unit or one time period, there are no
    // admissible matched DiD contrasts, so the natural helper object is a zero
    // matrix rather than a helper-level feasibility error.

    for (t = 2; t <= T_val; t++) {
        for (i = 1; i <= N; i++) {
            if (exist[t, i] == 0) continue
            if (exist[t - 1, i] == 0) continue
            if (same[t, i] == 1) continue

            k_it = obs_lookup[t, i]
            k_it1 = obs_lookup[t - 1, i]
            if (k_it == 0 | k_it1 == 0) continue

            c_it = weight_cit[k_it]
            t_it = treat_col[k_it]
            y_it1 = y_col[k_it1]
            target_treat = 1 - t_it
            W = J(T_val, N, 0)

            if (maxdev < 0) {
                has_match = 0
                min_diff = .
                match_unit = 0

                for (k = 1; k <= len_data; k++) {
                    if (unit_col[k] == i) continue
                    if (time_col[k] != (t - 1)) continue
                    if (treat_col[k] != target_treat) continue

                    i_prime = unit_col[k]
                    if (exist[t, i_prime] == 0) continue
                    if (same[t, i_prime] != 1) continue

                    diff = abs(y_it1 - y_col[k])
                    // R sorts by (unit_idx, time_idx) before entering the C
                    // generator, and the C scan replaces ties with the last
                    // admissible candidate. At the helper boundary we should
                    // therefore stabilize exact ties on the sorted unit order
                    // rather than on the caller's arbitrary row order.
                    if (!has_match | diff < min_diff | ///
                        (diff == min_diff & i_prime > match_unit)) {
                        min_diff = diff
                        match_unit = i_prime
                        has_match = 1
                    }
                }

                if (!has_match) continue

                W[t, i] = 1
                W[t - 1, i] = 1
                W[t, match_unit] = 1
                W[t - 1, match_unit] = -1
            }
            else if (t_it == 1) {
                match_count = 0

                for (k = 1; k <= len_data; k++) {
                    if (unit_col[k] == i) continue
                    if (time_col[k] != (t - 1)) continue
                    if (treat_col[k] != 0) continue

                    i_prime = unit_col[k]
                    if (exist[t, i_prime] == 0) continue
                    if (same[t, i_prime] != 1) continue

                    diff = abs(y_it1 - y_col[k])
                    if (diff <= maxdev) {
                        match_count = match_count + 1
                    }
                }

                if (match_count == 0) continue

                v_it = 1 / match_count
                W[t, i] = 1
                W[t - 1, i] = 1

                for (k = 1; k <= len_data; k++) {
                    if (unit_col[k] == i) continue
                    if (time_col[k] != (t - 1)) continue
                    if (treat_col[k] != 0) continue

                    i_prime = unit_col[k]
                    if (exist[t, i_prime] == 0) continue
                    if (same[t, i_prime] != 1) continue

                    diff = abs(y_it1 - y_col[k])
                    if (diff <= maxdev) {
                        W[t, i_prime] = v_it
                        W[t - 1, i_prime] = -v_it
                    }
                }
            }
            else if (t_it == 0) {
                match_count = 0

                for (k = 1; k <= len_data; k++) {
                    if (unit_col[k] == i) continue
                    if (time_col[k] != (t - 1)) continue
                    if (treat_col[k] != 1) continue

                    i_prime = unit_col[k]
                    if (exist[t, i_prime] == 0) continue
                    if (same[t, i_prime] != 1) continue

                    diff = abs(y_it1 - y_col[k])
                    if (diff <= maxdev) {
                        match_count = match_count + 1
                    }
                }

                if (match_count == 0) continue

                v_it = 1 / match_count
                W[t, i] = 1
                W[t - 1, i] = 1

                for (k = 1; k <= len_data; k++) {
                    if (unit_col[k] == i) continue
                    if (time_col[k] != (t - 1)) continue
                    if (treat_col[k] != 1) continue

                    i_prime = unit_col[k]
                    if (exist[t, i_prime] == 0) continue
                    if (same[t, i_prime] != 1) continue

                    diff = abs(y_it1 - y_col[k])
                    if (diff <= maxdev) {
                        W[t, i_prime] = v_it
                        W[t - 1, i_prime] = -v_it
                    }
                }
            }

            scale = c_it
            if (!is_ate) {
                scale = scale * t_it
            }
            Wdid = Wdid + scale * W
        }
    }

    return(Wdid)
}

end
