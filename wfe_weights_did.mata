// wfe_weights_did.mata — DiD weight matrix computation

version 16.0

mata:
mata set matastrict on

real matrix wfe_weights_did(
    real vector unit_idx,
    real vector time_idx,
    real vector treat,
    real vector C_it,
    real scalar N,
    real scalar T_val,
    real scalar len_data,
    string scalar qoi,
    real matrix exist,
    real matrix same
)
{
    real matrix Wdid, obs_lookup, cell_count
    real colvector cf_units, weight_cit
    real colvector unit_col, time_col, treat_col, cit_col
    real scalar t, i, i_prime, r, k, k_prime, ui, ti
    real scalar t_it, c_it, n_cf, v_it, scale
    real scalar is_ate

    if (qoi != "ate" & qoi != "att") {
        errprintf("wfe_weights_did: qoi must be 'ate' or 'att'\n")
        _error(3200)
    }
    if (missing(len_data) | len_data != floor(len_data) | len_data < 0) {
        errprintf("wfe_weights_did: len_data must be a nonnegative integer\n")
        _error(3200)
    }
    unit_col = _wfe_as_colvector(unit_idx)
    time_col = _wfe_as_colvector(time_idx)
    treat_col = _wfe_as_colvector(treat)
    cit_col = _wfe_as_colvector(C_it)
    if (rows(unit_col) != len_data | rows(time_col) != len_data | ///
        rows(treat_col) != len_data | rows(cit_col) != len_data) {
        errprintf("wfe_weights_did: input vector length mismatch\n")
        _error(3200)
    }
    if (any(treat_col :>= .)) {
        errprintf("wfe_weights_did: treat must not contain missing values\n")
        _error(3498)
    }
    if (any((treat_col :!= 0) :& (treat_col :!= 1))) {
        errprintf("wfe_weights_did: treat must contain only 0/1 values\n")
        _error(3498)
    }
    if (len_data == 0) {
        if (N != floor(N) | T_val != floor(T_val) | N < 0 | T_val < 0) {
            errprintf("wfe_weights_did: N and T must be nonnegative integers\n")
            _error(3200)
        }
        if ((N == 0 & T_val != 0) | (N != 0 & T_val == 0)) {
            errprintf("wfe_weights_did: N and T must both be zero or both be positive when len_data == 0\n")
            _error(3200)
        }
        weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_did")
        _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
            len_data, "wfe_weights_did")
        _wfe_validate_same_contract(same, unit_col, time_col, treat_col, N, T_val,
            len_data, "wfe_weights_did")
        return(J(T_val, N, 0))
    }
    if (N != floor(N) | T_val != floor(T_val) | N < 1 | T_val < 1) {
        errprintf("wfe_weights_did: N and T must be positive integers\n")
        _error(3200)
    }
    if (rows(exist) != T_val | cols(exist) != N) {
        errprintf("wfe_weights_did: exist matrix dimension mismatch\n")
        _error(3200)
    }
    if (rows(same) != T_val | cols(same) != N) {
        errprintf("wfe_weights_did: same matrix dimension mismatch\n")
        _error(3200)
    }
    if (any(unit_col :>= .) | any(time_col :>= .)) {
        errprintf("wfe_weights_did: unit_idx/time_idx must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :!= floor(unit_col)) | any(time_col :!= floor(time_col)) | ///
        any(unit_col :< 1) | any(unit_col :> N) | ///
        any(time_col :< 1) | any(time_col :> T_val)) {
        errprintf("wfe_weights_did: unit_idx/time_idx must be integer indices within 1..N and 1..T\n")
        _error(3200)
    }

    Wdid = J(T_val, N, 0)
    obs_lookup = J(T_val, N, 0)
    cell_count = J(T_val, N, 0)
    cf_units = J(N, 1, 0)
    weight_cit = _wfe_trunc_validate_cit(cit_col, "wfe_weights_did")
    is_ate = (qoi == "ate")

    for (k = 1; k <= len_data; k++) {
        ui = unit_col[k]
        ti = time_col[k]
        cell_count[ti, ui] = cell_count[ti, ui] + 1
        if (cell_count[ti, ui] > 1) {
            errprintf("wfe_weights_did: unit-time pair is not unique\n")
            _error(498)
        }
        obs_lookup[ti, ui] = k
    }
    _wfe_validate_exist_contract(exist, unit_col, time_col, N, T_val,
        len_data, "wfe_weights_did")
    _wfe_validate_same_contract(same, unit_col, time_col, treat_col, N, T_val,
        len_data, "wfe_weights_did")
    // Unlike the command-level compact indices, this public helper receives an
    // explicit declared grid (N by T) plus support masks. Missing whole
    // rows/columns are therefore legitimate and should remain zero in Wdid.
    // If the declared grid has only one unit or one time period, there are
    // simply no admissible DiD contrasts and the helper should return zeros;
    // command-level unidentified-support guards belong upstream.

    for (t = 2; t <= T_val; t++) {
        for (i = 1; i <= N; i++) {
            if (exist[t, i] == 0) continue
            if (exist[t - 1, i] == 0) continue
            if (same[t, i] == 1) continue

            k = obs_lookup[t, i]
            if (k == 0) continue

            t_it = treat_col[k]
            c_it = weight_cit[k]
            n_cf = 0

            for (i_prime = 1; i_prime <= N; i_prime++) {
                if (i_prime == i) continue
                if (exist[t, i_prime] == 0) continue

                k_prime = obs_lookup[t, i_prime]
                if (k_prime == 0) continue
                if (treat_col[k_prime] == t_it) continue
                if (exist[t - 1, i_prime] == 0) continue
                if (same[t, i_prime] != 1) continue

                n_cf = n_cf + 1
                cf_units[n_cf] = i_prime
            }

            if (n_cf == 0) continue

            v_it = 1 / n_cf
            scale = c_it
            if (!is_ate) {
                scale = scale * t_it
            }

            Wdid[t, i] = Wdid[t, i] + scale
            Wdid[t - 1, i] = Wdid[t - 1, i] + scale

            for (r = 1; r <= n_cf; r++) {
                i_prime = cf_units[r]
                Wdid[t, i_prime] = Wdid[t, i_prime] + scale * v_it
                Wdid[t - 1, i_prime] = Wdid[t - 1, i_prime] - scale * v_it
            }
        }
    }

    return(Wdid)
}

end
