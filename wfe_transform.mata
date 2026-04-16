// ============================================================
// wfe_transform.mata - Propensity-score weighted outcome transformation
//
// Implements IPW transformation: reweights outcomes by inverse propensity
// scores to adjust for selection bias in observational data.
// ============================================================

version 16.0
mata:

void _wfe_transform_require_vector(real matrix x, string scalar argname)
{
    if (rows(x) != 1 & cols(x) != 1) {
        errprintf("wfe_transform: %s must be a vector\n", argname)
        _error(3200)
    }
}

real colvector wfe_transform(real matrix Y, real matrix treat,
                             real matrix pscore)
{
    real scalar n, i, sumTreat, psDenom1, psDenom0
    real colvector Y_col, treat_col, pscore_col, Y_star

    _wfe_transform_require_vector(Y, "outcome")
    _wfe_transform_require_vector(treat, "treat")
    _wfe_transform_require_vector(pscore, "pscore")

    // Normalize inputs to column vectors
    Y_col = colshape(Y[., .], 1)
    treat_col = colshape(treat[., .], 1)
    pscore_col = colshape(pscore[., .], 1)

    n = rows(Y_col)
    if (rows(treat_col) != n | rows(pscore_col) != n) {
        errprintf("wfe_transform: outcome, treat, and pscore must have the same length\n")
        _error(3200)
    }
    Y_star = J(n, 1, 0)
    sumTreat = 0
    psDenom1 = 0
    psDenom0 = 0

    for (i = 1; i <= n; i++) {
        if (Y_col[i] >= .) {
            errprintf("wfe_transform: outcome must not contain missing values\n")
            _error(3498)
        }
        if (treat_col[i] >= .) {
            errprintf("wfe_transform: treat must not contain missing values\n")
            _error(3498)
        }
        if (treat_col[i] != 0 & treat_col[i] != 1) {
            errprintf("wfe_transform: treat must contain only 0/1 values\n")
            _error(3498)
        }
        if (pscore_col[i] >= .) {
            errprintf("wfe_transform: pscore must not contain missing values\n")
            _error(3498)
        }
        if (pscore_col[i] <= 0 | pscore_col[i] >= 1) {
            errprintf("wfe_transform: pscore must contain only values strictly between 0 and 1\n")
            _error(3498)
        }
        sumTreat = sumTreat + treat_col[i]
        if (treat_col[i] == 1) {
            psDenom1 = psDenom1 + 1 / pscore_col[i]
        }
        else {
            psDenom0 = psDenom0 + 1 / (1 - pscore_col[i])
        }
    }

    // All-same treat (sumTreat==0 or ==n) is allowed in within_unit mode:
    // no observations of the missing level exist, so the corresponding
    // psDenom division is never reached.

    for (i = 1; i <= n; i++) {
        if (treat_col[i] == 1) {
            Y_star[i] = Y_col[i] * sumTreat / (pscore_col[i] * psDenom1)
        }
        else {
            Y_star[i] = Y_col[i] * (n - sumTreat) / ((1 - pscore_col[i]) * psDenom0)
        }
    }

    return(Y_star)
}

end
