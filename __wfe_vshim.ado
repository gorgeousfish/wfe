*! __wfe_vshim.ado -- Shared Mata compatibility shim for command-level views

program define __wfe_vshim
    version 16
end

capture mata: mata drop _wfe_as_colvector()
capture mata: mata drop __wfe_vshim_version()
capture mata: mata drop _wfe_precompute_bread()
mata:
real colvector _wfe_as_colvector(real vector x)
{
    real matrix x_copy

    x_copy = x
    return(colshape(x_copy, 1))
}

real scalar __wfe_vshim_version()
{
    return(2)
}

void _wfe_precompute_bread(
    real matrix   X_tilde,
    real matrix   X_hat,
    real matrix   ginv_XX_tilde,
    real matrix   ginv_XX_hat,
    | real scalar NT
)
{
    if (rows(X_tilde) == 1 & cols(X_tilde) > 1) {
        X_tilde = X_tilde'
    }
    if (rows(X_hat) == 1 & cols(X_hat) > 1) {
        X_hat = X_hat'
    }

    if (args() == 5) {
        if (NT >= . | NT != floor(NT) | NT <= 0) {
            errprintf("_wfe_precompute_bread: NT must be a positive integer when provided\n")
            _error(3200)
        }
        if (rows(X_tilde) != NT | rows(X_hat) != NT) {
            errprintf("_wfe_precompute_bread: X_tilde/X_hat rows must equal NT when NT is provided\n")
            _error(3200)
        }
    }

    if (rows(X_tilde) != rows(X_hat)) {
        errprintf("_wfe_precompute_bread: X_tilde (%g rows) and X_hat (%g rows) mismatch\n",
                  rows(X_tilde), rows(X_hat))
        _error(3200)
    }
    if (cols(X_tilde) != cols(X_hat)) {
        errprintf("_wfe_precompute_bread: X_tilde (%g cols) and X_hat (%g cols) mismatch\n",
                  cols(X_tilde), cols(X_hat))
        _error(3200)
    }
    if (rows(X_tilde) == 0 | rows(X_hat) == 0) {
        errprintf("_wfe_precompute_bread: X_tilde and X_hat must each contain at least one observation\n")
        _error(3200)
    }
    if (cols(X_tilde) == 0 | cols(X_hat) == 0) {
        errprintf("_wfe_precompute_bread: X_tilde and X_hat must each contain at least one regressor\n")
        _error(3200)
    }
    if (any(X_tilde :>= .) | any(X_hat :>= .)) {
        errprintf("_wfe_precompute_bread: X_tilde and X_hat must not contain missing values\n")
        _error(3498)
    }

    ginv_XX_tilde = pinv(cross(X_tilde, X_tilde))
    ginv_XX_hat   = pinv(cross(X_hat,   X_hat))
}
end
