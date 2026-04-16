version 16.0
mata:
mata set matastrict on

struct wfe_weight_summary {
    real scalar T
    real scalar N
    real scalar total
    real scalar n_nonzero
    real scalar n_positive
    real scalar n_negative
    real scalar w_min
    real scalar w_max
    real scalar w_mean
    real scalar w_sd
    real scalar nz_min
    real scalar nz_max
    real scalar nz_mean
    real scalar nz_sd
    real scalar pos_min
    real scalar pos_max
    real scalar pos_mean
    real scalar pos_sd
    real scalar neg_min
    real scalar neg_max
    real scalar neg_mean
    real scalar neg_sd
    real scalar neg_ratio
}

real scalar wfe_post_sample_sd(real colvector values)
{
    real scalar n
    real scalar mu

    n = rows(values)
    if (n < 2) {
        return(.)
    }

    mu = mean(values)
    return(sqrt(quadcross(values :- mu, values :- mu) / (n - 1)))
}

struct wfe_weight_summary scalar wfe_summarize_weights(real matrix W)
{
    struct wfe_weight_summary scalar summary
    real colvector weights
    real colvector nz
    real colvector pos
    real colvector neg

    weights = colshape(W, 1)
    nz = select(weights, weights :!= 0)
    pos = select(weights, weights :> 0)
    neg = select(weights, weights :< 0)

    summary.T = rows(W)
    summary.N = cols(W)
    summary.total = rows(weights)
    summary.n_nonzero = rows(nz)
    summary.n_positive = rows(pos)
    summary.n_negative = rows(neg)
    summary.w_min = min(weights)
    summary.w_max = max(weights)
    summary.w_mean = mean(weights)
    summary.w_sd = wfe_post_sample_sd(weights)

    if (rows(nz) > 0) {
        summary.nz_min = min(nz)
        summary.nz_max = max(nz)
        summary.nz_mean = mean(nz)
        summary.nz_sd = wfe_post_sample_sd(nz)
    }
    else {
        summary.nz_min = .
        summary.nz_max = .
        summary.nz_mean = .
        summary.nz_sd = .
    }

    if (rows(pos) > 0) {
        summary.pos_min = min(pos)
        summary.pos_max = max(pos)
        summary.pos_mean = mean(pos)
        summary.pos_sd = wfe_post_sample_sd(pos)
    }
    else {
        summary.pos_min = .
        summary.pos_max = .
        summary.pos_mean = .
        summary.pos_sd = .
    }

    if (rows(neg) > 0) {
        summary.neg_min = min(neg)
        summary.neg_max = max(neg)
        summary.neg_mean = mean(neg)
        summary.neg_sd = wfe_post_sample_sd(neg)
        summary.neg_ratio = rows(neg) / rows(nz)
    }
    else {
        summary.neg_min = .
        summary.neg_max = .
        summary.neg_mean = .
        summary.neg_sd = .
        summary.neg_ratio = .
    }

    return(summary)
}

void wfe_store_weight_summary(real matrix W)
{
    struct wfe_weight_summary scalar summary

    summary = wfe_summarize_weights(W)
    st_numscalar("__wfe_sum_T", summary.T)
    st_numscalar("__wfe_sum_N", summary.N)
    st_numscalar("__wfe_sum_total", summary.total)
    st_numscalar("__wfe_sum_n_nonzero", summary.n_nonzero)
    st_numscalar("__wfe_sum_n_positive", summary.n_positive)
    st_numscalar("__wfe_sum_n_negative", summary.n_negative)
    st_numscalar("__wfe_sum_w_min", summary.w_min)
    st_numscalar("__wfe_sum_w_max", summary.w_max)
    st_numscalar("__wfe_sum_w_mean", summary.w_mean)
    st_numscalar("__wfe_sum_w_sd", summary.w_sd)
    st_numscalar("__wfe_sum_nz_min", summary.nz_min)
    st_numscalar("__wfe_sum_nz_max", summary.nz_max)
    st_numscalar("__wfe_sum_nz_mean", summary.nz_mean)
    st_numscalar("__wfe_sum_nz_sd", summary.nz_sd)
    st_numscalar("__wfe_sum_pos_min", summary.pos_min)
    st_numscalar("__wfe_sum_pos_max", summary.pos_max)
    st_numscalar("__wfe_sum_pos_mean", summary.pos_mean)
    st_numscalar("__wfe_sum_pos_sd", summary.pos_sd)
    st_numscalar("__wfe_sum_neg_min", summary.neg_min)
    st_numscalar("__wfe_sum_neg_max", summary.neg_max)
    st_numscalar("__wfe_sum_neg_mean", summary.neg_mean)
    st_numscalar("__wfe_sum_neg_sd", summary.neg_sd)
    st_numscalar("__wfe_sum_neg_ratio", summary.neg_ratio)
}

end
