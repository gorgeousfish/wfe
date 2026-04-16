// ============================================================
// wfe_se_hac.mata — HAC/HC 标准误计算 (单向 FE)
//
// 无外部 Mata 依赖
//
// 函数清单 (5):
//   _wfe_precompute_bread()    — ginv(X'X) 预计算
//   _wfe_omega_hac()           — Arellano CRVE Omega (按单位聚类)
//   _wfe_omega_hc()            — HC Omega (向量化)
//   _wfe_stockwatson_bias()    — Stock-Watson 偏差修正
//   _wfe_compute_se()          — 公共入口: 路由 + 组装 vcov
//
// 适用范围: method="unit", method="time", estimator="fd"
// ⚠️ 不适用于 estimator="did" / "Mdid" (使用 wfe_se_gmm.mata)
// ============================================================

version 16.0
mata:
mata set matastrict on

real colvector _wfe_se_as_colvector(real vector x)
{
    return(colshape(x[., .], 1))
}

real matrix _wfe_precompute_as_design(real matrix X, | real scalar n_obs)
{
    if (args() == 2) {
        if (n_obs >= . | n_obs != floor(n_obs) | n_obs <= 0) {
            errprintf("_wfe_precompute_as_design: n_obs must be a positive integer when provided\n")
            _error(3200)
        }
        if (rows(X) == 1 & cols(X) == n_obs & n_obs > 1) {
            return(X')
        }
        return(X)
    }

    return(X)
}

// ============================================================
// _wfe_precompute_bread() — 广义逆预计算
//
// @param X_tilde       real matrix [NT×p]  WFE 加权去均值后的 X
// @param X_hat         real matrix [NT×p]  标准 FE 去均值后的 X
// @param ginv_XX_tilde real matrix [out]   pinv(X_tilde'X_tilde)
// @param ginv_XX_hat   real matrix [out]   pinv(X_hat'X_hat)
// ============================================================
void _wfe_precompute_bread(
    real matrix   X_tilde,
    real matrix   X_hat,
    real matrix   ginv_XX_tilde,
    real matrix   ginv_XX_hat,
    | real scalar NT
)
{
    // cross(X, X) expects observation-by-regressor matrices. For a
    // one-regressor design, a 1 x NT rowvector is just an alternate storage of
    // the same 1D observation stream and should be canonicalized to NT x 1
    // before computing X'X.
    if (args() == 5) {
        if (NT >= . | NT != floor(NT) | NT <= 0) {
            errprintf("_wfe_precompute_bread: NT must be a positive integer when provided\n")
            _error(3200)
        }
        X_tilde = _wfe_precompute_as_design(X_tilde, NT)
        X_hat   = _wfe_precompute_as_design(X_hat, NT)
        if (rows(X_tilde) != NT | rows(X_hat) != NT) {
            errprintf("_wfe_precompute_bread: X_tilde/X_hat rows must equal NT when NT is provided\n")
            _error(3200)
        }
    }

    // 维度断言
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

    // ginv(X'X): cross(X,X) = X'X (BLAS 优化)，pinv = Moore-Penrose 广义逆
    ginv_XX_tilde = pinv(cross(X_tilde, X_tilde))
    ginv_XX_hat   = pinv(cross(X_hat,   X_hat))
}


real colvector _wfe_validate_se_unit_idx(
    real vector    unit_idx,
    real scalar    J_u,
    string scalar  caller
)
{
    real scalar i, n_obs, n_unique
    real colvector unit_col, compact_idx
    real colvector sort_order, unit_sorted
    real scalar prev_label

    if (J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("%s: J_u must be a positive integer\n", caller)
        _error(3200)
    }

    unit_col = _wfe_se_as_colvector(unit_idx)

    if (any(unit_idx :>= .)) {
        errprintf("%s: unit_idx must not contain missing values\n", caller)
        _error(3498)
    }

    for (i = 1; i <= rows(unit_col); i++) {
        if (unit_col[i] != floor(unit_col[i]) | unit_col[i] < 1) {
            errprintf("%s: unit_idx must contain positive integer labels\n", caller)
            _error(3200)
        }
    }

    n_obs = rows(unit_col)
    if (n_obs == 0) {
        errprintf("%s: unit_idx must not be empty when J_u is positive\n", caller)
        _error(3200)
    }

    sort_order = _wfe_se_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    compact_idx = J(n_obs, 1, .)
    n_unique = 0
    prev_label = .

    for (i = 1; i <= n_obs; i++) {
        if (i == 1 | unit_sorted[i] != prev_label) {
            n_unique++
            prev_label = unit_sorted[i]
        }
        compact_idx[sort_order[i]] = n_unique
    }

    if (n_unique != J_u) {
        errprintf("%s: unit_idx must identify exactly J_u unique units\n", caller)
        _error(3200)
    }

    return(compact_idx)
}

real colvector _wfe_se_unit_sort_order(real colvector unit_idx)
{
    real scalar n

    n = rows(unit_idx)
    if (n <= 1) {
        return(1::n)
    }

    // HAC/Stock-Watson aggregate by unit clusters only, so row order is not
    // part of the public contract. Sort local copies to satisfy panelsetup()
    // without changing the statistic.
    return(order((unit_idx, (1::n)), (1, 2)))
}


void _wfe_validate_sw_balance(
    real vector    unit_idx,
    real scalar    J_u
)
{
    real matrix info
    real colvector counts
    real colvector sort_order, unit_sorted, unit_col

    unit_col = _wfe_se_as_colvector(unit_idx)

    if (J_u >= . | J_u != floor(J_u) | J_u < 0) {
        errprintf("_wfe_validate_sw_balance: J_u must be a nonnegative integer\n")
        _error(3200)
    }

    if (J_u == 0) {
        if (rows(unit_col) != 0) {
            errprintf("_wfe_validate_sw_balance: J_u must be positive when unit_idx is nonempty\n")
            _error(3200)
        }
        return
    }

    unit_col = _wfe_validate_se_unit_idx(unit_col, J_u, "_wfe_validate_sw_balance")
    sort_order = _wfe_se_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    info = panelsetup(unit_sorted, 1)
    counts = info[., 2] - info[., 1] :+ 1

    if (rows(counts) != J_u | min(counts) != max(counts)) {
        errprintf("unbiased_se(on) is allowed only when panel is balanced\n")
        _error(498)
    }
}


// ============================================================
// _wfe_omega_hac() — Arellano (1987) CRVE Omega
//
// 计算: Ω_raw = Σᵢ Xᵢ'uᵢuᵢ'Xᵢ  (聚类稳健方差估计器核心)
//
// 优化: 使用秩一更新 Xuᵢ * Xuᵢ' 替代完整外积 Xᵢ'uᵢuᵢ'Xᵢ
//   数学等价: (Xᵢ'uᵢ)(uᵢ'Xᵢ) = Xᵢ'(uᵢuᵢ')Xᵢ
//   复杂度: O(Tᵢp + p²) vs O(Tᵢ²p + Tᵢp²)
//
// @param X         real matrix    [NT×p]  去均值后的 X
// @param u         real colvector [NT×1]  残差
// @param unit_idx  real colvector [NT×1]  单位索引 (连续整数)
// @param J_u       real scalar           单位数
// @return          real matrix    [p×p]   原始 Omega (未缩放)
// ============================================================
real matrix _wfe_omega_hac(
    real matrix    X,
    real vector    u,
    real vector    unit_idx,
    real scalar    J_u
)
{
    real scalar    p, i
    real matrix    Omega, info
    real matrix    Xi, X_sorted
    real colvector ui, Xui, sort_order, unit_sorted, u_sorted, u_col, unit_col

    u_col = _wfe_se_as_colvector(u)
    unit_col = _wfe_se_as_colvector(unit_idx)
    X = _wfe_precompute_as_design(X, rows(u_col))

    if (rows(X) != rows(u_col)) {
        errprintf("_wfe_omega_hac: X (%g rows) and u (%g rows) mismatch\n",
                  rows(X), rows(u_col))
        _error(3200)
    }
    if (rows(unit_col) != rows(X)) {
        errprintf("_wfe_omega_hac: unit_idx (%g rows) mismatch with X (%g rows)\n",
                  rows(unit_col), rows(X))
        _error(3200)
    }
    if (rows(X) < 1) {
        errprintf("_wfe_omega_hac: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (cols(X) < 1) {
        errprintf("_wfe_omega_hac: X must contain at least one regressor\n")
        _error(3200)
    }
    if (any(X :>= .) | any(u_col :>= .)) {
        errprintf("_wfe_omega_hac: X and u must not contain missing values\n")
        _error(3498)
    }
    if (any(unit_col :>= .)) {
        errprintf("_wfe_omega_hac: unit_idx must not contain missing values\n")
        _error(3498)
    }
    unit_col = _wfe_validate_se_unit_idx(unit_col, J_u, "_wfe_omega_hac")

    p     = cols(X)
    Omega = J(p, p, 0)
    sort_order = _wfe_se_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    X_sorted = X[sort_order, .]
    u_sorted = u_col[sort_order]

    // panelsetup 按 unit_idx 分组
    info = panelsetup(unit_sorted, 1)

    // 遍历每个单位，秩一更新累加
    for (i = 1; i <= rows(info); i++) {
        // 提取单位 i 的子矩阵/子向量
        Xi = panelsubmatrix(X_sorted, i, info)      // T_i × p
        ui = panelsubmatrix(u_sorted, i, info)      // T_i × 1

        // 秩一更新: Xui = Xi'ui (p×1), Omega += Xui * Xui'
        // 优化: 跳过完整外积 uu_i[T_i × T_i]，复杂度从 O(Tᵢ²p) 降至 O(Tᵢp)
        Xui   = cross(Xi, ui)               // p × 1
        Omega = Omega + Xui * Xui'           // p × p 秩一更新
    }

    // 返回未缩放的原始 Omega
    // ⚠️ df_adj 缩放在调用者 _wfe_compute_se() 中执行
    return(Omega)
}


// ============================================================
// _wfe_omega_hac_cross() — 不对称交叉 HAC meat (按单位聚类)
//
// 计算: Ω_cross = Σᵢ (Σ_{t:i} X_a[t] u_a[t]) × (Σ_{t:i} X_b[t] u_b[t])'
//
// 用于 White test 的 Lambda12 = B_wfe * (Ω_cross / J_u) * B_fe
// 保证 joint [[Ω_wfe, Ω_cross], [Ω_cross', Ω_fe]] 半正定
//
// @param X_a       real matrix    [NT×p]   第一组设计矩阵 (如 X_tilde)
// @param u_a       real vector    [NT×1]   第一组残差
// @param X_b       real matrix    [NT×p]   第二组设计矩阵 (如 X_hat)
// @param u_b       real vector    [NT×1]   第二组残差
// @param unit_idx  real vector    [NT×1]   单位标识
// @param J_u       real scalar             单位数
// @return          real matrix    [p×p]    原始 Omega_cross (未缩放)
// ============================================================
real matrix _wfe_omega_hac_cross(
    real matrix    X_a,
    real vector    u_a,
    real matrix    X_b,
    real vector    u_b,
    real vector    unit_idx,
    real scalar    J_u
)
{
    real scalar    p, i
    real matrix    Omega, info
    real matrix    X_a_sorted, X_b_sorted
    real colvector u_a_sorted, u_b_sorted, sort_order, unit_sorted
    real colvector u_a_col, u_b_col, unit_col
    real colvector score_a, score_b

    u_a_col = _wfe_se_as_colvector(u_a)
    u_b_col = _wfe_se_as_colvector(u_b)
    unit_col = _wfe_se_as_colvector(unit_idx)
    X_a = _wfe_precompute_as_design(X_a, rows(u_a_col))
    X_b = _wfe_precompute_as_design(X_b, rows(u_b_col))

    if (rows(X_a) != rows(u_a_col) | rows(X_b) != rows(u_b_col)) {
        errprintf("_wfe_omega_hac_cross: X and u row mismatch\n")
        _error(3200)
    }
    if (rows(X_a) != rows(X_b)) {
        errprintf("_wfe_omega_hac_cross: X_a (%g rows) and X_b (%g rows) mismatch\n",
                  rows(X_a), rows(X_b))
        _error(3200)
    }
    if (cols(X_a) != cols(X_b)) {
        errprintf("_wfe_omega_hac_cross: X_a (%g cols) and X_b (%g cols) mismatch\n",
                  cols(X_a), cols(X_b))
        _error(3200)
    }
    if (rows(unit_col) != rows(X_a)) {
        errprintf("_wfe_omega_hac_cross: unit_idx rows mismatch\n")
        _error(3200)
    }
    if (rows(X_a) < 1 | rows(X_b) < 1) {
        errprintf("_wfe_omega_hac_cross: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (cols(X_a) < 1 | cols(X_b) < 1) {
        errprintf("_wfe_omega_hac_cross: X_a and X_b must each contain at least one regressor\n")
        _error(3200)
    }
    if (any(X_a :>= .) | any(X_b :>= .) | ///
        any(u_a_col :>= .) | any(u_b_col :>= .)) {
        errprintf("_wfe_omega_hac_cross: inputs must not contain missing values\n")
        _error(3498)
    }
    unit_col = _wfe_validate_se_unit_idx(unit_col, J_u, "_wfe_omega_hac_cross")

    p     = cols(X_a)
    Omega = J(p, p, 0)
    sort_order = _wfe_se_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    X_a_sorted = X_a[sort_order, .]
    X_b_sorted = X_b[sort_order, .]
    u_a_sorted = u_a_col[sort_order]
    u_b_sorted = u_b_col[sort_order]

    info = panelsetup(unit_sorted, 1)

    for (i = 1; i <= rows(info); i++) {
        score_a = cross(panelsubmatrix(X_a_sorted, i, info),
                        panelsubmatrix(u_a_sorted, i, info))   // p × 1
        score_b = cross(panelsubmatrix(X_b_sorted, i, info),
                        panelsubmatrix(u_b_sorted, i, info))   // p × 1
        Omega = Omega + score_a * score_b'                     // p × p
    }

    return(Omega)
}


// ============================================================
// _wfe_omega_hc() — HC (异方差一致性) Omega
//
// 计算: Ω_HC = (1/J_u) * X' diag(u²) X
//            = (1/J_u) * Σₖ uₖ² xₖ xₖ'
//
// 缩放: 已含 (1/J_u) 因子，使 bread 可直接乘以原始 Omega。
//
// @param X       real matrix    [NT×p]  去均值后的 X
// @param diag_ee real vector    [NT×1 or 1×NT]  u² (残差平方)
// @param J_u     real scalar           单位数
// @return        real matrix    [p×p]   已含 1/J_u 缩放的 Omega
// ============================================================
real matrix _wfe_omega_hc(
    real matrix    X,
    real vector    diag_ee,
    real scalar    J_u
)
{
    real colvector diag_ee_col

    diag_ee_col = _wfe_se_as_colvector(diag_ee)
    X = _wfe_precompute_as_design(X, rows(diag_ee_col))

    // 维度断言
    if (rows(X) != rows(diag_ee_col)) {
        errprintf("_wfe_omega_hc: X (%g rows) and diag_ee (%g rows) mismatch\n",
                  rows(X), rows(diag_ee_col))
        _error(3200)
    }
    if (rows(X) < 1) {
        errprintf("_wfe_omega_hc: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (cols(X) < 1) {
        errprintf("_wfe_omega_hc: X must contain at least one regressor\n")
        _error(3200)
    }
    if (any(X :>= .)) {
        errprintf("_wfe_omega_hc: X must not contain missing values\n")
        _error(3498)
    }
    if (any(diag_ee_col :>= .)) {
        errprintf("_wfe_omega_hc: diag_ee must not contain missing values\n")
        _error(3498)
    }
    if (any(diag_ee_col :< 0)) {
        errprintf("_wfe_omega_hc: diag_ee must be non-negative\n")
        _error(3498)
    }
    if (J_u >= . | J_u != floor(J_u) | J_u <= 0) {
        errprintf("_wfe_omega_hc: J_u must be a positive integer\n")
        _error(3200)
    }
    if (J_u > rows(X)) {
        errprintf("_wfe_omega_hc: J_u must not exceed the number of observations\n")
        _error(3200)
    }

    // 向量化实现 (无需 panelsetup 循环)
    // X :* diag_ee: Mata 广播语义 — diag_ee[k] 乘以 X[k,.] 整行
    // cross(A, B) = A'B (BLAS 优化)
    return((1/J_u) * cross(X :* diag_ee_col, X))
}


void _wfe_sw_refchk(
    real matrix    X,
    real colvector u_col,
    real scalar    J_u,
    real matrix    Omega_HC,
    real matrix    ginv_XX
)
{
    real matrix omega_expected, ginv_expected
    real scalar rel_omega, rel_ginv
    real scalar omega_scale, ginv_scale

    omega_expected = _wfe_omega_hc(X, u_col :^ 2, J_u)
    ginv_expected  = pinv(cross(X, X))

    omega_scale = max((1, max(abs(vec(omega_expected)))))
    rel_omega = max(abs(vec(Omega_HC - omega_expected))) / omega_scale
    if (rel_omega > sqrt(epsilon(1))) {
        errprintf("_wfe_stockwatson_bias: Omega_HC must match _wfe_omega_hc(X, u:^2, J_u)\n")
        _error(498)
    }

    ginv_scale = max((1, max(abs(vec(ginv_expected)))))
    rel_ginv = max(abs(vec(ginv_XX - ginv_expected))) / ginv_scale
    if (rel_ginv > sqrt(epsilon(1))) {
        errprintf("_wfe_stockwatson_bias: ginv_XX must match pinv(X'X)\n")
        _error(498)
    }
}


// ============================================================
// _wfe_stockwatson_bias() — Stock-Watson 偏差修正
//
// 实现: Stock-Watson (Econometrica 2008, Eq. 6) 偏差调整
//
// @param X         real matrix    [NT×p]  去均值后的 X
// @param u         real colvector [NT×1]  残差
// @param unit_idx  real colvector [NT×1]  单位索引
// @param J_u       real scalar           单位数
// @param J_t       real scalar           时间数 (仅此函数使用)
// @param Omega_HC  real matrix    [p×p]   已缩放的 HC Omega (只读)
// @param ginv_XX   real matrix    [p×p]   pinv(X'X) (只读)
// @param Psi_sw    real matrix    [out]   Stock-Watson Psi (虽被覆盖)
// ============================================================
void _wfe_stockwatson_bias(
    real matrix    X,
    real vector    u,
    real vector    unit_idx,
    real scalar    J_u,
    real scalar    J_t,
    real matrix    Omega_HC,
    real matrix    ginv_XX,
    real matrix    Psi_sw
)
{
    real scalar    p, i, T_i, sum_u2_i
    real matrix    B_hat, info, Xi, XX_i, Sigma_HRFE, Bread, X_sorted
    real matrix    Omega_HC_use, ginv_XX_use
    real colvector ui, block_lengths, sort_order, unit_sorted, u_sorted
    real colvector u_col, unit_col
    real rowvector evals, omega_evals
    real scalar    asym, scale, omega_eig_tol, omega_eig_scale

    u_col = _wfe_se_as_colvector(u)
    unit_col = _wfe_se_as_colvector(unit_idx)
    X = _wfe_precompute_as_design(X, rows(u_col))

    if (rows(unit_col) != rows(X)) {
        errprintf("_wfe_stockwatson_bias: unit_idx (%g rows) mismatch with X (%g rows)\n",
                  rows(unit_col), rows(X))
        _error(3200)
    }
    if (rows(u_col) != rows(X)) {
        errprintf("_wfe_stockwatson_bias: X (%g rows) and u (%g rows) mismatch\n",
                  rows(X), rows(u_col))
        _error(3200)
    }
    if (rows(X) == 0) {
        errprintf("_wfe_stockwatson_bias: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (any(X :>= .) | any(u_col :>= .)) {
        errprintf("_wfe_stockwatson_bias: X and u must not contain missing values\n")
        _error(3498)
    }
    if (cols(X) < 1) {
        errprintf("_wfe_stockwatson_bias: X must contain at least one regressor\n")
        _error(3200)
    }

    p     = cols(X)
    B_hat = J(p, p, 0)

    if (rows(Omega_HC) != p | cols(Omega_HC) != p) {
        errprintf("_wfe_stockwatson_bias: Omega_HC must be %gx%g to match cols(X)\n",
                  p, p)
        _error(3200)
    }
    if (any(Omega_HC :>= .)) {
        errprintf("_wfe_stockwatson_bias: Omega_HC must not contain missing values\n")
        _error(3498)
    }
    asym = max(abs(vec(Omega_HC - Omega_HC')))
    scale = max((1, max(abs(vec(Omega_HC)))))
    if (asym > sqrt(epsilon(1)) * scale) {
        errprintf("_wfe_stockwatson_bias: Omega_HC must be symmetric\n")
        _error(3200)
    }
    // Stock-Watson only depends on the symmetric covariance part. Remove
    // harmless roundoff-level skew so numerically equivalent HC inputs follow
    // the same path as their exactly symmetric counterparts.
    Omega_HC_use = 0.5 :* (Omega_HC + Omega_HC')
    symeigensystem(Omega_HC_use, ., omega_evals)
    omega_eig_scale = max((1, max(abs(omega_evals))))
    omega_eig_tol = sqrt(epsilon(1)) * omega_eig_scale
    if (min(omega_evals) < -omega_eig_tol) {
        errprintf("_wfe_stockwatson_bias: Omega_HC must be positive semidefinite\n")
        _error(3200)
    }
    if (rows(ginv_XX) != p | cols(ginv_XX) != p) {
        errprintf("_wfe_stockwatson_bias: ginv_XX must be %gx%g to match cols(X)\n",
                  p, p)
        _error(3200)
    }
    if (any(ginv_XX :>= .)) {
        errprintf("_wfe_stockwatson_bias: ginv_XX must not contain missing values\n")
        _error(3498)
    }
    asym = max(abs(vec(ginv_XX - ginv_XX')))
    scale = max((1, max(abs(vec(ginv_XX)))))
    if (asym > sqrt(epsilon(1)) * scale) {
        errprintf("_wfe_stockwatson_bias: ginv_XX must be symmetric\n")
        _error(3200)
    }
    ginv_XX_use = 0.5 :* (ginv_XX + ginv_XX')
    symeigensystem(ginv_XX_use, ., evals)
    scale = max((1, max(abs(evals))))
    if (min(evals) < -sqrt(epsilon(1)) * scale) {
        errprintf("_wfe_stockwatson_bias: ginv_XX must be positive semidefinite\n")
        _error(3200)
    }
    if (J_u != floor(J_u) | J_u <= 0) {
        errprintf("_wfe_stockwatson_bias: J_u must be a positive integer\n")
        _error(3200)
    }
    if (J_u < 2) {
        errprintf("_wfe_stockwatson_bias: need at least 2 units for Stock-Watson standard errors\n")
        _error(498)
    }
    unit_col = _wfe_validate_se_unit_idx(unit_col, J_u, "_wfe_stockwatson_bias")

    if (J_t >= . | J_t != floor(J_t) | J_t <= 0) {
        errprintf("_wfe_stockwatson_bias: J_t must be a positive integer\n")
        _error(3200)
    }

    // J_t >= 3 检查 (公式中 (J_t-1)/(J_t-2) 在 J_t=2 时除零)
    if (J_t < 3) {
        errprintf("_wfe_stockwatson_bias: J_t=%g must be >= 3\n", J_t)
        _error(498)
    }

    _wfe_sw_refchk(X, u_col, J_u, Omega_HC_use, ginv_XX_use)

    sort_order = _wfe_se_unit_sort_order(unit_col)
    unit_sorted = unit_col[sort_order]
    X_sorted = X[sort_order, .]
    u_sorted = u_col[sort_order]

    info = panelsetup(unit_sorted, 1)
    block_lengths = info[, 2] - info[, 1] :+ 1

    if (any(block_lengths :!= J_t)) {
        errprintf("_wfe_stockwatson_bias: unit_idx must define a balanced panel with exactly J_t observations per unit\n")
        _error(498)
    }

    // 遍历每个单位，构造偏差修正矩阵 B_hat
    for (i = 1; i <= rows(info); i++) {
        Xi = panelsubmatrix(X_sorted, i, info)     // T_i × p
        ui = panelsubmatrix(u_sorted, i, info)     // T_i × 1
        T_i = rows(Xi)

        // 仅对观测数 > 1 的单位执行
        if (T_i > 1) {
            XX_i     = cross(Xi, Xi)           // p × p
            sum_u2_i = cross(ui, ui)           // 标量
            B_hat    = B_hat + (1/J_t) * XX_i * (1/(J_t-1)) * sum_u2_i
        }
    }

    // 缩放
    B_hat = B_hat / J_u

    // Sigma_HRFE: Stock-Watson (2008) Eq. 6 偏差修正后的 Sigma
    // ⚠️ Omega_HC 已包含 1/J_u 缩放
    Sigma_HRFE = ((J_t-1)/(J_t-2)) * (Omega_HC_use - (1/(J_t-1)) * B_hat)

    // 最终 Psi_sw
    Bread  = J_u * ginv_XX_use
    Psi_sw = Bread * Sigma_HRFE * Bread
}


// ============================================================
// _wfe_compute_se() — 公共入口函数
//
// 路由 HAC/HC/不支持路径，组装 Psi 和 vcov 矩阵。
//
// 适用范围: 单向 FE (method="unit"/"time", estimator="fd")
// ⚠️ 不适用于双向 FE (did/Mdid) → 使用 wfe_se_gmm
//
// @param X_tilde     real matrix   [NT×p]  WFE 加权去均值后的 X
// @param u_tilde     real vector           WFE residual stream (row/col accepted)
// @param X_hat       real matrix   [NT×p]  标准 FE 去均值后的 X
// @param u_hat       real vector           FE residual stream (row/col accepted)
// @param unit_idx    real vector           Unit-cluster index stream (row/col accepted)
// @param J_u         real scalar          单位数 (始终为 unit 数)
// @param p           real scalar          自变量数
// @param J_t         real scalar          时间数 (仅 SW 用)
// @param N_nonzero   real scalar          非零权重观测数 (dfHAC 用)
// @param df_adj_on   real scalar          1=使用 dfHAC 自由度修正, 0=仅 1/J_u 缩放
// @param hetero_se   string scalar        "on"/"off"
// @param auto_se     string scalar        "on"/"off"
// @param unbiased_se string scalar        "on"/"off"
// @param is_balanced real scalar          1=平衡, 0=非平衡
//
// [输出参数 — Mata 通过引用返回]:
// @param vcov_wfe      real matrix [p×p]  WFE 方差-协方差矩阵
// @param vcov_fe       real matrix [p×p]  FE 方差-协方差矩阵
// @param Psi_hat_wfe   real matrix [p×p]  WFE Sandwich (除 J_u 前)
// @param Psi_hat_fe    real matrix [p×p]  FE Sandwich (除 J_u 前)
// @param ginv_XX_tilde real matrix [p×p]  pinv(X_tilde'X_tilde)
// @param ginv_XX_hat   real matrix [p×p]  pinv(X_hat'X_hat)
// @param se_type       string scalar      标准误类型描述
// ============================================================
void _wfe_compute_se(
    real matrix    X_tilde,
    real vector    u_tilde,
    real matrix    X_hat,
    real vector    u_hat,
    real vector    unit_idx,
    real scalar    J_u,
    real scalar    p,
    real scalar    J_t,
    real scalar    N_nonzero,
    real scalar    df_adj_on,
    string scalar  hetero_se,
    string scalar  auto_se,
    string scalar  unbiased_se,
    real scalar    is_balanced,
    real matrix    vcov_wfe,
    real matrix    vcov_fe,
    real matrix    Psi_hat_wfe,
    real matrix    Psi_hat_fe,
    real matrix    ginv_XX_tilde,
    real matrix    ginv_XX_hat,
    string scalar  se_type
)
{
    real scalar    NT, df_adj, dfHAC
    real colvector u_tilde_col, u_hat_col, unit_col
    real colvector diag_ee_tilde, diag_ee_hat
    real matrix    Omega_wfe, Omega_fe, Omega_raw_wfe, Omega_raw_fe
    real matrix    Bread_wfe, Bread_fe
    real matrix    Psi_sw_wfe, Psi_sw_fe

    u_tilde_col = _wfe_se_as_colvector(u_tilde)
    u_hat_col   = _wfe_se_as_colvector(u_hat)
    unit_col    = _wfe_se_as_colvector(unit_idx)
    X_tilde = _wfe_precompute_as_design(X_tilde, rows(u_tilde_col))
    X_hat   = _wfe_precompute_as_design(X_hat, rows(u_hat_col))
    NT = rows(X_tilde)

    // ──────────────────────────────────────
    // 步骤 0: 输入验证
    // ──────────────────────────────────────
    if (rows(X_tilde) != rows(u_tilde_col)) {
        errprintf("_wfe_compute_se: X_tilde (%g rows) and u_tilde (%g rows) mismatch\n",
                  rows(X_tilde), rows(u_tilde_col))
        _error(3200)
    }
    if (rows(X_hat) != rows(u_hat_col)) {
        errprintf("_wfe_compute_se: X_hat (%g rows) and u_hat (%g rows) mismatch\n",
                  rows(X_hat), rows(u_hat_col))
        _error(3200)
    }
    if (rows(X_tilde) != rows(X_hat)) {
        errprintf("_wfe_compute_se: X_tilde (%g rows) and X_hat (%g rows) mismatch\n",
                  rows(X_tilde), rows(X_hat))
        _error(3200)
    }
    if (NT < 1) {
        errprintf("_wfe_compute_se: inputs must contain at least one observation\n")
        _error(3200)
    }
    if (p >= . | p != floor(p) | p <= 0) {
        errprintf("_wfe_compute_se: p must be a positive integer\n")
        _error(3200)
    }
    if (cols(X_tilde) != p | cols(X_hat) != p) {
        errprintf("_wfe_compute_se: X cols (%g, %g) mismatch with p=%g\n",
                  cols(X_tilde), cols(X_hat), p)
        _error(3200)
    }
    if (rows(unit_col) != NT) {
        errprintf("_wfe_compute_se: unit_idx (%g rows) mismatch with NT=%g\n",
                  rows(unit_col), NT)
        _error(3200)
    }
    if (any(X_tilde :>= .) | any(X_hat :>= .) | ///
        any(u_tilde_col :>= .) | any(u_hat_col :>= .)) {
        errprintf("_wfe_compute_se: inputs must not contain missing values\n")
        _error(3498)
    }
    if (hetero_se != "on" & hetero_se != "off") {
        errprintf("_wfe_compute_se: hetero_se must be on or off; got %s\n",
                  hetero_se)
        _error(3200)
    }
    if (auto_se != "on" & auto_se != "off") {
        errprintf("_wfe_compute_se: auto_se must be on or off; got %s\n",
                  auto_se)
        _error(3200)
    }
    if (unbiased_se != "on" & unbiased_se != "off") {
        errprintf("_wfe_compute_se: unbiased_se must be on or off; got %s\n",
                  unbiased_se)
        _error(3200)
    }
    if (unbiased_se == "on" & ///
        (hetero_se != "on" | auto_se != "off")) {
        errprintf("_wfe_compute_se: unbiased_se(on) requires hetero_se(on) and auto_se(off)\n")
        _error(198)
    }
    if (unbiased_se == "on" & ///
        (J_t >= . | J_t != floor(J_t) | J_t <= 0)) {
        errprintf("_wfe_compute_se: J_t must be a positive integer\n")
        _error(3200)
    }
    if (df_adj_on != 0 & df_adj_on != 1) {
        errprintf("_wfe_compute_se: df_adj_on must be 0 or 1; got %g\n",
                  df_adj_on)
        _error(3200)
    }
    if (is_balanced != 0 & is_balanced != 1) {
        errprintf("_wfe_compute_se: is_balanced must be 0 or 1; got %g\n",
                  is_balanced)
        _error(3200)
    }
    if (N_nonzero >= . | N_nonzero != floor(N_nonzero) | N_nonzero < 0) {
        errprintf("_wfe_compute_se: N_nonzero must be a nonnegative integer\n")
        _error(3200)
    }
    if (N_nonzero > NT) {
        errprintf("_wfe_compute_se: N_nonzero must not exceed NT\n")
        _error(3200)
    }

    unit_col = _wfe_validate_se_unit_idx(unit_col, J_u, "_wfe_compute_se")

    // Reject unsupported one-way SE combinations before any rank-specific
    // numerical feasibility checks. Collinearity must not mask the
    // "unsupported" boundary when the SE family itself is unavailable.
    if (hetero_se == "off" & auto_se == "off") {
        errprintf("_wfe_compute_se: standard errors with independence and homoskedasticity is not supported\n")
        _error(198)
    }
    if (hetero_se == "off" & auto_se == "on") {
        errprintf("_wfe_compute_se: robust standard errors with autocorrelation and homoskedasticity is not supported\n")
        _error(198)
    }

    // Reject rank-deficient designs before pinv(X'X)-based sandwich construction;
    // aliased regressors are not supported under the WFE no-intercept design.
    if (rank(X_tilde) < p) {
        errprintf("_wfe_compute_se: X_tilde regressors are collinear\n")
        _error(498)
    }
    if (rank(X_hat) < p) {
        errprintf("_wfe_compute_se: X_hat regressors are collinear\n")
        _error(498)
    }

    // ──────────────────────────────────────
    // 步骤 1: 预计算 ginv(X'X)
    // ──────────────────────────────────────
    _wfe_precompute_bread(X_tilde, X_hat, ginv_XX_tilde, ginv_XX_hat)

    // ──────────────────────────────────────
    // 步骤 2: 预计算 diag(u²) — 在分支判断之前
    // ──────────────────────────────────────
    diag_ee_tilde = u_tilde_col :^ 2        // NT × 1
    diag_ee_hat   = u_hat_col   :^ 2        // NT × 1

    // ──────────────────────────────────────
    // 步骤 3: 分支路由
    // ──────────────────────────────────────

    if (hetero_se == "on" & auto_se == "on") {
        // ── HAC 路径 (Arellano 聚类稳健, 默认) ──
        se_type = "Heteroscedastic / Autocorrelation Robust Standard Error"

        // 步骤 3a: 基本检查
        if (J_u < 2) {
            errprintf("_wfe_compute_se: J_u must be >= 2 for HAC standard errors\n")
            _error(498)
        }

        // 步骤 3b: Omega_raw 计算 + one-way HAC scaling
        //   Omega = (1/J_u) * Σᵢ Xᵢ'uᵢuᵢ'Xᵢ
        //   dfHAC = (J_u/(J_u-1)) * (N_nonzero/(N_nonzero-p))
        //   Psi   = dfHAC * (J_u * ginv_XX) * Omega * (J_u * ginv_XX)
        //   with df_adjustment(off): vcov = ginv_XX * Omega_raw * ginv_XX
        Omega_raw_wfe = _wfe_omega_hac(X_tilde, u_tilde_col, unit_col, J_u)
        Omega_raw_fe  = _wfe_omega_hac(X_hat,   u_hat_col,   unit_col, J_u)
        Omega_wfe = (1 / J_u) * Omega_raw_wfe
        Omega_fe  = (1 / J_u) * Omega_raw_fe
        dfHAC = 1
        if (df_adj_on) {
            if (N_nonzero - p <= 0) {
                errprintf("_wfe_compute_se: one-way HAC df_adjustment(on) requires N_nonzero - p > 0\n")
                _error(498)
            }
            dfHAC = (J_u / (J_u - 1)) * (N_nonzero / (N_nonzero - p))
        }

        // 步骤 3c: Psi (三明治矩阵)
        //   Psi = (J_u * ginv_XX) * Omega * (J_u * ginv_XX)
        Bread_wfe   = J_u * ginv_XX_tilde
        Bread_fe    = J_u * ginv_XX_hat
        Psi_hat_wfe = dfHAC * (Bread_wfe * Omega_wfe * Bread_wfe)
        Psi_hat_fe  = dfHAC * (Bread_fe  * Omega_fe  * Bread_fe)

    } else if (hetero_se == "on" & auto_se == "off") {
        // ── HC 路径 (异方差稳健, 无自相关) ──
        se_type = "Heteroscedastic Robust Standard Error"

        // 步骤 3a: Omega 计算
        //   已含 1/J_u 缩放
        Omega_wfe = _wfe_omega_hc(X_tilde, diag_ee_tilde, J_u)
        Omega_fe  = _wfe_omega_hc(X_hat,   diag_ee_hat,   J_u)

        // 步骤 3b: Stock-Watson 偏差修正 (Stock and Watson 2008, Eq. 6).
        // The plain HC sandwich is still well defined when J_u == 1 because it
        // only needs X' diag(u^2) X; the "at least 2 units" boundary belongs
        // specifically to the Stock-Watson correction.
        if (unbiased_se == "on") {
            // Validate using the actual unit_idx/J_t structure so direct
            // helper callers cannot be rejected by a stale parser-side summary.
            _wfe_validate_sw_balance(unit_col, J_u)
            se_type = "Heteroskedastic Standard Error (Stock-Watson Bias-Corrected)"

            // 计算 Stock-Watson 偏差修正 Psi
            _wfe_stockwatson_bias(X_tilde, u_tilde_col, unit_col, J_u, J_t,
                                  Omega_wfe, ginv_XX_tilde, Psi_sw_wfe)
            _wfe_stockwatson_bias(X_hat,   u_hat_col,   unit_col, J_u, J_t,
                                  Omega_fe,  ginv_XX_hat,   Psi_sw_fe)
            Psi_hat_wfe = Psi_sw_wfe
            Psi_hat_fe  = Psi_sw_fe
        }
        else {
            // 步骤 3c: 标准 HC Psi
            //   Psi = (J_u * ginv_XX) * Omega * (J_u * ginv_XX)
            Bread_wfe   = J_u * ginv_XX_tilde
            Bread_fe    = J_u * ginv_XX_hat
            Psi_hat_wfe = Bread_wfe * Omega_wfe * Bread_wfe
            Psi_hat_fe  = Bread_fe  * Omega_fe  * Bread_fe

            // 步骤 3d: HC 自由度修正 (与 HAC 同一公式)
            //   dfHC = (J_u/(J_u-1)) * (N_nonzero/(N_nonzero-p))
            if (df_adj_on) {
                if (N_nonzero - p <= 0) {
                    errprintf("_wfe_compute_se: HC df_adjustment(on) requires N_nonzero - p > 0\n")
                    _error(498)
                }
                dfHAC = (J_u / (J_u - 1)) * (N_nonzero / (N_nonzero - p))
                Psi_hat_wfe = dfHAC * Psi_hat_wfe
                Psi_hat_fe  = dfHAC * Psi_hat_fe
            }
        }

    }

    // ──────────────────────────────────────
    // 步骤 4: 最终 vcov 缩放
    //   vcov = Psi / J_u (⚠️ 除数是 J_u 单位数，不是 NT)
    // ──────────────────────────────────────
    vcov_wfe = Psi_hat_wfe / J_u
    vcov_fe  = Psi_hat_fe  / J_u
}

end
