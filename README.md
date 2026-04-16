# wfe

**Weighted Fixed Effects Estimators for Causal Inference with Panel Data**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](https://www.gnu.org/licenses/agpl-3.0.html)
[![Version: 1.0.0](https://img.shields.io/badge/Version-1.0.0-green.svg)]()

![wfe](image/image.jpg)

## Overview

`wfe` implements the **Weighted Fixed Effects** (WFE) and **Propensity-score Weighted Fixed Effects** (PWFE) estimators introduced by Imai and Kim (2021, *Political Analysis*) for causal inference with panel data.

Standard two-way fixed effects (TWFE) regressions implicitly assign regression weights to each observation, and these weights can be negative — leading to comparisons between observations with the *same* treatment status rather than *opposite* status. Imai and Kim (2021) show that this "mismatch" problem can bias causal estimates. The WFE approach derives observation-specific weights that target a well-defined causal estimand (ATE or ATT) and restricts matched comparison sets to observations with the opposite treatment status.

**Features:**

- Weighted fixed effects estimation targeting ATE or ATT
- One-way unit FE, one-way time FE, first-difference, multi-period DiD, and matched DiD estimators
- Propensity-score weighted fixed effects with internal penalized logit or user-supplied scores
- Cluster-robust HAC, HC, and Stock-Watson bias-corrected standard errors
- Built-in White (1980) misspecification test comparing WFE vs standard FE
- Postestimation: predicted values, residuals, and weight diagnostics

## Key Concepts

### Regression Weights and the Matching Representation

Imai and Kim (2021) analyze the standard two-way fixed effects (TWFE) estimator through a matching representation:

- **Proposition 1**: The standard TWFE estimator is equivalent to a two-way matching estimator. The counterfactual outcome for each observation is estimated from a within-unit matched set, a within-time matched set, and an adjustment set. Crucially, these matched sets may include "mismatches"—comparisons with observations of the *same* treatment status—which can bias causal estimates.
- **Proposition 2**: A weighted 2FE estimator (the WFE) restricts the within-unit matched set 𝓜\*ᵢₜ and within-time matched set 𝓝\*ᵢₜ to observations with the *opposite* treatment status, eliminating mismatches from these two sets. Some mismatches may remain in the adjustment set 𝓐\*ᵢₜ, which are corrected via a deflation factor Kᵢₜ. For one-way fixed effects (unit or time), no adjustment set arises and the resulting ATE and ATT weights are non-negative by construction.
- **Theorem 1**: Under the parallel trend assumption, the multi-period difference-in-differences estimator equals a weighted 2FE estimator, but some observations receive *negative* regression weights—specifically, control observations that frequently serve as adjustments for multiple treated observations.

### White Misspecification Test

The built-in White (1980) test compares the WFE specification against the standard (unweighted) FE specification:

- **H0**: No misspecification — the WFE and standard FE produce the same coefficients.
- **Reject H0**: Evidence that the standard FE regression weights produce a different estimate than the causally motivated WFE weights, suggesting potential bias from negative or misaligned weights.

### Estimator Families

| Estimator                  | Command / Option                        | Description                                             |
| -------------------------- | --------------------------------------- | ------------------------------------------------------- |
| **One-way unit FE**  | `wfe` with `method(unit)` (default) | Unit fixed effects with WFE weights                     |
| **One-way time FE**  | `wfe` with `method(time)`           | Time fixed effects with WFE weights                     |
| **First-difference** | `wfe` with `estimator(fd)`          | First-difference design with WFE weights                |
| **Multi-period DiD** | `wfe` with `estimator(did)`         | Two-way FE DiD with observation-specific weights        |
| **Matched DiD**      | `wfe` with `estimator(Mdid)`        | DiD with matching on pre-treatment outcomes             |
| **PWFE**             | `pwfe` (separate command)             | Propensity-score weighted FE via outcome transformation |

## Requirements

- Stata 16.0 or later (Mata complex type support)
- No additional dependencies

## Installation

### Install the package

```stata
net install wfe, from("https://raw.githubusercontent.com/gorgeousfish/wfe/main/")
```

This automatically installs:

- All commands and help files (`wfe`, `pwfe`, postestimation utilities)
- Compiled Mata library (`lwfe.mlib`) and Mata source files
- Example datasets (`dem.csv`, `castle.csv`, `acemoglu2019.csv`)

> **macOS note:** If you get `r(603)` during installation, fix PLUS directory permissions first:
> ```bash
> sudo chown -R $(whoami) /Applications/Stata/plus/
> ```

### Verify Installation

```stata
which wfe
which pwfe
help wfe
help pwfe
```

## Quick Start

### Loading Example Data

The bundled CSV datasets are installed with the package. Use `findfile` to locate and load them:

```stata
qui findfile castle.csv
import delimited using "`r(fn)'", clear
```

Or copy all datasets to your working directory once:

```stata
foreach ds in dem castle acemoglu2019 {
    qui findfile `ds'.csv
    copy "`r(fn)'" "`ds'.csv", replace
}
```

### Quick Example with Castle Doctrine Data

```stata
import delimited using "castle.csv", clear

* One-way unit FE, ATE
wfe l_homicide, treat(cdl_binary) unit(sid) time(year)

* Postestimation
predict double xbhat, xb
predict double resid, residuals
estat wfe_weights
```

Expected output (abbreviated):

```
Weighted Fixed Effects Estimation               Number of obs     =       550
  Method:      unit                             Number of units   =        50
  Quantity:    ate                              Time periods      =        11
                                                Non-zero wt       =       231
                                                Residual df       =       499
                                                Neg. weights      =         0
                                                Sigma             =  .1827925
------------------------------------------------------------------------------
             |  Heteroscedastic / Autocorrelation Robust Standard Error
  l_homicide | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
  cdl_binary |   .0242308   .0498756     0.49   0.627    -.0737613    .1222229
------------------------------------------------------------------------------
White (1980) Misspecification Test
  H0: No misspecification (WFE = Standard FE)
  Chi2(1)     =    0.1319
  P-value     =    0.7165
  -> Fail to reject H0 at alpha = 0.050
```

## Bundled Datasets

| Dataset              | Source                  | Units         | Periods    | Treatment                            | Outcome                            |
| -------------------- | ----------------------- | ------------- | ---------- | ------------------------------------ | ---------------------------------- |
| `dem.csv`          | Acemoglu et al. (2008)  | 184 countries | 1960–2010 | `dem` (democracy)                  | `y` (log GDP per capita)         |
| `castle.csv`       | Cheng & Hoekstra (2013) | 50 US states  | 2000–2010 | `cdl_binary` (castle doctrine law) | `l_homicide` (log homicide rate) |
| `acemoglu2019.csv` | Acemoglu et al. (2019)  | 175 countries | 1960–2010 | `dem` (democracy)                  | `ln_gdppc` (log GDP per capita)  |

**Notes on missing values:**

- `dem.csv` contains `NA` string values in columns `dem`, `y`, and `tradewb`. When importing, you must `destring` all three columns and drop missing cases before estimation:

```stata
import delimited using "dem.csv", clear
destring dem y tradewb, replace force
drop if missing(y, dem)
```

- `acemoglu2019.csv` is an unbalanced panel: not all countries have data for all years (e.g., AFG begins in 2000). It also contains `NA` string values in the `trade` column. If `trade` is used as a covariate, run `destring trade, replace force` after importing. The core variables (`dem`, `ln_gdppc`) have no missing values.
- `castle.csv` is a balanced panel with no missing values.

## Commands

| Command               | Description                                        |
| --------------------- | -------------------------------------------------- |
| `wfe`               | Weighted fixed effects estimation                  |
| `pwfe`              | Propensity-score weighted fixed effects estimation |
| `predict`           | Postestimation:`xb`, `fitted`, `residuals`   |
| `estat wfe_weights` | Weight matrix diagnostics                          |

## Options

### wfe Options

| Option                    | Description                                                                | Default                                      |
| ------------------------- | -------------------------------------------------------------------------- | -------------------------------------------- |
| **Required**        |                                                                            |                                              |
| `treat(varname)`        | Binary treatment indicator (0/1)                                           | Required                                     |
| `unit(varname)`         | Panel unit identifier                                                      | Required                                     |
| **Model**           |                                                                            |                                              |
| `time(varname)`         | Time identifier; required for `method(time)`, `estimator(fd/did/Mdid)` | Auto-generated for `method(unit)` baseline |
| `method(string)`        | `unit` or `time`                                                       | `unit`                                     |
| `qoi(string)`           | `ate` or `att`                                                         | `ate`                                      |
| `estimator(string)`     | `fd`, `did`, or `Mdid`                                               | Omit for baseline                            |
| `cit(varname)`          | Non-negative user-supplied C_it weights                                    | 1 for all                                    |
| `unweighted`            | Standard unweighted FE model                                               | —                                           |
| **Standard errors** |                                                                            |                                              |
| `hetero_se(on\|off)`     | Heteroskedasticity-robust SE                                               | `on`                                       |
| `auto_se(on\|off)`       | Autocorrelation-robust SE                                                  | `on`                                       |
| `df_adjustment(on\|off)` | Degrees-of-freedom correction                                              | `on`                                       |
| `unbiased_se`           | Stock-Watson bias-corrected SE (one-way only)                              | —                                           |
| **White test**      |                                                                            |                                              |
| `[no]white`             | White misspecification test                                                | `white`                                    |
| `white_alpha(#)`        | Significance level for White test                                          | `0.05`                                     |
| **DiD / Mdid**      |                                                                            |                                              |
| `maxdev_did(#)`         | Max deviation for matched DiD                                              | Nearest-neighbor                             |
| `tol(#)`                | Generalized inverse tolerance                                              | `sqrt(epsdouble)`                          |
| **Output**          |                                                                            |                                              |
| `[no]verbose`           | Progress messages                                                          | `verbose`                                  |
| `store_wdm`             | Store weighted demeaned data (one-way only)                                | —                                           |
| `diagnose`              | Print parsed state and exit                                                | —                                           |

### pwfe Options

| Option                     | Description                                                          | Default                                      |
| -------------------------- | -------------------------------------------------------------------- | -------------------------------------------- |
| **Required**         |                                                                      |                                              |
| `treat(varname)`         | Binary treatment indicator (0/1)                                     | Required                                     |
| `outcome(varname)`       | Outcome variable                                                     | Required                                     |
| `unit(varname)`          | Panel unit identifier                                                | Required                                     |
| **Model**            |                                                                      |                                              |
| `time(varname)`          | Time identifier; required for `method(time)` and `estimator(fd)` | Auto-generated for `method(unit)` baseline |
| `method(string)`         | `unit` or `time`                                                 | `unit`                                     |
| `qoi(string)`            | `ate` or `att`                                                   | `ate`                                      |
| `estimator(string)`      | `fd` or omit                                                       | Omit for baseline                            |
| `cit(varname)`           | Non-negative user-supplied C_it weights                              | 1 for all                                    |
| **Propensity score** |                                                                      |                                              |
| `pscore(varname)`        | Precomputed propensity score in (0,1)                                | —                                           |
| `nowithin_unit`          | Use pooled logit instead of split-by-unit/time                       | —                                           |
| **Standard errors**  |                                                                      |                                              |
| `hetero_se(on\|off)`      | Heteroskedasticity-robust SE                                         | `on`                                       |
| `auto_se(on\|off)`        | Autocorrelation-robust SE                                            | `on`                                       |
| `unbiased_se`            | Stock-Watson bias-corrected SE                                       | —                                           |
| **White test**       |                                                                      |                                              |
| `[no]white`              | White misspecification test                                          | `white`                                    |
| `white_alpha(#)`         | Significance level                                                   | `0.05`                                     |
| **Output**           |                                                                      |                                              |
| `[no]verbose`            | Progress messages                                                    | `verbose`                                  |
| `diagnose`               | Print parsed state and exit                                          | —                                           |

**Note on pwfe syntax:** The optional `varlist` in `pwfe` specifies propensity-score covariates, not outcome regressors. The outcome variable is supplied via `outcome()`. Factor-variable notation (e.g., `i.group`) is supported in the propensity-score formula.

## Examples

All examples below assume the bundled CSV files are in the current working directory. See [Loading Example Data](#loading-example-data) above to copy them.

### Example 1: One-Way Unit FE

```stata
import delimited using "castle.csv", clear

* Basic: unit FE, ATE (default)
wfe l_homicide, treat(cdl_binary) unit(sid) time(year)

* Variations
wfe l_homicide, treat(cdl_binary) unit(sid) time(year) qoi(att)          // ATT
wfe l_homicide, treat(cdl_binary) unit(sid) time(year) method(time)      // time FE
wfe l_homicide l_assault l_robbery, treat(cdl_binary) unit(sid) time(year) // covariates
```

Output (basic ATE):

```
Weighted Fixed Effects Estimation               Number of obs     =       550
  Method:      unit                             Number of units   =        50
  Quantity:    ate                              Time periods      =        11
                                                Non-zero wt       =       231
                                                Residual df       =       499
                                                Neg. weights      =         0
                                                Sigma             =  .1827925
------------------------------------------------------------------------------
             |  Heteroscedastic / Autocorrelation Robust Standard Error
  l_homicide | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
  cdl_binary |   .0242308   .0498756     0.49   0.627    -.0737613    .1222229
------------------------------------------------------------------------------
White (1980) Misspecification Test
  H0: No misspecification (WFE = Standard FE)
  Chi2(1)     =    0.1319
  P-value     =    0.7165
  -> Fail to reject H0 at alpha = 0.050
```

### Example 2: First-Difference and Two-Way DiD

```stata
* First-difference (Castle data)
import delimited using "castle.csv", clear
wfe l_homicide, treat(cdl_binary) unit(sid) time(year) estimator(fd)

* Two-way DiD (Democracy data)
import delimited using "dem.csv", clear
destring dem y tradewb, replace force
drop if missing(y, dem)

wfe y, treat(dem) unit(wbcode2) time(year) estimator(did)

* Matched DiD with nearest-neighbor matching
wfe y, treat(dem) unit(wbcode2) time(year) estimator(Mdid)

* Standard unweighted TWFE for comparison
wfe y, treat(dem) unit(wbcode2) time(year) estimator(did) unweighted
```

Output (DiD):

```
Weighted Fixed Effects Estimation               Number of obs     =      6934
  Method:      Weighted Two-way                 Number of units   =       175
  Quantity:    ate                              Time periods      =        51
  Estimator:   DiD                              Non-zero wt       =      6200
                                                Residual df       =      6200
                                                Neg. weights      =      2504
                                                Sigma             =  142.2282
------------------------------------------------------------------------------
             |  Heteroscedastic / Autocorrelation Robust Standard Error
           y | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
         dem |   .7661579     1.4875     0.52   0.607    -2.149857    3.682173
------------------------------------------------------------------------------
White (1980) Misspecification Test
  H0: No misspecification (WFE = Standard FE)
  Chi2(1)     =    0.0319
  P-value     =    0.8583
  -> Fail to reject H0 at alpha = 0.050

  Warning: 2504 observations have negative weights.
  Negative weights are admissible on the weighted two-way DiD path.
```

Note: The standard unweighted TWFE produces `dem = -10.11 (p=0.020)`, a significant *negative* effect — in sharp contrast to the WFE-DiD estimate of `dem = 0.77 (p=0.607)`. This illustrates how misaligned regression weights in standard TWFE can bias causal estimates.

### Example 3: SE Options

```stata
import delimited using "castle.csv", clear

* Skip White test for faster estimation
wfe l_homicide, treat(cdl_binary) unit(sid) time(year) nowhite

* HC standard errors only (no autocorrelation correction)
wfe l_homicide, treat(cdl_binary) unit(sid) time(year) ///
    hetero_se(on) auto_se(off)
```

### Example 4: Propensity-Score Weighted FE (pwfe)

```stata
import delimited using "castle.csv", clear

* PWFE with internal propensity-score estimation
* varlist = propensity-score covariates; outcome() = dependent variable
pwfe l_income l_prisoner, ///
    treat(cdl_binary) outcome(l_homicide) unit(sid) time(year)

* User-supplied propensity score
set seed 12345
gen double ps = runiform(0.05, 0.95)
pwfe, treat(cdl_binary) outcome(l_homicide) unit(sid) time(year) pscore(ps)
```

Output (within-unit):

```
Propensity-Score Weighted FE Estimation         Number of obs     =       550
  Method:      unit                             Number of units   =        50
  Quantity:    ate                              Time periods      =        11
  P-score:     estimated                        Non-zero wt       =       231
  Transform:   unit                             Residual df       =       499
                                                Sigma             =   .485687
------------------------------------------------------------------------------
             |  Heteroscedastic / Autocorrelation Robust Standard Error
  l_homicide | Coefficient  std. err.      t    P>|t|     [95% conf. interval]
-------------+----------------------------------------------------------------
  cdl_binary |   .0224899   .0655654     0.34   0.732    -.1063283    .1513081
------------------------------------------------------------------------------
```

### Example 5: Postestimation

```stata
import delimited using "castle.csv", clear
wfe l_homicide, treat(cdl_binary) unit(sid) time(year)

* Predicted values and residuals
predict double xbhat, xb
predict double resid, residuals

* Weight matrix diagnostics
estat wfe_weights
```

Output (`estat wfe_weights`):

```
Weight summary from e(W)
------------------------------------------------------------
Rows (T)                                         11
Columns (N)                                      50
Total elements                                  550
Nonzero weights                                 231
Positive weights                                231
Negative weights                                  0
------------------------------------------------------------
All weights:   min       0.000000      max      5.500000
All weights:  mean       0.840000      sd       1.052413
------------------------------------------------------------
Nonzero:      min        1.222222      max      5.500000
Nonzero:     mean        2.000000      sd       0.560036
------------------------------------------------------------
Positive:     min        1.222222      max      5.500000
Positive:    mean        2.000000      sd       0.560036
------------------------------------------------------------
Negative weight ratio                      0.000000
```

## Stored Results

### wfe Stored Results

**Scalars:**

| Result              | Description                                            |
| ------------------- | ------------------------------------------------------ |
| `e(N)`            | Number of observations                                 |
| `e(N_units)`      | Number of units                                        |
| `e(N_times)`      | Number of time periods                                 |
| `e(df_r)`         | Residual degrees of freedom                            |
| `e(sigma)`        | Root mean squared error                                |
| `e(sigma2)`       | Error variance                                         |
| `e(N_nonzero)`    | Number of nonzero regression weights                   |
| `e(N_negative)`   | Number of negative regression weights                  |
| `e(white_stat)`   | White test statistic (when enabled)                    |
| `e(white_pvalue)` | White test p-value (when enabled)                      |
| `e(white_alpha)`  | White test significance level (when enabled)           |
| `e(maxdev_did)`   | Matched DiD deviation bound (`estimator(Mdid)` only) |

**Matrices:**

| Result       | Description                                          |
| ------------ | ---------------------------------------------------- |
| `e(b)`     | Coefficient row vector                               |
| `e(V)`     | Variance-covariance matrix                           |
| `e(W)`     | Weight matrix (T × N)                               |
| `e(b_fe)`  | Standard FE coefficients (when White enabled)        |
| `e(V_fe)`  | Standard FE variance (when White enabled)            |
| `e(Y_wdm)` | Weighted demeaned response (with `store_wdm`)      |
| `e(X_wdm)` | Weighted demeaned design matrix (with `store_wdm`) |

**Macros:**

| Result                | Description                                 |
| --------------------- | ------------------------------------------- |
| `e(cmd)`            | `wfe`                                     |
| `e(depvar)`         | Dependent variable name                     |
| `e(method)`         | `unit`, `time`, or `Weighted Two-way` |
| `e(qoi)`            | `ate`, `att`, or `unweighted`         |
| `e(estimator)`      | `NULL`, `fd`, `did`, or `Mdid`      |
| `e(vcetype)`        | Variance estimator description              |
| `e(white_test)`     | `TRUE` or `FALSE` (when White enabled)  |
| `e(qoi_desc)`       | Human-readable QOI description              |
| `e(estimator_desc)` | Human-readable estimator description        |
| `e(predict)`        | `wfe_predict`                             |
| `e(estat_cmd)`      | `wfe_estat`                               |

**Functions:**

| Result        | Description                |
| ------------- | -------------------------- |
| `e(sample)` | Marks estimation sample    |

### pwfe Stored Results

**Scalars:**

| Result              | Description                                  |
| ------------------- | -------------------------------------------- |
| `e(N)`            | Number of observations                       |
| `e(N_units)`      | Number of units                              |
| `e(N_times)`      | Number of time periods                       |
| `e(df_r)`         | Residual degrees of freedom                  |
| `e(sigma)`        | Root mean squared error                      |
| `e(sigma2)`       | Error variance                               |
| `e(N_nonzero)`    | Number of nonzero regression weights         |
| `e(white_stat)`   | White test statistic (when enabled)          |
| `e(white_pvalue)` | White test p-value (when enabled)            |
| `e(white_alpha)`  | White test significance level (when enabled) |

**Matrices:**

| Result        | Description                 |
| ------------- | --------------------------- |
| `e(b)`      | PWFE coefficient row vector |
| `e(V)`      | PWFE covariance matrix      |
| `e(W)`      | Weight matrix (T × N)      |
| `e(b_fe)`   | FE benchmark coefficients   |
| `e(V_fe)`   | FE benchmark variance       |
| `e(y_star)` | Transformed outcome vector  |
| `e(pscore)` | Propensity scores used      |

**Macros:**

| Result                 | Description                                |
| ---------------------- | ------------------------------------------ |
| `e(cmd)`             | `pwfe`                                   |
| `e(depvar)`          | Outcome variable name                      |
| `e(method)`          | `unit` or `time`                       |
| `e(qoi)`             | `ate` or `att`                         |
| `e(estimator)`       | `NULL` or `fd`                         |
| `e(vcetype)`         | Variance estimator description             |
| `e(qoi_desc)`        | Human-readable QOI description             |
| `e(estimator_desc)`  | Human-readable estimator description       |
| `e(pscore_source)`   | `user` or `estimated`                  |
| `e(transform_scope)` | `global`, `unit`, or `time`          |
| `e(white_test)`      | `TRUE` or `FALSE` (when White enabled) |
| `e(predict)`         | `wfe_predict`                            |
| `e(estat_cmd)`       | `wfe_estat`                              |

**Functions:**

| Result        | Description                |
| ------------- | -------------------------- |
| `e(sample)` | Marks estimation sample    |

## Methodology

### Weighted Fixed Effects Estimation

Imai and Kim (2021) show that the standard TWFE estimator is equivalent to a matching estimator (Proposition 1). For each observation (i, t), the counterfactual outcome under the opposite treatment status is estimated from three components:

- The **within-unit matched set** 𝓜ᵢₜ: other time periods of the same unit
- The **within-time matched set** 𝓝ᵢₜ: other units at the same time period
- The **adjustment set** 𝓐ᵢₜ: observations correcting for double-counted unit and time effects

These sets may include "mismatches"—observations sharing the same treatment status as (i, t)—which can attenuate causal estimates. The WFE approach (Proposition 2) restricts 𝓜ᵢₜ and 𝓝ᵢₜ to observations with the opposite treatment status, and estimates β via weighted least squares:

```
β̂_WFE = argmin Σ_i Σ_t W_it · {(Y_it - Ȳ_i* - Ȳ_t* + Ȳ*) - β(X_it - X̄_i* - X̄_t* + X̄*)}²
```

where the asterisks denote weighted averages and the weights W_it are derived from the restricted matched sets. For **one-way** estimators (`method(unit)` or `method(time)`), the weights target the ATE or ATT and are non-negative. For **two-way DiD** estimators (`did`/`Mdid`), the weight construction follows Theorem 1 of the paper and some observations may receive negative weights.

### Standard Errors

For **one-way** estimators, the default reports cluster-robust HAC standard errors following Arellano (1987):

```
V̂ = dfHAC · (X'X)⁻ · Ω̂ · (X'X)⁻
```

where `Ω̂ = (1/J_u) Σ_i X_i' û_i û_i' X_i` is the unit-clustered meat (scaled by the number of units J_u) and `dfHAC = (J_u/(J_u-1)) · (N_nz/(N_nz-p))`.

For **two-way** estimators (`did`/`Mdid`), a GMM sandwich variance with unit clustering is used.

### White Misspecification Test

The White test compares the WFE coefficients against the standard FE coefficients:

```
H0: β_WFE = β_FE
Test statistic: W = δ' · Φ⁻¹ · δ  ~  χ²(p)
```

where `δ = β̂_WFE - β̂_FE` and `Φ = V̂_WFE + V̂_FE - Λ₁₂ - Λ₁₂'`. Rejection suggests that the implicit regression weights in the standard FE specification produce a meaningfully different estimate from the causally motivated WFE weights.

## References

**Methodology:**

Imai K, Kim IS. On the Use of Two-Way Fixed Effects Regression Models for Causal Inference with Panel Data. Political Analysis. 2021;29(3):405-415. doi:10.1017/pan.2020.33

Imai K, Kim IS. When Should We Use Unit Fixed Effects Regression Models for Causal Inference with Longitudinal Data? American Journal of Political Science. 2019;63(2):467-490. doi:10.1111/ajps.12417

**Standard errors:**

Arellano M. Computing Robust Standard Errors for Within-Groups Estimators. Oxford Bulletin of Economics and Statistics. 1987;49(4):431-434.

Stock JH, Watson MW. Heteroskedasticity-Robust Standard Errors for Fixed Effects Panel Data Regression. Econometrica. 2008;76(1):155-174.

White H. Using Least Squares to Approximate Unknown Regression Functions. International Economic Review. 1980;21(1):149-170.

**Bundled datasets:**

Acemoglu D, Johnson S, Robinson JA, Yared P. Income and Democracy. American Economic Review. 2008;98(3):808-842.

Acemoglu D, Naidu S, Restrepo P, Robinson JA. Democracy Does Cause Growth. Journal of Political Economy. 2019;127(1):47-100.

Cheng C, Hoekstra M. Does Strengthening Self-Defense Law Deter Crime or Escalate Violence? Evidence from Expansions to Castle Doctrine. Journal of Human Resources. 2013;48(3):821-854.

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Kosuke Imai**, Harvard University
- **In Song Kim**, MIT

## License

AGPL-3.0. See the [GNU AGPL-3.0 license](https://www.gnu.org/licenses/agpl-3.0.html) for details.

## Citation

If you use this package in your research, please cite both the methodology paper and the Stata implementation:

**APA Format:**

> Cai, X., & Xu, W. (2026). *wfe: Stata module for Weighted Fixed Effects estimation* (Version 1.0.0) [Computer software]. GitHub. https://github.com/gorgeousfish/wfe
>
> Imai, K., & Kim, I. S. (2021). On the Use of Two-Way Fixed Effects Regression Models for Causal Inference with Panel Data. *Political Analysis*, 29(3), 405–415. https://doi.org/10.1017/pan.2020.33

**BibTeX:**

```bibtex
@software{wfe2026stata,
  title={wfe: Stata module for Weighted Fixed Effects estimation},
  author={Xuanyu Cai and Wenli Xu},
  year={2026},
  version={1.0.0},
  url={https://github.com/gorgeousfish/wfe}
}

@article{imai2021twoway,
  title={On the Use of Two-Way Fixed Effects Regression Models for Causal Inference with Panel Data},
  author={Imai, Kosuke and Kim, In Song},
  journal={Political Analysis},
  volume={29},
  number={3},
  pages={405--415},
  year={2021},
  doi={10.1017/pan.2020.33}
}
```

## See Also

- R package by Imai and Kim: https://cran.r-project.org/package=wfe
- Paper: Imai, K., & Kim, I. S. (2021). https://doi.org/10.1017/pan.2020.33
