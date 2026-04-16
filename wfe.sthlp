{smcl}
{* *! wfe.sthlp}{...}
{vieweralsosee "pwfe" "help pwfe"}{...}
{vieweralsosee "xtreg" "help xtreg"}{...}
{vieweralsosee "xtdidregress" "help xtdidregress"}{...}
{vieweralsosee "predict" "help predict"}{...}
{vieweralsosee "estat" "help estat"}{...}
{viewerjumpto "Syntax" "wfe##syntax"}{...}
{viewerjumpto "Description" "wfe##description"}{...}
{viewerjumpto "Options" "wfe##options"}{...}
{viewerjumpto "Remarks" "wfe##remarks"}{...}
{viewerjumpto "Stored results" "wfe##results"}{...}
{viewerjumpto "Examples" "wfe##examples"}{...}
{viewerjumpto "References" "wfe##references"}{...}
{viewerjumpto "Also see" "wfe##alsosee"}{...}

{title:Title}

{p2colset 5 16 24 2}{...}
{p2col:{cmd:wfe} {hline 2}}Weighted Fixed Effects Estimator for Causal Inference with Panel Data{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:wfe} {it:depvar} [{it:indepvars}] {ifin} {cmd:,}
{opth treat(varname)} {opth unit(varname)} [{it:options}]
{p_end}

{pstd}
If {it:indepvars} is omitted, {cmd:wfe} uses the variable supplied in
{cmd:treat()} as the sole scored regressor.  This matches the
treatment-only WFE specification described in
{help wfe##IK2021:Imai and Kim (2021)}.
The current release requires Stata 16 or later and fails fast with
{cmd:wfe requires Stata 16 or later (Mata complex type support)} when
run under an earlier release.  Standard Stata
factor-variable notation in {it:indepvars}, such as {cmd:i.time_id} or
{cmd:c.x##i.group}, is expanded on the estimation sample before the Mata
bridge is called.
The dependent variable itself must be a raw numeric variable name; factor-variable notation is not allowed in {it:depvar} and is rejected with {cmd:dependent variable must be a numeric variable name; factor-variable notation is not allowed}.

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth treat(varname)}}binary treatment indicator coded as 0/1{p_end}
{synopt:{opth unit(varname)}}panel unit identifier{p_end}

{syntab:Model}
{synopt:{opth time(varname)}}time identifier; auto-generated if omitted except for {cmd:method(time)}, {cmd:estimator(fd)}, {cmd:estimator(did)}, and {cmd:estimator(Mdid)}{p_end}
{synopt:{opt method(string)}}{cmd:unit} (default) or {cmd:time}{p_end}
{synopt:{opt qoi(string)}}{cmd:ate} (default) or {cmd:att}{p_end}
{synopt:{opt estimator(string)}}{cmd:fd}, {cmd:did}, or {cmd:Mdid}{p_end}
{synopt:{opth cit(varname)}}nonnegative user-supplied {it:C_it} weights; fractional values are truncated toward zero{p_end}
{synopt:{opt unweighted}}fit the standard unweighted FE model{p_end}

{syntab:Standard errors}
{synopt:{opt hetero_se(on|off)}}heteroskedasticity-robust SE; default {cmd:on}{p_end}
{synopt:{opt auto_se(on|off)}}autocorrelation-robust SE; default {cmd:on}{p_end}
{synopt:{opt df_adjustment(on|off)}}degrees-of-freedom correction for one-way HAC and two-way GMM paths; default {cmd:on}{p_end}
{synopt:{opt unbiased_se}}post Stock-Watson bias-corrected covariance matrices; one-way only{p_end}

{syntab:White test}
{synopt:{opt [no]white}}run the White misspecification test; default {cmd:white}{p_end}
{synopt:{opt white_alpha(#)}}significance level in {cmd:(0,1)} for the White test; default {cmd:0.05}{p_end}

{syntab:DiD and Mdid}
{synopt:{opt maxdev_did(#)}}maximum deviation for matched DiD; {cmd:estimator(Mdid)} only{p_end}
{synopt:{opt tol(#)}}positive relative tolerance for generalized inverses{p_end}

{syntab:Output}
{synopt:{opt [no]verbose}}display progress messages; default {cmd:verbose}{p_end}
{synopt:{opt store_wdm}}store weighted demeaned data in {cmd:e()} for one-way estimators; not supported for {cmd:estimator(did)} / {cmd:estimator(Mdid)}{p_end}
{synopt:{opt diagnose}}prints the parsed execution state and exits without estimation{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:wfe} fits the weighted fixed effects (WFE) estimators introduced by
{help wfe##IK2021:Imai and Kim (2021)} for causal inference with panel
data.  It also derives the regression weights associated with the chosen
causal quantity of interest.

{pstd}
The key idea is to restrict matched comparison sets to observations with
the opposite treatment status.  In the two-way setting, this reduces the
"mismatch" problem emphasized in Proposition 2 of the paper: standard
two-way fixed effects comparisons may partially compare observations that
share the same treatment status, whereas the weighted estimator shrinks
that mismatch component through observation-specific weights.

{pstd}
The current Stata command exposes the following estimator families:
one-way WFE with unit fixed effects ({cmd:method(unit)}), one-way WFE with
time fixed effects ({cmd:method(time)}), first difference
({cmd:estimator(fd)}), multi-period DiD ({cmd:estimator(did)}), and
matched DiD ({cmd:estimator(Mdid)}).

{pstd}
{cmd:wfe} supports both the average treatment effect ({cmd:qoi(ate)}) and
the average treatment effect for the treated ({cmd:qoi(att)}).  A custom
nonnegative {it:C_it} variable may be supplied through {cmd:cit()} to
target a different population.  The White misspecification test is run by
default unless {cmd:nowhite} is specified.

{pstd}
The model specification should not include dummy variables for fixed
effects; fixed effects are handled internally.  If you prefer a
propensity-score-weighted alternative, see {help pwfe}.  The current
implementation requires Stata 16 or later because the Mata layer uses
complex-number support.

{pstd}
{cmd:wfe} sorts internally by unit-time for preprocessing and Mata
estimation, but it restores the caller's original observation order before
returning.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth treat(varname)} specifies the binary treatment indicator.  The
current implementation requires integer values 0 and 1 only; missing,
noninteger, and nonbinary values are not allowed.  If the estimation
sample contains more than two observed treatment values, {cmd:wfe} fails
fast with {cmd:'treat' must be a binary vector: there are more than two values of treatment}.
If exactly two observed treatment values remain but they are not {cmd:0}
and {cmd:1}, {cmd:wfe} fails fast with
{cmd:'treat' must be a either zero or one where one indicates treatment}.
If the estimation sample contains only one observed treatment level, {cmd:wfe} fails fast
with {cmd:'treat' must contain both 0 and 1 values}.
Nonnumeric {cmd:treat()} input is rejected with
{cmd:treat() must specify a numeric variable}.  Missing values are
rejected with {cmd:missing values in treat() are not allowed} before any
preprocessing or matched-set construction is attempted.

{phang}
{opth unit(varname)} specifies the panel unit identifier.  Because the
unit index defines the unit-time dyads of the estimator rather than an
ordinary regressand-side covariate, missing values in {cmd:unit()} are
not silently dropped; the command fails fast with
{cmd:missing values in unit() are not allowed}.  Both
{cmd:unit()} and {cmd:time()} may be numeric or string identifiers.

{phang}
{it:indepvars} may use standard Stata factor-variable notation.  The
current implementation expands factor terms to the same no-intercept
dummy-and-interaction columns that the equivalent hand-written regressor
list would produce on the estimation sample.  The treatment variable
supplied in {cmd:treat()} is always included in the scored design.  If
{it:indepvars} is omitted, it is the sole scored regressor.  If
additional {it:indepvars} are supplied, {cmd:wfe} prepends the
{cmd:treat()} main effect unless the right-hand side already spells that
same treatment main effect explicitly.  This keeps the interface
aligned with the paper's specification, where the treatment indicator
always appears on the right-hand side of the regression formula.
If {it:depvar}, {it:indepvars}, and any {cmd:if}/{cmd:in} restriction
leave no estimation-sample observations after the model-frame-style
missing-value screen, {cmd:wfe} stops with the standard Stata error
{cmd:no observations}.
The dependent variable itself must be a raw numeric variable name, so factor-variable notation is not allowed in {it:depvar}; any attempt to spell the response as a factor token is rejected with {cmd:dependent variable must be a numeric variable name; factor-variable notation is not allowed}.

{dlgtab:Model}

{phang}
{opth time(varname)} specifies the time identifier.  If omitted, the
command internally generates a consecutive time index when the selected
estimator allows it.  That internally generated clock follows the current
within-unit observation order before estimation.  The first-difference estimator
does not allow an omitted time variable and returns
{cmd:First Difference cannot calculate when 'time.index' is missing}.
The two-way DiD and matched-DiD estimators also require an explicit
{cmd:time()} variable and fail fast with
{cmd:'time.index' should be provided}.
When {cmd:time()} is supplied, missing values are rejected with
{cmd:missing values in time() are not allowed}.  When {cmd:time()}
is supplied, duplicate {cmd:(unit,time)} observations are rejected with
{cmd:duplicate observations found for unit-time combination}.  Duplicate
{cmd:(unit,time)} observations are checked only on the formula-valid
sample that survives missing-value removal.  Duplicate dyads that appear only in rows
dropped because {it:depvar} or {it:indepvars} are missing do not trigger
that parser error.  Current {cmd:wfe} factorizes either numeric or string panel identifiers to dense internal indices.

{phang}
{opt method(string)} chooses the one-way fixed effect dimension.  Allowed
values are {cmd:unit} (default) and {cmd:time}.  Any other value is
rejected with {cmd:method should be either unit or time}.  The current
command requires {cmd:time()} when {cmd:method(time)} is used; otherwise
it returns {cmd:'time.index' should be provided}.  The following
combinations are rejected with their exact messages:

{phang2}
{cmd:method(time)} + {cmd:estimator(fd)}:
{cmd:First Difference is not compatible with 'time' method: set method == 'unit'}

{phang2}
{cmd:method(time)} + {cmd:estimator(did)}:
{cmd:Difference-in-Differences is not compatible with 'time' method: set method == 'unit'}

{phang2}
{cmd:method(time)} + {cmd:estimator(Mdid)}:
{cmd:Match-Difference-in-Differences is not compatible with 'time' method: set method == 'unit'}

{phang2}
For matched-set feasibility, {cmd:method(time)} still requires at least
2 units in the estimation sample.  If only one unit remains after the
formula-valid missing-value screen, {cmd:wfe} stops during preprocessing
with {cmd:panel has only one unit; estimation not possible}.

{phang}
For two-way weighted fixed effects regressions, keep the default
{cmd:method(unit)}.

{phang}
{opt qoi(string)} selects the causal quantity of interest: {cmd:ate}
(default) for the average treatment effect or {cmd:att} for the average
treatment effect for the treated.  Under {cmd:att} with
{cmd:estimator(fd)}, {cmd:estimator(did)}, or {cmd:estimator(Mdid)}, the
target weight is multiplied by the current-period treatment indicator.
Only transitions ending the current period treated contribute to the ATT
target, so switch-off observations receive zero ATT weight and treated
observations are compared only against control observations from the
previous period, not the next period.  Any other value is rejected with
{cmd:qoi() must be 'ate' or 'att'}.

{phang}
{opt estimator(string)} requests an alternative estimator.  Omit the
option for the baseline one-way WFE estimator; use {cmd:fd} for the
first-difference design, {cmd:did} for the multi-period DiD estimator,
and {cmd:Mdid} for matched DiD.  The current Stata implementation is
case-sensitive and accepts only {cmd:fd}, {cmd:did}, or {cmd:Mdid}; any
other value returns {cmd:estimator() must be 'fd', 'did', or 'Mdid'}.
For {cmd:did} and {cmd:Mdid}, two-way fixed effect standard errors are
computed with a GMM sandwich estimator and only
{cmd:hetero_se(on) auto_se(on)} is supported.
The two-way DiD and matched-DiD paths also require at least 2 units and at least 2 time periods.
Theorem 1 also requires at least one treated switch with a valid
opposite-treatment comparison so that the public DiD/MDiD weight vector
contains at least one non-zero entry.  If no such matched contrast
exists, {cmd:wfe} fails fast with
{cmd:wfe_twoway: no non-zero DiD/MDiD weights; weighted two-way estimator is unidentified}.
For matched-set feasibility, the one-way {cmd:method(unit)} path and
{cmd:estimator(fd)} still require at least 2 time periods in the
estimation sample.  If only one period remains after preprocessing,
{cmd:wfe} stops with {cmd:panel has only one time period; estimation not possible}.
For one-way {cmd:method(unit)} paths, residual variance and covariance
also require {cmd:NT - p - J_u > 0} after preprocessing.  If this
denominator is nonpositive, {cmd:wfe} stops before entering the one-way
Mata backend with {cmd:wfe: insufficient degrees of freedom: NT=..., p=..., J_u=..., df_r=...}.

{phang}
{opth cit(varname)} supplies a user-defined nonnegative {it:C_it} weight
variable.  If omitted, the command sets {it:C_it} to one for all
observations.  Fractional {it:C_it} values are truncated toward zero before weights are constructed.  Nonnumeric input is rejected with {cmd:'C.it' must be a numeric vector with length equal to number of observations}.  Negative values are rejected with
{cmd:'C.it' must be a non-negative numeric vector}, and missing values are
rejected with {cmd:missing values in cit() are not allowed}.  These {cmd:cit()} guards are evaluated on the estimation
sample after missing-value removal for {it:depvar} and
{it:indepvars}, so rows dropped only because {it:y}/{it:X} are missing do
not trigger {cmd:cit()} errors.

{phang}
{opt unweighted} requests the standard unweighted fixed effect model,
setting all regression weights to one.  Under this option, the user does
not need to specify {cmd:qoi()}, and the White statistic should be close
to zero when the weighted and unweighted specifications coincide.  For
the standard two-way fixed effects regression, combine {cmd:unweighted}
with {cmd:estimator(did)} or {cmd:estimator(Mdid)}.  Both routes post
the same unweighted two-way FE coefficient and covariance objects; under
{cmd:estimator(Mdid)}, {cmd:maxdev_did()} only changes the stored
metadata {cmd:e(maxdev_did)}.

{phang}
{opt tol(#)} sets the relative cutoff used by the generalized inverse in
the two-way complex-projection path.  The option must be strictly
positive; {cmd:tol() <= 0} is rejected at parse time with
{cmd:tol() must be positive}.  If omitted, {cmd:wfe} uses
{cmd:sqrt(c(epsdouble))}.
Values outside {cmd:(0,+inf)} or nonnumeric input fail fast with
{cmd:tol() must be positive}.

{dlgtab:Standard errors}

{phang}
{opt hetero_se(on|off)} controls whether heteroskedasticity across units
is allowed in the reported standard errors.  The default is {cmd:on}.
Allowed values are {cmd:on} and {cmd:off}; any other value returns
{cmd:hetero_se() must be 'on' or 'off'}.

{phang}
{opt auto_se(on|off)} controls whether arbitrary serial correlation is
allowed in the reported standard errors.  The default is {cmd:on}.
Allowed values are {cmd:on} and {cmd:off}; any other value returns
{cmd:auto_se() must be 'on' or 'off'}.

{phang}
For one-way estimators, the current implementation rejects the following
exact combinations:

{phang2}
{cmd:hetero_se(off) auto_se(off)}:
{cmd:standard errors with independence and homoskedasticity is not supported}

{phang2}
{cmd:hetero_se(off) auto_se(on)}:
{cmd:robust standard errors with autocorrelation and homoskedasticity is not supported}

{phang}
For two-way estimators ({cmd:did} and {cmd:Mdid}), the following exact
messages apply:

{phang2}
{cmd:hetero_se(on) auto_se(off)} or
{cmd:hetero_se(off) auto_se(off)}:
{cmd:two-way FE requires hetero_se(on) and auto_se(on)}

{phang2}
{cmd:hetero_se(off) auto_se(on)}:
{cmd:Robust standard errors with autocorrelation and homoskedasticity is not supported}

{phang}
{opt df_adjustment(on|off)} toggles the small-sample degrees-of-freedom
correction used by the supported HAC/GMM variance estimators.  Under the
one-way default {cmd:hetero_se(on) auto_se(on)} path it rescales the HAC
Omega matrix by the finite-sample factor reported in the Remarks section.
Under {cmd:did} and {cmd:Mdid} it likewise toggles the two-way GMM
degrees-of-freedom correction.  The default is {cmd:on}.  Allowed values
are {cmd:on} and {cmd:off}; any other value returns
{cmd:df_adjustment() must be 'on' or 'off'}.

{phang}
{opt unbiased_se} requests the Stock-Watson (2008) bias-corrected covariance matrices
for the admissible one-way heteroskedasticity-robust path.  The correction
is intended only for balanced panels under {cmd:hetero_se(on) auto_se(off)}
and is only mathematically defined when the panel has at least 3 time periods
because the Stock-Watson factor uses {cmd:(J_t-1)/(J_t-2)}.  The current implementation posts the
Stock-Watson bias-corrected covariance matrices on the admissible one-way
path.  The current Stata parser rejects unsupported HAC and two-way uses with
{cmd:unbiased_se is allowed only for one-way hetero_se(on) auto_se(off)}.

{dlgtab:White test}

{phang}
{opt [no]white} controls whether the White misspecification test is run.
The default is {cmd:white}.  In Stata syntax, the disabled state is
written as {cmd:nowhite}.

{phang}
{opt white_alpha(#)} sets the significance level used when interpreting
the White test.  The default is {cmd:0.05}.  As a probability threshold,
Values outside this interval or nonnumeric input fail fast with
{cmd:white_alpha() must lie in (0,1)}.

{phang}
The current Stata White path requires the covariance-difference matrix
{cmd:Phi} to be positive semidefinite and full rank for chi-square
inference.  If {cmd:Phi} is not positive semidefinite or not full rank,
{cmd:wfe} aborts the White path without posting results.  If the
White statistic is numerically negative despite {cmd:Phi} passing those
checks, {cmd:wfe} issues a warning and posts {cmd:e(white_stat)} with
{cmd:e(white_pvalue) = 1} rather than aborting the entire command.

{dlgtab:DiD and Mdid}

{phang}
{opt maxdev_did(#)} sets the maximum allowed difference in pre-treatment
outcomes for matched DiD.  Outside {cmd:estimator(Mdid)}, the current
parser rejects the option with
{cmd:maxdev_did() is allowed only with estimator(Mdid)}.  When the
option is omitted, the implementation falls back to nearest-neighbor
matching.  A value of {cmd:0} requests exact matching.  Negative values are rejected with
{cmd:maxdev_did() must be non-negative}.  Nonnumeric
input is rejected with {cmd:maxdev_did() must be a numeric value}.  Only when
{cmd:estimator(Mdid)} is active does omitting {cmd:maxdev_did()} post {cmd:e(maxdev_did)} = {cmd:-1}
to denote nearest-neighbor matching; on every non-{cmd:Mdid} path, that
saved result is not posted.

{phang}
{opt tol(#)} sets the relative tolerance used when generalized inverses
detect zero singular values.  If omitted, the command uses the Stata
equivalent of {cmd:sqrt(epsdouble())}, approximately
{cmd:1.4901161193847656e-8}.

{dlgtab:Output}

{phang}
{opt [no]verbose} controls whether progress messages are printed.  The
default is {cmd:verbose}.  In Stata syntax, the disabled state is written
as {cmd:noverbose}.

{phang}
{opt store_wdm} stores the weighted demeaned response and design matrix in
{cmd:e(Y_wdm)} and {cmd:e(X_wdm)} for the one-way estimators.  It is not
supported for {cmd:estimator(did)} or {cmd:estimator(Mdid)} because the
two-way weighted-demeaned objects are complex-valued and therefore do
not fit the current public {cmd:e()} matrix contract.  The two-way path
therefore fails fast with
{cmd:store_wdm is not supported for two-way estimators (did/Mdid)}.
The default is not to store them.

{phang}
{opt diagnose} prints the parsed execution state and exits without
estimation.  The current implementation reports the parser-selected
values for {cmd:method}, {cmd:qoi}, {cmd:estimator}, {cmd:hetero_se},
{cmd:auto_se}, diagnostic output labels {cmd:df_adj} (reflecting the
state of {cmd:df_adjustment()}) and {cmd:unbiased} (reflecting the
state of {cmd:unbiased_se}), {cmd:white}, {cmd:white_alpha}, {cmd:tol},
{cmd:unweighted}, {cmd:verbose}, {cmd:store_wdm},
{cmd:maxdev_did}, the derived causal label, and
{cmd:touse N}.  Because the command exits after that diagnostics block,
{cmd:diagnose} does not post estimation results in {cmd:e()} and clears
any active estimation results, so stale {cmd:predict} or {cmd:estat}
output from a previous model cannot be replayed after the diagnostic run.

{marker remarks}{...}
{title:Remarks}

{marker compatibility}{...}
{bf:Compatibility restrictions}

{pstd}
The current Stata implementation follows the public parser contract in
{cmd:wfe.ado}.  The exact method and estimator incompatibilities are:

{phang2}
{cmd:method(time)} with {cmd:estimator(fd)}:
{cmd:First Difference is not compatible with 'time' method: set method == 'unit'}

{phang2}
{cmd:method(time)} with {cmd:estimator(did)}:
{cmd:Difference-in-Differences is not compatible with 'time' method: set method == 'unit'}

{phang2}
{cmd:method(time)} with {cmd:estimator(Mdid)}:
{cmd:Match-Difference-in-Differences is not compatible with 'time' method: set method == 'unit'}

{phang2}
omitting {cmd:time()} with {cmd:estimator(fd)}:
{cmd:First Difference cannot calculate when 'time.index' is missing}

{phang2}
omitting {cmd:time()} with {cmd:estimator(did)} or {cmd:estimator(Mdid)}:
{cmd:'time.index' should be provided}

{pstd}
Standard error restrictions are equally binding.  One-way paths reject
{cmd:hetero_se(off) auto_se(off)} and
{cmd:hetero_se(off) auto_se(on)}.  Two-way paths allow only
{cmd:hetero_se(on) auto_se(on)}.  The exact messages are reported in
{help wfe##options:Options}.

{marker semethods}{...}
{bf:Standard error methods}

{pstd}
For one-way estimators ({cmd:method(unit)}, {cmd:method(time)}, and
{cmd:estimator(fd)}), the command uses cluster-robust fixed effects
variance formulas following {help wfe##ARELLANO1987:Arellano (1987)}.  The
default setting {cmd:hetero_se(on) auto_se(on)} corresponds to the HAC
variance estimator
{cmd:vcov = dfHAC * ginv(X'X) * Omega_raw * ginv(X'X)}, where
{cmd:Omega_raw = sum_i X_i' u_i u_i' X_i} is the unscaled unit-cluster
meat and
{cmd:dfHAC = (J_u/(J_u-1)) * (N_nonzero/(N_nonzero-p))}.  Here
{it:J_u} is the number of units, {it:N_nonzero} is the number of
observed nonzero regression weights, and {it:p} is the number of
scored regressors.  When {cmd:df_adjustment(off)} is specified,
the command removes only that finite-sample multiplier and keeps the
unadjusted HAC covariance {cmd:ginv(X'X) * Omega_raw * ginv(X'X)}.
The adjusted path is defined only when {cmd:N_nonzero - p > 0}; if
that denominator is nonpositive, the helper fails fast with
{cmd:_wfe_compute_se: one-way HAC df_adjustment(on) requires N_nonzero - p > 0}.

{pstd}
With {cmd:hetero_se(on) auto_se(off)}, the one-way estimator switches to a
heteroskedasticity-robust HC variance.  The option {cmd:unbiased_se}
tracks the Stock-Watson idea for balanced panels with at least 3 time
periods.  On that admissible path, the current Stata implementation posts
the Stock-Watson bias-corrected covariance matrices, so the reported
variance can differ from the plain HC result.  The
{cmd:df_adjustment()} toggle does not rescale the Stock-Watson path;
however, on the plain HC path (without {cmd:unbiased_se}),
{cmd:df_adjustment(on)} applies the same
{cmd:dfHAC = (J_u/(J_u-1)) * (N_nonzero/(N_nonzero-p))} correction as the
HAC path.

{pstd}
For two-way estimators ({cmd:estimator(did)} and {cmd:estimator(Mdid)}),
the command uses a GMM sandwich variance with unit clustering.  When
{cmd:df_adjustment(on)} is active, the multiplier is
{cmd:(N*/(N*-1)) * ((N*-p)/(N*-n_u*-n_t*-p))}, where {it:N*} is the number
of nonzero-weight observations, {it:n_u*} is the number of units with at
least one nonzero weight, {it:n_t*} is the number of time periods with at
least one nonzero weight, and {it:p} is the number of scored regressors.
That correction is only defined when the denominator
{cmd:N*-n_u*-n_t*-p} is positive.  If the denominator is
nonpositive, the current implementation fails fast and advises
{cmd:df_adjustment(off)}.
If the weighted GMM moment matrix is singular after weighting and
two-way demeaning, the current two-way path fails fast with
{cmd:wfe_twoway: weighted GMM moment matrix is singular; two-way covariance is undefined}.
When {cmd:white} is requested, the FE-side HAC benchmark used for the
White comparison also requires {cmd:NT - J_u - N_times - p > 0}; if
that denominator is nonpositive, the current implementation fails fast
with {cmd:wfe_se_fe_twoway: FE-side HAC requires NT - J_u - N_times - p > 0}.
The two-way White covariance-difference cross term also requires
{cmd:Mstar - N_nonzero_units - N_nonzero_times - p > 0}; if that
denominator is nonpositive, the current implementation fails fast with
{cmd:wfe_white_test_twoway: df.white requires Mstar - N_nonzero_units - N_nonzero_times - p > 0}.

{marker weights}{...}
{bf:Weight interpretation}

{pstd}
The weight matrix {cmd:e(W)} is stored in {it:T x N} form, with rows
indexed by time and columns indexed by units.  Positive weights count an
observation as a control comparison for another observation; negative
weights arise naturally in the DiD representation discussed in
{help wfe##IK2021:Imai and Kim (2021)} and are summarized by
{cmd:estat wfe_weights}.

{pstd}
For sparse panels under {cmd:unweighted} in the one-way paths, {cmd:e(W)} is the observed-support mask on the declared time-by-unit grid, so {cmd:e(N_nonzero)} still counts only the observed nonzero regression weights.

{pstd}
For sparse panels under {cmd:unweighted} in the two-way
{cmd:did}/{cmd:Mdid} path, the public matrix {cmd:e(W)} is filled with ones over the declared
{it:T x N} grid even when some dyads are unobserved.  The reported
{cmd:e(N_nonzero)} still counts observed nonzero regression weights, so
it can be smaller than the number of positive cells in {cmd:e(W)} on
unbalanced panels.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:wfe} stores the following in {cmd:e()}:

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2:Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations in the estimation sample{p_end}
{synopt:{cmd:e(N_units)}}number of distinct units{p_end}
{synopt:{cmd:e(N_times)}}distinct time periods, or the omitted-time internal index upper bound{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(sigma)}}root mean squared error{p_end}
{synopt:{cmd:e(sigma2)}}error variance, when available{p_end}
{synopt:{cmd:e(N_nonzero)}}number of nonzero regression weights{p_end}
{synopt:{cmd:e(N_negative)}}number of negative regression weights{p_end}
{synopt:{cmd:e(maxdev_did)}}requested matched DiD deviation bound; {cmd:estimator(Mdid)} only, including {cmd:unweighted}{p_end}
{synopt:{cmd:e(white_stat)}}White test statistic; stored only when White runs{p_end}
{synopt:{cmd:e(white_pvalue)}}White test p-value; stored only when White runs{p_end}
{synopt:{cmd:e(white_alpha)}}White test alpha; stored only when White runs{p_end}

{p2col 5 22 26 2:Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient row vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix for {cmd:e(b)}{p_end}
{synopt:{cmd:e(W)}}weight matrix in {it:T x N} form{p_end}
{synopt:{cmd:e(b_fe)}}standard FE coefficients; stored only when White runs{p_end}
{synopt:{cmd:e(V_fe)}}standard FE variance; stored only when White runs{p_end}
{synopt:{cmd:e(Y_wdm)}}weighted demeaned response; stored only with {cmd:store_wdm} on one-way estimators{p_end}
{synopt:{cmd:e(X_wdm)}}weighted demeaned design matrix; stored only with {cmd:store_wdm} on one-way estimators{p_end}

{pstd}
When {cmd:store_wdm} is active on a one-way path, {cmd:e(Y_wdm)} and
{cmd:e(X_wdm)} are posted in the internal dense {cmd:(unit,time)}
estimation order after estimation-sample restriction, not in the caller's
current row order.

{pstd}
When {cmd:time()} is omitted, {cmd:e(N_times)} records the maximum
internally generated dense within-unit time index after estimation-sample
restriction, rather than a count taken from a user-supplied time
variable.

{pstd}
When {cmd:estimator(Mdid)} is active and {cmd:maxdev_did()} is omitted,
{cmd:e(maxdev_did)} is posted as {cmd:-1} to denote nearest-neighbor
matching.  On every non-{cmd:Mdid} path, {cmd:e(maxdev_did)} is not
posted.

{pstd}
For all estimation paths, {cmd:e(sigma)} and {cmd:e(sigma2)} come from the weighted WFE regression associated with {cmd:e(b)}.  In plain terms, e(sigma) and e(sigma2) come from the weighted WFE regression and do not describe the standard FE comparator used for White diagnostics.

{pstd}
For one-way paths, {cmd:e(df_r)} and {cmd:e(sigma2)} use the weighted-WFE denominator
{cmd:e(N) - colsof(e(b)) - G}.  Here {cmd:G = e(N_units)} for
{cmd:method(unit)}, {cmd:method(time)}, and {cmd:estimator(fd)}.

{pstd}
For two-way estimators, {cmd:e(df_r)} equals the number of rows retained by the weighted two-way WFE regression.  With {cmd:white} and {cmd:nowhite}, the posted {cmd:e(df_r)} remains {cmd:e(N_nonzero)}, because zero-weight rows may still be retained for White diagnostics but do not contribute to the weighted residual scale.

{pstd}
The shared token {cmd:Heteroscedastic / Autocorrelation Robust Standard Error} is posted both by the one-way HAC path and by the supported two-way GMM path.  To distinguish those cases programmatically, inspect {cmd:e(method)} and {cmd:e(estimator)} rather than relying on {cmd:e(vcetype)} alone.

{pstd}
No {cmd:wfe} path posts {cmd:Homoskedastic Standard Error}; {cmd:hetero_se(off) auto_se(off)} fails fast instead.

{p2col 5 22 26 2:Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:wfe}{p_end}
{synopt:{cmd:e(depvar)}}dependent variable name{p_end}
{synopt:{cmd:e(method)}}public method label: {cmd:unit}, {cmd:time}, or {cmd:Weighted Two-way} for {cmd:estimator(did)} / {cmd:estimator(Mdid)}{p_end}
{synopt:{cmd:e(qoi)}}quantity code; {cmd:ate} or {cmd:att}, and {cmd:unweighted} under {cmd:unweighted}{p_end}
{synopt:{cmd:e(qoi_desc)}}exact quantity label: {cmd:ATE (Average Treatment Effect)}, {cmd:ATT (Average Treatment Effect for the Treated)}, or {cmd:Unweighted (Standard) Fixed Effect} under {cmd:unweighted}{p_end}
{synopt:{cmd:e(estimator)}}estimator code: {cmd:NULL}, {cmd:fd}, {cmd:did}, or {cmd:Mdid}{p_end}
{synopt:{cmd:e(estimator_desc)}}exact estimator label: {cmd:NULL}, {cmd:FD (First-Difference)}, {cmd:DID (Difference-in-Differences)}, or {cmd:DID (Difference-in-Differences) with Matching on Pre-treatment Outcome}{p_end}
{synopt:{cmd:e(vcetype)}}variance-estimator description; exact values are {cmd:Heteroscedastic / Autocorrelation Robust Standard Error}, {cmd:Heteroscedastic Robust Standard Error}, and {cmd:Heteroskedastic Standard Error (Stock-Watson Bias-Corrected)}{p_end}
{synopt:{cmd:e(white_test)}}White decision label; {cmd:TRUE} means reject the null of no misspecification at {cmd:e(white_alpha)}, and {cmd:FALSE} means fail to reject; stored only when White runs{p_end}
{synopt:{cmd:e(predict)}}{cmd:wfe_predict}{p_end}
{synopt:{cmd:e(estat_cmd)}}{cmd:wfe_estat}{p_end}
{synopt:{cmd:e(_caller_n)}}hidden caller-side observation count captured at estimation time; informational only for shared {cmd:predict} replay{p_end}
{synopt:{cmd:e(_xb_varspec)}}hidden score specification replayed by {cmd:predict}'s score-design guard{p_end}
{synopt:{cmd:e(_xb_datasignature)}}hidden value-only data-signature fingerprint replayed by {cmd:predict}'s score-design guard{p_end}
{synopt:{cmd:e(datasignature)}}stored value-only data-signature fingerprint recorded for shared {cmd:predict} replay{p_end}
{synopt:{cmd:e(datasignaturevars)}}dependent variable plus scored regressors replayed by {cmd:predict}, as recorded by {cmd:signestimationsample}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synoptlines}
{p2colreset}{...}

{pstd}
The postestimation interface is intentionally narrow:
{cmd:predict} supports only {cmd:xb}, {cmd:fitted}, and {cmd:residuals},
and those statistics are mutually exclusive,
while {cmd:estat} currently implements only {cmd:estat wfe_weights}.
{cmd:estat wfe_weights} accepts no additional arguments or options.
A bare trailing comma, as in {cmd:estat wfe_weights,}, is accepted and treated exactly like {cmd:estat wfe_weights}.
If extra arguments or options are supplied, {cmd:estat wfe_weights}
stops with {cmd:estat wfe_weights does not accept arguments or options}.
If the stored weight matrix {cmd:e(W)} is unavailable, {cmd:estat wfe_weights}
stops with {cmd:e(W) is not available}.
Failed {cmd:estat wfe_weights} calls clear any previously returned {cmd:r()} weight-summary scalars rather than leaving a stale summary replayable after the error.

{pstd}
The hidden macros {cmd:e(_xb_varspec)} and {cmd:e(_xb_datasignature)}
store the exact score specification and its value-only data-signature fingerprint replayed by {cmd:predict}.
They are the score-design replay contract used by the shared postestimation
code before it recomputes {cmd:xb}, {cmd:fitted}, or raw residuals.
The visible macros {cmd:e(datasignature)} and {cmd:e(datasignaturevars)}
come from {cmd:signestimationsample} and document the caller-side
estimation sample in the standard Stata way, while the replay guard
stores only the normalized numeric values from that sample so harmless
byte/float/double storage recasts do not break {cmd:predict}.
The score-design guard in {cmd:wfe_predict} replays the hidden {_xb_*}
pair above.
If either hidden score-design replay key {cmd:e(_xb_varspec)} or
{cmd:e(_xb_datasignature)} is missing or blank, {cmd:predict, xb}
and {cmd:predict, fitted} stop with
{cmd:wfe_predict: stored score-design specification is unavailable}
or
{cmd:wfe_predict: stored score-design signature is unavailable},
rather than silently downgrading to a weaker replay check.
The same hidden replay-key failures also stop {cmd:predict, residuals}
with those exact messages, because raw residual replay depends on the
same stored score-design contract.
If the visible macro {cmd:e(depvar)} is missing or blank, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored dependent variable name is unavailable}.
If the visible macro {cmd:e(datasignature)} is missing or blank, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored data-signature is unavailable}.
The hidden macro {cmd:e(_caller_n)} stores the caller-side observation
count captured at estimation time, but it is now informational only.
Shared {cmd:predict} no longer treats {cmd:e(_caller_n)} as a replay
guard; valid {cmd:predict} calls continue to work when that hidden macro
is blank, missing, nonpositive, or noninteger, provided the stored
estimation sample and the statistic-specific replay fingerprints remain
reproducible.
When {cmd:predict} is called without an option, it defaults to
{cmd:xb}; {cmd:fitted} is an alias of {cmd:xb}.
If {cmd:xb} and {cmd:fitted} are requested together, {cmd:predict} stops with {cmd:xb and fitted may not be combined}.  If either {cmd:xb} or {cmd:fitted} is combined with {cmd:residuals}, {cmd:predict} stops with {cmd:xb/fitted and residuals may not be combined}.

{pstd}
All shared {cmd:predict} paths require the current data to retain the
same estimation-sample rows that were used at estimation time.  Added or
dropped observations outside {cmd:e(sample)} do not by themselves break
replay.  If the stored estimation sample or the variables needed for the
requested statistic have changed so that the sample can no longer be
replayed, {cmd:predict} stops immediately
with {cmd:wfe_predict: current data no longer match stored estimation sample}
before downstream replay-signature checks are evaluated.
Statistic-specific exact-name guards for required current variables may fail earlier.

{pstd}
When {cmd:estat wfe_weights} is run, the command stores the displayed
weight-summary scalars in {cmd:r()}, including the matrix dimensions,
{cmd:r(total)}, all-weight moments, nonzero-weight moments, and, when
defined, positive/negative-weight moments and {cmd:r(neg_ratio)}.

{pstd}
The exact scalar keys are {cmd:r(T)}, {cmd:r(N)}, {cmd:r(total)},
{cmd:r(n_nonzero)}, {cmd:r(n_positive)}, {cmd:r(n_negative)},
{cmd:r(w_min)}, {cmd:r(w_max)}, {cmd:r(w_mean)}, {cmd:r(w_sd)},
{cmd:r(nz_min)}, {cmd:r(nz_max)}, {cmd:r(nz_mean)}, {cmd:r(nz_sd)},
{cmd:r(pos_min)}, {cmd:r(pos_max)}, {cmd:r(pos_mean)}, {cmd:r(pos_sd)},
{cmd:r(neg_min)}, {cmd:r(neg_max)}, {cmd:r(neg_mean)},
{cmd:r(neg_sd)}, and {cmd:r(neg_ratio)}.  The positive-weight moment
entries may be missing when the stored weight matrix contains no
positive weights.  The negative-weight moment entries may be missing
when the stored weight matrix contains no negative weights.
If {cmd:r(n_nonzero) = 0}, {cmd:r(neg_ratio)} is reported as missing because the share of negative weights among nonzero weights is undefined on an empty nonzero-weight support.
If {cmd:r(n_negative) = 0} but {cmd:r(n_nonzero) > 0}, {cmd:r(neg_ratio)} is still reported and displayed as {cmd:0}.
If the stored weight matrix {cmd:e(W)} contains missing values,
{cmd:estat wfe_weights} stops with
{cmd:estat wfe_weights: e(W) must not contain missing values}
rather than reporting undefined moments.  Here {cmd:r(neg_ratio)} is the share of negative weights among the nonzero weights, that is,
{cmd:r(n_negative) / r(n_nonzero)}.

{pstd}
Under {cmd:wfe}, {cmd:predict, xb} and {cmd:predict, fitted}
depend only on the scored regressors and stored coefficients.
Changing only the current variable named in {cmd:e(depvar)}
does not alter {cmd:predict, xb}.  If the current dependent variable
named in {cmd:e(depvar)} is renamed, dropped, or otherwise unavailable
after estimation, {cmd:predict, xb} and {cmd:predict, fitted} remain
available, because the fitted-value object depends only on the scored
regressors and stored coefficients.  Because the variable supplied in
{cmd:treat()} is always part of the scored regressors, editing its
current values after estimation triggers the same score-variable
mismatch guard as any other regressor edit.  Renaming or dropping that
treatment variable after estimation likewise makes one of the
estimation-time regressors unavailable to {cmd:predict}.  If the scored regressors have
been edited after estimation, {cmd:predict}
stops with {cmd:wfe_predict: current score variables no longer match the stored estimation sample}
rather than silently
changing the fitted-value object.  If renaming or dropping a scored
regressor means the estimation-time regressor is no longer available in the
current data, {cmd:predict, xb} and {cmd:predict, fitted} fail immediately
with {cmd:wfe_predict: regressors from estimation are no longer available in current data}.
Changing only the numeric storage type of the current scored variables while leaving their {cmd:e(sample)} values unchanged does not alter {cmd:predict, xb}, {cmd:predict, fitted}, or {cmd:predict, residuals}.
Re-sorting the data does not change these rowwise
{cmd:predict, xb}/{cmd:fitted} or {cmd:predict, residuals} contracts,
because both statistics are recomputed on the current estimation-sample
rows rather than replayed in estimation-time order.

{pstd}
Under {cmd:wfe}, {cmd:predict, residuals} returns the raw residual
{cmd:y - xb}.  In the
Stata implementation, {cmd:y} is taken from the current variable named in
{cmd:e(depvar)} on the estimation sample.  If the visible macro
{cmd:e(depvar)} is missing or blank, {cmd:predict, residuals} stops with
{cmd:wfe_predict: stored dependent variable name is unavailable}.  If the current dependent
variable named in {cmd:e(depvar)} is renamed, dropped, or otherwise
unavailable after estimation, {cmd:predict, residuals} stops with
{cmd:wfe_predict: stored dependent variable is no longer available in current data}.
If that dependent variable or the scored regressors have instead been
edited in place after estimation, {cmd:predict} stops with
{cmd:wfe_predict: current score variables no longer match the stored estimation sample}
rather than silently changing the
residual object.  Because the variable supplied in {cmd:treat()} is
always part of the scored regressors, editing its current values after
estimation triggers the same score-variable mismatch guard before the
residual object is recomputed.  It does not add absorbed unit/time fixed effects back in, so
this postestimation residual is not a within-transformed FE
disturbance.

{pstd}
Changing {cmd:unit()} or {cmd:time()} values alone does not alter {cmd:predict, xb}, {cmd:predict, fitted}, or {cmd:predict, residuals}.
Renaming or dropping the current {cmd:unit()} or {cmd:time()} variables alone likewise does not affect these shared {cmd:predict} statistics. These postestimation replays depend on the current dependent variable and scored regressors rather than on the panel identifiers once estimation is complete.

{marker examples}{...}
{title:Examples}

{pstd}
The following examples generate a balanced panel in memory and then walk
through progressively richer uses of {cmd:wfe}.  All examples assume
Stata 16 or later.

{phang}
{cmd:. clear}
{cmd:. set seed 12345}
{cmd:. local N = 10}
{cmd:. local TT = 15}
{cmd:. set obs `=`N' * `TT''}
{cmd:. gen int id = ceil(_n / `TT')}
{cmd:. gen int t = mod(_n - 1, `TT') + 1}
{cmd:. gen byte tr = rbinomial(1, 0.25)}
{cmd:. bysort id: egen double __unit_mean = mean(tr)}
{cmd:. bysort t: egen double __time_mean = mean(tr)}
{cmd:. count if inlist(__unit_mean, 0, 1) | inlist(__time_mean, 0, 1)}
{cmd:. while r(N) > 0 {c -(}}
{cmd:>     replace tr = rbinomial(1, 0.25)}
{cmd:>     drop __unit_mean __time_mean}
{cmd:>     bysort id: egen double __unit_mean = mean(tr)}
{cmd:>     bysort t: egen double __time_mean = mean(tr)}
{cmd:>     count if inlist(__unit_mean, 0, 1) | inlist(__time_mean, 0, 1)}
{cmd:> {c )-}}
{cmd:. drop __unit_mean __time_mean}
{cmd:. bysort id: gen double alpha_i = cond(_n == 1, rnormal(), .)}
{cmd:. bysort id: replace alpha_i = alpha_i[1]}
{cmd:. gen double x1 = rnormal(0.5, 1)}
{cmd:. gen double x2 = rbeta(5, 1)}
{cmd:. gen double myweight = 1 + 0.05 * id + 0.01 * t}
{cmd:. gen double y = alpha_i + tr + x1 + x2 + rnormal()}

{pstd}
Example 1: basic one-way unit FE with ATE.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id)}

{pstd}
Example 2: supply an explicit time variable.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t)}

{pstd}
Example 3: one-way time fixed effects.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) method(time)}

{pstd}
Example 4: ATT instead of ATE.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) qoi(att)}

{pstd}
Example 5: first-difference estimation.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(fd)}

{pstd}
Example 6: multi-period DiD.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(did)}

{pstd}
Example 7: matched DiD with nearest-neighbor matching.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(Mdid)}

{pstd}
Example 8: matched DiD with a maximum deviation of 0.5.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(Mdid) maxdev_did(0.5)}

{pstd}
Example 9: matched DiD with a wider deviation bound.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(Mdid) maxdev_did(2)}

{pstd}
Example 10: disable the White test.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) nowhite}

{pstd}
Example 11: heteroskedasticity-robust but not autocorrelation-robust
standard errors for a one-way model.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) hetero_se(on) auto_se(off)}

{pstd}
Example 12: custom {it:C_it} weights.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) cit(myweight)}

{pstd}
Example 13: standard unweighted two-way fixed effects.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(did) unweighted}

{pstd}
Example 14: turn off the two-way FE degrees-of-freedom adjustment.

{phang}
{cmd:. wfe y x1 x2, treat(tr) unit(id) time(t) estimator(did) df_adjustment(off)}

{pstd}
Example 15: fit a richer model with multiple regressors and store the
weighted demeaned data.

{phang}
{cmd:. gen double x3 = x1 * x2}
{cmd:. wfe y x1 x2 x3, treat(tr) unit(id) time(t) store_wdm}

{pstd}
Example 16: current postestimation commands.

{phang}
{cmd:. predict double xbhat, xb}
{cmd:. predict double yhat, fitted}
{cmd:. predict double resid, residuals}
{cmd:. estat wfe_weights}

{marker references}{...}
{title:References}

{marker IK2021}{...}
{phang}
Imai, K., and I. S. Kim. 2021. "On the Use of Two-Way Fixed Effects
Regression Models for Causal Inference with Panel Data."
{it:Political Analysis} 29(3): 405-415.

{marker ARELLANO1987}{...}
{phang}
Arellano, M. 1987. "Computing Robust Standard Errors for Within-Groups
Estimators." {it:Oxford Bulletin of Economics and Statistics}
49(4): 431-434.

{marker STOCKWATSON2008}{...}
{phang}
Stock, J. H., and M. W. Watson. 2008. "Heteroskedasticity-Robust Standard
Errors for Fixed Effects Panel Data Regression." {it:Econometrica}
76(1): 155-174.

{marker WHITE1980}{...}
{phang}
White, H. 1980. "Using Least Squares to Approximate Unknown Regression
Functions." {it:International Economic Review} 21(1): 149-170.

{marker alsosee}{...}
{title:Also see}

{psee}
Manual:  {helpb xtreg}, {helpb xtdidregress}

{psee}
Postestimation:  {helpb predict}, {helpb estat}

{psee}
Related command:  {helpb pwfe}
