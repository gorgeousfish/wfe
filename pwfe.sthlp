{smcl}
{* *! pwfe.sthlp}{...}
{vieweralsosee "wfe" "help wfe"}{...}
{vieweralsosee "xtreg" "help xtreg"}{...}
{vieweralsosee "logit" "help logit"}{...}
{vieweralsosee "predict" "help predict"}{...}
{vieweralsosee "estat" "help estat"}{...}
{viewerjumpto "Syntax" "pwfe##syntax"}{...}
{viewerjumpto "Description" "pwfe##description"}{...}
{viewerjumpto "Options" "pwfe##options"}{...}
{viewerjumpto "Remarks" "pwfe##remarks"}{...}
{viewerjumpto "Stored results" "pwfe##results"}{...}
{viewerjumpto "Examples" "pwfe##examples"}{...}
{viewerjumpto "References" "pwfe##references"}{...}
{viewerjumpto "Also see" "pwfe##alsosee"}{...}

{title:Title}

{p2colset 5 16 24 2}{...}
{p2col:{cmd:pwfe} {hline 2}}Propensity-score weighted fixed effects estimator for panel causal inference{p_end}
{p2colreset}{...}

{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:pwfe} [{it:varlist}] {ifin} {cmd:,}
{opth treat(varname)} {opth outcome(varname)} {opth unit(varname)} [{it:options}]
{p_end}

{pstd}
{it:varlist} contains the propensity-score covariates used to estimate
the treatment propensity score.  It is not the outcome
equation.  The transformed and reported outcome variable is supplied by
{cmd:outcome()}.

{pstd}
Standard Stata factor-variable notation is allowed in {it:varlist}.  For
example, terms such as {cmd:i.group} and {cmd:c.x##i.group} are expanded
before the propensity-score bridge is called.

{pstd}
If {cmd:pscore()} is not specified, {it:varlist} is optional.
Omitting it requests an intercept-only propensity-score specification.
If {cmd:pscore()} is
specified, {it:varlist} must be omitted.  The current implementation
requires Stata 16 or later.

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth treat(varname)}}binary treatment indicator coded as 0/1{p_end}
{synopt:{opth outcome(varname)}}outcome variable to be transformed and regressed{p_end}
{synopt:{opth unit(varname)}}panel unit identifier{p_end}

{syntab:Model}
{synopt:{opth time(varname)}}time identifier; required by {cmd:method(time)} and {cmd:estimator(fd)}{p_end}
{synopt:{opt method(string)}}{cmd:unit} (default) or {cmd:time}{p_end}
{synopt:{opt qoi(string)}}{cmd:ate} (default) or {cmd:att}{p_end}
{synopt:{opt estimator(string)}}baseline PWFE or {cmd:fd}{p_end}
{synopt:{opth cit(varname)}}nonnegative user-supplied {it:C_it} weights; fractional values are truncated toward zero{p_end}

{syntab:Propensity score}
{synopt:{opth pscore(varname)}}precomputed propensity score in (0,1){p_end}
{synopt:{opt nowithin_unit}}use pooled logit and pooled transform{p_end}

{syntab:Standard errors}
{synopt:{opt hetero_se(on|off)}}heteroskedasticity option; default {cmd:on}{p_end}
{synopt:{opt auto_se(on|off)}}autocorrelation option; default {cmd:on}{p_end}
{synopt:{opt unbiased_se}}post Stock-Watson bias-corrected covariance matrices{p_end}

{syntab:White test}
{synopt:{opt [no]white}}run the White misspecification test; default {cmd:white}{p_end}
{synopt:{opt white_alpha(#)}}significance level in {cmd:(0,1)} for the White test; default {cmd:0.05}{p_end}

{syntab:Output}
{synopt:{opt [no]verbose}}display progress messages; default {cmd:verbose}{p_end}
{synopt:{opt diagnose}}prints the parsed execution state and exits without estimation{p_end}
{synoptline}

{marker description}{...}
{title:Description}

{pstd}
{cmd:pwfe} first transforms the variable named in {cmd:outcome()} with a
propensity score and then estimates a fixed-effects or first-difference
regression on the transformed outcome.  The treatment indicator is given
by {cmd:treat()} and must be binary with values 0 and 1.

{pstd}
The command supports three propensity-score paths.  If {cmd:pscore()} is
provided, {cmd:pwfe} skips internal logit estimation and transforms the
full estimation sample directly.  Without {cmd:pscore()}, the default
{cmd:method(unit)} estimates a separate logit within each unit and
transforms within each unit.  With {cmd:method(time)}, the default logic
estimates a separate logit within each time period and transforms within
each time period.  {cmd:nowithin_unit} replaces these split paths with a
single pooled logit and pooled transform.

{pstd}
Whenever {cmd:pwfe} estimates propensity scores internally with a
nonempty {it:varlist}, it uses a no-constant logit, that is,
{cmd:logit, noconstant}, rather than a logit with an added intercept.
When
{it:varlist} is omitted, {cmd:pwfe} instead fits an
intercept-only logit.

{pstd}
Whenever {cmd:pwfe} estimates propensity scores internally, each split
or pooled propensity-score fit is handled directly by the penalized
binomial IRLS path under the
same no-constant formula described above.  There is no preliminary plain
{helpb logit} stage in the current implementation.  The command stops
with {cmd:pwfe propensity-score logit failed in ...; consider providing pscore()}
only if that weakly regularized fit still fails to deliver usable
propensity scores for the relevant pool.  The exact pool-scoped runtime
messages are {cmd:pwfe propensity-score logit failed in unit pool; consider providing pscore()},
{cmd:pwfe propensity-score logit failed in time pool; consider providing pscore()},
and {cmd:pwfe propensity-score logit failed in global sample; consider providing pscore()}.
Under this weakly regularized fit, a split unit/time pool may still be admissible when it contains only treated or only control observations.
Equivalently, split unit/time pools may still be single-support under the penalized propensity-score fit, so successful {cmd:method(unit)} and
{cmd:method(time)} runs do not require every pool to contain both
treatment arms.

{pstd}
Every internal
propensity-score fit is routed through a penalized binomial IRLS
routine with an adaptive Cauchy prior on the coefficients.  This applies to the default
split-by-unit path, {cmd:method(time)}, and the pooled
{cmd:nowithin_unit} path.  The resulting fitted propensity scores stay
strictly inside {cmd:(0,1)} instead of preserving successful
plain-logit fits as-is.

{pstd}
If {it:if}/{it:in} restrictions and required-variable or propensity-score
complete-case screening leave no estimation-sample observations,
{cmd:pwfe} stops with the standard Stata error {cmd:no observations}
before any transformed-outcome or weighting work begins.

{pstd}
{cmd:pwfe} sorts internally by unit-time for index construction and Mata
estimation, but it restores the caller's original observation order before
returning.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth treat(varname)} specifies the binary treatment indicator.  The
current implementation rejects missing values with
{cmd:missing values in treat() are not allowed}.  Nonbinary data are
rejected with {cmd:'treat' must be a binary vector} or
{cmd:'treat' must be either 0 or 1 where 1 indicates treated}.  If the
caller-selected sample contains only one observed treatment level,
{cmd:pwfe} fails fast with {cmd:'treat' must contain both 0 and 1 values}.
Because {cmd:treat()} is a structural estimator input rather than a
propensity-score covariate, nonnumeric input is rejected with
{cmd:treat() must specify a numeric variable}.  Missing values and nonbinary treatment codings are rejected on the caller-selected sample even when the same rows would later be dropped because {it:varlist} is missing.

{phang}
{opth outcome(varname)} specifies the observed outcome variable that is
transformed before estimation.  {cmd:e(depvar)} records this variable
name, not a name from {it:varlist}.  Nonnumeric input is rejected with
{cmd:outcome() must specify a numeric variable}.  Missing values are
rejected with {cmd:missing values in outcome() are not allowed} because
the transformed outcome is undefined without an observed response.

{phang}
{opth unit(varname)} specifies the panel unit identifier.  Missing values
are rejected with {cmd:missing values in unit() are not allowed} because
the unit-time panel dyad is undefined without them.  Both {cmd:unit()}
and {cmd:time()} may be numeric or string identifiers.

{dlgtab:Model}

{phang}
{opth time(varname)} specifies the time identifier.  If omitted, {cmd:pwfe}
generates an internal time index for the admissible one-way path.  When
{cmd:time()} is omitted, that internal index follows the current
within-unit observation order before estimation.  The option is required
by {cmd:method(time)} and by {cmd:estimator(fd)}.  When {cmd:time()} is
supplied, missing values are rejected with
{cmd:missing values in time() are not allowed}, and duplicate
{cmd:(unit,time)} observations are rejected with
{cmd:duplicate observations found for unit-time combination}.  Duplicate
{cmd:(unit,time)} observations are checked only on the formula-valid
sample that survives propensity-score covariate missing-value removal.
Duplicate dyads that appear only in rows dropped because {it:varlist}
covariates are missing do not trigger that parser error.  Current
{cmd:pwfe} factorizes either numeric or string panel identifiers to the
dense internal indices.

{phang}
{opt method(string)} chooses the fixed-effect dimension.  Allowed values
are {cmd:unit} (default) and {cmd:time}.  Any other value is rejected
with {cmd:method should be either unit or time}.  Exact parsing errors
are:

{phang2}
missing {cmd:time()} under {cmd:method(time)}:
{cmd:'time.index' should be provided}

{phang2}
{cmd:method(time)} with {cmd:estimator(fd)}:
{cmd:First Difference is not compatible with 'time' method}

{phang2}
For matched-set feasibility, {cmd:method(time)} requires at least 2 units
in the estimation sample.  If only one unit remains, the runtime path
fails fast from the time-weight helper with
{cmd:wfe_weights_time: time FE requires at least 2 units}.

{phang}
{opt qoi(string)} selects {cmd:ate} or {cmd:att}.  Any other value is
rejected with {cmd:qoi should be either ate or att}.

{phang}
{opt estimator(string)} leaves the baseline one-way PWFE path in place
when omitted and allows {cmd:fd} as the alternative successful path.
The following inputs are rejected and are documented here because they are
part of the current public contract:

{phang2}
missing {cmd:time()} with {cmd:estimator(fd)}:
{cmd:First Difference cannot calculate when 'time.index' is missing}

{phang2}
{cmd:estimator(did)}:
{cmd:Difference-in-Differences is not compatible with pwfe}

{phang2}
{cmd:estimator(Mdid)}:
{cmd:pwfe does not support estimator('Mdid')}

{phang2}
any other nonempty {cmd:estimator()}:
{cmd:estimator() must be 'fd'}

{phang2}
For matched-set feasibility, {cmd:method(unit)} and {cmd:estimator(fd)}
each require at least 2 time periods in the estimation sample.  If only
one period remains, the {cmd:method(unit)} runtime path fails fast from
{cmd:wfe_weights_unit: unit FE requires at least 2 time periods}.  If
only one period remains, the {cmd:estimator(fd)} runtime path fails fast
from {cmd:wfe_weights_fd: FD requires at least 2 time periods}.

{phang}
{opth cit(varname)} supplies a nonnegative {it:C_it} weight variable.
Fractional {it:C_it} values are truncated toward zero before weights are constructed.  Nonnumeric input is rejected with {cmd:'C.it' must be a numeric vector with length equal to number of observations}, missing values are rejected with {cmd:missing values in cit() are not allowed}, and negative values are rejected with {cmd:'C.it' must be a non-negative numeric vector}.
After truncation, {it:C_it} must also remain within the signed 32-bit integer
range required by the underlying one-way weight generators; otherwise
{cmd:pwfe} stops with {cmd:_wfe_pwfe_estimate: C_it must lie within the signed 32-bit integer range after truncation}.

{dlgtab:Propensity score}

{phang}
{opth pscore(varname)} supplies a precomputed propensity score.  The
variable must be numeric, nonmissing on the estimation sample, and
strictly between 0 and 1.  When this option is used, {it:varlist}
must be omitted; otherwise {cmd:pwfe} returns {cmd:'formula' should not be provided when pscore is specified}.  Nonnumeric input is rejected with {cmd:'pscore' must be a numeric vector with length equal to number of observations}.  Missing values in {cmd:pscore()} are rejected with {cmd:missing values in pscore() are not allowed}.  Out-of-range values are rejected with {cmd:'pscore' must be strictly between 0 and 1}.
When {cmd:pscore()} is supplied, {cmd:nowithin_unit} has no effect because {cmd:pwfe} skips internal propensity-score estimation and always applies the user-supplied scores globally.

{phang}
{opt nowithin_unit} replaces split-by-unit or split-by-time estimation
with a pooled logit and pooled transform over the full estimation sample.
If {it:varlist} is omitted, this pooled logit is intercept-only and
estimates a single propensity score over the full estimation sample.

{dlgtab:Standard errors}

{phang}
{opt hetero_se(on|off)} and {opt auto_se(on|off)} determine which of the
three currently supported covariance paths is used: HAC when both are
{cmd:on}, HC when {cmd:hetero_se(on) auto_se(off)}, and homoskedastic
when {cmd:hetero_se(off) auto_se(off)}.  HAC still requires at least 2 units, whereas plain HC does not.
The Stock-Watson correction keeps the 2-unit requirement described below.
Any other value is rejected at
the parser with {cmd:hetero_se() must be 'on' or 'off'} or
{cmd:auto_se() must be 'on' or 'off'}.

{phang}
{cmd:hetero_se(off) auto_se(on)} is rejected with the exact error
{cmd:robust standard errors with autocorrelation and homoskedasticity is not supported}.

{phang}
{opt unbiased_se} posts the Stock-Watson bias-corrected PWFE and FE covariance matrices.
It is allowed only with {cmd:hetero_se(on) auto_se(off)}; any HAC or
homoskedastic combination returns the exact parser message
{cmd:unbiased_se is allowed only for hetero_se(on) auto_se(off)}.  The correction is only defined for
balanced PWFE and FE unit panels with at least 3 time periods, and it
still requires at least 2 units before the Stock-Watson branch is
admissible.  Plain HC remains available with a single unit as long as
the HC denominator {cmd:NT - J_u - p} is positive.  These constraints match the cluster count and
{cmd:(J_t-1)/(J_t-2)} factors used by the implementation.  When the
balanced panel has fewer than 3 time periods, {cmd:pwfe} stops with
{cmd:unbiased_se requires at least 3 time periods}.
In the current implementation that branch replaces the final reported PWFE and FE covariance matrices with the Stock-Watson bias-corrected versions.

{dlgtab:White test}

{phang}
{opt [no]white} turns the White misspecification test on or off.  When
enabled, the command stores {cmd:e(white_stat)}, {cmd:e(white_pvalue)},
{cmd:e(white_alpha)}, and {cmd:e(white_test)}.  The saved
{cmd:e(white_test)} macro takes the string value {cmd:TRUE} or
{cmd:FALSE}.  {cmd:TRUE} means reject the null of no misspecification at
{cmd:e(white_alpha)}, while {cmd:FALSE} means fail to reject.

{phang}
The current Stata PWFE White path requires the covariance-difference
matrix {cmd:Phi} to be positive semidefinite and full rank for
chi-square inference.  If that condition fails, {cmd:pwfe} aborts the
White path instead of posting fallback pseudo-statistics.

{phang}
{opt white_alpha(#)} sets the White-test significance level.
The option is a probability threshold, so white_alpha() must lie in
{cmd:(0,1)}.  Values outside this interval or nonnumeric input fail fast
with the exact message {cmd:white_alpha() must lie in (0,1)}.

{dlgtab:Output}

{phang}
{opt [no]verbose} controls progress output from the internal logit and
estimation steps.

{phang}
{opt diagnose} prints the parsed execution state and exits without
estimation.  The current implementation reports the parser-selected
values for {cmd:method}, {cmd:qoi}, {cmd:estimator}, {cmd:within_unit},
{cmd:hetero_se}, {cmd:auto_se}, the diagnostic output label
{cmd:unbiased} (which reflects the state of the {cmd:unbiased_se}
option), {cmd:white},
{cmd:white_alpha}, {cmd:verbose}, {cmd:pscore}, {cmd:cit}, and
{cmd:N_parse}.  Because the command exits after that diagnostics block,
{cmd:diagnose} does not post estimation results in {cmd:e()} and clears
any active estimation results so stale {cmd:predict} or {cmd:estat}
output cannot be replayed after the diagnostic run.

{marker remarks}{...}
{title:Remarks}

{phang}
{bf:Differences from wfe.}  Relative to {helpb wfe}, {cmd:pwfe} adds
{cmd:outcome()}, {cmd:pscore()}, and {cmd:nowithin_unit}.  It does not
document or expose {cmd:wfe}-specific options such as the alternative
baseline-weight path or weighted-demeaned-data storage.

{phang}
{bf:Standard-error and White semantics.}  HAC, HC, and homoskedastic
paths are all available in {cmd:pwfe}.  HAC clustering is always tied to
the unit index.  The White statistic uses the cross residual product
{it:u_hat * u_tilde}.

{phang}
{bf:Method(time) FE benchmark.}  Under {cmd:method(time)}, {cmd:e(b_fe)} and
{cmd:e(V_fe)} come from the time fixed-effects benchmark.  The White
comparison uses that same time-demeaned FE side.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:pwfe} stores the following in {cmd:e()}:

{synoptset 28 tabbed}{...}
{p2coldent:* Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations posted by {cmd:ereturn post}{p_end}
{synopt:{cmd:e(N_units)}}number of units in the estimation sample{p_end}
{synopt:{cmd:e(N_times)}}number of estimation-sample time groups; with omitted {cmd:time()}, the maximum internally generated dense within-unit time index{p_end}
{synopt:{cmd:e(df_r)}}residual degrees of freedom{p_end}
{synopt:{cmd:e(sigma)}}root mean squared error scale{p_end}
{synopt:{cmd:e(sigma2)}}variance scale used for homoskedastic formulas{p_end}
{synopt:{cmd:e(N_nonzero)}}count of nonzero entries in the stored weight matrix {cmd:e(W)}{p_end}
{synopt:{cmd:e(white_stat)}}White test statistic when White is enabled{p_end}
{synopt:{cmd:e(white_pvalue)}}White test p-value when White is enabled{p_end}
{synopt:{cmd:e(white_alpha)}}White test significance level when White is enabled{p_end}

{pstd}
In the current implementation, {cmd:e(sigma)} and {cmd:e(sigma2)} come
from the PWFE/WFE weighted-demeaned regression.  They do not describe the FE benchmark
reported in {cmd:e(b_fe)} and {cmd:e(V_fe)}.

{p2coldent:* Matrices}{p_end}
{synopt:{cmd:e(b)}}PWFE coefficient row vector{p_end}
{synopt:{cmd:e(V)}}PWFE covariance matrix{p_end}
{synopt:{cmd:e(W)}}{it:T x N} PWFE regression weight matrix; rows are estimation-sample time groups and columns are units{p_end}
{synopt:{cmd:e(b_fe)}}fixed-effects benchmark coefficient row vector{p_end}
{synopt:{cmd:e(V_fe)}}fixed-effects benchmark covariance matrix{p_end}
{synopt:{cmd:e(y_star)}}transformed outcome vector used by {cmd:pwfe}{p_end}
{synopt:{cmd:e(pscore)}}user or fitted propensity scores used for transform{p_end}
{synopt:{cmd:e(_y_star_unit_idx)}}rowwise dense unit-index keys aligned with {cmd:e(y_star)}{p_end}
{synopt:{cmd:e(_y_star_time_idx)}}rowwise dense time-index keys aligned with {cmd:e(y_star)}{p_end}

{pstd}
The stored column vectors {cmd:e(y_star)} and {cmd:e(pscore)} follow the
estimation-sample order rather than the internal unit-time sort used
during Mata computation, so they map directly back to the preserved data
rows after {cmd:pwfe} returns.

{pstd}
The companion keyed matrices {cmd:e(_y_star_unit_idx)} and
{cmd:e(_y_star_time_idx)} use the same caller-order row alignment, but
their entries are dense internal indices rather than the raw
{cmd:unit()} / {cmd:time()} values.  They therefore identify the stored
objects by the internal {cmd:1..e(N_units)} and {cmd:1..e(N_times)}
clock that {cmd:pwfe} actually estimates on.

{pstd}
When {cmd:time()} is omitted, {cmd:e(N_times)} records the maximum
internally generated dense within-unit time index after estimation-sample
restriction, rather than a count taken from a user-supplied time
variable.

{pstd}
When {cmd:time()} is supplied, {cmd:e(N_times)} records the number
of distinct estimation-sample time groups induced by the supplied
{cmd:time()} variable.

{pstd}
The stored matrix {cmd:e(W)} uses {it:T x N} orientation: rows index
estimation-sample time groups and columns index units.

{pstd}
When {cmd:time()} is omitted, the rows of {cmd:e(W)} are indexed by that
same internal dense within-unit time clock rather than a user-supplied
time variable.

{pstd}
When {cmd:time()} is omitted, {cmd:e(_y_star_time_idx)} stores the same
dense within-unit observation-order key used internally after
estimation-sample restriction.  In that omitted-time path, the saved key
documents how {cmd:e(y_star)} and {cmd:e(pscore)} line up with the
preserved caller data while still matching the dense {cmd:1..e(N_times)}
clock used for weighting and demeaning.

{pstd}
When {cmd:time()} was supplied at estimation, {cmd:e(_y_star_unit_idx)}
and {cmd:e(_y_star_time_idx)} store rowwise dense keys aligned with
{cmd:e(y_star)} and {cmd:e(pscore)}.  These keyed matrices document the
stored unit-time identity of the transformed-outcome objects on the
internal dense support, not in the raw caller value space, and can be
used to re-align {cmd:e(y_star)} and {cmd:e(pscore)} programmatically
after a data re-sort.  The current {cmd:predict, residuals} contract
does not replay those stored transformed outcomes.

{pstd}
The macros {cmd:e(_unitvar)} and {cmd:e(_timevar)} identify which caller
variables those saved rowwise keys refer to.  {cmd:e(_unitvar)} stores
the name supplied to {cmd:unit()}, and {cmd:e(_timevar)} stores the name
supplied to {cmd:time()} or is empty when {cmd:time()} is omitted.
The companion macro {cmd:e(_y_star_var)} is currently the empty string,
which documents that {cmd:pwfe} leaves no caller-side transformed-outcome
cache variable behind and instead stores the transformed outcome only in
{cmd:e(y_star)}.

{p2coldent:* Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:pwfe}{p_end}
{synopt:{cmd:e(depvar)}}name supplied to {cmd:outcome()}{p_end}
{synopt:{cmd:e(_treatvar)}}hidden treatment variable name used by shared postestimation guards{p_end}
{synopt:{cmd:e(_unitvar)}}name supplied to {cmd:unit()}{p_end}
{synopt:{cmd:e(_timevar)}}name supplied to {cmd:time()}; empty when {cmd:time()} is omitted{p_end}
{synopt:{cmd:e(_y_star_var)}}empty string; no caller-side transformed-outcome cache variable{p_end}
{synopt:{cmd:e(_caller_n)}}hidden caller-side observation count captured at estimation time; informational only for shared {cmd:predict} replay{p_end}
{synopt:{cmd:e(_xb_datasignature)}}hidden treatment-only data-signature fingerprint used by {cmd:predict, xb} and {cmd:predict, fitted}{p_end}
{synopt:{cmd:e(method)}}exact method label: {cmd:unit} or {cmd:time}{p_end}
{synopt:{cmd:e(qoi)}}exact qoi code: {cmd:ate} or {cmd:att}{p_end}
{synopt:{cmd:e(qoi_desc)}}exact causal estimand label: {cmd:ATE (Average Treatment Effect)} or {cmd:ATT (Average Treatment Effect for the Treated)}{p_end}
{synopt:{cmd:e(estimator)}}exact estimator code: {cmd:NULL} under the default baseline path, or {cmd:fd} under {cmd:estimator(fd)}{p_end}
{synopt:{cmd:e(estimator_desc)}}exact estimator label: {cmd:NULL} under the default baseline path, or {cmd:FD (First-Difference)} under {cmd:estimator(fd)}{p_end}
{synopt:{cmd:e(datasignature)}}stored data-signature fingerprint used by shared {cmd:predict} guards{p_end}
{synopt:{cmd:e(datasignaturevars)}}outcome and treatment variables replayed by {cmd:predict}{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(vcetype)}}one of {cmd:Heteroscedastic / Autocorrelation Robust Standard Error}, {cmd:Heteroscedastic Robust Standard Error}, {cmd:Heteroskedastic Standard Error (Stock-Watson Bias-Corrected)}, or {cmd:Homoskedastic Standard Error}{p_end}
{synopt:{cmd:e(pscore_source)}}{cmd:user} or {cmd:estimated}{p_end}
{synopt:{cmd:e(transform_scope)}}{cmd:global}, {cmd:unit}, or {cmd:time}{p_end}
{synopt:{cmd:e(predict)}}{cmd:wfe_predict}{p_end}
{synopt:{cmd:e(estat_cmd)}}{cmd:wfe_estat}{p_end}
{synopt:{cmd:e(white_test)}}stored White decision label; {cmd:TRUE} means reject the null of no misspecification at {cmd:e(white_alpha)}, and {cmd:FALSE} means fail to reject; stored only when White runs{p_end}
{synoptlines}
{p2colreset}{...}

{pstd}
When {cmd:pscore()} is supplied, {cmd:e(pscore_source)} is {cmd:user} and {cmd:e(transform_scope)} is {cmd:global}.
This same {cmd:user/global} mapping is posted whether or not {cmd:nowithin_unit} is spelled out, because that option only changes the internal propensity-score estimation path that {cmd:pscore()} bypasses.

{pstd}
When {cmd:pwfe} estimates propensity scores internally, {cmd:e(pscore_source)} is {cmd:estimated}; {cmd:e(transform_scope)} is {cmd:unit} for the default split-by-unit path, {cmd:time} for {cmd:method(time)}, and {cmd:global} for {cmd:nowithin_unit}.

{pstd}
{cmd:e(vcetype)} records the exact covariance branch selected by the command:
HAC posts {cmd:Heteroscedastic / Autocorrelation Robust Standard Error}, HC
posts {cmd:Heteroscedastic Robust Standard Error}, the Stock-Watson corrected
branch posts {cmd:Heteroskedastic Standard Error (Stock-Watson Bias-Corrected)},
and {cmd:hetero_se(off) auto_se(off)} posts {cmd:Homoskedastic Standard Error}.

{pstd}
The visible macros {cmd:e(datasignature)} and {cmd:e(datasignaturevars)}
record the stored data-signature fingerprint for the current outcome and
treatment variables replayed by {cmd:predict}.  They document the sample
identity contract enforced before the current outcome/treat columns are
reused for {cmd:xb}, {cmd:fitted}, or {cmd:residuals}.

{pstd}
The hidden macro {cmd:e(_xb_datasignature)} stores the treatment-only
value-only data-signature fingerprint used by {cmd:predict, xb} and
{cmd:predict, fitted}.  It is intentionally narrower than the visible
{cmd:e(datasignature)} fingerprint because the treatment-only linear
predictor does not depend on the current outcome variable.
If the hidden macro {cmd:e(_xb_datasignature)} is missing or blank,
{cmd:predict, xb} and {cmd:predict, fitted} stop with
{cmd:wfe_predict: stored treatment-only signature is unavailable}
rather than silently replaying the treatment-only linear predictor
without its fingerprint guard.

{pstd}
If the hidden macro {cmd:e(_treatvar)} is missing or blank, both {cmd:predict, xb}/{cmd:fitted} and {cmd:predict, residuals} stop with {cmd:wfe_predict: stored treatment variable name is unavailable}.

{pstd}
If the visible macro {cmd:e(depvar)} is missing or blank, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored dependent variable name is unavailable}.

{pstd}
If the visible macro {cmd:e(datasignature)} is missing or blank, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored data-signature is unavailable}.

{pstd}
The hidden macro {cmd:e(_caller_n)} stores the caller-side observation count captured at estimation time, but it is now informational only.  Shared {cmd:predict} no longer treats {cmd:e(_caller_n)} as a replay guard; valid {cmd:predict} calls continue to work when that hidden macro is blank, missing, nonpositive, or noninteger, provided the stored estimation sample and the statistic-specific replay fingerprints remain reproducible.

{pstd}
The current postestimation surface is intentionally narrow:
{cmd:predict} supports the shared {cmd:wfe_predict} options, but only one
statistic may be requested at a time, and
{cmd:estat} currently supports only {cmd:estat wfe_weights}.
{cmd:estat wfe_weights} accepts no additional arguments or options.
A bare trailing comma, as in {cmd:estat wfe_weights,}, is accepted and treated exactly like {cmd:estat wfe_weights}.
If extra arguments or options are supplied, {cmd:estat wfe_weights}
stops with {cmd:estat wfe_weights does not accept arguments or options}.
If the stored weight matrix {cmd:e(W)} is unavailable, {cmd:estat wfe_weights}
stops with {cmd:e(W) is not available}.
Failed {cmd:estat wfe_weights} calls clear any previously returned {cmd:r()} weight-summary scalars rather than leaving a stale summary replayable after the error.
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
Under {cmd:pwfe}, {cmd:predict, xb} and {cmd:predict, fitted} return the treatment-only linear predictor implied by {cmd:e(b)}.  They do not add back absorbed fixed effects or reconstruct an outcome-scale fitted value.

{pstd}
Changing the current outcome values alone does not alter {cmd:predict, xb},
because that predictor depends only on the current treatment values and
the stored coefficients.  Editing the current {cmd:treat()} values
causes {cmd:predict, xb} and {cmd:predict, fitted} to stop with
{cmd:wfe_predict: current treatment variable no longer matches the stored estimation sample}.
Changing only the numeric storage type of the current treatment variable while leaving its {cmd:e(sample)} values unchanged does not alter {cmd:predict, xb} or {cmd:predict, fitted}.
If the current outcome variable named in {cmd:e(depvar)} is renamed, dropped, or otherwise unavailable after estimation, {cmd:predict, xb} and {cmd:predict, fitted} remain available, because the treatment-only linear predictor does not depend on the current outcome variable.

{pstd}
Changing {cmd:unit()} or {cmd:time()} values alone does not alter {cmd:predict, xb} or {cmd:predict, fitted}, because the treatment-only linear predictor does not use the stored unit-time keys published with {cmd:e(y_star)} and {cmd:e(pscore)}.

{pstd}
If the treatment variable named in {cmd:e(_treatvar)} is renamed, dropped, or otherwise unavailable after estimation, both {cmd:predict, xb}/{cmd:fitted} and {cmd:predict, residuals} stop with {cmd:wfe_predict: stored treatment variable is no longer available in current data}.
If the original {cmd:treat()} variable is renamed or dropped but a new
variable with the stored name is recreated with identical values on
{cmd:e(sample)}, {cmd:predict, xb} and {cmd:predict, fitted} remain available,
because this treatment-only linear predictor depends on the current
values of the stored-name treatment regressor rather than on variable
identity.

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
{cmd:r(neg_sd)}, and {cmd:r(neg_ratio)}.  The positive-weight moment entries may be missing when the stored weight matrix contains no positive weights.  The negative-weight moment entries may be missing when the stored weight matrix contains no negative weights.  If {cmd:r(n_nonzero) = 0}, {cmd:r(neg_ratio)} is reported as missing because the share of negative weights among nonzero weights is undefined on an empty nonzero-weight support.  If {cmd:r(n_negative) = 0} but {cmd:r(n_nonzero) > 0}, {cmd:r(neg_ratio)} is still reported and displayed as {cmd:0}.  Here {cmd:r(neg_ratio)} is the share of negative weights among the nonzero weights, that is,
{cmd:r(n_negative) / r(n_nonzero)}.  If the stored weight matrix
{cmd:e(W)} contains missing values, {cmd:estat wfe_weights} stops with
{cmd:estat wfe_weights: e(W) must not contain missing values}
rather than reporting undefined moments.

{pstd}
Under {cmd:pwfe}, {cmd:predict, residuals} returns residuals from the
observed-outcome scale, namely {cmd:outcome - treat * coef}.  Therefore
{cmd:xb + residuals} recovers the current variable named in
{cmd:outcome()} on the estimation sample rather than the stored
{cmd:e(y_star)} vector.  Re-sorting the data does not change this
residual identity.  Changing {cmd:unit()} or {cmd:time()} values alone
does not alter {cmd:predict, residuals}; the current implementation
recomputes outcome-scale residuals directly from the current data rather
than replaying stored transformed outcomes.  Renaming or dropping the
current {cmd:unit()} or {cmd:time()} variables alone does not alter
{cmd:predict, residuals}.  Only the current outcome
and treatment values themselves must remain unchanged.  If the current outcome values are edited in place after estimation, {cmd:predict, residuals} stops with {cmd:wfe_predict: current outcome/treat variables no longer match the stored estimation sample}.  If the current treatment values are edited in place after estimation, {cmd:predict, residuals} stops with {cmd:wfe_predict: current outcome/treat variables no longer match the stored estimation sample}.
Changing only the numeric storage type of the current outcome or treatment variable while leaving their {cmd:e(sample)} values unchanged does not alter {cmd:predict, residuals}.
If the visible macro {cmd:e(depvar)} is missing or blank, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored dependent variable name is unavailable}.
If the treatment variable named in {cmd:e(_treatvar)} is renamed, dropped, or otherwise unavailable after estimation, {cmd:predict, residuals} stops with {cmd:wfe_predict: stored treatment variable is no longer available in current data}.
This exact-name guard is narrower than the estimation-sample mismatch
check because the residual path must first verify that the stored
{cmd:treat()} variable itself still exists in the current data.
If the current outcome variable named in {cmd:e(depvar)} is renamed, dropped, or otherwise unavailable after estimation, {cmd:predict, residuals} stops with
{cmd:wfe_predict: stored dependent variable is no longer available in current data}
rather than replaying a stored residual vector.

{marker examples}{...}
{title:Examples}

{pstd}
The examples below first generate a small balanced panel in memory
and then walk through progressively richer uses of {cmd:pwfe}.
In logit-based examples, treatment variation only needs to
appear in the overall estimation sample.  Under the current
penalized propensity-score fit, split unit/time pools may
still be single-support under the penalized fit, so
{cmd:method(unit)} and {cmd:method(time)} do not require every pool to
contain both treated and control observations.

{phang}
{cmd:. clear}
{cmd:. set seed 12345}
{cmd:. local N = 10}
{cmd:. local TT = 15}
{cmd:. set obs `=`N' * `TT''}
{cmd:. gen int id = ceil(_n / `TT')}
{cmd:. gen int t = mod(_n - 1, `TT') + 1}
{cmd:. gen byte d = rbinomial(1, 0.25)}
{cmd:. bysort id: egen double __unit_mean = mean(d)}
{cmd:. bysort t: egen double __time_mean = mean(d)}
{cmd:. count if inlist(__unit_mean, 0, 1) | inlist(__time_mean, 0, 1)}
{cmd:. while r(N) > 0 {c -(}}
{cmd:>     replace d = rbinomial(1, 0.25)}
{cmd:>     drop __unit_mean __time_mean}
{cmd:>     bysort id: egen double __unit_mean = mean(d)}
{cmd:>     bysort t: egen double __time_mean = mean(d)}
{cmd:>     count if inlist(__unit_mean, 0, 1) | inlist(__time_mean, 0, 1)}
{cmd:> {c )-}}
{cmd:. drop __unit_mean __time_mean}
{cmd:. gen double x1 = rnormal(0.5, 1)}
{cmd:. gen double x2 = rbeta(5, 1)}
{cmd:. bysort id: gen double alpha_i = cond(_n == 1, rnormal(), .)}
{cmd:. bysort id: replace alpha_i = alpha_i[1]}
{cmd:. gen double y = alpha_i + d + x1 + x2 + rnormal()}
{cmd:. gen double ps = runiform(0.05, 0.95)}

{pstd}
Example 1: basic propensity-score weighted FE with unit effects.

{phang}
{cmd:. pwfe x1 x2, treat(d) outcome(y) unit(id)}

{pstd}
Example 2: supply an explicit time variable.

{phang}
{cmd:. pwfe x1 x2, treat(d) outcome(y) unit(id) time(t)}

{pstd}
Example 3: time-level propensity-score estimation.

{phang}
{cmd:. pwfe x1 x2, treat(d) outcome(y) unit(id) time(t) method(time)}

{pstd}
Example 4: first-difference estimation.

{phang}
{cmd:. pwfe x1 x2, treat(d) outcome(y) unit(id) time(t) estimator(fd)}

{pstd}
Example 5: user-supplied propensity scores.

{phang}
{cmd:. pwfe, treat(d) outcome(y) unit(id) time(t) pscore(ps)}

{pstd}
Example 6: homoskedastic standard errors.

{phang}
{cmd:. pwfe x1 x2, treat(d) outcome(y) unit(id) hetero_se(off) auto_se(off)}

{pstd}
Example 7: postestimation.

{phang}
{cmd:. predict double xbhat, xb}

{phang}
{cmd:. estat wfe_weights}

{marker references}{...}
{title:References}

{marker IK2021}{...}
{phang}
Imai, K., and I. S. Kim. 2021. "On the Use of Two-Way Fixed Effects
Regression Models for Causal Inference with Panel Data."
{it:Political Analysis} 29(3): 405-415.

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
Manual:  {helpb wfe}, {helpb xtreg}, {helpb logit}

{psee}
Postestimation:  {helpb predict}, {helpb estat}
