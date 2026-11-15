{smcl}
{* *! version 1.0.0 15Nov2026}{...}
{vieweralsosee "[R] reghdfe" "help reghdfe"}{...}
{vieweralsosee "[R] ivreg2" "help ivreg2"}{...}
{vieweralsosee "[R] ivreghdfe" "help ivreghdfe"}{...}
{vieweralsosee "" "--"}{...}
{viewerjumpto "Syntax" "ivreghdfe2##syntax"}{...}
{viewerjumpto "Description" "ivreghdfe2##description"}{...}
{viewerjumpto "Options" "ivreghdfe2##options"}{...}
{viewerjumpto "How it works" "ivreghdfe2##method"}{...}
{viewerjumpto "The constant term" "ivreghdfe2##cons"}{...}
{viewerjumpto "Limitations" "ivreghdfe2##limits"}{...}
{viewerjumpto "Examples" "ivreghdfe2##examples"}{...}
{viewerjumpto "Stored results" "ivreghdfe2##results"}{...}
{title:Title}

{p2colset 5 22 24 2}{...}
{p2col :{hi:ivreghdfe2} {hline 2}}ivreghdfe with a corrected second-stage VCE, a reported constant term, extra IV diagnostics, and both stages saved{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 16 2}
{cmd:ivreghdfe2} {depvar} [{it:exogvars}] {cmd:(}{it:endog} {cmd:=} {it:instruments}{cmd:)}
{ifin} {weight} {cmd:,} {opth a:bsorb(reghdfe##absvar:absvars)}
[{it:options}]

{pstd}
The model is written the same way as for {help ivreghdfe}. Exogenous control
variables may be listed before or after the parentheses. {it:weight} may be
{opt aweight}, {opt fweight} or {opt pweight}.

{synoptset 30 tabbed}{...}
{synopthdr}
{synoptline}
{syntab :Model}
{synopt :{opth a:bsorb(reghdfe##absvar:absvars)}}categorical variables to absorb
as fixed effects. May be omitted, in which case {cmd:ivreghdfe} already reports
a constant and results are passed through unchanged.{p_end}

{syntab :Reported output}
{synopt :{opt first}}also display the first stage regression table{p_end}
{synopt :{opt sto:re(name)}}save the estimates: {it:name}{cmd:_second}, plus
{it:name}{cmd:_first} (one endogenous regressor) or {it:name}{cmd:_first1},
{it:name}{cmd:_first2}, ... (several){p_end}
{synopt :{opt f:ormat(%fmt)}}display format for the summary block; default
{cmd:%10.4f}{p_end}
{synopt :{opt nosum:mary}}suppress the summary block{p_end}

{syntab :Everything else}
{synopt :{it:options}}any other {help ivreghdfe} or {help reghdfe} option, for
instance {opt cluster()}, {opt robust}, {opt small}, {opt gmm2s}, {opt liml},
{opt bw()}, {opt kernel()}, {opt tolerance()}. These are passed straight
through. {opt vce(robust)} and {opt vce(cluster }{it:varlist}{opt )} are
accepted and translated.{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
{cmd:ivreghdfe2} estimates a 2SLS model with high-dimensional fixed effects and
adds three things that {help ivreghdfe} does not provide:

{phang2}1. a reported {cmd:_cons} row, with a valid standard error;{p_end}
{phang2}2. both the first and the second stage saved as estimation results, so
that they can go straight into {help esttab} or {cmd:reg2docx};{p_end}
{phang2}3. a summary block with the Cragg-Donald Wald F, Kleibergen-Paap rk Wald
F, Kleibergen-Paap rk LM and Hansen J statistics.{p_end}

{pstd}
{bf:Why not just build the two stages by hand?} The usual way to get those three
things is to run {cmd:reghdfe} for the first stage, {cmd:predict} the fitted
values, then {cmd:reghdfe} the outcome on those fitted values. That recovers the
2SLS coefficients exactly, but {bf:its standard errors are wrong}. A manual second
stage forms its residuals as {it:y} - {it:xhat}*b, whereas 2SLS requires
{it:y} - {it:x}*b. Under homoskedasticity the two differ only by a scalar, but
under {opt robust} or {opt cluster()} the meat of the sandwich is structurally
different and no rescaling can repair it. On the worked example below, the manual
route understates the standard error on the endogenous regressor by about 11%,
which moves its p-value from 0.006 to 0.014.

{pstd}
{cmd:ivreghdfe2} instead takes its slopes and its variance matrix from
{cmd:ivreghdfe} itself, so they are correct by construction, and derives the
constant term without disturbing them.


{marker method}{...}
{title:How it works}

{pstd}
{cmd:absorb()} makes {cmd:ivreghdfe} partial the constant out (it forces
{opt noconstant} internally), which is why no {cmd:_cons} row is available to
begin with. {cmd:ivreghdfe2} recovers one as follows.

{phang2}
{bf:1. Reference run.} {cmd:ivreghdfe} is run on the model as typed. Its
{cmd:e(sample)}, its diagnostics and its slope estimates and variance matrix are
kept as the authority for everything that follows.{p_end}

{phang2}
{bf:2. First stage.} For each endogenous regressor, {cmd:reghdfe} is run on the
instruments and the exogenous controls. This is ordinary OLS, so its standard
errors need no correction, and it reports a constant of its own.{p_end}

{phang2}
{bf:3. Partialling out.} Each variable in the model is residualised on the
absorbed fixed effects with {cmd:reghdfe}, and its sample mean is then added
back. Adding the means back is what makes the intercept identified again; by the
Frisch-Waugh-Lovell theorem it leaves the slopes and the residuals untouched.{p_end}

{phang2}
{bf:4. Augmented second stage.} {cmd:ivreg2} is run on the transformed variables
{it:with} a constant. All of the variance machinery -- {opt robust}, one- and
two-way {opt cluster()}, HAC kernels, GMM -- is therefore {cmd:ivreg2}'s own,
including the {cmd:_cons} row. No variance formula is reimplemented here.{p_end}

{phang2}
{bf:5. Calibration and verification.} The augmented run is deliberately made
{it:without} {opt small}, so {cmd:ivreg2} leaves its variance matrix unscaled.
The scalar small-sample factor is then recovered as the ratio of the reference
variances to the augmented ones. Because {cmd:ivreg2}'s small-sample correction
multiplies the whole matrix by one number, that ratio must be the same for every
coefficient; {cmd:ivreghdfe2} checks that it is and aborts if it is not. This
sidesteps having to reproduce {cmd:ivreghdfe}'s internal degrees-of-freedom
bookkeeping, and doubles as a check that the augmented run really reproduced the
reference fit. The outcome is left in {cmd:e(chk_b)}, {cmd:e(chk_v)} and
{cmd:e(smallfac)}.{p_end}

{phang2}
{bf:6. Reporting.} The slope block is taken verbatim from the reference run, so
reported slopes and their standard errors are identical to {cmd:ivreghdfe}'s to
the last digit. Only the {cmd:_cons} row is new.{p_end}


{marker cons}{...}
{title:The constant term}

{pstd}
{cmd:_cons} is reported on {help reghdfe}'s convention: it equals
{it:ybar} - {it:Xbar}'{it:b} over the estimation sample, that is, the mean of the
absorbed fixed effects. Its variance comes from the same augmented regression as
the slopes, so it is internally consistent with them.

{pstd}
In a model with fixed effects the intercept is not separately identified from the
fixed effects, so {cmd:_cons} is {bf:not a structural parameter}. It is reported
because tables usually want a row for it, not because it should be interpreted.


{marker limits}{...}
{title:Limitations}

{phang2}
{bf:No factor variables in the variable list.} {cmd:i.}{it:var},
{cmd:i.}{it:a}{cmd:#i.}{it:b} and {cmd:c.}{it:var} are refused. Expanding them
would require tracking base and omitted levels through the augmented regression
in order to keep coefficient names aligned with the design matrix. Put
categorical variables in {opt absorb()} instead. Time-series operators
({cmd:L.}{it:x}, {cmd:D.}{it:x}) {it:are} supported.{p_end}

{phang2}
{bf:No} {opt partial()}{bf:.} Combined with {opt absorb()} its small-sample
degrees-of-freedom arithmetic depends on a quantity internal to {cmd:ivreg2}
that cannot be read from outside.{p_end}

{phang2}
{bf:Slower than} {cmd:ivreghdfe}. One reference run plus one partialling-out pass
per variable, so expect roughly two to three times the run time on a large
panel.{p_end}

{phang2}
{bf:First stage under} {opt bw()}{bf:/}{opt kernel()}{bf:.} {cmd:reghdfe} cannot
compute HAC standard errors, so the reported first stage falls back to
{opt vce(robust)} and says so. The second stage is unaffected.{p_end}

{phang2}
{bf:No} {opt iweight}{bf:.} {cmd:reghdfe} does not accept it.{p_end}

{phang2}
{bf:t versus z.} Inference follows {cmd:ivreghdfe}. With {opt absorb()},
{cmd:ivreghdfe} always applies the small-sample correction, so t statistics are
reported against {cmd:e(df_r)} whether or not {opt small} is typed. Without
{opt absorb()}, it behaves like {cmd:ivreg2}: t statistics with {opt small},
z statistics otherwise.{p_end}


{marker examples}{...}
{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse nlswork}{p_end}

{pstd}Basic use{p_end}
{phang2}{cmd:. ivreghdfe2 ln_wage age tenure not_smsa south (ttl_exp = hours union), absorb(idcode year) cluster(idcode)}{p_end}

{pstd}Show the first stage and save both stages{p_end}
{phang2}{cmd:. ivreghdfe2 ln_wage age tenure not_smsa south (ttl_exp = hours union), absorb(idcode year) cluster(idcode) first store(IV)}{p_end}

{pstd}Export a two-stage table{p_end}
{phang2}{cmd:. esttab IV_first IV_second, se scalars(N r2_a cdf rkf rklm hansenj)}{p_end}

{pstd}Or with {cmd:reg2docx}{p_end}
{phang2}{cmd:. reg2docx IV_first IV_second using "IV.docx", replace}{p_end}
{phang2}{cmd:      b(%20.4f) t(%20.4f)}{p_end}
{phang2}{cmd:      scalars(N(%20.0fc) r2_a(%20.4f) cdf(%20.4f) rkf(%20.4f) rklm(%20.4f) hansenj(%20.4f))}{p_end}
{phang2}{cmd:      addfe("Id=Yes" "Year=Yes") mtitles("First stage" "Second stage")}{p_end}

{pstd}Several endogenous regressors (first stages become {cmd:IV2_first1}, {cmd:IV2_first2}){p_end}
{phang2}{cmd:. ivreghdfe2 ln_wage not_smsa south (ttl_exp tenure = hours union msp nev_mar), absorb(idcode year) cluster(idcode) store(IV2)}{p_end}

{pstd}Set the display format of the summary block{p_end}
{phang2}{cmd:. ivreghdfe2 ln_wage age tenure (ttl_exp = hours union), absorb(idcode year) cluster(idcode) f(%9.3f)}{p_end}

{pstd}Confirm against {cmd:ivreghdfe} -- the slopes and standard errors match exactly{p_end}
{phang2}{cmd:. ivreghdfe ln_wage age tenure not_smsa south (ttl_exp = hours union), absorb(idcode year) cluster(idcode)}{p_end}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:ivreghdfe2} keeps everything {cmd:ivreghdfe} stores in {cmd:e()} and adds
the following.

{synoptset 22 tabbed}{...}
{p2col 5 22 26 2: Scalars}{p_end}
{synopt:{cmd:e(cdf)}}Cragg-Donald Wald F statistic{p_end}
{synopt:{cmd:e(rkf)}}Kleibergen-Paap rk Wald F statistic{p_end}
{synopt:{cmd:e(rklm)}}Kleibergen-Paap rk LM statistic (= {cmd:e(idstat)}){p_end}
{synopt:{cmd:e(rklmp)}}its p-value (= {cmd:e(idp)}){p_end}
{synopt:{cmd:e(hansenj)}}Hansen J statistic (= {cmd:e(j)}){p_end}
{synopt:{cmd:e(hansenjp)}}its p-value (= {cmd:e(jp)}){p_end}
{synopt:{cmd:e(smallfac)}}recovered scalar small-sample factor{p_end}
{synopt:{cmd:e(chk_b)}}worst relative coefficient difference against the reference run{p_end}
{synopt:{cmd:e(chk_v)}}spread of the variance ratio across coefficients{p_end}
{synopt:{cmd:e(df_r_ref)}}residual degrees of freedom from the reference run{p_end}

{p2col 5 22 26 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:ivreghdfe2}{p_end}
{synopt:{cmd:e(endovars)}}endogenous regressors{p_end}
{synopt:{cmd:e(exexvars)}}excluded instruments{p_end}
{synopt:{cmd:e(inexogvars)}}included exogenous regressors{p_end}
{synopt:{cmd:e(absvarlist)}}contents of {opt absorb()}{p_end}
{synopt:{cmd:e(fsstored)}}names of the stored first-stage estimates{p_end}
{synopt:{cmd:e(ssstored)}}name of the stored second-stage estimate{p_end}

{p2col 5 22 26 2: Matrices}{p_end}
{synopt:{cmd:e(b)}, {cmd:e(V)}}coefficients and variance matrix, including {cmd:_cons}{p_end}
{synopt:{cmd:e(fsb)}}first-stage coefficients, endogenous regressors by excluded instruments{p_end}
{synopt:{cmd:e(fsp)}}the corresponding p-values{p_end}
{p2colreset}{...}

{pstd}
{cmd:e(rklm)}, {cmd:e(rklmp)}, {cmd:e(hansenj)} and {cmd:e(hansenjp)} are short
aliases for {cmd:e(idstat)}, {cmd:e(idp)}, {cmd:e(j)} and {cmd:e(jp)}, so that
the four diagnostics can be named directly in a {cmd:scalars()} list. The two
p-values are taken from {cmd:ivreg2} rather than recomputed, which keeps the
degrees of freedom right when there is more than one endogenous regressor: the
Kleibergen-Paap rk LM statistic has L-K+1 of them and the Hansen J statistic
L-K, where L counts the excluded instruments and K the endogenous regressors.


{title:Author}

{pstd}公众号：凯恩斯学计量{p_end}


{title:Also see}

{psee}
Online: {help ivreghdfe}, {help ivreg2}, {help reghdfe}
{p_end}
