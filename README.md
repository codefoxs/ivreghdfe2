# ivreghdfe2

`ivreghdfe` with a corrected second-stage VCE, a reported constant term, extra IV
diagnostics, and both stages saved as estimation results.

## What this fixes

[`ivreghdfef`](https://github.com/codefoxs/ivreghdfef) produced a constant term
and two stored stages by running 2SLS by hand:

```stata
reghdfe x  iv controls, absorb(...)   // first stage
predict xhat
reghdfe y  xhat controls, absorb(...) // second stage
```

The coefficients that returns are correct — 2SLS slopes really do equal an OLS
regression on the fitted values. **The second-stage standard errors are not.** A
manual second stage builds its residuals as `y - xhat*b`, while 2SLS requires
`y - x*b`. Under homoskedasticity the two differ by a scalar; under `robust` or
`cluster()` the meat of the sandwich is structurally different and no rescaling
repairs it.

Measured on the example below, `ivreghdfef` understates the standard error on the
endogenous regressor by **10.7%**, which turns a *p* = 0.006 into a *p* = 0.014.

`ivreghdfe2` takes its slopes and variance matrix from `ivreghdfe` itself, so
they are correct by construction, and derives the constant term without
disturbing them.

## Install

```stata
* Latest version
cap ado uninstall ivreghdfe2
net install ivreghdfe2, from("https://raw.githubusercontent.com/codefoxs/ivreghdfe2/main/") replace
```

Requires [`ivreghdfe`](https://github.com/sergiocorreia/ivreghdfe),
[`ivreg2`](http://fmwww.bc.edu/RePEc/bocode/i/) and
[`reghdfe`](https://github.com/sergiocorreia/reghdfe).

## Syntax

The model is written the way `ivreghdfe` wants it, not the way `ivreghdfef` did:

```stata
ivreghdfe2 depvar [exogvars] (endog = instruments) [if] [in] [weight] ,
    absorb(absvars) [first store(name) format(%fmt) nosummary  ...ivreghdfe options]
```

Options such as `cluster()`, `robust`, `small`, `gmm2s`, `liml`, `bw()`,
`kernel()` and `tolerance()` are written directly and passed straight through —
no more `iv()` / `opt()` / `regopt()` wrappers.

## Example

```stata
webuse nlswork

ivreghdfe2 ln_wage age tenure not_smsa south (ttl_exp = hours union), ///
    absorb(idcode year) cluster(idcode) first store(IV)

esttab IV_first IV_second, se scalars(N r2_a cdf rkf rklm hansenj)
```

`store(IV)` leaves behind `IV_first` and `IV_second`, both carrying a `_cons`
row. With several endogenous regressors the first stages become `IV_first1`,
`IV_first2`, and so on.

## How it works

`absorb()` makes `ivreghdfe` partial the constant out (`ivreghdfe.ado` lines
674–681 force `noconstant`), which is why there is no `_cons` row to report in
the first place. `ivreghdfe2` recovers one like this:

1. **Reference run.** `ivreghdfe` on the model as typed. Its `e(sample)`,
   diagnostics, slopes and variance matrix are the authority for everything else.
2. **First stage.** `reghdfe` of each endogenous regressor on the instruments and
   controls. Ordinary OLS, so its standard errors need no correction, and it
   reports a constant of its own.
3. **Partialling out.** Every variable is residualised on the absorbed fixed
   effects, then has its sample mean added back. Adding the means back is what
   re-identifies the intercept; by Frisch–Waugh–Lovell it leaves the slopes and
   residuals untouched. This is the same "regression through the means"
   convention `reghdfe` uses (`reghdfe.mata` lines 3830–3870).
4. **Augmented second stage.** `ivreg2` on the transformed variables *with* a
   constant. Every variance variant — `robust`, one- and two-way `cluster()`,
   HAC kernels, GMM — is therefore `ivreg2`'s own code, including the `_cons`
   row. No variance formula is reimplemented.
5. **Calibration and verification.** The augmented run is made *without* `small`,
   so `ivreg2` leaves its variance matrix unscaled; the scalar small-sample
   factor is then recovered as the ratio of reference variances to augmented
   ones. Since `ivreg2`'s small-sample correction multiplies the whole matrix by
   one number, that ratio must be identical for every coefficient — the command
   checks this and aborts if it is not. This avoids having to reproduce
   `ivreghdfe`'s internal `sdofminus` bookkeeping, and doubles as a check that
   the augmented run really reproduced the reference fit. Results are left in
   `e(chk_b)`, `e(chk_v)` and `e(smallfac)`.
6. **Reporting.** The slope block is copied verbatim from the reference run, so
   reported slopes and their standard errors equal `ivreghdfe`'s to the last
   digit. Only the `_cons` row is new.

## Correctness evidence

`tests/test_ivreghdfe2.do` runs 68 checks on `nlswork` (68/68 passing). The ones
that matter:

- **Against `ivreghdfe`** across nine variance specifications (unadjusted,
  `robust`, one- and two-way `cluster()`, each with and without `small`): the
  command's own internal comparison of the *raw* augmented estimates against the
  reference agrees to `e(chk_b)` ≈ 4e-13 on coefficients, and the variance ratio
  is constant to between 1e-11 and 7e-9.
- **The constant term's value** equals `ybar - Xbar'b` to 6e-12.
- **The constant term's variance.** Using an endogenous regressor as its own
  instrument makes 2SLS numerically identical to OLS, so the full `e(b)` and
  `e(V)` — `_cons` row included — must equal `reghdfe`'s. They match to 4e-13
  (unadjusted), 1.7e-12 (`robust`) and 1.5e-12 (`cluster()`).
- **Nested clustering.** Clustering on a panel id that is also absorbed is the
  one case where `reghdfe` and `ivreghdfe` genuinely disagree: `reghdfe` subtracts
  one extra degree of freedom, so the two differ by a single scalar of order
  1/N. All six coefficients, `_cons` included, share that same scalar
  (1.00002808633 vs 1.00002808633) to 12 digits — evidence that the `_cons` row
  tracks the slopes rather than being computed independently.
- Weights (`aweight`, `fweight`, `pweight`), `if`/`in`, time-series operators,
  several endogenous regressors, `absorb()` omitted, `esttab` and `reg2docx`
  export, and replay.

## Limitations

- **No factor variables in the variable list** (`i.x`, `i.a#i.b`, `c.x`). Keeping
  coefficient names aligned with the design matrix through the augmented
  regression would require tracking base and omitted levels. Put categorical
  variables in `absorb()`. Time-series operators (`L.x`, `D.x`) are supported.
- **No `partial()`.** Combined with `absorb()`, its small-sample degrees-of-freedom
  arithmetic depends on a quantity internal to `ivreg2`.
- **Slower than `ivreghdfe`**: one reference run plus one partialling-out pass per
  variable, so roughly 2–3× the run time on a large panel.
- **First stage under `bw()`/`kernel()`** falls back to `vce(robust)`, because
  `reghdfe` cannot do HAC. The second stage is unaffected.
- **No `iweight`** — `reghdfe` does not accept it.
- `_cons` follows `reghdfe`'s convention (the mean of the absorbed fixed
  effects). With fixed effects it is not a structural parameter; it is reported
  because tables want a row for it.

## Reference

- https://github.com/sergiocorreia/ivreghdfe
- https://github.com/sergiocorreia/reghdfe
- https://github.com/codefoxs/ivreghdfef
