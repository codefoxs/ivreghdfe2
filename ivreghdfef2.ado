*! Version 0.1.0 30Jul2026
*! ivreghdfef2: ivreghdfe with a corrected second-stage VCE, a reported constant
*!              term, extra IV diagnostics, and both stages saved as estimates.

program define ivreghdfef2, eclass
    version 14

    if replay() {
        if `"`e(cmd)'"' != "ivreghdfef2" error 301
        Display
        exit
    }

    Estimate `0'
    Display
end


* ---------------------------------------------------------------------------
* Estimate
* ---------------------------------------------------------------------------
program define Estimate, eclass
    local cmdline `"ivreghdfef2 `0'"'

    syntax anything(equalok) [if] [in] [aw fw pw/] [,   ///
        Absorb(string)                                  ///
        FIRST                                           ///
        STOre(name)                                     ///
        Format(string)                                  ///
        NOSUMmary                                       ///
        CLuster(string)                                 ///
        Robust                                          ///
        SMall                                           ///
        VCE(string)                                     ///
        BW(string)                                      ///
        KERnel(string)                                  ///
        PARTIAL(string)                                 ///
        NOCONstant                                      ///
        TOLerance(string)                               ///
        ITERate(string)                                 ///
        TECHnique(string)                               ///
        TRAnsform(string)                               ///
        ACCELeration(string)                            ///
        PREConditioner(string)                          ///
        POOLsize(string)                                ///
        DOFadjustments(string)                          ///
        GROUPVar(string)                                ///
        RESiduals(string)                               ///
        VERBose(string)                                 ///
        COMPACT                                         ///
        PARallel(string)                                ///
        noWARN                                          ///
        * ]

    * ---------------------------------------------------------------- guards
    if `"`partial'"' != "" {
        di as err "option partial() is not supported by ivreghdfef2"
        di as err "  the small-sample dof arithmetic for partial() combined with"
        di as err "  absorb() cannot be reproduced from outside ivreg2"
        exit 198
    }

    if `"`vce'"' != "" {
        if `"`cluster'`robust'"' != "" {
            di as err "vce() may not be combined with robust or cluster()"
            exit 198
        }
        gettoken vcetype vcerest : vce
        local vcetype = lower(`"`vcetype'"')
        if "`vcetype'" == "robust" {
            local robust robust
        }
        else if inlist("`vcetype'", "cl", "clu", "clus", "clust", "cluste", "cluster") {
            local cluster `vcerest'
        }
        else if inlist("`vcetype'", "unadjusted", "ols") {
            * nothing to do
        }
        else {
            di as err "vce(`vce') is not supported; use robust, cluster(), or unadjusted"
            exit 198
        }
    }

    if `"`format'"' == "" local format "%10.4f"
    confirm format `format'

    * ------------------------------------------------------- parse the model
    local input `"`anything'"'
    local p1 = strpos(`"`input'"', "(")
    local p2 = strpos(`"`input'"', ")")
    if `p1' == 0 | `p2' == 0 | `p2' < `p1' {
        di as err "you must specify the endogenous regressors and their instruments"
        di as err "  as (endog = instruments)"
        exit 198
    }

    local pre  = trim(substr(`"`input'"', 1, `p1' - 1))
    local ivin = trim(substr(`"`input'"', `p1' + 1, `p2' - `p1' - 1))
    local post = trim(substr(`"`input'"', `p2' + 1, .))

    local peq = strpos(`"`ivin'"', "=")
    if `peq' == 0 {
        di as err "(endog = instruments) requires an equals sign"
        exit 198
    }
    local endo = trim(substr(`"`ivin'"', 1, `peq' - 1))
    local exex = trim(substr(`"`ivin'"', `peq' + 1, .))

    gettoken depvar preexog : pre
    local inexog = trim(itrim(`"`preexog' `post'"'))

    if "`depvar'" == "" {
        di as err "no dependent variable specified"
        exit 198
    }
    if "`endo'" == "" {
        di as err "no endogenous regressor specified"
        exit 198
    }
    if "`exex'" == "" {
        di as err "no excluded instrument specified"
        exit 198
    }

    * factor variables are not supported -- see the help file for why
    foreach v in `depvar' `endo' `exex' `inexog' {
        local lv = lower("`v'")
        if strpos("`v'", "#")                       ///
        |  regexm("`lv'", "^i[bn0-9()/., ]*\.")     ///
        |  regexm("`lv'", "^[0-9]+[bno]*\.")        ///
        |  regexm("`lv'", "^c\.") {
            di as err "ivreghdfef2 does not support factor variables ({bf:`v'})"
            di as err "  put categorical variables in absorb() instead"
            exit 198
        }
    }

    * resolve time-series operators into plain variables (1:1, order preserved)
    foreach grp in depvar endo exex inexog {
        local orig `"``grp''"'
        if `"`orig'"' == "" continue
        tsrevar `orig'
        local `grp'_r `r(varlist)'
    }

    marksample touse
    markout `touse' `depvar_r' `endo_r' `exex_r' `inexog_r'

    * -------------------------------------------------------- option bundles
    local wt
    local swt
    if "`weight'" != "" {
        local wt "[`weight'=`exp']"
        local sweight = cond("`weight'" == "pweight", "aweight", "`weight'")
        local swt "[`sweight'=`exp']"
    }

    local absopt
    if `"`absorb'"' != "" local absopt "absorb(`absorb')"

    * vce options in ivreg2 spelling ...
    local vceopt
    if "`robust'"    != "" local vceopt `vceopt' robust
    if `"`cluster'"' != "" local vceopt `vceopt' cluster(`cluster')

    * ... and the same thing in reghdfe spelling, which has no bare -robust-
    local hdfevce
    if `"`cluster'"' != ""      local hdfevce vce(cluster `cluster')
    else if "`robust'" != ""    local hdfevce vce(robust)
    else if `"`bw'`kernel'"' != "" {
        * reghdfe cannot do HAC; the first stage falls back to robust
        local hdfevce vce(robust)
        di as txt "note: reghdfe does not support bw()/kernel(); the reported first stage"
        di as txt "      uses vce(robust). The second stage is unaffected."
    }

    * hdfe options that are safe to repeat on every partialling-out call.
    * a tighter default tolerance keeps the augmented run and the reference run
    * in agreement well below display precision.
    if `"`tolerance'"' == "" local tolerance 1e-10
    local hdfeopts tolerance(`tolerance')
    foreach o in iterate technique transform acceleration preconditioner poolsize verbose {
        if `"``o''"' != "" local hdfeopts `hdfeopts' `o'(``o'')
    }
    if "`compact'" != "" local hdfeopts `hdfeopts' compact
    if "`warn'"    != "" local hdfeopts `hdfeopts' nowarn

    * what the reference ivreghdfe call sees
    local refopts `vceopt' `hdfeopts' `options'
    if "`small'"      != "" local refopts `refopts' small
    if "`noconstant'" != "" local refopts `refopts' noconstant
    if `"`bw'"'       != "" local refopts `refopts' bw(`bw')
    if `"`kernel'"'   != "" local refopts `refopts' kernel(`kernel')
    if `"`dofadjustments'"' != "" local refopts `refopts' dofadjustments(`dofadjustments')
    if `"`groupvar'"'  != "" local refopts `refopts' groupvar(`groupvar')
    if `"`residuals'"' != "" local refopts `refopts' residuals(`residuals')

    * what the augmented ivreg2 call sees -- deliberately WITHOUT small, so that
    * ivreg2 leaves V unscaled.  The scalar small-sample factor is then recovered
    * from the reference run, which avoids having to reproduce ivreghdfe's
    * internal sdofminus bookkeeping.
    local augopts `vceopt' `options'
    if "`noconstant'" != "" local augopts `augopts' noconstant
    if `"`bw'"'       != "" local augopts `augopts' bw(`bw')
    if `"`kernel'"'   != "" local augopts `augopts' kernel(`kernel')

    * ----------------------------------------------------- reference IV run
    tempvar esample
    qui ivreghdfe `input' if `touse' `wt', `absopt' `refopts' nooutput
    qui gen byte `esample' = e(sample)

    tempname bref Vref
    matrix `bref' = e(b)
    matrix `Vref' = e(V)

    local ref_df_r = e(df_r)
    local ref_N    = e(N)

    * harvest every e() result so it can be restored after ereturn post
    local drop_s "df_r"
    local drop_m "cmd cmdline depvar properties predict"
    local drop_x "b V"
    local eslist : e(scalars)
    local emlist : e(macros)
    local exlist : e(matrices)
    local eslist : list eslist - drop_s
    local emlist : list emlist - drop_m
    local exlist : list exlist - drop_x

    foreach s of local eslist {
        tempname S_`s'
        scalar `S_`s'' = e(`s')
    }
    foreach m of local emlist {
        local M_`m' `"`e(`m')'"'
    }
    foreach x of local exlist {
        tempname X_`x'
        matrix `X_`x'' = e(`x')
    }

    * ------------------------------------------------------- first stage(s)
    local qui1 = cond("`first'" == "", "qui", "")
    local nendo : word count `endo'
    local nexex : word count `exex'

    tempname fsb fsp rt
    matrix `fsb' = J(`nendo', `nexex', .)
    matrix `fsp' = J(`nendo', `nexex', .)

    local fsstored
    forvalues j = 1/`nendo' {
        local ev : word `j' of `endo'
        if "`first'" != "" {
            di as txt _n "{hline 78}"
            di as res "First stage: `ev'"
            di as txt "{hline 78}"
        }
        `qui1' reghdfe `ev' `exex' `inexog' if `esample' `wt', `absopt' `hdfevce' `hdfeopts'

        matrix `rt' = r(table)
        forvalues i = 1/`nexex' {
            local zn : word `i' of `exex'
            local c = colnumb(`rt', "`zn'")
            if `c' < . {
                matrix `fsb'[`j', `i'] = `rt'[1, `c']
                matrix `fsp'[`j', `i'] = `rt'[4, `c']
            }
        }

        if "`store'" != "" {
            local fsname = cond(`nendo' == 1, "`store'_first", "`store'_first`j'")
            capture estimates drop `fsname'
            estimates store `fsname'
            local fsstored `fsstored' `fsname'
        }
    }
    * cosmetic only -- names carrying time-series operators can be rejected here
    capture matrix rownames `fsb' = `endo'
    capture matrix colnames `fsb' = `exex'
    capture matrix rownames `fsp' = `endo'
    capture matrix colnames `fsp' = `exex'

    * ------------------------------------------- augmented second stage run
    tempname b V
    local smallfac 1
    local worstb 0
    local spread 0

    if `"`absorb'"' == "" {
        * without fixed effects ivreghdfe already reports the constant term
        matrix `b' = `bref'
        matrix `V' = `Vref'
    }
    else {
        * Partial out the fixed effects, then add the sample means back so that
        * the intercept is identified again.  This is exactly reghdfe's
        * "regression through the means" convention, and by FWL it leaves the
        * slopes and the residuals untouched.
        local dm_depvar
        local dm_endo
        local dm_inexog
        local dm_exex
        local ii 0
        local tmpall
        local dispall
        foreach grp in depvar endo inexog exex {
            local nv : word count ``grp'_r'
            forvalues t = 1/`nv' {
                local v  : word `t' of ``grp'_r'
                local vn : word `t' of ``grp''
                local ++ii
                tempvar d`ii'
                qui reghdfe `v' if `esample' `wt', `absopt' residuals(`d`ii'') `hdfeopts'
                if "`noconstant'" == "" {
                    qui summarize `v' if `esample' `swt', meanonly
                    qui replace `d`ii'' = `d`ii'' + r(mean) if `esample'
                }
                qui count if `esample' & missing(`d`ii'')
                if r(N) > 0 {
                    di as err "partialling out `vn' left `r(N)' missing values inside e(sample)"
                    di as err "  the estimation sample could not be reproduced"
                    exit 498
                }
                local dm_`grp' `dm_`grp'' `d`ii''
                local tmpall  `tmpall'  `d`ii''
                local dispall `dispall' `vn'
            }
        }

        * nooutput skips ivreg2's whole display block.  That matters beyond
        * tidiness: when an absorbed fixed effect is nested inside the cluster
        * variable, the residuals sum to zero within every cluster, so the added
        * constant column contributes a zero score and ivreg2 warns that the
        * moment covariance matrix is rank deficient.  That warning is expected
        * here and irrelevant -- every statistic reported below comes from the
        * reference run, and the warning only concerns ivreg2's own J statistic.
        qui ivreg2 `dm_depvar' `dm_inexog' (`dm_endo' = `dm_exex') ///
                   if `esample' `wt', `augopts' nooutput

        matrix `b' = e(b)
        matrix `V' = e(V)

        * map the temporary variables back onto the names the user typed
        local cn : colnames `b'
        local newnames
        local refidx
        foreach nm of local cn {
            if "`nm'" == "_cons" {
                local newnames `newnames' _cons
                local refidx `refidx' 0
            }
            else {
                local ix : list posof "`nm'" in tmpall
                if `ix' == 0 {
                    di as err "internal error: cannot map coefficient `nm' back to a variable"
                    exit 498
                }
                local dn : word `ix' of `dispall'
                local jx = colnumb(`bref', "`dn'")
                if `jx' >= . {
                    di as err "internal error: `dn' is absent from the reference estimates"
                    exit 498
                }
                local newnames `newnames' `dn'
                local refidx `refidx' `jx'
            }
        }
        matrix colnames `b' = `newnames'
        matrix colnames `V' = `newnames'
        matrix rownames `V' = `newnames'

        * ------------------------------------------------------- verification
        * The augmented run must reproduce the reference slopes, and its unscaled
        * variances must differ from the reference ones by a single common
        * scalar -- that scalar is ivreg2's small-sample correction.
        local k = colsof(`b')
        local fmin = .
        local fmax = 0
        local nchk 0
        local badb ""
        forvalues i = 1/`k' {
            local j : word `i' of `refidx'
            if `j' == 0 continue
            local nm : word `i' of `newnames'

            local a1 = `b'[1, `i']
            local a2 = `bref'[1, `j']
            local rd = abs(`a1' - `a2') / max(1, abs(`a2'))
            if `rd' > `worstb' local worstb = `rd'
            if `rd' > 1e-5 local badb `badb' `nm'

            local v1 = `V'[`i', `i']
            local v2 = `Vref'[`j', `j']
            if `v1' > 0 & `v1' < . & `v2' > 0 & `v2' < . {
                local r = `v2' / `v1'
                if `r' < `fmin' local fmin = `r'
                if `r' > `fmax' local fmax = `r'
                local ++nchk
            }
        }

        if "`badb'" != "" {
            di as err "verification failed: the augmented second stage does not reproduce"
            di as err "  the ivreghdfe point estimates for:`badb'"
            di as err "  worst relative difference = " %12.4e `worstb'
            di as err "  try a tighter tolerance(), e.g. tolerance(1e-12)"
            exit 498
        }
        if `worstb' > 1e-7 {
            di as txt "note: point estimates match ivreghdfe only to " %7.1e `worstb'
            di as txt "      (limited by the HDFE convergence tolerance)"
        }
        if `nchk' == 0 {
            di as err "verification failed: no usable variance to calibrate against"
            exit 498
        }

        local smallfac = (`fmin' + `fmax') / 2
        if `smallfac' <= 0 | `smallfac' >= . {
            di as err "verification failed: implausible small-sample factor `smallfac'"
            exit 498
        }
        * The spread is pure floating-point noise from the conditioning of the
        * meat matrix, so it grows as the number of clusters falls (measured:
        * ~5e-11 with 3410 clusters, ~7e-9 with 12) and does not shrink with a
        * tighter tolerance().  A structurally wrong VCE would show a spread of
        * order 0.1-1, hence the generous error threshold.
        local spread = (`fmax' - `fmin') / `smallfac'
        if `spread' > 1e-4 {
            di as err "verification failed: the variance ratio against ivreghdfe is not"
            di as err "  constant across coefficients (min = " %18.12g `fmin' ///
                      ", max = " %18.12g `fmax' ")"
            di as err "  the augmented run did not reproduce the reference VCE"
            exit 498
        }
        if `spread' > 1e-6 {
            di as txt "note: the variance ratio against ivreghdfe varies by " %7.1e `spread' _n ///
                      "      (ill-conditioned VCE; affects only the _cons row)"
        }

        matrix `V' = `V' * `smallfac'

        * adopt the reference slope block verbatim, so that the reported slopes
        * and their standard errors are identical to ivreghdfe's
        forvalues i = 1/`k' {
            local ji : word `i' of `refidx'
            if `ji' == 0 continue
            matrix `b'[1, `i'] = `bref'[1, `ji']
            forvalues c = 1/`k' {
                local jc : word `c' of `refidx'
                if `jc' == 0 continue
                matrix `V'[`i', `c'] = `Vref'[`ji', `jc']
            }
        }
    }

    * ------------------------------------------------------------ post it all
    local dofopt
    if "`small'" != "" & `ref_df_r' < . local dofopt "dof(`ref_df_r')"

    ereturn post `b' `V', esample(`esample') depname(`depvar') obs(`ref_N') `dofopt'

    foreach s of local eslist {
        ereturn scalar `s' = `S_`s''
    }
    foreach m of local emlist {
        ereturn local `m' `"`M_`m''"'
    }
    foreach x of local exlist {
        ereturn matrix `x' = `X_`x''
    }

    * diagnostics under the names ivreghdfef used, so that existing
    * reg2docx / esttab scalars() lists keep working unchanged
    ereturn scalar rklm     = e(idstat)
    ereturn scalar rklmp    = e(idp)
    ereturn scalar hansenj  = e(j)
    ereturn scalar hansenjp = e(jp)

    * outcome of the internal verification against the reference run, so that it
    * can be inspected (and tested) from outside
    ereturn scalar smallfac = `smallfac'
    ereturn scalar chk_b    = `worstb'
    ereturn scalar chk_v    = `spread'
    ereturn scalar df_r_ref = `ref_df_r'
    ereturn matrix fsb = `fsb'
    ereturn matrix fsp = `fsp'

    ereturn local endovars   "`endo'"
    ereturn local exexvars   "`exex'"
    ereturn local inexogvars "`inexog'"
    ereturn local absvarlist `"`absorb'"'
    ereturn local fsstored   "`fsstored'"
    ereturn local dispfmt    "`format'"
    ereturn local nosummary  "`nosummary'"
    ereturn local cmdline    `"`cmdline'"'
    ereturn local cmd        "ivreghdfef2"

    if "`store'" != "" {
        ereturn local ssstored "`store'_second"
        capture estimates drop `store'_second
        estimates store `store'_second
    }
end


* ---------------------------------------------------------------------------
* Display
* ---------------------------------------------------------------------------
program define Display
    local fmt `e(dispfmt)'
    if "`fmt'" == "" local fmt "%10.4f"

    di as txt _n "IV (2SLS) estimation with high-dimensional fixed effects"
    di as txt "{hline 78}"
    di as txt "Number of obs" _col(28) "=" as res %14.0fc e(N)
    if e(df_a) < .     di as txt "Absorbed FE dof" _col(28) "=" as res %14.0fc e(df_a)
    if e(N_clust) < .  di as txt "Number of clusters" _col(28) "=" as res %14.0fc e(N_clust)
    if e(F) < .        di as txt "Model F statistic" _col(28) "=" as res %14.4f e(F)
    if e(chi2) < .     di as txt "Model chi2" _col(28) "=" as res %14.4f e(chi2)
    if e(r2) < .       di as txt "R-squared" _col(28) "=" as res %14.4f e(r2)
    if e(r2_a) < .     di as txt "Adj R-squared" _col(28) "=" as res %14.4f e(r2_a)
    if e(rmse) < .     di as txt "Root MSE" _col(28) "=" as res %14.4f e(rmse)
    di as txt "{hline 78}"

    ereturn display

    if `"`e(absvarlist)'"' != "" {
        di as txt "Absorbed fixed effects: " as res `"`e(absvarlist)'"'
        di as txt "_cons follows reghdfe's convention (the mean of the absorbed fixed"
        di as txt "effects); with fixed effects it is not a structural parameter."
    }

    if "`e(nosummary)'" != "" exit

    tempname b V fsb fsp
    matrix `b'   = e(b)
    matrix `V'   = e(V)
    matrix `fsb' = e(fsb)
    matrix `fsp' = e(fsp)

    di as res _n "{hline 24} Results Summary {hline 37}"

    local endo `e(endovars)'
    local exex `e(exexvars)'

    di as res "==> First stage coefficients"
    local jj 0
    foreach ev of local endo {
        local ++jj
        local ii 0
        foreach zv of local exex {
            local ++ii
            local bb = `fsb'[`jj', `ii']
            local pp = `fsp'[`jj', `ii']
            if `bb' >= . continue
            Stars `pp'
            di as txt "beta for `ev' <- `zv': " _c
            di as res `fmt' `bb' as res "`r(stars)'"
        }
    }

    di as res "==> Second stage coefficients"
    foreach ev of local endo {
        local i = colnumb(`b', "`ev'")
        if `i' >= . continue
        local bb = `b'[1, `i']
        local se = sqrt(`V'[`i', `i'])
        local tt = `bb' / `se'
        if e(df_r) < . local pp = 2 * ttail(e(df_r), abs(`tt'))
        else           local pp = 2 * normal(-abs(`tt'))
        Stars `pp'
        di as txt "beta for `e(depvar)' <- `ev'_hat: " _c
        di as res `fmt' `bb' as res "`r(stars)'"
    }

    di as res "==> Weak identification"
    di as txt "Cragg-Donald Wald F statistic: " _c
    di as res `fmt' e(cdf)
    di as txt "Kleibergen-Paap rk Wald F statistic: " _c
    di as res `fmt' e(rkf)

    di as res "==> Underidentification"
    Stars `=e(rklmp)'
    di as txt "Kleibergen-Paap rk LM statistic: " _c
    di as res `fmt' e(rklm) as res "`r(stars)'"

    di as res "==> Overidentification"
    Stars `=e(hansenjp)'
    di as txt "Hansen J statistic: " _c
    di as res `fmt' e(hansenj) as res "`r(stars)'"

    if "`e(fsstored)'" != "" {
        di as txt _n "Stored estimates: " as res "`e(fsstored)' `e(ssstored)'"
    }
end


* ---------------------------------------------------------------------------
* Stars
* ---------------------------------------------------------------------------
program define Stars, rclass
    args p
    local s ""
    if `p' < 0.01      local s "***"
    else if `p' < 0.05 local s "**"
    else if `p' < 0.1  local s "*"
    return local stars "`s'"
end
