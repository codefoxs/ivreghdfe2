*! test_ivreghdfef2.do -- numerical correctness tests for ivreghdfef2
*  Run with:  do tests/test_ivreghdfef2.do
*  Results are written to tests/test_results.log

clear all
set more off
set linesize 160

capture log close _all
log using "D:/Code/Github/ivreghdfef2/tests/test_results.log", replace text

adopath ++ "D:/Code/Github/ivreghdfef2"
discard
which ivreghdfef2

global NFAIL 0
global NPASS 0

* ---------------------------------------------------------------- helpers
capture program drop sevec
program define sevec
    args out
    tempname VV bb ss
    matrix `VV' = e(V)
    matrix `bb' = e(b)
    local k = colsof(`bb')
    matrix `ss' = J(1, `k', .)
    forvalues i = 1/`k' {
        matrix `ss'[1, `i'] = sqrt(`VV'[`i', `i'])
    }
    matrix colnames `ss' = `: colnames `bb''
    matrix `out' = `ss'
end

* mchk <label> <ref matrix> <new matrix> [tol]
* iterates over the columns of the REFERENCE matrix, so extra columns in the
* new matrix (namely _cons) are simply ignored
capture program drop mchk
program define mchk
    args label m1 m2 tol
    if "`tol'" == "" local tol 1e-8
    local nm : colnames `m1'
    local bad 0
    local worst 0
    foreach v of local nm {
        local i = colnumb(`m1', "`v'")
        local j = colnumb(`m2', "`v'")
        if `j' >= . {
            di as err "    absent from new matrix: `v'"
            local bad 1
            continue
        }
        local a = `m1'[1, `i']
        local b = `m2'[1, `j']
        local d = abs(`a' - `b') / max(1, abs(`b'))
        if `d' > `worst' local worst = `d'
        if `d' > `tol' {
            local bad 1
            di as err "    `v': " %20.12g `a' " vs " %20.12g `b'
        }
    }
    if `bad' {
        di as err  "FAIL  `label'  worst rel diff " %8.2e `worst'
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt  "pass  `label'  worst rel diff " %8.2e `worst'
        global NPASS = $NPASS + 1
    }
end

* schk <label> <value> <expected> [tol]
capture program drop schk
program define schk
    args label a b tol
    if "`tol'" == "" local tol 1e-8
    local d = abs(`a' - `b') / max(1, abs(`b'))
    if `d' > `tol' {
        di as err "FAIL  `label'  " %20.12g `a' " vs " %20.12g `b' "  rel " %8.2e `d'
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt "pass  `label'  rel " %8.2e `d'
        global NPASS = $NPASS + 1
    }
end

* ------------------------------------------------------------------- data
webuse nlswork, clear
qui keep if !missing(ln_wage, ttl_exp, age, tenure, not_smsa, south, ///
                     hours, union, idcode, year, race, ind_code, msp, nev_mar, wks_work)
qui gen double ttl_exp_z = ttl_exp
di as res _n "N = " _N

local Y   ln_wage
local X   ttl_exp
local C   age tenure not_smsa south
local Z   hours union
local FE  absorb(idcode year)
local TOL tolerance(1e-10)

* =======================================================================
* TEST 1  slopes and standard errors must equal ivreghdfe
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 1  ivreghdfef2 vs ivreghdfe (slopes and SEs)"
di as res "{hline 70}"

local vcelist `" " " "robust" "cluster(idcode)" "cluster(ind_code)" "cluster(idcode year)" "small" "small robust" "small cluster(idcode)" "small cluster(ind_code)" "'

local i 0
foreach o of local vcelist {
    local ++i
    qui ivreghdfe `Y' `C' (`X' = `Z'), `FE' `TOL' `o'
    tempname b1 s1
    matrix `b1' = e(b)
    sevec `s1'

    qui ivreghdfef2 `Y' `C' (`X' = `Z'), `FE' `TOL' `o'
    tempname b2 s2
    matrix `b2' = e(b)
    sevec `s2'

    di as txt "  spec `i': [`o']   smallfac = " %14.10f e(smallfac)
    mchk b_s`i'  `b1' `b2' 1e-8
    mchk se_s`i' `s1' `s2' 1e-8

    * the slope block is copied from the reference run, so the comparison above
    * only checks the plumbing.  e(chk_b) / e(chk_v) hold the command's own
    * internal comparison of the raw augmented estimates against the reference:
    * chk_b is the worst relative coefficient difference, chk_v the spread of
    * the variance ratio across coefficients.  Both must be ~0.
    * chk_v is floating-point noise from the conditioning of the meat matrix and
    * is tolerance-independent (verified: 6.503e-09 at tolerance 1e-10, 1e-12 and
    * 1e-13 alike for cluster(ind_code), which has only 12 clusters; ~5e-11 for
    * cluster(idcode), which has 3410).  It enters the results only through
    * smallfac, which scales the _cons row.
    schk chk_b_s`i' `=e(chk_b)' 0 1e-10
    schk chk_v_s`i' `=e(chk_v)' 0 1e-7
}

* =======================================================================
* TEST 2  _cons equals ybar - Xbar'b  (reghdfe's convention)
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 2  constant term equals ybar - Xbar'b"
di as res "{hline 70}"

qui ivreghdfef2 `Y' `C' (`X' = `Z'), `FE' `TOL' cluster(idcode)
tempname bc
matrix `bc' = e(b)
tempvar smp
qui gen byte `smp' = e(sample)

qui summarize `Y' if `smp', meanonly
local cman = r(mean)
foreach v of local X {
    qui summarize `v' if `smp', meanonly
    local cman = `cman' - r(mean) * `bc'[1, colnumb(`bc', "`v'")]
}
foreach v of local C {
    qui summarize `v' if `smp', meanonly
    local cman = `cman' - r(mean) * `bc'[1, colnumb(`bc', "`v'")]
}
schk cons_formula `=`bc'[1, colnumb(`bc', "_cons")]' `cman' 1e-7

* =======================================================================
* TEST 3  degenerate IV (instrument = the endogenous regressor itself)
*         2SLS is then numerically identical to OLS, so the FULL b and V --
*         including the _cons row -- must equal reghdfe's.  This is the only
*         hard test of the constant term's variance.
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 3  degenerate IV must reproduce reghdfe exactly (incl. _cons)"
di as res "{hline 70}"

* reghdfe spells the vce differently from ivreg2, so translate inside the loop
local ivvce `" " " "robust" "cluster(ind_code)" "'
local i 0
foreach oiv of local ivvce {
    local ++i
    local o
    if strpos("`oiv'", "cluster(") {
        local cv = substr("`oiv'", strpos("`oiv'", "(") + 1, ///
                          strpos("`oiv'", ")") - strpos("`oiv'", "(") - 1)
        local o vce(cluster `cv')
    }
    else if strpos("`oiv'", "robust") {
        local o vce(robust)
    }
    qui reghdfe `Y' `X' `C', `FE' `TOL' `o'
    tempname bo so
    matrix `bo' = e(b)
    sevec `so'

    qui ivreghdfef2 `Y' `C' (`X' = ttl_exp_z), `FE' `TOL' small `oiv'
    tempname bi si
    matrix `bi' = e(b)
    sevec `si'

    di as txt "  spec `i': reghdfe [`o']  vs  ivreghdfef2 [`oiv']"
    mchk degen_b`i'  `bo' `bi' 1e-7
    mchk degen_se`i' `so' `si' 1e-6
}

* Nested case: clustering on the panel id that is also absorbed.  reghdfe
* subtracts one extra dof here (its nested_adj) which ivreghdfe does not, so the
* two differ by a single scalar of order 1/N.  Asserting that the SE ratio is the
* SAME for every coefficient -- _cons included -- is therefore still a real test
* of the constant term's variance: a broken _cons row would not track the slopes.
di as res _n "  nested cluster case: SE ratio must be one common scalar"
qui reghdfe `Y' `X' `C', `FE' `TOL' vce(cluster idcode)
tempname bo4 so4
matrix `bo4' = e(b)
sevec `so4'
qui ivreghdfef2 `Y' `C' (`X' = ttl_exp_z), `FE' `TOL' small cluster(idcode)
tempname bi4 si4
matrix `bi4' = e(b)
sevec `si4'

local rmin = .
local rmax = 0
local nms : colnames `so4'
foreach v of local nms {
    local r = `so4'[1, colnumb(`so4', "`v'")] / `si4'[1, colnumb(`si4', "`v'")]
    if `r' < `rmin' local rmin = `r'
    if `r' > `rmax' local rmax = `r'
    di as txt "    `v'" _col(16) "SE ratio = " %18.14f `r'
}
di as txt "  common scalar = " %18.14f `rmin' "  (so reghdfe and ivreghdfe differ" ///
          " by " %8.2e (`rmin' - 1) ")"
schk nested_ratio_uniform `=(`rmax' - `rmin') / `rmin'' 0 1e-9

* =======================================================================
* TEST 4  the bug being fixed: same point estimates as ivreghdfef,
*         different (correct) standard errors
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 4  ivreghdfef comparison (same b, different SE)"
di as res "{hline 70}"

qui ivreghdfef2 `Y' `C' (`X' = `Z'), `FE' `TOL' cluster(idcode)
tempname bnew snew
matrix `bnew' = e(b)
sevec `snew'

preserve
capture noisily {
    qui ivreghdfef `Y' `X' `C', iv(`Z') o(`FE' cluster(idcode))
    tempname bold sold
    matrix `bold' = e(b)
    sevec `sold'
}
local rc = _rc
restore

if `rc' == 0 {
    mchk old_b `bold' `bnew' 1e-6
    local seo = `sold'[1, colnumb(`sold', "`X'")]
    local sen = `snew'[1, colnumb(`snew', "`X'")]
    di as txt "  SE for `X':  ivreghdfef = " %12.8f `seo' "   ivreghdfef2 = " %12.8f `sen'
    di as txt "  ivreghdfef understates the SE by " %6.2f (1 - `seo'/`sen')*100 "%"
    if abs(`seo'/`sen' - 1) < 1e-6 {
        di as err "FAIL  old_se_differs  (the old SE was already correct?)"
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt "pass  old_se_differs"
        global NPASS = $NPASS + 1
    }
}
else {
    di as txt "  skipped (ivreghdfef not available, rc=`rc')"
}

* =======================================================================
* TEST 5  multiple endogenous regressors
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 5  multiple endogenous regressors"
di as res "{hline 70}"

qui ivreghdfe `Y' not_smsa south (ttl_exp tenure = hours union nev_mar msp), `FE' `TOL' cluster(idcode)
tempname bm1 sm1
matrix `bm1' = e(b)
sevec `sm1'

qui ivreghdfef2 `Y' not_smsa south (ttl_exp tenure = hours union nev_mar msp), ///
    `FE' `TOL' cluster(idcode) store(MULTI)
tempname bm2 sm2
matrix `bm2' = e(b)
sevec `sm2'

mchk multi_b  `bm1' `bm2' 1e-8
mchk multi_se `sm1' `sm2' 1e-8

foreach nm in MULTI_first1 MULTI_first2 MULTI_second {
    capture estimates describe `nm'
    if _rc {
        di as err "FAIL  stored_`nm'"
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt "pass  stored_`nm'"
        global NPASS = $NPASS + 1
    }
}

* =======================================================================
* TEST 6  weights, if/in, no absorb()
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 6  weights, subsamples, no absorb()"
di as res "{hline 70}"

foreach w in aweight fweight pweight {
    qui ivreghdfe `Y' `C' (`X' = `Z') [`w'=wks_work], `FE' `TOL' cluster(idcode)
    tempname bw1 sw1
    matrix `bw1' = e(b)
    sevec `sw1'
    qui ivreghdfef2 `Y' `C' (`X' = `Z') [`w'=wks_work], `FE' `TOL' cluster(idcode)
    tempname bw2 sw2
    matrix `bw2' = e(b)
    sevec `sw2'
    mchk w_b_`w'  `bw1' `bw2' 1e-8
    mchk w_se_`w' `sw1' `sw2' 1e-8
}

qui ivreghdfe `Y' `C' (`X' = `Z') if race == 1 & year > 75, `FE' `TOL' cluster(idcode)
tempname bif1 sif1
matrix `bif1' = e(b)
sevec `sif1'
qui ivreghdfef2 `Y' `C' (`X' = `Z') if race == 1 & year > 75, `FE' `TOL' cluster(idcode)
tempname bif2 sif2
matrix `bif2' = e(b)
sevec `sif2'
mchk if_b  `bif1' `bif2' 1e-8
mchk if_se `sif1' `sif2' 1e-8

* no absorb(): ivreghdfe already reports _cons, so results must be identical
qui ivreghdfe `Y' `C' (`X' = `Z'), cluster(idcode) small
tempname bn1 sn1
matrix `bn1' = e(b)
sevec `sn1'
qui ivreghdfef2 `Y' `C' (`X' = `Z'), cluster(idcode) small
tempname bn2 sn2
matrix `bn2' = e(b)
sevec `sn2'
mchk noabs_b  `bn1' `bn2' 1e-12
mchk noabs_se `sn1' `sn2' 1e-12

* =======================================================================
* TEST 7  time-series operators
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 7  time-series operators in the variable list"
di as res "{hline 70}"

qui xtset idcode year
capture noisily {
    qui ivreghdfe `Y' `C' L.hours (`X' = `Z'), `FE' `TOL' cluster(idcode)
    tempname bt1 st1
    matrix `bt1' = e(b)
    sevec `st1'
    qui ivreghdfef2 `Y' `C' L.hours (`X' = `Z'), `FE' `TOL' cluster(idcode)
    tempname bt2 st2
    matrix `bt2' = e(b)
    sevec `st2'
}
if _rc == 0 {
    mchk ts_b  `bt1' `bt2' 1e-8
    mchk ts_se `st1' `st2' 1e-8
}
else {
    di as err "FAIL  ts_operators (rc=" _rc ")"
    global NFAIL = $NFAIL + 1
}

* =======================================================================
* TEST 8  factor variables and partial() must be refused
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 8  unsupported syntax is refused"
di as res "{hline 70}"

capture ivreghdfef2 `Y' `C' i.race (`X' = `Z'), `FE'
if _rc == 198 {
    di as txt "pass  reject_factor_var"
    global NPASS = $NPASS + 1
}
else {
    di as err "FAIL  reject_factor_var (rc=" _rc ")"
    global NFAIL = $NFAIL + 1
}

capture ivreghdfef2 `Y' `C' (`X' = `Z'), `FE' partial(age)
if _rc == 198 {
    di as txt "pass  reject_partial"
    global NPASS = $NPASS + 1
}
else {
    di as err "FAIL  reject_partial (rc=" _rc ")"
    global NFAIL = $NFAIL + 1
}

* =======================================================================
* TEST 9  table output and full display
* =======================================================================
di as res _n "{hline 70}"
di as res "TEST 9  display and table output"
di as res "{hline 70}"

ivreghdfef2 `Y' `C' (`X' = `Z'), `FE' `TOL' cluster(idcode) first store(IV)

capture which esttab
if !_rc {
    capture noisily esttab IV_first IV_second, se scalars(N r2_a cdf rkf rklm hansenj)
    if _rc {
        di as err "FAIL  esttab (rc=" _rc ")"
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt "pass  esttab"
        global NPASS = $NPASS + 1
    }
}
else {
    di as txt "  esttab not installed, skipped"
}

capture which reg2docx
if !_rc {
    capture reg2docx IV_first IV_second using "D:/Code/Github/ivreghdfef2/tests/_tmp_iv.docx", ///
        replace b(%12.4f) t(%12.4f) scalars(N(%12.0fc) r2_a(%12.4f) cdf(%12.4f) rkf(%12.4f))
    if _rc {
        di as err "FAIL  reg2docx (rc=" _rc ")"
        global NFAIL = $NFAIL + 1
    }
    else {
        di as txt "pass  reg2docx"
        global NPASS = $NPASS + 1
        capture erase "D:/Code/Github/ivreghdfef2/tests/_tmp_iv.docx"
    }
}
else {
    di as txt "  reg2docx not installed, skipped"
}

* replay
capture noisily ivreghdfef2
if _rc {
    di as err "FAIL  replay (rc=" _rc ")"
    global NFAIL = $NFAIL + 1
}
else {
    di as txt "pass  replay"
    global NPASS = $NPASS + 1
}

* =======================================================================
di as res _n "{hline 70}"
di as res "TEST SUMMARY   passed = $NPASS   failed = $NFAIL"
di as res "{hline 70}"

log close
