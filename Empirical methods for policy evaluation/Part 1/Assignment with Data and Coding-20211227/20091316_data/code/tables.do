/**********************************************************
 *
 * TABLES.DO: MAKES TABLES FOR VOTING PAPER
 *
 *
 **********************************************************/

cap log close
set linesize 255

**********************************************************
* PRELIMINARIES
**********************************************************
version 10
clear all
set mem 1g
set matsize 5000
set more off
adopath + "..\external\"
loadglob using "input_param.txt"
global maxwindow = 1

cap erase ..\output\tables.txt
cap erase ..\output\appendixtable.txt
cap erase ..\output\onlinetables.txt

*****************************************************************
* TABLE <tab:transition>
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

gen numcomp_trunc = numdailies
replace numcomp_trunc = 4 if numdailies>=4
gen Lnumcomp_trunc = L.numcomp_trunc
drop if Lnumcomp_trunc==.
tab numcomp_trunc, gen(znumcomp_trunc)
collapse (sum) znumcomp_trunc*, by(Lnumcomp_trunc)
foreach level of numlist 1(1)5{
	local level1 = `level'-1
	rename znumcomp_trunc`level' znumcomp_trunc`level1'
	replace znumcomp_trunc`level1' = . if Lnumcomp_trunc==`level1'
	}
mkmat znumcomp*, matrix(TABLE)

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.0f) title(<tab:transition>) append

*****************************************************************
* TABLE <tab:nparchive>
*****************************************************************
use "../temp/nparchive.dta", clear
cap matrix drop TABLE
keep if tot_count_nparchive>0&tot_count_nparchive<.&presrepshare~=.

areg rep_share_nparchive polaff_rep, absorb(year) cluster(permid)
	matrix TABLE = (nullmat(TABLE), (_b[polaff_rep] \ _se[polaff_rep] \ . \ . \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ .))
areg rep_share_nparchive polaff_rep presrepshare, absorb(year) cluster(permid)
	matrix TABLE = (nullmat(TABLE), (_b[polaff_rep] \ _se[polaff_rep] \ _b[presrepshare] \ _se[presrepshare] \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ .))
areg rep_share_nparchive polaff_rep if polaff_now~="I", absorb(year) cluster(permid)
	matrix TABLE = (nullmat(TABLE), (_b[polaff_rep] \ _se[polaff_rep] \ . \ . \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ .))
areg rep_share_nparchive polaff_rep if polaff_now=="I", absorb(year) cluster(permid)
	matrix TABLE = (nullmat(TABLE), (_b[polaff_rep] \ _se[polaff_rep] \ . \ . \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ .))


areg rep_share_nparchive polaff_rep polaff_rep_now_ind polaff_now_ind, absorb(year) cluster(permid)
	test polaff_rep_now_ind=0
	matrix TABLE = (nullmat(TABLE), (. \ . \ . \ . \ . \ . \ . \ . \ round(r(p), 0.0001)))

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:nparchive>) append

*****************************************************************
* TABLE <tab:Turnout_main>
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))
end

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x, changein(numdailies) maxchange($maxchange) window(1)

areg D.readshare_hhld x_0 $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	addtotable

areg D.prestout x_0 if mainsample, absorb(styr) cluster(cnty90)
	addtotable
areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	addtotable

areg D.congtout x_0 $demolist $misdemolist if mainsample&abs(D.congtout)<1, absorb(styr) cluster(cnty90)
	addtotable

areg D_congtout_long x_0 $demolist $misdemolist if mainsample&abs(D_congtout_long)<1, absorb(styr) cluster(cnty90)
	addtotable

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main>) append

*****************************************************************
* TABLE <tab:Turnout_comp>
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x1, changein(numdailies>=1) $ifmulti window(1)
define_event x2, changein(numdailies>=2) $ifmulti window(1)
define_event x3, changein(numdailies>=3) $ifmulti window(1)
define_event xdiv, changein(num_polaff_R>0&num_polaff_D>0) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.readshare_hhld x1_0 x2_0 x3_0 $demolist $misdemolist $multicontrol if mainsample_circ, absorb(styr) cluster(cnty90)
	lincom x1_0+x2_0
	lincom x1_0+x2_0+x3_0
	test x1_0=x2_0=x3_0
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
	estimates save "..\temp\est_read", replace
areg D.prestout x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	lincom x1_0+x2_0
	lincom x1_0+x2_0+x3_0
	test x1_0=x2_0=x3_0
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
	estimates save "..\temp\est_pout", replace
areg D.prestout x1_0 x2_0 x3_0 xdiv_0 $multicontrol $demolist $misdemolist if mainsample&num_polaff_I==0&L.num_polaff_I==0&num_polaff_none==0&L.num_polaff_none==0, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ _b[xdiv_0] \ _se[xdiv_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
areg D.congshare x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample&abs(D.congshare)<1, absorb(styr) cluster(cnty90)
	lincom x1_0+x2_0
	lincom x1_0+x2_0+x3_0
	test x1_0=x2_0=x3_0
	local p = round(r(p), 0.0001)
*	outreg2 x1_0 x2_0 x3_0 using "..\output\tables\tab_turnout_comp.tex", se noaster noobs nonotes nocons bdec(4) label rdec(4) nor2 addstat(F-test of equality of coefficients, r(F), p-value, `p', R2, e(r2), Number of counties, e(N_clust), Number of county-years, e(N)) append comma

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_comp>) append
	
*****************************************************************
* TABLE <tab:Turnout_time>
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[x1_0_newspaper] \ _se[x1_0_newspaper] \ _b[x1_0_radio] \ _se[x1_0_radio] \ _b[x1_0_tv] \ _se[x1_0_tv] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
end

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x1, changein(numdailies>=1) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)
gen x1_0_newspaper = x1_0*newspaper
gen x1_0_radio = x1_0*radio
gen x1_0_tv = x1_0*tv

areg D.readshare_hhld x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time_circ, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	addtotable
areg D.prestout x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	addtotable
areg D.congtout x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time&abs(D.congtout)<1, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	addtotable
areg D_congtout_long x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time&abs(D_congtout_long)<1, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	addtotable

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_time>) append

*****************************************************************
* TABLE <tab:Turnout_placebo>
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[xpol_0] \ _se[xpol_0] \ _b[xnonpol_0] \ _se[xnonpol_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
end

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

gen num_pol = numdailies-num_nonpol
define_event xpol, changein(num_pol) $ifmulti maxchange($maxchange) window(1)
define_event xnonpol, changein(num_nonpol) $ifmulti maxchange($maxchange) window(1)
define_event xmulti, eventis(multi) window(1)

areg D.lpreseligible xpol_0 xnonpol_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test xpol_0 = xnonpol_0
	addtotable
areg D.prestout xpol_0 xnonpol_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test xpol_0 = xnonpol_0
	addtotable

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_placebo>) append

*****************************************************************
* TABLE <tab:Repshare_main>
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.circrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.presrepshare xRDdiff_0 $multicontrol if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))
	estimates save "..\temp\est_ciu", replace
areg D.congrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Repshare_main>) append

*****************************************************************
* TABLE <tab:Repshare_comp>
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff1_0] \ _se[xRDdiff1_0] \ _b[xRDdiff2_0] \ _se[xRDdiff2_0] \ _b[xRDdiff3_0] \ _se[xRDdiff3_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
end

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x1, changein(numdailies==1) $ifmulti window(1)
define_event x2, changein(numdailies==2) $ifmulti window(1)
define_event x3, changein(numdailies>=3) $ifmulti window(1)
define_event xRDdiff1, changein((num_polaff_R-num_polaff_D)*(numdailies==1)) maxchange($maxchange) $ifmulti window(1)
define_event xRDdiff2, changein((num_polaff_R-num_polaff_D)*(numdailies==2)) maxchange($maxchange) $ifmulti window(1)
define_event xRDdiff3, changein((num_polaff_R-num_polaff_D)*(numdailies>=3)) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.circrepshare xRDdiff1_0 xRDdiff2_0 xRDdiff3_0 x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	test xRDdiff1_0=xRDdiff2_0=xRDdiff3_0
	addtotable
areg D.presrepshare xRDdiff1_0 xRDdiff2_0 xRDdiff3_0 x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test xRDdiff1_0=xRDdiff2_0=xRDdiff3_0
	addtotable
areg D.congrepshare xRDdiff1_0 xRDdiff2_0 xRDdiff3_0 x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test xRDdiff1_0=xRDdiff2_0=xRDdiff3_0
	addtotable

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Repshare_comp>) append

*****************************************************************
* TABLE {tab_incumb_comp}
*****************************************************************
use "../temp/voting_district_clean.dta", clear
cap matrix drop TABLE
global panelvar "stdist"
define_event x1, changein(numdailies>=1) $ifmulti window(1)
define_event x2, changein(numdailies>=2) $ifmulti window(1)
define_event x3, changein(numdailies>=3) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

xi: areg D.congshare_inccand x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(stdist)
	test x1_0=x2_0=x3_0
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ e(r2) \ e(N_clust) \ e(N)))
xi: areg D.uncontested_inc x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample_unc, absorb(styr) cluster(stdist)
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ e(r2) \ e(N_clust) \ e(N)))

global panelvar "cnty90"

matrix_to_txt, saving(..\output\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Incumb_comp>) append

*****************************************************************
* APPENDIXTABLE <tab:Robustness>
*****************************************************************
* MAIN SPECIFICATIONS (FOR COMPARISON)
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE1
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE1 = (nullmat(TABLE1), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE1 = (nullmat(TABLE1), (_b[xRDdiff_0] \ _se[xRDdiff_0]))

* TRUNCATION
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE2
define_event x, changein(numdailies) maxchange(1) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange(1) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE2 = (nullmat(TABLE2), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE2 = (nullmat(TABLE2), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
 
* ISOLATED MARKETS
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE3
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample&pmsa90==.&nummarkets==1, absorb(styr) cluster(cnty90)
	matrix TABLE3 = (nullmat(TABLE3), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample&pmsa90==.&nummarkets==1, absorb(styr) cluster(cnty90)
	matrix TABLE3 = (nullmat(TABLE3), (_b[xRDdiff_0] \ _se[xRDdiff_0]))

* CUMULATIVE EFFECTS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE4
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

gen cumul = 1
areg D.prestout x_0 x_p1 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[x_0]+_b[x_p1], post
	matrix TABLE4 = (nullmat(TABLE4), (_b[cumul] \ _se[cumul]))
areg D.presrepshare xRDdiff_0 xRDdiff_p1 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_p1], post
	matrix TABLE4 = (nullmat(TABLE4), (_b[cumul] \ _se[cumul]))
global maxwindow = 1

* CONTROLLING FOR POLYNOMIAL TRENDS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE5
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

areg D.prestout x_0 `x_poly$maxpolyorder' if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE5 = (nullmat(TABLE5), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 `xRDdiff_poly$maxpolyorder' if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE5 = (nullmat(TABLE5), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
global maxwindow = 1
	
* CONTROLLING FOR RESTRICTED TRENDS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE6
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

quietly areg D.lpreseligible `x_all' if mainsample&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
predict pred_Dlpreselig_xall
quietly areg D.lpreseligible `xRDdiff_all' if mainsample&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
quietly predict pred_Dlpreselig_xall_RDdiff

areg D.prestout x_0 pred_Dlpreselig_xall if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE6 = (nullmat(TABLE6), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 pred_Dlpreselig_xall_RDdiff if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE6 = (nullmat(TABLE6), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
global maxwindow = 1

* SMALL VS LARGE
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE7
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)
gen x_0_small = x_0*(abs(D.circ_elig)<0.1)
gen x_0_large = x_0*(abs(D.circ_elig)>=0.1)
gen xRDdiff_0_small = xRDdiff_0*(abs(D.circrepshare)<0.2)
gen xRDdiff_0_large = xRDdiff_0*(abs(D.circrepshare)>=0.2)

areg D.prestout x_0_small x_0_large $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE7 = (nullmat(TABLE7), (_b[x_0_small] \ _se[x_0_small] \ _b[x_0_large] \ _se[x_0_large]))
areg D.presrepshare xRDdiff_0_small xRDdiff_0_large $multicontrol $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE7 = (nullmat(TABLE7), (_b[xRDdiff_0_small] \ _se[xRDdiff_0_small] \ _b[xRDdiff_0_large] \ _se[xRDdiff_0_large]))

* EXCLUDE COUNTIES WITH BORDER CHANGES
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE8
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample & border_changes == 0, absorb(styr) cluster(cnty90)
	matrix TABLE8 = (nullmat(TABLE8), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample & border_changes == 0, absorb(styr) cluster(cnty90)
	matrix TABLE8 = (nullmat(TABLE8), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* EXCLUDE SHORT-LIVED PAPERS
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE9
define_event x, changein(numdailies_y4) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R_y4-num_polaff_D_y4) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE9 = (nullmat(TABLE9), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE9 = (nullmat(TABLE9), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* CONTROL FOR LITERACY
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE10
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist D_ishare_lit mis_D_ishare_lit if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE10 = (nullmat(TABLE10), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist D_ishare_lit mis_D_ishare_lit if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE10 = (nullmat(TABLE10), (_b[xRDdiff_0] \ _se[xRDdiff_0]))

* BLOCK BOOTSTRAP
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE11
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)
gen D_prestout = D.prestout
gen D_presrepshare = D.presrepshare
xtset, clear

bootstrap, reps(100) cluster(cnty90) seed(120879): areg D_prestout x_0 $demolist $misdemolist if mainsample, absorb(styr)
	matrix TABLE11 = (nullmat(TABLE11), (_b[x_0] \ _se[x_0]))
bootstrap, reps(100) cluster(cnty90) seed(120879): areg D_presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr)
	matrix TABLE11 = (nullmat(TABLE11), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* ALLOWING FOR AR1
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE12
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

xtset cnty90 year, delta(4)
xtreg D.prestout x_0 $demolist $misdemolist i.styr if mainsample, pa corr(ar1)
	matrix TABLE12 = (nullmat(TABLE12), (_b[x_0] \ _se[x_0]))
xtset cnty90 year, delta(4)
xtreg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist  i.styr if mainsample, pa corr(ar1)
	matrix TABLE12 = (nullmat(TABLE12), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* JOIN ALL TABLES TO CREATE APPENDIX TABLE
cap matrix drop TABLE
matrix TABLE = (TABLE1 \ TABLE2 \ TABLE3 \ TABLE4 \ TABLE5 \ TABLE6 \ TABLE7 \ TABLE8 \ TABLE9 \ TABLE10 \ TABLE11 \ TABLE12)

matrix_to_txt, saving(..\output\appendixtable.txt) mat(TABLE) format(%20.4f) title(<tab:Robustness>) append

*****************************************************************
* ONLINE APPENDIX TABLE <tab:Serialcorr>: Alternative corrections for serial correlation
*****************************************************************

* MAIN SPECIFICATIONS (FOR COMPARISON)
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE1
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE1 = (nullmat(TABLE1), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE1 = (nullmat(TABLE1), (_b[xRDdiff_0] \ _se[xRDdiff_0]))

* BLOCK BOOTSTRAP
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE2
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)
gen D_prestout = D.prestout
gen D_presrepshare = D.presrepshare
xtset, clear

bootstrap, reps(100) cluster(cnty90) seed(120879): areg D_prestout x_0 $demolist $misdemolist if mainsample, absorb(styr)
	matrix TABLE2 = (nullmat(TABLE2), (_b[x_0] \ _se[x_0]))
bootstrap, reps(100) cluster(cnty90) seed(120879): areg D_presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr)
	matrix TABLE2 = (nullmat(TABLE2), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
 
 * CLUSTER BY STATE DECADE (BESTER CONLEY HANSEN)
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE3
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)
gen decade = int(year/10)
egen state_decade = group(state decade)

areg D.prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(state_decade)
	matrix TABLE3 = (nullmat(TABLE3), (_b[x_0] \ _se[x_0]))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(state_decade)
	matrix TABLE3 = (nullmat(TABLE3), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* RANDOM EFFECTS
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE4
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

xtset cnty90 year, delta(4)
xtreg D.prestout x_0 $demolist $misdemolist i.styr if mainsample, re
	matrix TABLE4 = (nullmat(TABLE4), (_b[x_0] \ _se[x_0]))
xtset cnty90 year, delta(4)
xtreg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist  i.styr if mainsample, re
	matrix TABLE4 = (nullmat(TABLE4), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* ALLOWING FOR AR1
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE5
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

xtset cnty90 year, delta(4)
xtreg D.prestout x_0 $demolist $misdemolist i.styr if mainsample, pa corr(ar1)
	matrix TABLE5 = (nullmat(TABLE5), (_b[x_0] \ _se[x_0]))
xtset cnty90 year, delta(4)
xtreg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist  i.styr if mainsample, pa corr(ar1)
	matrix TABLE5 = (nullmat(TABLE5), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* ALLOWING FOR AR2
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE6
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

xtset cnty90 year, delta(4)
xtreg D.prestout x_0 $demolist $misdemolist i.styr if mainsample, pa corr(ar2)
	matrix TABLE6 = (nullmat(TABLE6), (_b[x_0] \ _se[x_0]))
xtset cnty90 year, delta(4)
xtreg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist  i.styr if mainsample, pa corr(ar2)
	matrix TABLE6 = (nullmat(TABLE6), (_b[xRDdiff_0] \ _se[xRDdiff_0]))
	
* JOIN ALL TABLES TO CREATE TABLE
cap matrix drop TABLE
matrix TABLE = (TABLE1 \ TABLE2 \ TABLE3 \ TABLE4 \ TABLE5 \ TABLE6)

matrix_to_txt, saving(..\output\onlinetables.txt) mat(TABLE) format(%20.4f) title(<tab:text_Serialcorr>) append

*****************************************************************
* ONLINE APPENDIX TABLE <tab:Governor>: GUBERNATORIAL ELECTIONS
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
define_event x1, changein(numdailies>=1) $ifmulti window(1)
gen x1_0_newspaper = x1_0*newspaper
gen x1_0_radio = x1_0*radio
gen x1_0_tv = x1_0*tv
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

cap matrix drop TABLE

areg D.govtout x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	matrix TABLE = (nullmat(TABLE), ( _b[x1_0_newspaper] \ _se[x1_0_newspaper] \ _b[x1_0_radio] \ _se[x1_0_radio] \ _b[x1_0_tv] \ _se[x1_0_tv] \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

areg D.govrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ _b[xRDdiff_0] \ _se[xRDdiff_0] \ . \ . \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\onlinetables.txt) mat(TABLE) format(%20.4f) title(<tab:Governor>) append

*****************************************************************
* ONLINE APPENDIX TABLE <tab:Senator>: SENATORIAL ELECTIONS
*****************************************************************
* note: direct election of senators began in 1913 so newspaper period is poorly identified
use "../temp/voting_cnty_clean.dta", clear
define_event x1, changein(numdailies>=1) $ifmulti window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

gen x1_0_newspaper = x1_0*newspaper
gen x1_0_radio = x1_0*radio
gen x1_0_tv = x1_0*tv

gen xRDdiff_0_newspaper = xRDdiff_0*newspaper
gen xRDdiff_0_radio = xRDdiff_0*radio
gen xRDdiff_0_tv = xRDdiff_0*tv

cap matrix drop TABLE
areg D.sentout x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test x1_0_radio=x1_0_tv
	matrix TABLE = (nullmat(TABLE), ( _b[x1_0_radio] \ _se[x1_0_radio] \ _b[x1_0_tv] \ _se[x1_0_tv] \ . \ . \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
areg D.senrepshare xRDdiff_0_newspaper xRDdiff_0_radio xRDdiff_0_tv  $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test xRDdiff_0_radio=xRDdiff_0_tv
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[xRDdiff_0_radio] \ _se[xRDdiff_0_radio] \ _b[xRDdiff_0_tv]	\ _se[xRDdiff_0_tv] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
matrix_to_txt, saving(..\output\onlinetables.txt) mat(TABLE) format(%20.4f) title(<tab:Senator>) append

*****************************************************************
* ONLINE APPENDIX TABLE <tab:Repalways>: PERMANENT AFFILIATES
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

define_event xRDdiff, changein(num_always_R-num_always_D) maxchange($maxchange)

cap matrix drop TABLE
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.congrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0] \ _se[xRDdiff_0] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\onlinetables.txt) mat(TABLE) format(%20.4f) title(<tab:Repalways>) append

*****************************************************************
* TABLE <tab:Updown>: TURNOUT EFFECT BY DIRECTION OF CHANGE
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
	gen x_0_up = x_0*(x_0>0)
	gen x_0_down = x_0*(x_0<0)
define_event xmulti, eventis(multi)

areg D.prestout x_0_up x_0_down $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test x_0_up = x_0_down
	matrix TABLE = (nullmat(TABLE), (_b[x_0_up] \ _se[x_0_up] \ _b[x_0_down] \ _se[x_0_down] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\onlinetables.txt) mat(TABLE) format(%20.4f) title(<tab:Updown>) append

