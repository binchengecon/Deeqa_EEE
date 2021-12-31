/**********************************************************
 *
 * TEXT.DO: REGRESSIONS AND OTHER CALCULATIONS 
 *               TO SUPPORT CLAIMS IN TEXT OF VOTING PAPER
 *
 **********************************************************/

cap log close
log using "../output/text/text.log", text replace
set linesize 255

**********************************************************
* PRELIMINARIES
**********************************************************
version 11
clear all
set mem 1g
set matsize 5000
set more off
adopath + "..\external\"
loadglob using "input_param.txt"

cap erase ..\output\texttables.txt

*****************************************************************
* BASIC SUMMARY STATISTICS
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
define_event x, changein(numdailies) maxchange(1)

* number of years with entry/exit
tab x_0

cap matrix drop TABLE
foreach B of numlist -1/1 {
	cap matrix drop ROW
	quietly count if x_0 == `B'
	matrix TABLE = (nullmat(TABLE) \ r(N))
}

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_entryexit>) append

* circulation changes in event years
gen lcirc = log(circ)
gen absDlcirc = abs(D.lcirc)

sum absDlcirc if miss_circ==0&L.miss_circ==0&x_0~=0, det

* turnout in counties with no newspapers
sum prestout if mainsample&numdailies==0

* size of excluded vs included counties; LR test of equal state-year effects between these two groups
cap matrix drop TABLE
sum preseligible if year>=$samplestart&year<=$sampleend&num_changes==0
	local numerator = r(mean)
	matrix TABLE = (nullmat(TABLE) \ r(mean), .)
sum preseligible if year>=$samplestart&year<=$sampleend&num_changes>0
	matrix TABLE = (nullmat(TABLE) \ r(mean), .)
	local denominator = r(mean)
	local ratio = `numerator'/`denominator'
	matrix TABLE = (nullmat(TABLE) \ `ratio', .)

gen insample = num_changes>0
xtmixed D.prestout insample || styr: insample if year>=$samplestart&year<=$sampleend, covariance(unstructured)
	estimates store unrestricted
xtmixed D.prestout insample || styr:          if year>=$samplestart&year<=$sampleend, covariance(unstructured)
	estimates store restricted
lrtest unrestricted restricted
	matrix TABLE = (nullmat(TABLE) \ r(chi2), r(p))
	
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_excluded>) append
	
use "../external/newspapers_constant.dta", clear
cap matrix drop TABLE
tab polaff

quietly count if polaff=="D"
	matrix TABLE = (nullmat(TABLE) \ r(N), .)
quietly count if polaff=="I"
	matrix TABLE = (nullmat(TABLE) \ r(N), .)
quietly count if polaff=="R"
	matrix TABLE = (nullmat(TABLE) \ r(N), .)
count if polaff==""
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* number of newspapers
unique permid
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* number of papers with an affiliation in 1952
use "../external/newspapers_yearly.dta", clear
count if year==1952
	matrix TABLE[5,2] = r(N)
count if year==1952&polaff~="I"&polaff~=""
	matrix TABLE = (nullmat(TABLE) \ . , r(N))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_polaff>) append

*****************************************************************
* CONTENT ANALYSIS
***************************************************************** 

* NPARCHIVE
use "../temp/nparchive.dta", clear
cap matrix drop TABLE
sort permid polaff
keep if tot_count_nparchive>0&tot_count_nparchive<.

* unique newspapers with nparchive hits
unique permid
quietly count if permid!=permid[_n+1]
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* #R papers with nparchive hits
codebook permid if polaff=="R"
quietly count if polaff=="R" & permid!=permid[_n+1]
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* #D papers with nparchive hits
codebook permid if polaff=="D"
quietly count if polaff=="D" & permid!=permid[_n+1]
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* #I papers with nparchive hits
codebook permid if polaff=="I"
quietly count if polaff=="I" & permid!=permid[_n+1]
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* #no-aff papers with nparchive hits
codebook permid if polaff==""
quietly count if polaff=="" & permid!=permid[_n+1]
	matrix TABLE = (nullmat(TABLE) \ r(N), .)

* total number of hits
sum tot_count_nparchive
display "total hits:" r(N)*r(mean)
	matrix TABLE[1,2] = r(N)*r(mean)

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_nparchive_polaff>) append

* ENDORSEMENTS
use "../temp/endorse.dta", clear
tab polaff party if party==-1|party==1, row
quietly gen rep = party==1 if party==-1|party==1

cap matrix drop TABLE
foreach A in D I R {
	cap matrix drop ROW
	quietly sum rep if polaff == "`A'"
	matrix TABLE = (nullmat(TABLE) \ r(mean))
}

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_endorse>) append

* NON-POLITICAL PAPERS
use "../external/newspapers_constant.dta", clear
gen nopolaff = polaff==""
tab nonpol nopolaff, row

mean nopolaff, over(nonpol)
matrix define TABLE = e(b)
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_nonpol>) append

*****************************************************************
* JOURNALISTS
***************************************************************** 
use "../temp/journalists.dta", clear
cap matrix drop TABLE

* 1880 journalists v. dailies
* graph
scatter journalists_1880 numdailies1880, xtitle("Number of Dailies in 1880") ytitle("Number of Journalists in 1880") mlabel(city)
	graph export "../output/text/text_dailies_journalists_1880.eps", replace
* regression
reg journalists_1880 numdailies1880
	matrix TABLE = (nullmat(TABLE), (_b[numdailies1880] \ _se[numdailies1880] \ . \ . \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ . \ . \ .))

* change 1869-1880
* graph
scatter chnjournalists chndailies,  ytitle("Change in no. of reporters" "1870-1880") xtitle("Change in no. of dailies 1869-1880") mlabel(city)
	graph export "../output/text/text_dailies_journalists_1869_1880.eps", replace
* regression
reg chnjournalists chndailies
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[chndailies] \ _se[chndailies] \ _b[_cons] \ _se[_cons] \ e(N) \ e(r2) \ . \ . \ .))
sum journalists_1870 if e(sample)
	matrix TABLE[9,2] = r(mean)

* summary statistics
sum chnjournalists if chndailies>0, det
	matrix TABLE[10,2] = r(mean)
sum chnjournalists if chndailies<=0, det
	matrix TABLE[11,2] = r(mean)

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_journalists>) append

*****************************************************************
* ANALYZE CORRESPONDENCE BETWEEN INCOME CONCEPTS
*****************************************************************
use "../temp/compare_income.dta", clear
cap matrix drop TABLE

scatter lmfg lincome, xtitle("Log(Income per Capita)") ytitle("Log(Manufacturing Output per Capita)") mlabel(state)
graph export "../output/text/text_mfg_income.eps", replace

reg lmfg lincome
	matrix TABLE = (nullmat(TABLE), (_b[lincome] \ _se[lincome] \ e(N) \ e(r2) \ sqrt(e(r2))))
reg lmfg lincome if state~="ID"&state~="WY"&state~="NV"
	matrix TABLE = (nullmat(TABLE), (_b[lincome] \ _se[lincome] \ e(N) \ e(r2) \ sqrt(e(r2))))
	
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_mfg_income>) append

*****************************************************************
* EFFECTS ON CIRCULATION
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)

areg D.circ_elig x_0 if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0] \ . \ . \ .\ . \ e(r2) \ e(N_clust) \ e(N)))
areg D.circ_elig x_0 $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0] \ . \ . \ .\ . \ e(r2) \ e(N_clust) \ e(N)))
areg D.circ_elig x_n2 x_n1 x_0 x_p1 x_p2 if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_n2] \ _se[x_n2] \ _b[x_n1] \ _se[x_n1] \ _b[x_0] \ _se[x_0] \ _b[x_p1] \ _se[x_p1] \ _b[x_p2] \ _se[x_p2] \ e(r2) \ e(N_clust) \ e(N)))
areg D.circ_elig x_n2 x_n1 x_0 x_p1 x_p2 $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_n2] \ _se[x_n2] \ _b[x_n1] \ _se[x_n1] \ _b[x_0] \ _se[x_0] \ _b[x_p1] \ _se[x_p1] \ _b[x_p2] \ _se[x_p2] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_circ>) append

*****************************************************************
* LIFECYCLE ANALYSIS
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[zsince_entry2] \ _se[zsince_entry2] \ _b[zsince_entry3] \ _se[zsince_entry3] \ _b[zto_exit2] \ _se[zto_exit2] \ _b[zto_exit1] \ _se[zto_exit1] \ e(N) \ e(r2)))
end

use "../temp/lifecycle.dta", clear
cap matrix drop TABLE
* CORE PROPERTIES OF CIRCULATION

** PERSISTENT SHOCKS **
reg lcirc L.lcirc
	matrix TABLE = (nullmat(TABLE), (_b[L.lcirc] \ _se[L.lcirc] \ e(r2) \ e(N)))
areg lcirc L.lcirc, absorb(year)
	matrix TABLE = (nullmat(TABLE), (_b[L.lcirc] \ _se[L.lcirc] \ e(r2) \ e(N)))

** STRONG TIME TRENDS **
** (note both constant and R2) **
areg D.lcirc, absorb(year)

** SMALL SHARE OF VARIANCE EXPLAINED BY CITY POPULATION **
areg D.lcirc D.licitypop, absorb(year)
reg D.lcirc D.licitypop

** TYPE OF CIRCULATION MATTERS BUT DOES NOT EXPLAIN A LOT OF VARIANCE **
xi: areg D.lcirc i.L_circtype_group i.circtype_group, absorb(year)
testparm _I*
drop _I*

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_lifecirc>) append

* BASIC LIFECYCLE FACTS
cap matrix drop TABLE

sum D.lcirc if mainsample, det

foreach type in max median mean{
	** relative circulation variables
	egen `type'_circ_norm = `type'(circ_norm), by(permid)
	gen rel_`type'_circ_norm = circ_norm/`type'_circ_norm
	** balancing variables
	egen num_yeartype_n1_`type' = sum(yeartype==-1&rel_`type'_circ_norm~=.), by(permid)
	egen num_yeartype_0_`type' = sum(yeartype== 0&rel_`type'_circ_norm~=.), by(permid)
	egen num_yeartype_1_`type'  = sum(yeartype== 1&rel_`type'_circ_norm~=.), by(permid)
	oo CIRCULATION RELATIVE TO `type' NORMALIZED FOR TIME
	sum rel_`type'_circ_norm, det
	tab lifespan yeartype if num_yeartype_n1_`type'==1&num_yeartype_1_`type'==1&mainsample, sum(rel_`type'_circ_norm)
	}
mean rel_mean_circ_norm if mainsample, over(yeartype)
	matrix TABLE = e(b)
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_lifemeans>) append

* REGRESSIONS
cap matrix drop TABLE

** MAIN MODELS **
areg D.lcirc zsince_entry2 zsince_entry3 zto_exit2 zto_exit1 if mainsample, absorb(year)
	addtotable
areg D.lcirc zsince_entry2 zsince_entry3 zto_exit2 zto_exit1 D.licitypop if mainsample, absorb(year)
	addtotable
areg D.lcirc zsince_entry2 zsince_entry3 zto_exit2 zto_exit1 D.licitypop if mainsample&circtype_group==L.circtype_group, absorb(year)
	addtotable
areg D.lcirc zsince_entry2 zsince_entry3 zto_exit2 zto_exit1 D.licitypop if mainsample&circtype_group==L.circtype_group&abs(D.lcirc)<1, absorb(year)
	addtotable

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_lifecycle>) append
	
** SEPARATELY BY LIFESPAN TO ADJUST FOR COMPOSITION **
foreach span of numlist 20(4)40{
	areg D.lcirc zsince_entry2 zsince_entry3 zto_exit2 zto_exit1 if mainsample&lifespan==`span', absorb(year)
}
	
*****************************************************************
* CROSS-SECTIONAL DETERMINANTS OF NUMBER OF NEWSPAPERS
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

reg samplemean_numdailies samplemean_preseligible if year==1900&mainsample
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ _b[samplemean_preseligible] \ _se[samplemean_preseligible] \ . \ . \ e(N) \ e(r2)))
reg samplemean_numdailies samplemean_outputpercap if year==1900&mainsample
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ _b[samplemean_outputpercap] \ _se[samplemean_outputpercap] \ e(N) \ e(r2)))
reg samplemean_numdailies samplemean_outputpercap samplemean_preseligible if year==1900&mainsample
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ _b[samplemean_preseligible] \ _se[samplemean_preseligible] \ _b[samplemean_outputpercap] \ _se[samplemean_outputpercap] \ e(N) \ e(r2)))
reg samplemean_numdailies $samplemean_demolevel if year==1900&mainsample
	matrix TABLE = (nullmat(TABLE), (_b[samplemean_ishare_foreign] \ _se[samplemean_ishare_foreign] \ _b[samplemean_ishare_manuf] \ _se[samplemean_ishare_manuf] \ _b[samplemean_ishare_male] \ _se[samplemean_ishare_male] \ _b[samplemean_ishare_urb] \ _se[samplemean_ishare_urb] \ _b[samplemean_ishare_town] \ _se[samplemean_ishare_town] \ _b[samplemean_ishare_white] \ _se[samplemean_ishare_white] \ . \ . \ . \ . \ e(N) \ e(r2)))
reg samplemean_numdailies samplemean_preseligible $samplemean_demolevel if year==1900&mainsample
	gen demo_sample = e(sample)
	matrix TABLE = (nullmat(TABLE), (_b[samplemean_ishare_foreign] \ _se[samplemean_ishare_foreign] \ _b[samplemean_ishare_manuf] \ _se[samplemean_ishare_manuf] \ _b[samplemean_ishare_male] \ _se[samplemean_ishare_male] \ _b[samplemean_ishare_urb] \ _se[samplemean_ishare_urb] \ _b[samplemean_ishare_town] \ _se[samplemean_ishare_town] \ _b[samplemean_ishare_white] \ _se[samplemean_ishare_white] \ _b[samplemean_preseligible] \ _se[samplemean_preseligible] \ . \ . \ e(N) \ e(r2)))
reg samplemean_numdailies samplemean_preseligible if year==1900&mainsample&demo_sample
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ . \ _b[samplemean_preseligible] \ _se[samplemean_preseligible] \ . \ . \ e(N) \ e(r2)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_xsection>) append

*****************************************************************
* CROSS-SECTIONAL DETERMINANTS OF AFFILIATION
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
gen entry_Rshare = entry_polaff_R/(entry_polaff_R+entry_polaff_D)
gen entrynum = entry_polaff_R+entry_polaff_D
gen correct_pred = (L.presrepshare>=0.5)*entry_Rshare + (L.presrepshare<0.5)*(1-entry_Rshare)

sum entry_Rshare if L.presrepshare>=0.5 [fw=entrynum]
	matrix TABLE = (nullmat(TABLE) \ r(mean))
sum entry_Rshare if L.presrepshare<0.5  [fw=entrynum]
	matrix TABLE = (nullmat(TABLE) \ r(mean))
sum correct_pred [fw=entrynum]
	matrix TABLE = (nullmat(TABLE) \ r(mean))
reg entry_Rshare L.presrepshare [fw=entrynum]

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_xaffil>) append

*****************************************************************
* EXPECTED DIRECTION OF CONFOUNDS
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)

areg D.lpreseligible x_0 if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ . \ e(r2) \ e(N_clust) \ e(N)))
areg D.prestout D.lpreseligible if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[D.lpreseligible] \ _se[D.lpreseligible] \ . \ . \ e(r2) \ e(N_clust) \ e(N)))

areg D_ilog_manufoutputpercap x_0 if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ . \ e(r2) \ e(N_clust) \ e(N)))
areg D.prestout D_ilog_manufoutputpercap if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[D_ilog_manufoutputpercap] \ _se[D_ilog_manufoutputpercap] \ e(r2) \ e(N_clust) \ e(N)))
	
areg D_ilog_manufoutputpercap x_0 D.lpreseligible if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ _b[D.lpreseligible] \ _se[D.lpreseligible] \ . \ . \ e(r2) \ e(N_clust) \ e(N)))
areg D.prestout D_ilog_manufoutputpercap D.lpreseligible if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[D.lpreseligible] \ _se[D.lpreseligible] \ _b[D_ilog_manufoutputpercap] \ _se[D_ilog_manufoutputpercap] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_popinc>) append

*****************************************************************
* FALSIFICATION TEST: SHIFTERS OF POLITICAL INTEREST
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
egen mean_x_0 = mean(x_0), by(styr)

* Confirm effect of state-level competitiveness on turnout
* Note: changing fixed effects to year and cluster level to state
areg D.prestout x_0 D.state_presevenness if mainsample, absorb(year) cluster(state)
	gen D_prestout_pred = _b[D.state_presevenness]*D.state_presevenness
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ . \ . \ . \ e(N) \ e(r2)))

* Check "effect" of newspapers on evenness
areg D.state_presevenness x_0 if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ . \ . \ . \ e(N) \ e(r2)))

* Scaled measure of bias
areg D_prestout_pred x_0 if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ . \ . \ . \ e(N) \ e(r2)))

* Scaled measure of bias controlling for evenness in this county
areg D_prestout_pred x_0 D.presevenness if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ _b[D.presevenness] \ _se[D.presevenness] \ . \ . \ . \ . \ e(N) \ e(r2)))

* Scaled measure of bias: using state mean event
areg D_prestout_pred mean_x_0 if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[mean_x_0] \ _se[mean_x_0] \ . \ . \ e(N) \ e(r2)))
areg D_prestout_pred mean_x_0 D.presevenness if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[D.presevenness] \ _se[D.presevenness] \ _b[mean_x_0] \ _se[mean_x_0] \ . \ . \ e(N) \ e(r2)))

* Check effect of competitiveness on circulation
areg D.circ_elig D.state_presevenness if mainsample_circ, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ _b[D.state_presevenness] \ _se[D.state_presevenness] \ e(N) \ e(r2)))
areg D_lcirc_continue D.state_presevenness if mainsample_circ, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ . \ . \ _b[D.state_presevenness] \ _se[D.state_presevenness] \ e(N) \ e(r2)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_competitive>) append

*****************************************************************
* FALSIFICATION TEST: ELSEWHERE IN STATE
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)

areg D_prestout_state x_0 $demolist $misdemolist if mainsample, absorb(year) cluster(state)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_elsewhere>) append

*****************************************************************
* FIGURE NOT USING MISSING DATA
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

define_event x, changein(numdailies) maxchange($maxchange)
local D_prestout_title "Change in Turnout per Eligible Voter"

* TURNOUT CONTROLLING FOR DEMOS: LAG TERMS ONLY
areg D.prestout x_n10-x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
testparm x_n10-x_n1
plotcoeffs `x_all', label(-40 -36 -32 -28 -24 -20 -16 -14 -8 -4 0) ytitle(`D_prestout_title') xtitle($xtitle) ylabel(-.005(.005).01)
graph export ..\output\text\text_turnout_lagonly.eps, as(eps) replace

* TURNOUT CONTROLLING FOR DEMOS: BRINGING IN PRE-SAMPLE DATA
areg D.prestout x_n10-x_0 $demolist $misdemolist if num_changes>0&year<=$sampleend&abs(D.prestout)<1, absorb(styr) cluster(cnty90)
testparm x_n10-x_n1
plotcoeffs `x_all', label(-40 -36 -32 -28 -24 -20 -16 -14 -8 -4 0) ytitle(`D_prestout_title') xtitle($xtitle) ylabel(-.005(.005).01)
graph export ..\output\text\text_turnout_presample.eps, as(eps) replace

*****************************************************************
* TEST OF READERSHIP MODEL
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x1, changein(numdailies>=1) $ifmulti
define_event x2, changein(numdailies>=2) $ifmulti
define_event x3, changein(numdailies>=3) $ifmulti
define_event xmulti, eventis(multi)

* READERSHIP *
areg D.readshare_hhld x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	predict pred_readshare_hhld
	test x1_0=x2_0=x3_0
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ . \ . \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N) \ . \ .))
areg D.prestout pred_readshare_hhld x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test x2_0 x3_0
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ _b[pred_readshare_hhld] \ _se[pred_readshare_hhld] \ . \ . \ . \ . \ e(r2) \ e(N_clust) \ e(N) \ r(F) \ round(r(p), 0.0001)))
	drop pred_readshare_hhld
	
* CIRCULATION *
areg D.circ_elig x1_0 x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample_circ, absorb(styr) cluster(cnty90)
	predict pred_circ_elig
	test x1_0=x2_0=x3_0
	matrix TABLE = (nullmat(TABLE), (_b[x1_0] \ _se[x1_0] \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ . \ . \ . \ . \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N) \ . \ .))
areg D.prestout pred_circ_elig x2_0 x3_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test x2_0 x3_0
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[x2_0] \ _se[x2_0] \ _b[x3_0] \ _se[x3_0] \ . \ . \ _b[pred_circ_elig] \ _se[pred_circ_elig] \ . \ . \ e(r2) \ e(N_clust) \ e(N) \ r(F) \ round(r(p), 0.0001)))
	drop pred_circ_elig

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_readership>) append
	
*****************************************************************
* ROBUSTNESS TO ADDING EXCLUDED COUNTIES
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

* ADDING EXCLUDED COUNTIES BACK
areg D.prestout x_0 $demolist $misdemolist if num_changes>=0&year>=$samplestart&year<=$sampleend&abs(D.prestout)<1, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ .))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if num_changes>=0&year>=$samplestart&year<=$sampleend&abs(D.prestout)<1, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[xRDdiff_0] \ _se[xRDdiff_0] \ . \ .))
areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if num_changes>=0&year>=$samplestart&year<=$sampleend&abs(D.prestout)<1 & inc_unkwn_h==0 & inc_run_unkwn_h==0 & L.inc_unkwn_h==0 & L.inc_run_unkwn_h==0 & redistricted==0 & inc_cand_h_index==L.inc_cand_h_index & inc_cand_h_index~=0, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0]))

* ADDING OUTLIERS BACK
areg D.prestout x_0 $demolist $misdemolist if num_changes>0&year>=$samplestart&year<=$sampleend, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ .))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if num_changes>0&year>=$samplestart&year<=$sampleend, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[xRDdiff_0] \ _se[xRDdiff_0] \ . \ .))
areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if num_changes>0&year>=$samplestart&year<=$sampleend & inc_unkwn_h==0 & inc_run_unkwn_h==0 & L.inc_unkwn_h==0 & L.inc_run_unkwn_h==0 & redistricted==0 & inc_cand_h_index==L.inc_cand_h_index & inc_cand_h_index~=0, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0]))

* DROPPING 1872
areg D.prestout x_0 $demolist $misdemolist if mainsample & year>1872, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ .))
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample & year>1872, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[xRDdiff_0] \ _se[xRDdiff_0] \ . \ .))
areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc & year>1872, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0]))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_sample>) append

*****************************************************************
* REPSHARE EFFECTS OVER TIME
*****************************************************************
* program to add results to table
cap program drop addtotable
program addtotable
	matrix TABLE = (nullmat(TABLE), (_b[xRDdiff_0_newspaper] \ _se[xRDdiff_0_newspaper] \ _b[xRDdiff_0_radio] \ _se[xRDdiff_0_radio] \ _b[xRDdiff_0_tv] \ _se[xRDdiff_0_tv] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
end

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)
gen xRDdiff_0_newspaper = xRDdiff_0*newspaper
gen xRDdiff_0_radio = xRDdiff_0*radio
gen xRDdiff_0_tv = xRDdiff_0*tv

areg D.circrepshare xRDdiff_0_newspaper xRDdiff_0_radio xRDdiff_0_tv $multicontrol $demolist $misdemolist if mainsample_time_circ, absorb(styr) cluster(cnty90)
	test xRDdiff_0_newspaper=xRDdiff_0_radio=xRDdiff_0_tv
	addtotable
areg D.presrepshare xRDdiff_0_newspaper xRDdiff_0_radio xRDdiff_0_tv  $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test xRDdiff_0_newspaper=xRDdiff_0_radio=xRDdiff_0_tv
	addtotable
areg D.congrepshare xRDdiff_0_newspaper xRDdiff_0_radio xRDdiff_0_tv $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	test xRDdiff_0_newspaper=xRDdiff_0_radio=xRDdiff_0_tv
	addtotable

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_repshare_time>) append

*****************************************************************
* INCUMBENCY EFFECTS OVER TIME
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x1, changein(numdailies>=1) $ifmulti
define_event xmulti, eventis(multi)
gen x1_0_newspaper = x1_0*newspaper
gen x1_0_radio = x1_0*radio
gen x1_0_tv = x1_0*tv

areg D.congshare_inccand x1_0_newspaper x1_0_radio x1_0_tv $multicontrol $demolist $misdemolist if mainsample_time_inc, absorb(styr) cluster(cnty90)
	test x1_0_newspaper=x1_0_radio=x1_0_tv
	matrix TABLE = (nullmat(TABLE), (_b[x1_0_newspaper] \ _se[x1_0_newspaper] \ _b[x1_0_radio] \ _se[x1_0_radio] \ _b[x1_0_tv] \ _se[x1_0_tv] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_incumb_time>) append

******************************************************************
* DFBETA: CHECK FOR INFLUENTIAL OBSERVATIONS
******************************************************************
global maxwindow = 1

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange) window(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
define_event xmulti, eventis(multi) window(1)

* prestout
quietly: reg D.prestout x_0 $demolist $misdemolist i.styr if mainsample
local upper = 2/sqrt(e(N))
local lower = -2/sqrt(e(N))
dfbeta x_0
sum _dfbeta_1, det
histogram _dfbeta_1, xline(`upper' `lower') scheme(s1mono)
graph export ..\output\text\text_dfbeta_prestout.eps, as(eps) replace
areg D.prestout x_0 $demolist $misdemolist if mainsample & abs(_dfbeta_1)<=`upper', absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ . \ . \ . \ .))

* voteshare
quietly: reg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist i.styr if mainsample
local upper = 2/sqrt(e(N))
local lower = -2/sqrt(e(N))
dfbeta xRDdiff_0
sum _dfbeta_2, det
histogram _dfbeta_2, xline(`upper' `lower') scheme(s1mono)
graph export ..\output\text\text_dfbeta_voteshare.eps, as(eps) replace
areg D.presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample & abs(_dfbeta_2)<=`upper', absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ _b[xRDdiff_0] \ _se[xRDdiff_0] \ . \ .))

* incumb
quietly: reg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist i.styr if mainsample_inc
local upper = 2/sqrt(e(N))
local lower = -2/sqrt(e(N))
dfbeta x_0
sum _dfbeta_3, det
histogram _dfbeta_3, xline(`upper' `lower') scheme(s1mono)
graph export ..\output\text\text_dfbeta_incumb.eps, as(eps) replace
areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc & abs(_dfbeta_3)<=`upper', absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), ( . \ . \ . \ . \ _b[x_0] \ _se[x_0]))
	
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_dfbeta>) append

global maxwindow = 10

*****************************************************************
* CALCULATE MAGNITUDES FOR DISCUSSION
*****************************************************************
global maxwindow = 1
cap matrix drop TABLE

* *Numbers Referenced in Discussion of Effects of Newspapers Political Participation* *
use "../temp/voting_cnty_clean.dta", clear
estimates use "..\temp\est_read"
* Percentage point increase in people reading at least one paper
scalar read = _b[x1_0]

estimates use "..\temp\est_pout"
* Percentage point increase in turnout
scalar ptout = _b[x1_0]

* Intent to Treat Calculation
scalar itt55=  scalar(ptout) / scalar(read)

* Average turnout in counties with no newspaper
sum prestout if numdailies==0&mainsample
scalar mean_ptout = r(mean)
scalar mean_p_not_out = (1-r(mean))

* Percentage point increase in readership among non-voters
scalar novote_read= scalar(mean_p_not_out) * scalar(read)

* Persuation rate
scalar persrate = scalar(ptout) / scalar(novote_read)

* Calculate stdev
sum D.prestout if mainsample
scalar stdev_tout = r(sd)
scalar stdev_tout_share = scalar(ptout)/scalar(stdev_tout)

* Calculate RMSE
gen D_prestout = D.prestout
areg D_prestout $multicontrol $demolist $misdemolist if mainsample, absorb(styr)
predict e, resid
sum e
scalar rmse_tout = r(sd)
scalar rmse_tout_share = scalar(ptout)/scalar(rmse_tout)

matrix TABLE = (nullmat(TABLE) \ (scalar(ptout) , scalar(read) , scalar(itt55) , . ) \ (scalar(mean_ptout) , scalar(mean_p_not_out) , scalar(novote_read) , scalar(persrate)) \ (scalar(stdev_tout) , scalar(rmse_tout), scalar(stdev_tout_share),  scalar(rmse_tout_share)) )

* *Numbers Referenced in Discussion of Effects of Newspapers on Party Vote Shares* *
use "../temp/voting_cnty_clean.dta", clear
estimates use "..\temp\est_ciu"

* Calculate upper bound of confidence interval
scalar ciu= _b[xRDdiff] + invttail(e(df_r),.025)*_se[xRDdiff] /*Note invtail is a one-tail test*/

* Unit change when going from one Rep paper to one Dem paper
scalar unitch = 2

* Upper bound of change in republican share of vote
scalar chrepshare = scalar(ciu) * scalar(unitch)

* Intent to treat calculation
scalar itt64 = scalar(chrepshare) /  scalar(read)

* % Dem in population (assuming evenly split market)
scalar pctdem = .5

* Percentage point increase in readership among democrats
scalar dem_read = scalar(pctdem) * scalar(read)

* Persuasion rate
scalar persrate_repshare = scalar(chrepshare) / scalar(dem_read)

* Persuasion rate using average repshare in places Republican papers enter
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti window(1)
sum presrepshare if mainsample & xRDdiff_0<0
scalar pctrep_repexit=r(mean)
sum presrepshare if mainsample & xRDdiff_0>0
scalar pctrep_repentry=r(mean)
scalar pctdem_repentry=1-scalar(pctrep_repentry)
scalar dem_read_repentry=scalar(pctdem_repentry)*scalar(read)
scalar persrate_pctrep_repentry=scalar(chrepshare) / ((scalar(pctdem_repentry))*scalar(read))
display scalar(persrate_pctrep_repentry)

* Persuasion rate using average repshare in places Republican papers are in
sum presrepshare if mainsample & num_polaff_R==1&numdailies==1
scalar pctrep_repin=r(mean)
scalar pctdem_repin=1-scalar(pctrep_repin)
scalar dem_read_repin=scalar(pctdem_repin)*scalar(read)
scalar persrate_pctrep_repin=scalar(chrepshare) / ((scalar(pctdem_repin))*scalar(read))
display scalar(persrate_pctrep_repin)

* Dem readership rate, given various assumptions (see internal_appendix.pdf)
scalar cd = scalar(read) - ((2 * scalar(chrepshare) * ( scalar(mean_ptout) + scalar(ptout)))  /  scalar(itt55)) 

* Republican readership rate, given various assumptions (see internal_appendix.pdf)
scalar cr = scalar(read) + ((2 * scalar(chrepshare) * ( scalar(mean_ptout) + scalar(ptout)))  /  scalar(itt55))

* Rep readership rate divided by Dem readership rate
scalar multiple = scalar(cr) / scalar(cd)

* Calculate stdev
sum D.presrepshare if mainsample
scalar stdev_repshare = r(sd)
scalar stdev_repshare_share = scalar(ciu)/scalar(stdev_repshare)

* Calculate RMSE
gen D_presrepshare = D.presrepshare
areg D_presrepshare $multicontrol $demolist $misdemolist if mainsample, absorb(styr)
predict e, resid
sum e
scalar rmse_repshare = r(sd)
scalar rmse_repshare_share = scalar(ciu)/scalar(rmse_repshare)

matrix TABLE = (nullmat(TABLE) \ (scalar(ciu) , scalar(unitch) , scalar(itt64) , .) \ (scalar(pctdem) , scalar(dem_read) , scalar(persrate_repshare) , .) \ (scalar(cd) , scalar(cr) , scalar(multiple) , .) \ (scalar(stdev_repshare) , scalar(rmse_repshare), scalar(stdev_repshare_share),  scalar(rmse_repshare_share)) )

* Export Table - main values
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:discussion_numbers>) append



* Export Table - robustness of values
matrix TABLE2 = (nullmat(TABLE2) \ (scalar(ciu) , scalar(unitch) , scalar(read)) \ (scalar(pctrep_repin) , scalar(pctrep_repentry) , scalar(pctrep_repexit)) \ (scalar(pctdem_repentry) , scalar(dem_read_repentry) , scalar(persrate_pctrep_repentry)) \ (scalar(pctdem_repin) , scalar(dem_read_repin) , scalar(persrate_pctrep_repin)) )
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE2) format(%20.4f) title(<tab:discussion_numbers_rob>) append
cap matrix drop TABLE2
global maxwindow = 10

*****************************************************************
* TABLE <tab:text_affil>: TURNOUT EFFECT BY AFFILIATION
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event xaffil, changein(num_polaff_R+num_polaff_D) maxchange($maxchange) $ifmulti
define_event xI, changein(num_polaff_I) maxchange($maxchange) $ifmulti
define_event xnone, changein(num_polaff_none) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

areg D.prestout xaffil_0 xI_0 xnone_0 $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)
	test xaffil_0=xI_0
	matrix TABLE = (nullmat(TABLE), (_b[xaffil_0] \ _se[xaffil_0] \ _b[xI_0] \ _se[xI_0]  \ _b[xnone_0] \ _se[xnone_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_affil>) append

*****************************************************************
* TABLE <tab:text_confirmed>: TURNOUT EFFECT BY CONFIRMED ENTRY
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
	gen x_0_up = x_0*(x_0>0)
	gen x_0_down = x_0*(x_0<0)
	* separates cases where number of newspaper goes up according to whether there is an entry "confirmed" by date established during the period
	gen x_0_nentry = x_0_up*(entry_conf>0)
	gen x_0_wentry = x_0_up*(entry_conf==0)
define_event xmulti, eventis(multi)

areg D.prestout x_0_nentry x_0_wentry x_0_down $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	test x_0_nentry = x_0_wentry
	matrix TABLE = (nullmat(TABLE), (_b[x_0_nentry] \ _se[x_0_nentry] \ _b[x_0_wentry] \ _se[x_0_wentry] \ _b[x_0_down] \ _se[x_0_down] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_confirmed>) append

*****************************************************************
* TABLE <tab:text_symmetry>: SYMMETRY OF REPSHARE EFFECT
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event xR, changein(num_polaff_R) maxchange($maxchange) $ifmulti
define_event xD, changein(num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

areg D.presrepshare xR_0 xD_0 $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)
	test xR_0 = -xD_0
	matrix TABLE = (nullmat(TABLE), (_b[xR_0] \ _se[xR_0] \ _b[xD_0] \ _se[xD_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_symmetry>) append

*****************************************************************
* TABLE <tab:text_cumul>: CUMULATIVE REPSHARE EFFECTS
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmain, changein(numdailies) $ifmulti
define_event xmulti, eventis(multi)

gen cumul = 1
areg D.presrepshare xRDdiff_0 xRDdiff_p1 xRDdiff_p2 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_p1]+_b[xRDdiff_p2], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))
areg D.presrepshare xRDdiff_0 xRDdiff_p1 xRDdiff_p2 $multicontrol $demolist $misdemolist if mainsample&pmsa90==.&nummarkets==1, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_p1]+_b[xRDdiff_p2], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))
areg D.congrepshare xRDdiff_0 xRDdiff_p1 xRDdiff_p2 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_p1]+_b[xRDdiff_p2], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))
areg D.congrepshare xRDdiff_0 xRDdiff_p1 xRDdiff_p2 $multicontrol $demolist $misdemolist if mainsample&pmsa90==.&nummarkets==1, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_p1]+_b[xRDdiff_p2], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_cumul>) append

*****************************************************************
* FIGURE {text_incumb}
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

define_event x, changein(numdailies) maxchange($maxchange)
local D_incumb_title "Effect on Incumbent Vote Share"

areg D.congshare_inccand `x_all' if mainsample_inc, absorb(styr) cluster(cnty90)
estimates store rd
areg D.congshare_inccand `x_all' $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
quietly predict_list $demolist, gen(pred_D_incumb_demo)
areg D.pred_D_incumb_demo `x_all' $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
estimates store pl
plotcoeffs `x_all', estimates(rd pl) graphs(err linenose) legend(order(1 3) label(1 actual) label(3 predicted) symxsize(*.2)) label("$xlabel") ytitle(`D_incumb_title') xtitle($xtitle) ylabel(-.02(.01).04)
graph export ..\output\text\text_incumb_a.eps, as(eps) replace

areg D.congshare_inccand `x_all' $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
testparm x_n10-x_n1
plotcoeffs `x_all', label("$xlabel") ytitle(`D_incumb_title') xtitle($xtitle) ylabel(-.02(.01).04)
graph export ..\output\text\text_incumb_b.eps, as(eps) replace

*****************************************************************
* TABLE <tab:text_Incumb_main>
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x, changein(numdailies) maxchange($maxchange)

areg D.congshare_inccand x_0 $multicontrol if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))
	estimates save "..\temp\est_ci_inc", replace
use "../temp/voting_district_clean.dta", clear

global panelvar "stdist"
define_event x, changein(numdailies) maxchange($maxchange)

areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(stdist)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.uncontested_inc x_0 $demolist $misdemolist if mainsample_unc, absorb(styr) cluster(stdist)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))

global panelvar "cnty90"

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_Incumb_main>) append

*****************************************************************
* TABLE <tab:text_incumb_robustness>
*****************************************************************
* MAIN SPECIFICATIONS (FOR COMPARISON)
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE1
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE1 = (nullmat(TABLE1), (_b[x_0] \ _se[x_0]))

* TRUNCATION
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE2
define_event x, changein(numdailies) maxchange(1)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange(1)

areg D.congshare_inccand x_0 $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE2 = (nullmat(TABLE2), (_b[x_0] \ _se[x_0]))
 
* ISOLATED MARKETS
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE3
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc&pmsa90==.&nummarkets==1, absorb(styr) cluster(cnty90)
	matrix TABLE3 = (nullmat(TABLE3), (_b[x_0] \ _se[x_0]))

* CUMULATIVE EFFECTS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE4
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

gen cumul = 1
xi: areg D.congshare_inccand x_0 x_p1 x_p2 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[x_0]+_b[x_p1]+_b[x_p2], post
	matrix TABLE4 = (nullmat(TABLE4), (_b[cumul] \ _se[cumul]))
global maxwindow = 1

* CONTROLLING FOR POLYNOMIAL TRENDS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE5
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

areg D.congshare_inccand x_0 `x_poly$maxpolyorder' $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE5 = (nullmat(TABLE5), (_b[x_0] \ _se[x_0]))
global maxwindow = 1
	
* CONTROLLING FOR RESTRICTED TRENDS
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE6
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

quietly areg D.lpreseligible `x_all' if mainsample_inc&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
quietly predict pred_Dlpreselig_xall_xinc

areg D.congshare_inccand x_0 pred_Dlpreselig_xall_xinc $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE6 = (nullmat(TABLE6), (_b[x_0] \ _se[x_0]))
global maxwindow = 1

* SMALL VS LARGE
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE7
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

foreach indepvar in x{
	gen `indepvar'_0_small = `indepvar'_0*(abs(D.circ_elig)<0.1)
	gen `indepvar'_0_large = `indepvar'_0*(abs(D.circ_elig)>=0.1)
}

areg D.congshare_inccand x_0_small x_0_large $multicontrol $demolist $misdemolist if mainsample_inc&mainsample_circ, absorb(styr) cluster(cnty90)
	matrix TABLE7 = (nullmat(TABLE7), (_b[x_0_small] \ _se[x_0_small] \ _b[x_0_large] \ _se[x_0_large]))

* EXCLUDE COUNTIES WITH BORDER CHANGES
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE8
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc & border_changes == 0, absorb(styr) cluster(cnty90)
	matrix TABLE8 = (nullmat(TABLE8), (_b[x_0] \ _se[x_0]))
	
* EXCLUDE SHORT-LIVED PAPERS
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE9
define_event x, changein(numdailies_y4) maxchange($maxchange)
define_event xmulti, eventis(multi)

areg D.congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)
	matrix TABLE9 = (nullmat(TABLE9), (_b[x_0] \ _se[x_0]))

* JOIN ALL TABLES TO CREATE APPENDIX TABLE
cap matrix drop TABLE
matrix TABLE = (TABLE1 \ TABLE2 \ TABLE3 \ TABLE4 \ TABLE5 \ TABLE6 \ TABLE7 \ TABLE8 \ TABLE9)

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_incumb_robustness>) append

global maxwindow = 10

*****************************************************************
* TABLE <tab:text_governor>: GUBERNATORIAL ELECTIONS
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
define_event x, changein(numdailies) maxchange($maxchange)
define_event xmulti, eventis(multi)

cap matrix drop TABLE
areg D.govtout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))
areg D.govshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc_gov, absorb(styr) cluster(cnty90)
	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_governor>) append

*****************************************************************
* TABLE <tab:text_repnum>: EFFECT OF NUMBER OF NEWSPAPERS ON REPSHARE
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x, changein(numdailies) maxchange($maxchange)

areg D.presrepshare x_0 $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)

	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_repnum>) append

*****************************************************************
* TABLE <tab:text_logspec>: LOG-LOG SPECIFICATION
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE

define_event x, changein(numdailies) maxchange($maxchange)

areg D.lpresttlvote x_0 D.lpreseligible $demolist $misdemolist $multicontrol if mainsample & abs(D.lpresttlvote)<1, absorb(styr) cluster(cnty90)

	matrix TABLE = (nullmat(TABLE), (_b[x_0] \ _se[x_0] \ _b[D.lpreseligible] \ _se[D.lpreseligible] \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_logspec>) append

*****************************************************************
* TABLE <tab:text_presinc>: SPLIT BY INCUMBENCY STATUS OF PRESIDENT
*****************************************************************
/* Cases where the sitting president had completed two terms and therefore
by convention could not run again:
1876 (Grant could not run) 
1896 (Cleveland could not run) 
1904 (McKinley could not run) 
1920 (Wilson could not run)
*/ 

use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
gen twoterm = (year==1876|year==1896|year==1904|year==1920)
define_event x, changein(numdailies) maxchange($maxchange)
gen x_two_0 = x_0*twoterm
gen x_notwo_0 = x_0*(1-twoterm)

areg D.prestout x_two_0 x_notwo_0 $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)
test x_two_0 = x_notwo_0

	matrix TABLE = (nullmat(TABLE), (_b[x_two_0] \ _se[x_two_0] \ _b[x_notwo_0] \ _se[x_notwo_0] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_presinc>) append

* Note: structurally the model calls for an additional control. It does not make much difference in practice but is worth checking.

gen DzN = D.twoterm*L.numdailies
areg D.prestout x_two_0 x_notwo_0 DzN $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)

*****************************************************************
* FIGURE {text_manufout}: Trends in manufacturing output per capita before and after
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
define_event x, changein(numdailies) maxchange($maxchange)

areg D_ilog_manufoutputpercap `x_all' if mainsample&abs(D_ilog_manufoutputpercap)<1, absorb(styr) cluster(cnty90)
plotcoeffs `x_all', label("$xlabel") ytitle(Change in Log(Manufacturing Output per Capita)) xtitle(Years Relative to Change in Number of Newspapers)
graph export ..\output\text\text_manufout.eps, as(eps) replace


*****************************************************************
* TABLE <tab:text_rural>: SPLIT BY RURAL/NOT RURAL
*****************************************************************

use "../temp/voting_cnty_clean.dta", clear

define_event x, changein(numdailies) maxchange($maxchange)
egen mean_ishare_urb = mean(ishare_urb) if mainsample, by(cnty90)
gen rural = mean_ishare_urb==0
gen x_0_rural = x_0*rural
gen x_0_notrural = x_0*(1-rural)

cap matrix drop TABLE
areg D.prestout x_0_notrural x_0_rural $demolist $misdemolist $multicontrol if mainsample, absorb(styr) cluster(cnty90)
test x_0_notrural = x_0_rural

	matrix TABLE = (nullmat(TABLE), (_b[x_0_notrural] \ _se[x_0_notrural] \ _b[x_0_rural] \ _se[x_0_rural] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_rural>) append
lincom x_0_rural-x_0_notrural


*****************************************************************
* TABLE <tab:text_updown_time>: ENTRIES VS EXITS AND CHANGES OVER TIME
*****************************************************************

global maxwindow = 1
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE TABLE2

define_event x1, changein(numdailies>=1) $ifmulti window(1)
gen x1_0_up = x1_0*(x1_0>0)
gen x1_0_down = x1_0*(x1_0<0)
define_event xmulti, eventis(multi) window(1)
foreach type in up down {
	gen x1_0_`type'_newspaper = x1_0_`type'*newspaper
	gen x1_0_`type'_radio = x1_0_`type'*radio
	gen x1_0_`type'_tv = x1_0_`type'*tv
}

areg D.prestout x1_0_up_newspaper x1_0_down_newspaper x1_0_up_radio x1_0_down_radio x1_0_up_tv x1_0_down_tv $multicontrol $demolist $misdemolist if mainsample_time, absorb(styr) cluster(cnty90)
	foreach type in up down {
		test x1_0_`type'_newspaper=x1_0_`type'_radio=x1_0_`type'_tv
		matrix TABLE = (nullmat(TABLE), (_b[x1_0_`type'_newspaper] \ _se[x1_0_`type'_newspaper] \ _b[x1_0_`type'_radio] \ _se[x1_0_`type'_radio] \ _b[x1_0_`type'_tv] \ _se[x1_0_`type'_tv] \ r(F) \ round(r(p), 0.0001) \ e(r2) \ e(N_clust) \ e(N)))
	}

	foreach period in newspaper radio tv {
		test x1_0_up_`period' = x1_0_down_`period'
		matrix TABLE2 = (nullmat(TABLE2) \ r(F) \ round(r(p), 0.0001))
	}

	matrix TABLE = (nullmat(TABLE), (TABLE2 \ . \ . \ . \ . \. ))

matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_updown_time>) append

***********************************************************************
* TABLE <tab:text_before>: CUMULATIVE EFFECTS: BEFORE INSTEAD OF AFTER
***********************************************************************
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)

gen cumul = 1
areg D.prestout x_0 x_n1 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[x_0]+_b[x_n1], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))
areg D.presrepshare xRDdiff_0 xRDdiff_n1 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
	nlcom cumul : _b[xRDdiff_0]+_b[xRDdiff_n1], post
	matrix TABLE = (nullmat(TABLE), (_b[cumul] \ _se[cumul]))
global maxwindow = 1
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_before>) append

***********************************************************************
* TABLE <tab:text_suffrage>: WOMEN'S SUFFRAGE AND NEWSPAPER ENTRY/EXIT
***********************************************************************
global maxwindow = 10
global panelvar "st"
 
use "../temp/suffrage.dta", clear
cap matrix drop TABLE
areg D.numdailies D.suf if sufsample, absorb(year) cluster(sta) 
	matrix TABLE = (nullmat(TABLE), (_b[D.suf] \ _se[D.suf] \ e(r2) \ e(N_clust) \ e(N)))
 
global panelvar "cnty90"
global maxwindow = 1
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_suffrage>) append

***********************************************************************
* TABLE <tab:text_serial>: FORMAL TESTS FOR SERIAL CORRELATION
***********************************************************************
global maxwindow = 10
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)
gen D_prestout = D.prestout
gen D_presrepshare = D.presrepshare
gen D_congshare_inccand = D.congshare_inccand

areg D_prestout x_0 $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
predict eD_prestout, resid
reg eD_prestout L.eD_prestout
	test L.eD_prestout = -0.5
	matrix TABLE = (nullmat(TABLE), (_b[L.eD_prestout] \ _se[L.eD_prestout]))
	
areg D_presrepshare xRDdiff_0 $multicontrol $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
predict eD_presrepshare, resid
reg eD_presrepshare L.eD_presrepshare
	test L.eD_presrepshare = -0.5
	matrix TABLE = (nullmat(TABLE), (_b[L.eD_presrepshare] \ _se[L.eD_presrepshare]))
	
areg D_congshare_inccand x_0 $multicontrol $demolist $misdemolist if mainsample_inc, absorb(styr) cluster(cnty90)	
predict eD_congshare_inccand, resid
reg eD_congshare_inccand L.eD_congshare_inccand
	test L.eD_congshare_inccand = -0.5
	matrix TABLE = (nullmat(TABLE), (_b[L.eD_congshare_inccand] \ _se[L.eD_congshare_inccand]))
		
global maxwindow = 1
matrix_to_txt, saving(..\output\texttables.txt) mat(TABLE) format(%20.4f) title(<tab:text_serial>) append

cap log close
