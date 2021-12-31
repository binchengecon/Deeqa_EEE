/**********************************************************
 *
 * FIGURES.DO: MAKES FIGURES FOR VOTING PAPER
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

*****************************************************************
* FIGURE {fig_sum}
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

* NUMBER OF NEWSPAPERS BY YEAR
graph bar (sum) numdailies if year >= $samplestart, over(year, label(angle(vertical))) ytitle( "Number of Daily Newspapers") bar(1,color(black)) scheme(s1mono)
graph export ..\output\figures\fig_sum_a.eps, as(eps) replace

* NUMBER OF COMPETITORS BY YEAR
tab numcomp, gen(znumcomp)
graph bar (sum) znumcomp4 znumcomp3 znumcomp2 if year>=$samplestart, over(year, label(angle(vertical))) stack ytitle( "Number of Counties") legend(symxsize(*0.25) order( 1 "3+ newspapers" 2 "2 newspapers" 3 "1 newspaper")) scheme(s1mono) bar(1, color(black)) bar(2, color(gs6)) bar(3, color(gs12))
graph export ..\output\figures\fig_sum_b.eps, as(eps) replace

* CIRCULATION PER ELIGIBLE VOTER
graph bar (mean) circ_elig [w=preseligible] if year >= $samplestart & miss_circ==0, over(year, label(angle(vertical))) ytitle( "Circulation per Eligible Voter") bar(1,color(black)) scheme(s1mono)
graph export ..\output\figures\fig_sum_c.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_events}: NUMBER OF EVENTS BY YEAR
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

* hard-coding maxchange because we only want to compute number gaining and number losing
define_event x, changein(numdailies) maxchange(1)
tab x_0, gen(zx_0)
graph bar (sum) zx_01 zx_03 if year >= $samplestart, over(year, label(angle(vertical))) ytitle( "Number of Counties") legend(symxsize(*0.25) order( 1 "Losing Papers" 2 "Gaining Papers")) bar(1, color(gs8)) bar(2, color(black)) scheme(s1mono)
graph export ..\output\figures\fig_events.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_pop}
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
define_event x, changein(numdailies) maxchange($maxchange)

* RAW
areg D.lpreseligible `x_all' if mainsample&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
plotcoeffs `x_all', label("$xlabel") ytitle(Change in Log(Voting-Eligible Population)) xtitle(Years Relative to Change in Number of Newspapers)
graph export ..\output\figures\fig_pop_a.eps, as(eps) replace

* WHITENED
areg lpreseligible L(1/10).lpreseligible x_p1-x_p10 if mainsample&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
predict pred_lpreseligible_t1
gen innov_lpreseligible = lpreseligible-pred_lpreseligible_t1
areg innov_lpreseligible `x_all' if mainsample&abs(D.lpreseligible)<1, absorb(styr) cluster(cnty90)
plotcoeffs `x_all', label("$xlabel") ytitle(Innovation to Log(Voting-Eligible Population)) xtitle(Years Relative to Change in Number of Newspapers) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_pop_b.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_circ} EFFECTS ON READERSHIP PER ELIGIBLE
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event x, changein(numdailies) maxchange($maxchange)

areg D.readshare_hhld `x_all' if mainsample_circ, absorb(styr) cluster(cnty90)
	plotcoeffs `x_all', label("$xlabel") ytitle(Change in Readership per Eligible Voter) xtitle(Years Relative to Change in Number of Newspapers)
	graph export ..\output\figures\fig_circ.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_turnout}
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

define_event x, changein(numdailies) maxchange($maxchange)
local D_prestout_title "Change in Turnout per Eligible Voter"

* TURNOUT & TURNOUT PREDICTED FROM DEMOS
areg D.prestout `x_all' if mainsample, absorb(styr) cluster(cnty90)
estimates store rd
areg D.prestout `x_all' $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
quietly predict_list $demolist, gen(pred_D_prestout_demo)
areg pred_D_prestout_demo `x_all' $misdemolist if mainsample, absorb(styr) cluster(cnty90)
estimates store pl
plotcoeffs `x_all', estimates(rd pl) graphs(err linenose) legend(order(1 3) label(1 actual) label(3 predicted) symxsize(*.2)) label("$xlabel") ytitle(`D_prestout_title') xtitle(Years Relative to Change in Number of Newspapers) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_turnout_a.eps, as(eps) replace

* TURNOUT CONTROLLING FOR DEMOS
areg D.prestout `x_all' $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
testparm x_n10-x_n1
plotcoeffs `x_all', label("$xlabel") ytitle(`D_prestout_title') xtitle(Years Relative to Change in Number of Newspapers) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_turnout_b.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_circ_repshare} EFFECTS ON REPUBLICAN READERSHIP SHARE
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti

areg D.circrepshare `xRDdiff_all' if mainsample_circ, absorb(styr) cluster(cnty90)
	plotcoeffs `xRDdiff_all', label("$xlabel") ytitle(Change in Rep - Dem Readership Share) xtitle(Years Relative to Change in (#Rep-#Dem))
	graph export ..\output\figures\fig_circ_repshare.eps, as(eps) replace

*****************************************************************
* FIGURE {fig_repshare}
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear

define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti
define_event xmulti, eventis(multi)
local D_presrepshare_title "Change in Republican Vote Share"

* REPSHARE & REPSHARE PREDICTED FROM STATE-YEAR FIXED EFFECTS
areg D.presrepshare `xRDdiff_all' $multicontrol if mainsample, absorb(year) cluster(cnty90)
estimates store rd
gen D_presrepshare = D.presrepshare
areg D_presrepshare `xRDdiff_all' $multicontrol if mainsample, absorb(styr) cluster(cnty90)
quietly predict pred_styr, d
areg pred_styr `xRDdiff_all' $multicontrol if mainsample, absorb(year) cluster(cnty90)
estimates store pl
plotcoeffs `xRDdiff_all', estimates(rd pl) graphs(err linenose) legend(order(1 3) label(1 actual) label(3 predicted) symxsize(*.2)) label("$xlabel") ytitle(`D_presrepshare_title') xtitle(Years Relative to Change in (#Rep-#Dem)) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_repshare_a.eps, as(eps) replace

* REPSHARE CONTROLLING FOR STATE-YEAR FIXED EFFECTS
areg D.presrepshare `xRDdiff_all' $multicontrol if mainsample, absorb(styr) cluster(cnty90)
testparm xRDdiff_n10-xRDdiff_n1
plotcoeffs `xRDdiff_all', label("$xlabel") ytitle(`D_presrepshare_title') xtitle(Years Relative to Change in (#Rep-#Dem)) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_repshare_b.eps, as(eps) replace

