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
set maxvar 32767
set more off
adopath + "..\external\"
loadglob using "input_param.txt"



*****************************************************************
* FIGURE {fig_circ_repshare} EFFECTS ON REPUBLICAN READERSHIP SHARE
*****************************************************************
use "../temp/voting_cnty_clean.dta", clear
cap matrix drop TABLE
define_event xRDdiff, changein(num_polaff_R-num_polaff_D) maxchange($maxchange) $ifmulti

areg D.circrepshare `xRDdiff_all' if mainsample_circ, absorb(styr) cluster(cnty90)
	plotcoeffs `xRDdiff_all', label("$xlabel") ytitle(Change in Rep - Dem Readership Share) xtitle(Years Relative to Change in (#Rep-#Dem))
	graph export ..\output\figures\fig_circ_repshare_BinCHENG.eps, as(eps) replace

	
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
graph export ..\output\figures\fig_turnout_a_BinCHENG.eps, as(eps) replace

* TURNOUT CONTROLLING FOR DEMOS
areg D.prestout `x_all' $demolist $misdemolist if mainsample, absorb(styr) cluster(cnty90)
testparm x_n10-x_n1
plotcoeffs `x_all', label("$xlabel") ytitle(`D_prestout_title') xtitle(Years Relative to Change in Number of Newspapers) ylabel(-.005(.005).01)
graph export ..\output\figures\fig_turnout_b_BinCHENG.eps, as(eps) replace