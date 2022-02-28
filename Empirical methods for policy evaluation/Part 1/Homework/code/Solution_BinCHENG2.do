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

	log using "Log_BinCHENG6",replace
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

matrix_to_txt, saving(..\output2\tables.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main>) append



*****************************************************************
* TABLE <tab:Turnout_main_replication>
*****************************************************************

use "http://fmwww.bc.edu/repec/bocode/t/turnout_dailies_1868-1928.dta", clear

* program to add results to table
cap program drop addtotable_LATE
program addtotable_LATE
	matrix TABLE = (nullmat(TABLE), ( e(b_LATE) \ e(se_LATE) \ e(N)))
end

cap matrix drop TABLE

sum pres_turnout numdailies // check that we have the same output as on the paper

* 1868 and 1872 elections
gen G1872=(fd_numdailies>0) if (year==1872) & fd_numdailies!=. & fd_numdailies>= 0 & sample==1
sort cnty90 year
replace G1872=G1872[_n+1] if cnty90==cnty90[_n+1] & year==1868

***DID TC CIC without controls
fuzzydid pres_turnout G1872 year numdailies, did tc cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
***LQTE
gen numdailies_bin = (numdailies >= 1)
fuzzydid pres_turnout G1872 year numdailies_bin, lqte breps(200) cluster(cnty90)

*****************************************************************
* TABLE <tab:Turnout_main_replication_DID_1868+1872>
*****************************************************************

fuzzydid readshare_hhld G1872 year numdailies, did newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G1872 year numdailies, did newcateg(0 1 2 45)  qualitative(st1-st48 ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G1872 year numdailies, did newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G1872 year numdailies, did newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G1872 year numdailies, did newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE

matrix_to_txt, saving(..\output2\tables_did6872.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_did6872>) append

*****************************************************************
* TABLE <tab:Turnout_main_replication_TC_1868+1872>
*****************************************************************

cap matrix drop TABLE

fuzzydid readshare_hhld G1872 year numdailies, tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G1872 year numdailies, tc newcateg(0 1 2 45)  qualitative(st1-st48 ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G1872 year numdailies, tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G1872 year numdailies, tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G1872 year numdailies, tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
matrix_to_txt, saving(..\output2\tables_tc6872.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_tc6872>) append

*****************************************************************
* TABLE <tab:Turnout_main_replication_CIC_1868+1872>
*****************************************************************

cap matrix drop TABLE

fuzzydid readshare_hhld G1872 year numdailies, cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G1872 year numdailies,  cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G1872 year numdailies,  cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G1872 year numdailies,  cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE

matrix_to_txt, saving(..\output2\tables_cic6872.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_cic6872>) append

*****************************************************************
* TABLE <tab:Turnout_main_replication_DID_Full>
*****************************************************************

cap matrix drop TABLE

* Full sample of elections

sort cnty90 year
by cnty90 year: egen mean_D = mean(numdailies)
by cnty90: g lag_mean_D = mean_D[_n-1] if cnty90==cnty90[_n-1]&year-4==year[_n-1]
gen G_T = sign(mean_D - lag_mean_D) if sample==1
gen G_Tplus1 = G_T[_n+1] if cnty90==cnty90[_n+1]&year+4==year[_n+1]

***DID TC CIC without controls
fuzzydid pres_turnout G_T G_Tplus1 year numdailies, did tc cic newcateg(0 1 2 45)  breps(100) cluster(cnty90)



***DID TC with demolist and misdemolist controls and state specific trends

fuzzydid readshare_hhld G_T G_Tplus1 year numdailies, did  newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
***DID without controls but with state specific trends
fuzzydid pres_turnout G_T G_Tplus1 year numdailies, did newcateg(0 1 2 45) qualitative(st1-st48) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G_T G_Tplus1 year numdailies, did  newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G_T G_Tplus1 year numdailies, did  newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G_T G_Tplus1 year numdailies, did  newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE


matrix_to_txt, saving(..\output2\tables_didfull.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_did_full>) append

*****************************************************************
* TABLE <tab:Turnout_main_replication_TC_Full>
*****************************************************************

cap matrix drop TABLE

fuzzydid readshare_hhld G_T G_Tplus1 year numdailies,  tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
***TC without controls but with state specific trends
fuzzydid pres_turnout G_T G_Tplus1 year numdailies, tc newcateg(0 1 2 45) qualitative(st1-st48) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid pres_turnout G_T G_Tplus1 year numdailies,  tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G_T G_Tplus1 year numdailies,  tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G_T G_Tplus1 year numdailies,  tc newcateg(0 1 2 45) continuous($demolist) qualitative(st1-st48  $misdemolist ) breps(100) cluster(cnty90)
addtotable_LATE

matrix_to_txt, saving(..\output2\tables_tcfull.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_tc_full>) append

*****************************************************************
* TABLE <tab:Turnout_main_replication_CIC_Full>
*****************************************************************

cap matrix drop TABLE

fuzzydid readshare_hhld G_T G_Tplus1 year numdailies,  cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE
***TC without controls but with state specific trends
fuzzydid pres_turnout G_T G_Tplus1 year numdailies, cic newcateg(0 1 2 45) breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid congtout G_T G_Tplus1 year numdailies,  cic newcateg(0 1 2 45)  breps(100) cluster(cnty90)
addtotable_LATE
fuzzydid D_congtout_long G_T G_Tplus1 year numdailies,  cic newcateg(0 1 2 45)  breps(100) cluster(cnty90)
addtotable_LATE

matrix_to_txt, saving(..\output2\tables_cicfull.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main_cic_full>) append

*** Placebo Wald-DID to assess Ass 4 and 5

xtset cnty90 year
gen fd_numdailies_l1=l4.fd_numdailies
gen pres_turnout_l1=l4.pres_turnout
sort cnty90 year
g G_T_placebo = sign(mean_D - lag_mean_D) if sample==1&fd_numdailies_l1==0
g G_Tplus1_placebo = G_T_placebo[_n+1] if cnty90==cnty90[_n+1]&year+4==year[_n+1]

fuzzydid pres_turnout_l1 G_T_placebo G_Tplus1_placebo year numdailies, did tc newcateg(0 1 2 45) qualitative(st1-st48) breps(100) cluster(cnty90)


log close