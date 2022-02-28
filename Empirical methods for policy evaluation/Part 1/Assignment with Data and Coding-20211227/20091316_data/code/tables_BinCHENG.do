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

matrix_to_txt, saving(..\output\tables_BinCHENG.txt) mat(TABLE) format(%20.4f) title(<tab:Turnout_main>) append