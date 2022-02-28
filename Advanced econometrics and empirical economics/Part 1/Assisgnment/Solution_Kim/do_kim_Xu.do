cd "C:\Users\alienware\Desktop\Deeqa_EEE\Advanced econometrics and empirical economics\Part 1\Assisgnment\Solution_Kim"
use "data_assignment.dta", clear
log using "log_kim_Xu"

*** labelling
label define aa 1 "Male" 0 "Female"
label define bb 1 "c parent" 0 "no c parent"
label define dd 1 "attend" 0 "drop"
label define ee 1 "h sch grad" 0 "no h sch grad"

label values male aa
label values parentcollege bb
label values school dd
label values degree ee

label variable male "male"
label variable parentcollege ""
label variable ability_math "math"
label variable ability_language "language"





************************ Q1

** period 1
tab male school if period==1 
corr school male if period==1
corr school parentcollege if period==1
corr school ability_math if period==1
corr school ability_language if period==1


** corrtex package
//corrtex school male if period==1, file(cor1) replace
//corrtex school parentcollege if period==1, file(cor2) replace
//corrtex school ability_math if period==1, file(cor3) replace
//corrtex school ability_language if period==1, file(cor4) replace
//corrtex school dist if period==1, file(cor5) replace



** period 2
tab male school if period==2 & degree==1
corr school male if period==2 & degree==1 
corr school parentcollege if period==2 & degree==1
corr school ability_math if period==2 & degree==1
corr school ability_language if period==2 & degree==1
corr school dist if period==2 & degree==1

//corrtex school male if period==2 & degree==1, file(cor7) replace
//corrtex school parentcollege if period==2 & degree==1, file(cor8) replace
//corrtex school ability_math if period==2 & degree==1, file(cor9) replace
//corrtex school ability_language if period==2 & degree==1, file(cor10) replace
//corrtex schoo dist if period==2 & degree==1, file(cor11) replace






************ Q2 c

** period 1
mlogit school male parentcollege ability_math ability_language if period==1, base(0)
//eststo: mlogit school male parentcollege ability_math ability_language if period==1, base(0)


** period 2
mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 
//eststo: mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 

** probability getting a high school degree.
mlogit degree male parentcollege ability_math ability_language if period==1 & school==1, base(0) 
gen xb_degree = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math +  [1]_b[ability_language]*ability_language if period==1
gen p_degree = invlogit(xb_degree)
//eststo: mlogit degree male parentcollege ability_math ability_language if period==1 & school==1, base(0) 

*OUTPUT Q2C
//esttab using Q2c.tex, se replace
//eststo clear

*proba going to school in period 1 by static model
quietly mlogit school male parentcollege ability_math ability_language if period==1, base(0)
gen xb_period1 = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language if period==1
gen p_period1 = invlogit(xb_period1)

*proba going to school in period 2 by static model
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0)  
gen xb_period2 = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist
gen p_period2 = invlogit(xb_period2)




************ Q2 e Counterfactual 1, distance zero
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0)   
gen xb_cf1 = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*0 if degree==1
gen p_cf1= invlogit(xb_cf1)
gen xb_cf1_p2_unc = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*0 if period==1
gen p_cf1_p2_unc= invlogit(xb_cf1_p2_unc)

************ Q2 e Counterfactual 2, half utility
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 
gen xb_cf2 = 0.5*([1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist) if degree==1
gen p_cf2 = invlogit(xb_cf2)
gen xb_cf2_p2_unc = 0.5*([1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist) if period==1
gen p_cf2_p2_unc= invlogit(xb_cf2_p2_unc)

**** result 
sum school p_period2 p_cf1 p_cf2 if period==2 & degree==1

gen p_unc = p_cf1*p_degree*p_period2
gen p_cf1_unc = p_cf1*p_degree*p_cf1_p2_unc
gen p_cf2_unc = p_cf1*p_degree*p_cf2_p2_unc


* static for Table 12. static, static c1, static c2
sum p_unc p_cf1_unc p_cf2_unc if period==1













****** Dynamic model

************ Q3 c
*decision period 1
gen deci_1 = 1  if school==1 & period==1
replace deci_1 = 0 if school==0 & period==1
replace deci_1 = -1 if deci_1==.


*********** try 2 step ********************
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 
gen xb_period2_dym = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist if period==1


capture program drop dymodel2
 program dymodel2
  version 14
  args lnf xb1
  local y1 "$ML_y1"
  tempvar c 
  quietly gen double `c' = `xb1'+0.95*p_degree*(log(1+exp(xb_period2_dym)))
  
  // contribution to log-likelihood for each outcome y 
  quietly replace `lnf' = ln(  invlogit(`c')) if `y1'==1 
  quietly replace `lnf' = ln(1-invlogit(`c')) if `y1'==0 
  quietly replace `lnf' = 0 if `y1'==-1 
end
ml model lf dymodel2 (xb1: deci_1= male parentcollege ability_math ability_language) if period==1
ml search
eststo: ml maximize
//esttab using Q3c.tex, se replace
//eststo clear
 
*proba going to school in period 1 by dym model
gen xb_period1_dym = [xb1]_b[_cons] + [xb1]_b[male]*male + [xb1]_b[parentcollege]*parentcollege + [xb1]_b[ability_math]*ability_math + [xb1]_b[ability_language]*ability_language if period==1 
gen xb_period1_dym_sum = xb_period1_dym + 0.95*p_degree*xb_period2_dym
gen p_period1_dym = invlogit(xb_period1_dym_sum)

*proba going to school in period 2 by dym model = use static model, p_period2
gen p_period2_dym = invlogit(xb_period2_dym)

gen p_dym_current = p_period1_dym*p_degree*p_period2_dym

************ Q3 e Counterfactual 1, distance zero
*period 2, SAME AS Q2


*period 1	
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0)
gen xb_cf1_dym_2x = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*0 if period==1
gen xb_cf1_dym_1 = xb_period1_dym +0.95*p_degree*xb_cf1_dym_2x
gen p_cf1_dym_1 = invlogit(xb_cf1_dym_1)
gen p_cf1_dym_2 = invlogit(xb_cf1_dym_2x)

gen p_dym_dist0 = p_cf1_dym_1*p_degree*p_cf1_dym_2 

* table 12 dynamic predicted and C1
sum p_dym_current p_dym_dist0 

************ Q3 f Counterfactual 2, half utility
*period 2 SAME AS Q2


*period 1
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0)
gen xb_cf2_dym_2x = 0.5*([1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist) if period==1
gen xb_cf2_dym_1 = xb_period1_dym + 0.95*p_degree*xb_cf2_dym_2x
gen p_cf2_dym_2 = invlogit(xb_cf2_dym_2x)
gen p_cf2_dym_1 = invlogit(xb_cf2_dym_1)

gen p_dym_50 = p_cf2_dym_1*p_degree*p_cf2_dym_2 


*predicted proba going to school in period 1, counter 1,2 by dym model





gen ph_current = p_period1_dym*p_degree
gen ph_dist0 = p_cf1_dym_1*p_degre
gen ph_50 = p_cf2_dym_1*p_degree

*table 10 predicted counter 1 counter 2
sum p_period1_dym p_cf1_dym_1 p_cf2_dym_1 if period==1

*table11 predicted counter 1 counter 2
sum ph_current ph_dist0 ph_50 if period==1

*table12 predicted c1 c2
sum p_dym_current p_dym_dist0 p_dym_50 if period==1

*Actual proportions table 10, 11, 12
sum school if period ==1
sum degree if period ==1
sum school if period ==2

gen phs_static = p_period1
gen pgrad_static = p_period1*p_degree

mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 
predict p_period2_static
gen pcollege_static = p_period1*p_degree*p_period2_static

*static estimations table 10, 11
sum phs_static pgrad_static pcollege_static



************ Q4. Q5


** CCP
quietly mlogit school male parentcollege ability_math ability_language dist if period==2 & degree==1, base(0) 
gen ccpfit = [1]_b[_cons] + [1]_b[male]*male + [1]_b[parentcollege]*parentcollege + [1]_b[ability_math]*ability_math + [1]_b[ability_language]*ability_language + [1]_b[dist]*dist if period==1
gen ccp = invlogit(ccpfit) // proba going to college

* set j*=0, hence uij*2=0, and corresponding ccp for j*=0 is 1-ccp


capture program drop ccpmodel
 program ccpmodel
  version 14
  args lnf xb
  local y1 "$ML_y1"
  tempvar c 
  quietly gen double `c' = `xb' + 0.95*p_degree*(-log(1-ccp))
  
  // contribution to log-likelihood for each outcome y 
  quietly replace `lnf' = ln(  invlogit(`c')) if `y1'==1 
  quietly replace `lnf' = ln(1-invlogit(`c')) if `y1'==0 
  quietly replace `lnf' = 0 if `y1'==-1 
end
 ml model lf ccpmodel (xb: deci_1= male parentcollege ability_math ability_language) if period==1
 ml search
eststo: ml maximize
//esttab using Q45.tex, se replace
//eststo clear

log close
