/// solution of assignment DEEQA 2020-2021
cd "C:\Users\33678\Example"

cap log close
log using assignment,replace

/*
Author: Bin Cheng

Terms of use data, code and log files:
do not distribute and do not use the data for purposes outside this course
*/


use Data,clear

xtset ID period

*1)	Some descriptives

	sum male ability_math ability_language parentcollege if period==1
		//ability by construction normalized 0-1 (in initial data), balanced in terms of gender, 30% has parent who went to college

	sum school degree 	if period==1
		//91% reaches final stage of high school, 78% finishes high school successfully

	sum school 	dist	if period==2
		//58% goes to college, distance is 8.5km on average

	sum school 			if period==2 & L.degree==1
		//74.8% of high school graduates go to college

	corr male ability_math ability_language parentcollege if period==1
		//ability measures highly correlated, also a bit with parental education
		
	reg school male ability_math ability_language parentcollege if period==1, robust
		//males, low ability and parents who didn't go to college -> more likely to drop out from high school

	reg degree male ability_math ability_language parentcollege if period==1 & school==1, robust
		//similar effects on failing in high school
		
	reg school male ability_math ability_language parentcollege if period==2 & L.degree==1, robust
		//similar effects but larger for college entry
	
	
*2) A static discrete choice model
	
	*high school stage (period 1)
	logit school male ability_math ability_language parentcollege if period==1 
	est store hs_static
	predict hs_static_pr
	
	*state transition: degree
	logit degree male ability_math ability_language parentcollege if period==1 & school==1
	est store degree
	predict degree_pr
	
	*college stage (period 2)	
	logit school male ability_math ability_language parentcollege dist if period==2 & L.degree==1
	est store college
	predict college_pr
	
	*joint probability that specifies college enrollment number
	gen college_pr_static_total=hs_static_pr*degree_pr*college_pr
	
	*Impact on the cost of going to college of not having a parent who went to college in kilometers: 
	di -_b[parentcollege]/_b[dist]

	*counterfactual 1: distance set to 0
	est restore college

	gen college_pr_static_counter1=invlogit(_b[_cons]+_b[male]*male ///
				+_b[ability_math]*ability_math+_b[ability_language]*ability_language ///
				+_b[parentcollege]*parentcollege +_b[dist]*(dist*0)) 						
		/*
		*Alternative (can be useful for bigger models):
		clonevar dist_clone=dist			//copy the variable you want to change
		replace dist=0						//set to counterfactual value
		predict college_pr_static_counter1	//use predict command
		replace dist=dist_clone				//set back to actual value
		drop dist_clone						//remove the copy
		*/
		
	gen college_pr_static_total_counter1=hs_static_pr*degree_pr*college_pr_static_counter1
		
	gen stat_effect1=college_pr_static_total_counter1 ///
					- college_pr_static_total
	

	sum college_pr_static_total college_pr_static_total_counter1 stat_effect1  if period==1
		//decrease in the college enrollment rate of 0.979 %points (or a relative change of 1.682%)

	*counterfactual 2: 50% allowed in college
	gen college_pr_static_total_counter2=hs_static_pr*degree_pr*0.5*college_pr
	
	gen stat_effect2=college_pr_static_total_counter2 ///
					- college_pr_static_total
							
	
	sum college_pr_static_total college_pr_static_total_counter2 stat_effect2 if period==1
		//decrease in the college enrollment rate of 29.098 %points (relative = 50%)
	
		
*3) A dynamic discrete choice model
	*calculate the expected value of behaving optimally in the period 2
	est restore college //results in college period remain valid
	
	gen u_college=_b[_cons]+_b[male]*male ///
				+_b[ability_math]*ability_math+_b[ability_language]*ability_language ///
				+_b[parentcollege]*parentcollege +_b[dist]*dist
		/*
		*Alternative:
		predict u_college, xb
		*/
	
	gen emax=0.5772156649 ///
			+degree_pr*ln(exp(0)+exp(u_college)) ///
			+(1-degree_pr)*ln(exp(0)) //this term is actually 0 so can leave it out	
		/*
		Note: 0.5772156649 is the Euler constant. You can also leave it out, which 
		changes your interpretation of the constant in the utility function (but
		doesn't impact counterfactual simulations). However, do not only include it 
		after "degree_pr" as you also receive it if you do not get the degree.
		*/
	
	*high school stage
	constraint 1 emax=0.95
	logit school male ability_math ability_language parentcollege emax if period==1 , constraint(1)
	est store hs_dynamic
	predict hs_dynamic_pr
	
	gen college_pr_dynamic_total=hs_dynamic_pr*degree_pr*college_pr if period==2
	sum school college_pr_dynamic_total  if period==2
	

	*counterfactual: distance decreases by 100%
		*predictions of college and the emax (=expected value of behaving optimally)
		est restore college
		gen u_college_counter1=_b[_cons]+_b[male]*male ///
					+_b[ability_math]*ability_math+_b[ability_language]*ability_language ///
					+_b[parentcollege]*parentcollege +_b[dist]*(dist*0)
		gen emax_counter1=0.5772156649 ///
			+degree_pr*ln(exp(0)+exp(u_college_counter1)) ///
			+(1-degree_pr)*ln(exp(0))	
		
		gen college_pr_dynamic_counter1=invlogit(u_college_counter1)
		
		est restore hs_dynamic
		gen v_counter1=_b[_cons]+_b[male]*male ///
					+_b[ability_math]*ability_math+_b[ability_language]*ability_language ///
					+_b[parentcollege]*parentcollege +_b[emax]*emax_counter1
				
		gen hs_dynamic_pr_counter1=invlogit(v_counter)
		
		gen college_pr_dyn_total_counter1=hs_dynamic_pr_counter1*degree_pr*college_pr_dynamic_counter1 if period==2
						
		gen dynamic_effect1=college_pr_dyn_total_counter1-college_pr_dynamic_total
		
		sum college_pr_dynamic_total college_pr_dyn_total_counter1 dynamic_effect1	
			//increase in the college enrollment rate of 1.080 %points (1.854% relative)
	
		gen dynamic_effect_hs1=hs_dynamic_pr_counter1-hs_dynamic_pr
		
		sum hs_dynamic_pr hs_dynamic_pr_counter1 dynamic_effect_hs1 if period==1
			//increase high school of 0.205 %points (relative 0.227%)
			
	*counterfactual: 50% allowed in college
		*choice probs in counterfactuals	
		est restore college
		gen emax_counter2=0.5772156649 ///
			+0.5*degree_pr*ln(exp(0)+exp(u_college)) ///
			+(1-0.5*degree_pr)*ln(exp(0))	
		
		est restore hs_dynamic
		gen v_counter2=_b[_cons]+_b[male]*male ///
					+_b[ability_math]*ability_math+_b[ability_language]*ability_language ///
					+_b[parentcollege]*parentcollege +_b[emax]*emax_counter2
					
		gen hs_dynamic_pr_counter2=invlogit(v_counter2)
		
		gen college_pr_dyn_total_counter2=hs_dynamic_pr_counter2*degree_pr*0.5*college_pr if period==2
						
		gen dynamic_effect2=college_pr_dyn_total_counter2-college_pr_dynamic_total
		
		sum college_pr_dynamic_total college_pr_dyn_total_counter2 dynamic_effect2
				//decrease in the college enrollment rate of 29.902 %points (relative = 51.331% > 50%)
				
		gen dynamic_effect_hs2=hs_dynamic_pr_counter2-hs_dynamic_pr
		
		sum hs_dynamic_pr hs_dynamic_pr_counter2 dynamic_effect_hs2 if period==1
			//decrease of 2.470 %points (relative 2.728%)
		
*Extra (not asked on assignment): the difference in counterfactual predictions 
	*counterfactual 1
	gen bias1=stat_effect1-dynamic_effect1
		scatter bias1 college_pr_dynamic_total 
	
	*counterfactual 2
	gen bias2=stat_effect2-dynamic_effect2
		scatter bias2 college_pr_dynamic_total
	
	*some explanation
	scatter dynamic_effect_hs1 hs_dynamic_pr 	
	scatter dynamic_effect_hs2 hs_dynamic_pr 	
			
		/*At the extremes, dynamic incentives are not important (e.g. because high ability students always go to 
		school, low ability students never go). However, there is a group of more "marginal" students in high school that 
		would see their willingness to go to high school changed by the policy which has repercussions on the total 
		college enrollment rate. 
		Note that for a small fraction of students the "biases" seem to go in the opposite direction but this can be 
		explained by the fact that we are using different parametric forms (a linear utility function in a dynamic model 
		creates nonlinear conditional value functions).
		*/
	
			
*4) CCP estimator
	*recover the CCP 
	est restore college
	predict ccp_college
	gen ccp_dropout=1-ccp_college

	gen emax_ccp=0.5772156649 ///
			+degree_pr*(-ln(ccp_dropout)) ///
			+(1-degree_pr)*ln(exp(0)) 
	
	*high school stage
	constraint 1 emax_ccp=0.95
	logit school male ability_math ability_language parentcollege emax_ccp if period==1 , constraint(1)
	est store hs_dynamic_ccp

	est table hs_dynamic hs_dynamic_ccp
		
	/*
	Note that "emax" and "max_ccp" are identical here so estimating the dynamic model using
	this will give you identical estimates.
	In general, the model for the second period is more complicated and you
	don't want to solve it (that's the purpose of the CCP estimator). But
	note that you don't really need to know the model, you only need to know
	"ccp_dropout". Let's be agnostic about the model and get an estimate of
	the CCP from a more "flexible" logit:
	*/
	
	logit school male ability_math ability_language parentcollege dist ///
					 c.male#c.(ability_math ability_language parentcollege dist) ///
					 c.ability_math#c.(ability_math ability_language parentcollege dist) ///
					 c.ability_language#c.(ability_math ability_language parentcollege dist) ///
					 c.parentcollege#c.(ability_math ability_language dist) ///
					 c.dist#c.(ability_math ability_language parentcollege dist) ///
					if period==2 & L.degree==1
	est store ccp_flex
	predict ccp_college_flex
	gen ccp_dropout_flex=1-ccp_college_flex
	
	gen emax_ccp_flex=0.5772156649 ///
			+degree_pr*(-ln(ccp_dropout_flex)) ///
			+(1-degree_pr)*ln(exp(0)) 
			
	scatter emax emax_ccp_flex 
		//more flexible model but actually differences are small
	
	*now try to estimate the dynamic model again
		*high school stage
	constraint 1 emax_ccp_flex=0.95
	logit school male ability_math ability_language parentcollege emax_ccp_flex if period==1 , constraint(1)
	est store hs_dynamic_ccp_flex

	est table hs_dynamic hs_dynamic_ccp hs_dynamic_ccp_flex
		
	/*
	Conclusion: almost but not exactly identical results because we are approximating the CCP. In general
	estimates will be consistent if the prediction for "ccp_college_flex" is consistent.
	*/
	
*5) what if more periods?

	/*
	The estimation results we got from "ccp_flex" should not be interpreted as utility estimates, rather
	they are predictors for "ccp_college_flex", the only component we needed. This still holds if the model 
	continuous after period 2.
	*/
	

log close
		

		
		