	// Environment

	cd "C:\Users\33678\Desktop\Deeqa_EEE\Advanced econometrics and empirical economics\Part 1\Assisgnment\Solution"
	use "data_assignment.dta", clear

	log using "Log_BinCHENG"


	xtset period
	// Question 1

	***************Period 1******************
		sum male ability_math ability_language parentcollege school degree  if period==1

		corr male ability_math ability_language parentcollege school degree if period==1
			
		reg school male ability_math ability_language parentcollege if period==1

		reg degree male ability_math ability_language parentcollege if period==1 & school==1

	***************Period 2******************

		sum male ability_math ability_language parentcollege school degree 	dist	if period==2
		corr male ability_math ability_language parentcollege school degree  dist	if period==2

		reg school male ability_math ability_language parentcollege dist if period==2 & L.degree==1

	// Question 2
	
	

			
			
