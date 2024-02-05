/*Effect modification heterogeneity tests and supplemental tables
Manuscript "Simulating the Impact of Green Space Exposure on Cardiometabolic Biomarkers in a
Diverse Population Living in San Diego, California: A G-Computation Application"
Version 1.31.2024
Contact: Anais Teyton, ateyton@ucsd.edu


Example provided is for Blood Glucose Level Outcome
Can be substituted for any of the other 6 biomarkers
*/

//BLOOD GLUCOSE LEVEL 
import excel "ENTER DATASET HERE", sheet("Sheet 1") firstrow clear

keep if substr(V1, 1, 2) == "RD"
destring se, replace

gen result = string(mean, "%9.2f") + " (" + string(ll, "%9.2f") + ", " + string(ul, "%9.2f") + ")"
export excel using "SAVE SUPPLEMENTAL MATERIAL TABLE HERE", firstrow(variables) replace

//Sex
levelsof V1, local(levels) 
 foreach l of local levels {
 	preserve
keep if name=="bootstrap_gluc_male" | name=="bootstrap_gluc_female"
keep if V1=="`l'"
ta V1
meta set mean se, civartolerance(.05)
meta summarize
restore
}

//Hispanic
levelsof V1, local(levels) 
 foreach l of local levels {
 	preserve
keep if name=="bootstrap_gluc_hisp" | name=="bootstrap_gluc_nothisp"
keep if V1=="`l'"
ta V1
meta set mean se, civartolerance(.05)
meta summarize
restore
}

//Income
levelsof V1, local(levels) 
 foreach l of local levels {
 	preserve
keep if name=="bootstrap_gluc_inc1" | name=="bootstrap_gluc_inc2" | name=="bootstrap_gluc_inc3"
keep if V1=="`l'"
ta V1
meta set mean se, civartolerance(.05)
meta summarize
restore
}

//Age
levelsof V1, local(levels) 
 foreach l of local levels {
 	preserve
keep if name=="bootstrap_gluc_ageover65" | name=="bootstrap_gluc_ageunder65"
keep if V1=="`l'"
ta V1
meta set mean se, civartolerance(.05)
meta summarize
restore
}
