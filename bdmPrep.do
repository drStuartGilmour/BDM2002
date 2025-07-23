********************
*
* Code to reproduce and analyze data for BDM2002
*	By Stuart Gilmour, 2025/7/16
*******************

****************** Import the data **********************
* be warned this takes a long time!
use "http://www.nber.org/morg/annual/morg79.dta", clear
forvalues i=80/99 {
  display `i'
    quietly append using "http://www.nber.org/morg/annual/morg`i'.dta"
}
* file paths have been omitted for "security reasons"
cd ""
save fullMorg,replace
//We extract data on women in their fourth interview month in
//the Merged Outgoing Rotation Group of the CPS
keep if sex==2&minsamp==4

//We focus on all women between the ages 25 and 50.
keep if age > 24 & age < 51 

//We extract information on weekly earnings, employment status,
//education, age, and state of residence.
keep age uearnwk earnwke state esr grade92 hhid hhnum hurespli year

describe
gen lEarn=ln(earnwke)
* save this (remember I have removed file paths - you might not need this cd anyway!)
cd ""
save morgDat,replace

***************************
*
* Plot some arbitrary data
* pick a state and an age group and plot it. We produce 4 panels of a plot
******************************
* I set a file path here for saving images - you might not need it
cd ""

* local macros for state choice, an arbitrary age, that age +10, and the label associated
* with the state number
local statePick=11
local agePick=25
local agePick2=`agePick'+10
local lab1: value label state
local legLab: label `lab1' `statePick'
twoway (scatter lEarn year if state==`statePick'&age==`agePick',color($sp4)) ///
	(scatter lEarn year if state==`statePick'&age==`agePick2',color($sp7)), ///
	subtitle("State `legLab'") legend(label(1 "Age `agePick'" )  label(2 "Age `agePick2'")) ///
	xtitle("year") ytitle("Log(earnings)")
graph save fig1Panel1,replace

local statePick=54
local agePick=30
local agePick2=`agePick'+10
local lab1: value label state
local legLab: label `lab1' `statePick'
twoway (scatter lEarn year if state==`statePick'&age==`agePick',color($sp4)) ///
	(scatter lEarn year if state==`statePick'&age==`agePick2',color($sp7)), ///
	subtitle("State `legLab'") legend(label(1 "Age `agePick'" )  label(2 "Age `agePick2'")) ///
	xtitle("year") ytitle("Log(earnings)")
graph save fig1Panel2,replace

local statePick=74
local agePick=29
local agePick2=`agePick'+10
local lab1: value label state
local legLab: label `lab1' `statePick'
twoway (scatter lEarn year if state==`statePick'&age==`agePick',color($sp4)) ///
	(scatter lEarn year if state==`statePick'&age==`agePick2',color($sp7)), ///
	subtitle("State `legLab'") legend(label(1 "Age `agePick'" )  label(2 "Age `agePick2'")) ///
	xtitle("year") ytitle("Log(earnings)")
graph save fig1Panel3,replace


local statePick=93
local agePick=32
local agePick2=`agePick'+10
local lab1: value label state
local legLab: label `lab1' `statePick'
twoway (scatter lEarn year if state==`statePick'&age==`agePick',color($sp4)) ///
	(scatter lEarn year if state==`statePick'&age==`agePick2',color($sp7)), ///
	subtitle("State `legLab'") legend(label(1 "Age `agePick'" )  label(2 "Age `agePick2'")) ///
	xtitle("year") ytitle("Log(earnings)")
graph save fig1Panel4,replace

* combine and produce final figure
graph combine fig1Panel1.gph fig1Panel2.gph fig1Panel3.gph fig1Panel4.gph
graph save BDMPostFig1,replace
graph export BDMPostFig1.png,replace


******************
*
* Non-DiD regression model to produce figure 2 and an estimate of the annual change
*
******************
regress lEarn year i.state i.age i.age#c.year i.state#i.age
predict pEarn
twoway  (scatter lEarn year if state==11&age==30,color($sp4)) ///
	(line pEarn year if state==11&age==30,lcolor($sp4)), ///
	legend(label(1 "Observed") label(2 "Predicted"))
graph save bdmFig2,replace
graph export bdmFig2.png,replace
* get teh slope for 30 yearolds
lincom year+30.age#year
	
**********************
*
* Brief digression to check if the data is time series
*
**********************	
* this data set is huge, but I'm not convinced that there are repeat observations of any woman
* across multiple years. Let's check
bysort hurespli year: egen minYear=min(year)
bysort hurespli year: egen maxYear=max(year)
gen yearCheck=(minYear==maxYear)
tab yearCheck
* only 1s! So no woman has observations in multiple years!
	


*******************
*
* producing figure 3: regression results from an arbitrary intervention
*
*******************	
	
* let's try with an intervention
* set it arbitrarily to 1990
gen stepVal=(year>1990)
* need to assign treatment groups, let's do it just based on state numbers
gen treatVal=(state<55)

* also generate age groups for simplicity
gen age5=5*floor(age/5)
tab treatVal
tab age5

* generate hte twfe term
gen twfeTerm=stepVal*treatVal

* simple, unadjusted twfe model
regress lEarn i.twfeTerm i.year i.state



*******************************
*
* Building a better model
*
*******************************
* we try and build a more sophisticated model based on 5 year age groups

/* here we include:
- a continuous trend
- a step for treatment
- a step for period
- a change in level in all age groups at the time of the intervention, that is
common across control and intervention groups
- different trends by age group
- different trends in the treatment and control groups
- we also include the step/period interaction (which measures the DiD effect)
*/



* now build a model
drop pEarn
gen year2=year-1990
regress lEarn year2 i.age5 i.stepVal i.treatVal i.age5#c.year2 i.age5#i.treatVal i.treatVal#c.year2 i.stepVal#i.treatVal
predict pEarn
twoway  (scatter lEarn year if state==59&age5==25,color($sp4)) ///
	(line pEarn year if state==59&age5==25,color($sp4)) ///
	(scatter lEarn year if state==59&age5==35,color($sp7)) ///
	(line pEarn year if state==59&age5==35,color($sp7)), ///
	legend(order(1 "Age 25" 3 "Age 35"))


* now export for use in R
export delimited using "morgDat.csv",replace
