cap program drop synth2
program define synth2, eclass
	version 10
	
	cap findfile "synthopt.plugin"
	if _rc {
		di as err `"{cmd:synth2} requires {cmd:synth}. To install type or click {stata "ssc install synth"}."'
		error 601
	}

	syntax varlist(ts fv numeric) [if] [in], TRUnit(numlist max=1 int) TRPeriod(numlist max=1 int) ///
		[NESTed MSPEPERiod(numlist int) customV(string) margin(real 0.005) maxiter(integer 1000) sigf(integer 12) bound(integer 10)]
	
	if `"`customV'"' != "" confirm matrix `customV'

	local tvar: char _dta[tis]
	local ivar: char _dta[iis]
	if "`tvar'"=="" | "`ivar'"=="" {
		di as err "You must {help xtset} the data to specify the panel and time variables."
		error 459
	}
	
	fvexpand `varlist'
	local varlist `r(varlist)'
	tokenize `r(varlist)'
	macro shift
	local predictornames `*'

	fvrevar `varlist'
	tokenize `r(varlist)'
	local depvar `1'
	macro shift
	local predictors `*'
	
	tempname A l u bslack wsol H c boundmat marginmat maxitermat sigfmat touse
	mat `bslack' = 1
	mat `boundmat' = `bound'
	mat `marginmat' = `margin'
	mat `maxitermat' = `maxiter'
	mat `sigfmat' = `sigf'
	local plugincmd qui plugin call synthopt, `c' `H' `A' `bslack' `l' `u' `boundmat' `marginmat' `maxitermat' `sigfmat' `wsol'

	preserve
	keep `ivar' `tvar' `r(varlist)'
	cap keep `in' `if'

	if "`mspeperiod'"=="" {
		local delta = `:char _dta[_TSdelta]'
		sum `tvar' if `depvar'<., mean
		numlist "`r(min)'(`delta')`=`trperiod' - `delta''"
		local mspeperiod `r(numlist)'
	}
	
	tsfill, full // give data canonical structure

	ereturn clear
	mata synth2("`ivar'","`tvar'",`trunit',`trperiod',"`depvar'","`predictors'","`predictornames'","`mspeperiod'","`c'","`H'","`A'","`l'","`u'","`wsol'","`plugincmd'")
end

cap program synthopt, plugin
