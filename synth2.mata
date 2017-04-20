mata
mata clear
mata set matastrict on
mata set mataoptimize on
mata set matalnum off

void synth2_d0(todo, real rowvector b, real colvector Xtr, real matrix Xco, real matrix XcoT, real colvector Ztr, real matrix Zco, string scalar cname, string scalar Hname, string scalar wsolname, string scalar pluginline, real scalar v, g, H) {
	real matrix t
	pragma unset todo; pragma unset g; pragma unset H
	t = abs(b) 
	t = XcoT :* (t/sum(t))
	st_matrix(Hname,  t * Xco)
	st_matrix(cname, (-t) * Xtr)
	stata(pluginline) // solve the quadratic programming problem
	t = Ztr - Zco * st_matrix(wsolname)
	v = t't
}

void synth2(string scalar ivar, string scalar tvar, real scalar treati, real scalar treatt, string scalar depvar, string scalar predictors, string scalar predictornames, 
		string scalar mspeperiod, string scalar cname, string scalar Hname, string scalar Aname, string scalar lname, string scalar uname, string scalar wsolname, string scalar plugincmd) {
	real scalar j, std, Nco
	real colvector i, controli, t, mspesample, _Xtr, Ztr, wsol
	real rowvector _mspeperiod, XtrT, _XtrT
	real matrix X, XcoT, _Xco, _XcoT, Zco, ZcoT, V, beta
	string scalar label
	string colvector ivarlabels
	string matrix names
	transmorphic scalar S

	t = st_data(., tvar)
	i = st_data(., ivar)
	controli = uniqrows(i); controli = select(controli, controli :!= treati)
	Nco = rows(controli)

	st_matrix(Aname, J(1, Nco, 1))
	st_matrix(lname, J(1, Nco, 0))
	st_matrix(uname, J(Nco, 1, 1))
	st_matrix(wsolname, J(Nco, 1, 1))

	_mspeperiod = strtoreal(tokens(mspeperiod))

	mspesample = J(st_nobs(), 1, 0)
	for (j=cols(_mspeperiod);j;j--) mspesample = mspesample :| (t:==_mspeperiod[j])
	
	XtrT =          st_data(selectindex(i:==treati :& t:==treatt), predictors)
	XcoT =          st_data(selectindex(i:!=treati :& t:==treatt), predictors)
	Ztr  =          st_data(selectindex(i:==treati :& mspesample), depvar    ) // colvector, row for each time in mspeperiod
	ZcoT = colshape(st_data(selectindex(i:!=treati :& mspesample), depvar    ), cols(_mspeperiod))

	if (hasmissing(XtrT) | hasmissing(XcoT) | hasmissing (Ztr) | hasmissing(ZcoT)) return

	// normalize X vars
	X = XtrT \ XcoT
	std = sqrt(colsum(X:*X) - mean(X):^2)
	_XtrT = XtrT :/ std
	_XcoT = XcoT :/ std

	_Xtr = _XtrT' // colvector, row for each predictor for treatment unit
	_Xco = _XcoT' // maxtrix, row for each predicor, col for each control unit
	Zco  = ZcoT'  // matrix, row for each time in mspeperiod, col for each control unit
	
	if (st_local("customV")=="") { // synth's regression-based method
		X = _XtrT \ _XcoT
		X = X :- mean(X) // partial out constant, which is part of this regression
		beta = invsym(cross(X, X)) * cross(X, (Ztr' \ ZcoT))
		V = rowsum(beta :* beta)'
	} else {
		V = st_matrix(st_local("customV"))
		if (rows(V) > 1)
			V = diagonal(V)'
	}
	V = V / rowsum(V)
	
	if (st_local("nested") != "") {
		S = optimize_init()
		optimize_init_evaluator(S, &synth2_d0())
		optimize_init_which(S, "min")
		optimize_init_params(S, V)
		optimize_init_technique(S, "nm")
		optimize_init_nmsimplexdeltas(S, .1)
		optimize_init_tracelevel(S, "none")
		optimize_init_argument(S, 1, _Xtr)
		optimize_init_argument(S, 2, _Xco)
		optimize_init_argument(S, 3, _XcoT)
		optimize_init_argument(S, 4, Ztr)
		optimize_init_argument(S, 5, Zco)
		optimize_init_argument(S, 6, cname)
		optimize_init_argument(S, 7, Hname)
		optimize_init_argument(S, 8, wsolname)
		optimize_init_argument(S, 9, plugincmd)

		(void) _optimize(S)
		V = optimize_result_converged(S)? optimize_result_params(S) : . // or otherwise handle error
	}

	t = _XcoT :* V // compute weight matrix one last time, at optimum (if nested), or only time (otherwise)
	st_matrix(Hname, t * _Xco)
	st_matrix(cname, -(t * _Xtr))
	stata(plugincmd)
		
	wsol = round(st_matrix(wsolname), .001)
	t = Ztr - Zco * wsol; st_numscalar("e(RMSPE)", t't)

	if ((label=st_varvaluelabel(ivar)) != "") // label matrices with labels, if any
		ivarlabels = st_vlmap(label, controli)
	else
		ivarlabels = strofreal(controli)

	st_matrix("e(W)", wsol)
	st_matrixrowstripe("e(W)", (J(rows(wsol),1,""),ivarlabels))
	
	st_matrix("e(Y_treated)", XtrT')
	st_matrix("e(Y_synthetic)", XcoT ' wsol)
	names = J(cols(XtrT),1,""), tokens(predictornames)'
	st_matrixrowstripe("e(Y_treated)"  , names)
	st_matrixrowstripe("e(Y_synthetic)", names)
	
	st_matrix("e(V_matrix)", diag(V))
	names = J(cols(V),1,""), tokens(mspeperiod)'
	st_matrixrowstripe("e(V_matrix)", names)
	st_matrixcolstripe("e(V_matrix)", names)
}

mata mlib create lsynth2, dir("`c(sysdir_plus)'l") replace
mata mlib add lsynth2 *(), dir("`c(sysdir_plus)'l")
mata mlib index
end
