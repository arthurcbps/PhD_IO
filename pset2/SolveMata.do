
global b0=_b[/b0]
global bl=_b[/bl]
global bk=_b[/bk]

clear mata

mata

	b0=$b0
	bl=exp($bl)
	bk=exp($bk)

	void AckerbergSolve(todo, p, lnf,  S, H){
		
		al=p[1]
		ak=p[2]
		theta=p[3]
		
		lnf = (al*theta)/(1-theta*(1-al-ak))-.3965 \   
		(ak*theta)/(1-theta*(1-al-ak))-.2650 \
		theta*((1/(1-theta*(1-al-ak))) *(ln(theta*(1-al-ak))+theta^2/10)) +1.136
		
		lnf=lnf'*lnf
		
	}

	S = optimize_init()

	optimize_init_evaluator(S, &AckerbergSolve())

	optimize_init_evaluatortype(S, "v0")

	optimize_init_params(S, (.4,.4,.1))

	optimize_init_which(S,  "min" )

	optimize_init_tracelevel(S,"none")

	optimize_init_conv_ptol(S, 1e-16)

	optimize_init_conv_vtol(S, 1e-16)

	p = optimize(S)

	p 
end

