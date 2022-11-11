
clear

set obs 1000
set seed 2022

gen Firm=_n

global Rho=0.9

gen Omega_0=rnormal(0,1)
forvalues t=1/100{
	local tt=`t'-1
	gen Omega_`t'=$Rho * Omega_`tt'+rnormal(0,1)
}

forvalues t=0/100{
    if `t'<90{
		drop Omega_`t'
	}
	else{
	    local tt=`t'-90
		rename Omega_`t' Omega_`tt'
	}
}

gen Epsilon_0=rnormal(0,3)
forvalues t=0/9{
    
    local tt=`t'+1
	gen Epsilon_`tt'=rnormal(0,3)
	gen U_`tt'=rnormal(0,.2)
    gen L_`tt'=(0.108)^3*(0.73)^0.8*exp(3.69+3.6*Omega_`t')
	gen K_`tt'=0.73*L_`tt'
	gen M_`tt'=(2/5)^(5/3)*exp(8/75)*(L_`tt'^(0.3)*K_`tt'^(0.2)*exp(0.3*Epsilon_`tt'+Omega_`tt'))^(4/3)
    
	gen Y_`tt'=(exp(U_`tt'+Omega_`tt'+0.3*Epsilon_`t')*L_`tt'^(0.3)*K_`tt'^(0.2)*M_`tt'^(0.5))^(0.8)
	replace Y_`tt'=Y_`tt'^((1-5)/(-5))
	
}

reshape long Y_ L_ K_ M_ Omega_ U_ Epsilon_, i(Firm) j(t)

drop if t==0

foreach var in Y K L M {
    gen log_`var'_=log(`var'_)
}

reg  log_Y_ log_K_ log_L_ log_M_

predict log_Y_Hat

sort Firm t
bysort Firm: gen log_Y_Hat_lag=log_Y_Hat[_n-1]
bysort Firm: gen log_K_lag=log_K_[_n-1]
bysort Firm: gen log_L_lag=log_L_[_n-1]
bysort Firm: gen log_M_lag=log_M_[_n-1]


gmm (log_Y_-{b0}-{bk}*log_K_-{bl}*log_L_-(1-{bl}-{bk})*log_M_-{RhoHat}*(log_Y_Hat_lag-{b0}-{bk}*log_K_lag-{bl}*log_L_lag-(1-{bl}-{bk})*log_M_lag)), instruments(log_K_ log_L_lag log_Y_Hat_lag) 

*gmm (log_Y_-{b0}-{bk}*log_K_-{bl}*log_L_-{RhoHat}*(log_Y_Hat_lag-{b0}-{bk}*log_K_lag-{bl}*log_L_lag)), instruments(log_K_ log_L_lag log_Y_Hat_lag) 



