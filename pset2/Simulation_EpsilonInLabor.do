
*Initialize simulation
clear
set obs 1000
set seed 2000

*Create data
gen Firm=_n

global Rho=0.9

*Simulate 100 periods and keep the last 10 to make sure that the series is close to a stationary state
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

*Use analytical relations to construct inputs and revenue (in the process we also create some unobserved variables but we will not use those for estimation)
gen Epsilon_0=rnormal(0,3)
forvalues t=0/9{
    
    local tt=`t'+1
	gen Epsilon_`tt'=rnormal(0,3)
	gen U_`tt'=rnormal(0,.2)
    gen L_`tt'=(0.108)^3*(0.73)^0.8*exp(3.69+3.6*Omega_`t')
	gen K_`tt'=0.73*L_`tt'
	replace L_`tt'=L_`tt'*exp(Epsilon_`tt')
	gen M_`tt'=(2/5)^(5/3)*exp(8/75)*(L_`tt'^(0.3)*K_`tt'^(0.2)*exp(Omega_`tt'))^(4/3)
    
	gen Y_`tt'=(exp(U_`tt'+Omega_`tt')*L_`tt'^(0.3)*K_`tt'^(0.2)*M_`tt'^(0.5))^(0.8)
	
}

*Reshape data and get log-variables
reshape long Y_ L_ K_ M_ Omega_ U_ Epsilon_, i(Firm) j(t)

drop if t==0

foreach var in Y K L M {
    gen log_`var'_=log(`var'_)
}

*First stage predictions
reg log_Y_ log_K_ log_L_ log_M_ 
predict log_Y_Hat

*Get lags
sort Firm t
bysort Firm: gen log_Y_Hat_lag=log_Y_Hat[_n-1]
bysort Firm: gen log_K_lag=log_K_[_n-1]
bysort Firm: gen log_L_lag=log_L_[_n-1]
bysort Firm: gen log_M_lag=log_M_[_n-1]
bysort Firm: gen log_Y_lag=log_Y_[_n-1]

*************************************************************
*GMM- Ackerberg et Al
*************************************************************
gmm (log_Y_-({b0=-.1}+{bk}*log_K_+{bl}*log_L_)-{RhoHat=.9}*(log_Y_Hat_lag-({b0}+{bk}*log_K_lag+{bl}*log_L_lag))), instruments(log_K_ log_L_lag log_Y_Hat_lag ) 

global b0=_b[/b0]
global bl=_b[/bl]
global bk=_b[/bk]
global RhoHat=_b[/RhoHat]

do D:\DDC\PhD_IO\pset2\SolveMata

putexcel set D:\DDC\PhD_IO\pset2\Results, modify

putexcel c4=$bl
putexcel d4=$bk
putexcel e4=$b0
putexcel f4=$RhoHat

putexcel c9=matrix(r)
putexcel f9=$RhoHat

*************************************************************
*GMM- Blundell-Bond
*************************************************************
gmm (log_Y_-{RhoHat=.8}*log_Y_Hat_lag-{b0}*(1-{RhoHat})-{bk}*(log_K_-{RhoHat}*log_K_lag)-{bl}*(log_L_-{RhoHat}*log_L_lag)), instruments(log_K_ log_K_lag log_L_ log_L_lag ) 

global b0=_b[/b0]
global bl=_b[/bl]
global bk=_b[/bk]
global RhoHat=_b[/RhoHat]

do D:\DDC\PhD_IO\pset2\SolveMata

putexcel i4=$bl
putexcel j4=$bk
putexcel k4=$b0
putexcel l4=$RhoHat


putexcel i9=matrix(r)
putexcel l9=$RhoHat
