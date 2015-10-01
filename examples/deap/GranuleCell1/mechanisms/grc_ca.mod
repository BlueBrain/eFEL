TITLE Cerebellum Granule Cell Model

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_CaHVA 
	USEION ca READ eca WRITE ica 
	RANGE gcabar, ica, g, alpha_s, beta_s, alpha_u, beta_u 
	RANGE Aalpha_s, Kalpha_s, V0alpha_s
	RANGE Abeta_s, Kbeta_s, V0beta_s
	RANGE Aalpha_u, Kalpha_u, V0alpha_u
	RANGE Abeta_u, Kbeta_u, V0beta_u
	RANGE s_inf, tau_s, u_inf, tau_u 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 


	Aalpha_s = 0.04944 (/ms)
	Kalpha_s =  15.87301587302 (mV)
	V0alpha_s = -29.06 (mV)
	
	Abeta_s = 0.08298 (/ms)
	Kbeta_s =  -25.641 (mV)
	V0beta_s = -18.66 (mV)

	Aalpha_u = 0.0013 (/ms)
	Kalpha_u =  -18.183 (mV)
	V0alpha_u = -48 (mV)
		
	Abeta_u = 0.0013 (/ms)
	Kbeta_u =   83.33 (mV)
	V0beta_u = -48 (mV)

	v (mV) 
	gcabar= 0.00046 (mho/cm2) 
	eca = 129.33 (mV) 
	celsius = 30 (degC) 
} 

STATE { 
	s 
	u 
} 

ASSIGNED { 
	ica (mA/cm2) 
	s_inf 
	u_inf 
	tau_s (ms) 
	tau_u (ms) 
	g (mho/cm2) 
	alpha_s (/ms)
	beta_s (/ms)
	alpha_u (/ms)
	beta_u (/ms)
} 
 
INITIAL { 
	rate(v) 
	s = s_inf 
	u = u_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gcabar*s*s*u 
	ica = g*(v - eca) 
	alpha_s = alp_s(v)
	beta_s = bet_s(v)
	alpha_u = alp_u(v)
	beta_u = bet_u(v)
}
 
DERIVATIVE states { 
	rate(v) 
	s' =(s_inf - s)/tau_s 
	u' =(u_inf - u)/tau_u 
} 
 
FUNCTION alp_s(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	alp_s = Q10*Aalpha_s*exp((v-V0alpha_s)/Kalpha_s) 
} 
 
FUNCTION bet_s(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	bet_s = Q10*Abeta_s*exp((v-V0beta_s)/Kbeta_s) 
} 
 
FUNCTION alp_u(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	alp_u = Q10*Aalpha_u*exp((v-V0alpha_u)/Kalpha_u) 
} 
 
FUNCTION bet_u(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	bet_u = Q10*Abeta_u*exp((v-V0beta_u)/Kbeta_u) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_s, b_s, a_u, b_u 
	TABLE s_inf, tau_s, u_inf, tau_u 
	DEPEND Aalpha_s, Kalpha_s, V0alpha_s, 
	       Abeta_s, Kbeta_s, V0beta_s,
               Aalpha_u, Kalpha_u, V0alpha_u,
               Abeta_u, Kbeta_u, V0beta_u, celsius FROM -100 TO 100 WITH 200 
	a_s = alp_s(v)  
	b_s = bet_s(v) 
	a_u = alp_u(v)  
	b_u = bet_u(v) 
	s_inf = a_s/(a_s + b_s) 
	tau_s = 1/(a_s + b_s) 
	u_inf = a_u/(a_u + b_u) 
	tau_u = 1/(a_u + b_u) 
}

