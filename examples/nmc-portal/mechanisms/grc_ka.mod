TITLE Cerebellum Granule Cell Model

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
ENDCOMMENT

NEURON { 
	SUFFIX GrC_KA
	USEION k READ ek WRITE ik 
	RANGE gkbar, ik, g, alpha_a, beta_a, alpha_b, beta_b
	RANGE Aalpha_a, Kalpha_a, V0alpha_a
	RANGE Abeta_a, Kbeta_a, V0beta_a
	RANGE Aalpha_b, Kalpha_b, V0alpha_b
	RANGE Abeta_b, Kbeta_b, V0beta_b
	RANGE V0_ainf, K_ainf, V0_binf, K_binf
	RANGE a_inf, tau_a, b_inf, tau_b 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_a = 4.88826 (/ms)
	Kalpha_a = -23.32708 (mV)
	V0alpha_a = -9.17203 (mV)
	Abeta_a = 0.99285 (/ms)
	Kbeta_a = 19.47175 (mV)
	V0beta_a = -18.27914 (mV)

	Aalpha_b = 0.11042 (/ms)
	Kalpha_b = 12.8433 (mV)
	V0alpha_b = -111.33209 (mV)
	Abeta_b = 0.10353 (/ms)
	Kbeta_b = -8.90123 (mV)
	V0beta_b = -49.9537 (mV)

	V0_ainf = -46.7 (mV)
	K_ainf = -19.8 (mV)

	V0_binf = -78.8 (mV)
	K_binf = 8.4 (mV)
	v (mV) 
	gkbar= 0.004 (mho/cm2) 
	ek = -84.69 (mV) 
	celsius = 30 (degC) 
} 

STATE { 
	a
	b 
} 

ASSIGNED { 
	ik (mA/cm2) 
	a_inf 
	b_inf 
	tau_a (ms) 
	tau_b (ms) 
	g (mho/cm2) 
	alpha_a (/ms)
	beta_a (/ms)
	alpha_b (/ms)
	beta_b (/ms)
} 
 
INITIAL { 
	rate(v) 
	a = a_inf 
	b = b_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gkbar*a*a*a*b 
	ik = g*(v - ek)
	alpha_a = alp_a(v)
	beta_a = bet_a(v) 
	alpha_b = alp_b(v)
	beta_b = bet_b(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	a' =(a_inf - a)/tau_a 
	b' =(b_inf - b)/tau_b 
} 
 
FUNCTION alp_a(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	alp_a = Q10*Aalpha_a*sigm(v-V0alpha_a,Kalpha_a)
} 
 
FUNCTION bet_a(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	bet_a = Q10*Abeta_a/(exp((v-V0beta_a)/Kbeta_a))
} 
 
FUNCTION alp_b(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	alp_b = Q10*Aalpha_b*sigm(v-V0alpha_b,Kalpha_b)
} 
 
FUNCTION bet_b(v(mV))(/ms) { LOCAL Q10
	Q10 = 3^((celsius-20(degC))/10(degC))
	bet_b = Q10*Abeta_b*sigm(v-V0beta_b,Kbeta_b)
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_a, b_a, a_b, b_b 
	TABLE a_inf, tau_a, b_inf, tau_b 
	DEPEND Aalpha_a, Kalpha_a, V0alpha_a, 
	       Abeta_a, Kbeta_a, V0beta_a,
               Aalpha_b, Kalpha_b, V0alpha_b,
               Abeta_b, Kbeta_b, V0beta_b, celsius FROM -100 TO 100 WITH 200 
	a_a = alp_a(v)  
	b_a = bet_a(v) 
	a_b = alp_b(v)  
	b_b = bet_b(v) 
	a_inf = 1/(1+exp((v-V0_ainf)/K_ainf)) 
	tau_a = 1/(a_a + b_a) 
	b_inf = 1/(1+exp((v-V0_binf)/K_binf))
	tau_b = 1/(a_b + b_b) 
}

FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
                linoid = x/(exp(x/y) - 1)
        }
} 

FUNCTION sigm(x (mV),y (mV)) {
                sigm = 1/(exp(x/y) + 1)
}
