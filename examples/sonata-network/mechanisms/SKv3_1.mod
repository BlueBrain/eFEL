:Comment :
:Reference : :	Rettig et.al (1992) EMBO J 11, no. 7: 2473-86.
:								Methods: Grupe et al. (1990) EMBO J 9, 1749-1756.
: LJP: OK, no LJP, "Patch pipetes were filed with the normal bathing solution in al experiments.""

NEURON	{
	SUFFIX SKv3_1
	USEION k READ ek WRITE ik
	RANGE gSKv3_1bar, gSKv3_1, ik
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gSKv3_1bar = 0.00001 (S/cm2)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gSKv3_1	(S/cm2)
	mInf
	mTau
}

STATE	{
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gSKv3_1 = gSKv3_1bar*m
	ik = gSKv3_1*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
	rates()
	m = mInf
}

PROCEDURE rates(){
	UNITSOFF
		mInf =  1/(1+exp(((v -(18.700))/(-9.700))))
		mTau =  0.2*20.000/(1+exp(((v -(-46.560))/(-44.140))))
	UNITSON
}
