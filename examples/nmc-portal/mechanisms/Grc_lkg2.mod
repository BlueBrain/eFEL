TITLE Cerebellum Granule Cell Model

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_Lkg2 
	NONSPECIFIC_CURRENT il
	RANGE egaba, ggaba , i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	ggaba = 2.17e-5 (mho/cm2)
	celsius = 30 (degC)
	egaba = -65 (mV)
} 

ASSIGNED { 
	il (mA/cm2) 
	i (mA/cm2) 
}

BREAKPOINT { 
	il = ggaba*(v - egaba) 
	i =il
} 
