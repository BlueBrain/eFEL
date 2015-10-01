TITLE Cerebellum Granule Cell Model

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
ENDCOMMENT
 
NEURON { 
	SUFFIX GrC_Lkg1 
	NONSPECIFIC_CURRENT il
	RANGE el, gl,i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gl = 5.68e-5 (mho/cm2)
	celsius = 30 (degC)
	el = -58 (mV)
} 

ASSIGNED { 
	il (mA/cm2) 
	i (mA/cm2) 
}
  
BREAKPOINT { 
	il = gl*(v - el) 
	i = il
} 
