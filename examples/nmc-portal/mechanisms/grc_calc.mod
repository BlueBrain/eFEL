TITLE Cerebellum Granule Cell Model

COMMENT
	Reference: Theta-Frequency Bursting and Resonance in Cerebellar Granule Cells:Experimental
	Evidence and Modeling of a Slow K+-Dependent Mechanism
	Egidio D'Angelo,Thierry Nieus,Arianna Maffei,Simona Armano,Paola Rossi,Vanni Taglietti,
	Andrea Fontana and Giovanni Naldi
ENDCOMMENT

NEURON {
        SUFFIX Calc
        USEION ca READ ica, cao WRITE cai
        RANGE d, beta, cai0
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	F      = (faraday) (coulomb)
}

PARAMETER {
        ica             (mA/cm2)
        celsius = 30    (degC)
        d = .2          (um)
        cao = 2.        (mM)         
        cai0 = 1e-4     (mM)         
        beta = 1.5        (/ms)
}

STATE {
	cai (mM)
}

INITIAL {
        cai = cai0 
}

BREAKPOINT {
       SOLVE conc METHOD derivimplicit
}

DERIVATIVE conc {    
	cai' = -ica/(2*F*d)*(1e4) - beta*(cai-cai0)
}

