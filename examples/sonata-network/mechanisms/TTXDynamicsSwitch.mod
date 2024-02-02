: Simple switching TTX dynamics
: Outside TTX concentration level stays fixed at ttxo_level, which can be
: set as a range variable
: Implemented by Werner Van Geit @ BlueBrain Project, Jan 2015

NEURON	{
	SUFFIX TTXDynamicsSwitch
	USEION ttx WRITE ttxo, ttxi VALENCE 1
    RANGE ttxo_level
}

PARAMETER	{
    : 1e-12 represent 'no TTX', using 0.0 could generate problems during the 
    : possible calculation of ettx 
    ttxo_level = 1e-12 (mM)

    : Check at the end of the mod file for the reasoning behind this sentinel
    : value
    ttxi_sentinel = 0.015625 (mM) : 1.0/64.0
}

ASSIGNED {
    ttxo (mM)
    ttxi (mM)
}

BREAKPOINT	{
    ttxo = ttxo_level
}

INITIAL {
    ttxo = ttxo_level
    ttxi = ttxi_sentinel
}

: WVG @ BBP, Jan 2015
: The internal ttx concentration is a sentinel value, this value should
: not be used for any calculation
: The reason it is here is that Neuron's default value for a concentration
: is 1.0. In case a sodium channel uses the TTX concentration to 
: enable/disable itself, it shouldn't use the default value. It should only
: use the ttxo value when it has been initialised by this TTXDynamicsSwitch
: mechanism.
: (or any other value that manages the outside ttx concentration)
: The channel should read the ttxi value, and check for this sentinel value
: If it matches, it means this mechanism is control the ttxo concentration
: otherwise ttxo should be ignored
: Chose 1.0/64.0 as sentinel because it can be exactly represented in binary
: floating-point representation
