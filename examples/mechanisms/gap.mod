NEURON {
POINT_PROCESS Gap
RANGE vgap
RANGE g, i
NONSPECIFIC_CURRENT i
RANGE synapseID, selected_for_report
}
PARAMETER {
    g = 1 (nanosiemens)
    selected_for_report = 0
}
ASSIGNED {
v (millivolt)
vgap (millivolt)
i (nanoamp)
synapseID
}
BREAKPOINT { i = (v - vgap)*(g*1e-3) }
