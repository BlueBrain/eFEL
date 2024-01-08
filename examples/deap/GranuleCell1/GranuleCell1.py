# tracer for mono compartment granule cell
# Stefano "Bremen" Masoli - 21 june 2014

from neuron import h
import numpy as np


class GrCmono:

    def __init__(self, conductances):

        # Properties of the cell soma
        self.soma = h.Section(name='soma')
        self.soma.nseg = 1
        self.soma.diam = 9.76
        self.soma.L = 9.76
        self.soma.Ra = 100
        self.soma.cm = 1

        # Soma channels
        self.soma.insert('GrC_CaHVA')
        self.soma.gcabar_GrC_CaHVA = conductances['GrC_CaHVA']

        self.soma.insert('GrC_KA')
        self.soma.gkbar_GrC_KA = conductances['GrC_KA']  # for the 5 + 0.009

        self.soma.insert('GrC_Lkg1')
        self.soma.gl_GrC_Lkg1 = conductances['GrC_Lkg1']

        self.soma.insert('GrG_KM')
        self.soma.gkbar_GrG_KM = conductances['GrG_KM']

        self.soma.insert('Calc')
        self.soma.beta_Calc = conductances['Calc']

        self.soma.insert('GrC_KCa')
        self.soma.gkbar_GrC_KCa = conductances['GrC_KCa']

        self.soma.insert('GrC_Kir')
        self.soma.gkbar_GrC_Kir = conductances['GrC_Kir']

        self.soma.insert('GrC_Lkg2')
        self.soma.ggaba_GrC_Lkg2 = conductances['GrC_Lkg2']

        self.soma.insert('GrC_pNa')
        self.soma.gnabar_GrC_pNa = conductances['GrC_pNa']

        self.soma.insert('GrG_KV')
        self.soma.gkbar_GrG_KV = conductances['GrG_KV']

        self.soma.insert('GrG_Na')
        self.soma.gnabar_GrG_Na = conductances['GrG_Na']

        self.soma.insert('GrG_Nar')
        self.soma.gnabar_GrG_Nar = conductances['GrG_Nar']

        # Specific Reversal potential. The one for the two leakage are taken
        # directly from the mod files.
        self.soma.ena = 87.39
        self.soma.ek = -84.69
        self.soma.eca = 129.33


def steps(step_number, parameters):

    # instantiate the class
    cell = GrCmono(parameters)

    # no table are in use from the mod files
    h.usetable_GrG_Na = 0
    h.usetable_GrC_pNa = 0
    h.usetable_GrC_CaHVA = 0
    h.usetable_GrG_KV = 0
    h.usetable_GrC_KA = 0
    h.usetable_GrC_Kir = 0
    h.usetable_GrC_KCa = 0
    h.usetable_GrG_KM = 0

    # load graph for the membrane voltage in mV
    # h('load_file("vm.ses")')
    # h.nrncontrolmenu()
    # Stimulation data. Del = delay in ms, dur = duration in ms, amp =
    # amplitude in nA
    stimdata = dict()

    if step_number == 0:
        stimdata['stim0del'] = 100
        stimdata['stim0dur'] = 500
        stimdata['stim0amp'] = 0.010
    elif step_number == 1:
        stimdata['stim1del'] = 100
        stimdata['stim1dur'] = 500
        stimdata['stim1amp'] = 0.016
    elif step_number == 2:
        stimdata['stim2del'] = 100
        stimdata['stim2dur'] = 500
        stimdata['stim2amp'] = 0.022

    # totale time of the simulation
    stimdata['timeglobal'] = 700

    # integration step
    h.dt = 0.025

    h.cvode_active(1)
    # temperature
    h.celsius = 30
    # time at which the simulation have to end
    h.tstop = stimdata['timeglobal']
    # initial voltage
    h.v_init = -80

    if step_number == 0:
        # actual current injections, one for each current step.
        stim = [h.IClamp(0.5, sec=cell.soma)]

        stim[0].delay = stimdata['stim0del']
        stim[0].dur = stimdata['stim0dur']
        stim[0].amp = stimdata['stim0amp']
    elif step_number == 1:
        stim2 = [h.IClamp(0.5, sec=cell.soma)]

        stim2[0].delay = stimdata['stim1del']
        stim2[0].dur = stimdata['stim1dur']
        stim2[0].amp = stimdata['stim1amp']
    elif step_number == 2:
        stim3 = [h.IClamp(0.5, sec=cell.soma)]

        stim3[0].delay = stimdata['stim2del']
        stim3[0].dur = stimdata['stim2dur']
        stim3[0].amp = stimdata['stim2amp']

    # Code to record everything.
    time = h.Vector()
    time.record(h._ref_t)

    vm = h.Vector()
    vm.record(cell.soma(0.5)._ref_v)

    # function to initialize the simulation.
    def initialize():
        h.finitialize()
        h.run()

    initialize()

    #np.savetxt('trace_%d.txt' % step_number, np.column_stack((np.array(time),np.array(vm))), delimiter = ' ')

    return np.array(time), np.array(vm)
