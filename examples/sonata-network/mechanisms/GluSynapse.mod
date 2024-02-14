COMMENT
/**
 * @file GluSynapse.mod
 * @brief Probabilistic synapse featuring long-term plasticity
 * @author king, chindemi, rossert
 * @date 2021-05-19
 * @version 1.0.1
 * @remark Copyright 2005-2023 Blue Brain Project / EPFL
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * 
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
 Glutamatergic synapse model featuring:
1) AMPA receptor with a dual-exponential conductance profile.
2) NMDA receptor  with a dual-exponential conductance profile and magnesium
   block as described in Jahr and Stevens 1990.
3) Tsodyks-Markram presynaptic short-term plasticity as Barros et al. 2019.
   Implementation based on the work of Eilif Muller, Michael Reimann and
   Srikanth Ramaswamy (Blue Brain Project, August 2011), who introduced the
   2-state Markov model of vesicle release. The new model is an extension of
   Fuhrmann et al. 2002, motivated by the following constraints:
        a) No consumption on failure
        b) No release until recovery
        c) Same ensemble averaged trace as canonical Tsodyks-Markram using same
           parameters determined from experiment.
   For a pre-synaptic spike or external spontaneous release trigger event, the
   synapse will only release if it is in the recovered state, and with
   probability u (which follows facilitation dynamics). If it releases, it will
   transition to the unrecovered state. Recovery is as a Poisson process with
   rate 1/Dep.
   John Rahmon and Giuseppe Chindemi introduced multi-vesicular release as an
   extension of the 2-state Markov model of vesicle release described above
   (Blue Brain Project, February 2017).
4) NMDAR-mediated calcium current. Fractional calcium current Pf_NMDA from
   Schneggenburger et al. 1993. Fractional NMDAR conductance treated as a
   calcium-only permeable channel with Erev = 40 mV independent of extracellular
   calcium concentration (see Jahr and Stevens 1993). Implemented by Christian
   Rossert and Giuseppe Chindemi (Blue Brain Project, 2016).
5) Spine volume.
6) VDCC.
7) Postsynaptic calcium dynamics.
8) Long-term synaptic plasticity. Calcium-based STDP model based on Graupner and
   Brunel 2012.
Model implementation, optimization and simulation curated by James King (Blue
Brain Project, 2021).
ENDCOMMENT


TITLE Glutamatergic synapse

NEURON {
    THREADSAFE
    POINT_PROCESS GluSynapse
    : AMPA Receptor
    GLOBAL tau_r_AMPA, E_AMPA
    RANGE tau_d_AMPA, gmax0_AMPA, gmax_d_AMPA, gmax_p_AMPA, g_AMPA
    : NMDA Receptor
    GLOBAL scale_NMDA, slope_NMDA
    GLOBAL tau_r_NMDA, tau_d_NMDA, E_NMDA
    RANGE gmax_NMDA, g_NMDA
    GLOBAL mg
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    RANGE Use, Dep, Fac, Nrrp
    RANGE Use_d, Use_p
    BBCOREPOINTER rng
    : NMDAR-mediated calcium current
    RANGE ica_NMDA
    : Spine
    RANGE volume_CR
    : VDCC (R-type)
    GLOBAL ljp_VDCC, vhm_VDCC, km_VDCC, mtau_VDCC, vhh_VDCC, kh_VDCC, htau_VDCC, gca_bar_VDCC
    RANGE ica_VDCC
    : Postsynaptic Ca2+ dynamics
    GLOBAL gamma_ca_CR, tau_ca_CR, min_ca_CR, cao_CR
    : Long-term synaptic plasticity
    GLOBAL rho_star_GB, tau_ind_GB, tau_exp_GB, tau_effca_GB
    GLOBAL gamma_d_GB, gamma_p_GB
    RANGE theta_d_GB, theta_p_GB, rho0_GB, dep_GB, pot_GB
    : Misc
    RANGE vsyn, synapseID, selected_for_report, verbose
    NONSPECIFIC_CURRENT i
    RANGE conductance
    RANGE next_delay
    BBCOREPOINTER delay_times, delay_weights
    GLOBAL nc_type_param
    GLOBAL minis_single_vesicle
    GLOBAL init_depleted
    : For debugging
    :RANGE sgid, tgid
}


UNITS {
    (nA)    = (nanoamp)
    (mV)    = (millivolt)
    (uS)    = (microsiemens)
    (nS)    = (nanosiemens)
    (pS)    = (picosiemens)
    (umho)  = (micromho)
    (um)    = (micrometers)
    (mM)    = (milli/liter)
    (uM)    = (micro/liter)
    FARADAY = (faraday) (coulomb)
    PI      = (pi)      (1)
    R       = (k-mole)  (joule/degC)
}


PARAMETER {
    celsius                     (degC)
    : AMPA Receptor
    tau_r_AMPA      = 0.2       (ms)        : Tau rise, dual-exponential conductance profile
    tau_d_AMPA      = 1.7       (ms)        : Tau decay, IMPORTANT: tau_r < tau_d
    E_AMPA          = 0         (mV)        : Reversal potential
    gmax0_AMPA      = 1.0       (nS)        : Initial peak conductance
    gmax_d_AMPA     = 1.0       (nS)        : Peak conductance in the depressed state
    gmax_p_AMPA     = 2.0       (nS)        : Peak conductance in the potentitated state
    : NMDA Receptor
    mg              = 1         (mM)        : Extracellular magnesium concentration
    scale_NMDA      = 2.552     (mM)        : Scale of the mg block (Vargas-Caballero and Robinson 2003)
    slope_NMDA      = 0.072     (/mV)       : Slope of the ma block (Vargas-Caballero and Robinson 2003)
    tau_r_NMDA      = 0.29      (ms)        : Tau rise, dual-exponential conductance profile
    tau_d_NMDA      = 70        (ms)        : Tau decay, IMPORTANT: tau_r < tau_d
    E_NMDA          = -3        (mV)        : Reversal potential (Vargas-Caballero and Robinson 2003)
    gmax_NMDA       = 0.55      (nS)        : Peak conductance
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use             = 0.5       (1)         : Initial utilization of synaptic efficacy
    Dep             = 100       (ms)        : Relaxation time constant from depression
    Fac             = 10        (ms)        : Relaxation time constant from facilitation
    Nrrp            = 1         (1)         : Number of release sites for given contact
    Use_d           = 0.2       (1)         : Depressed Use
    Use_p           = 0.8       (1)         : Potentiated Use
    : Spine
    volume_CR       = 0.087     (um3)       : From spine data by Ruth Benavides-Piccione (unpublished)
    : VDCC (R-type)
    gca_bar_VDCC    = 0.0744    (nS/um2)    : Density spines: 20 um-2 (Sabatini 2000), unitary conductance VGCC 3.72 pS (Bartol 2015)
    ljp_VDCC        = 0         (mV)
    vhm_VDCC        = -5.9      (mV)        : v 1/2 for act, Magee and Johnston 1995 (corrected for m*m)
    km_VDCC         = 9.5       (mV)        : act slope, Magee and Johnston 1995 (corrected for m*m)
    vhh_VDCC        = -39       (mV)        : v 1/2 for inact, Magee and Johnston 1995
    kh_VDCC         = -9.2      (mV)        : inact, Magee and Johnston 1995
    mtau_VDCC       = 1         (ms)        : max time constant (guess)
    htau_VDCC       = 27        (ms)        : max time constant 100*0.27
    : Postsynaptic Ca2+ dynamics
    gamma_ca_CR     = 0.04      (1)         : Percent of free calcium (not buffered), Sabatini et al 2002: kappa_e = 24+-11 (also 14 (2-31) or 22 (18-33))
    tau_ca_CR       = 12        (ms)        : Rate of removal of calcium, Sabatini et al 2002: 14ms (12-20ms)
    min_ca_CR       = 70e-6     (mM)        : Sabatini et al 2002: 70+-29 nM, per AP: 1.1 (0.6-8.2) uM = 1100 e-6 mM = 1100 nM
    cao_CR          = 2.0       (mM)        : Extracellular calcium concentration in slices
    : Long-term synaptic plasticity
    rho_star_GB     = 0.5       (1)
    tau_ind_GB      = 70        (s)
    tau_exp_GB      = 100       (s)
    tau_effca_GB    = 200       (ms)
    gamma_d_GB      = 100       (1)
    gamma_p_GB      = 450       (1)
    theta_d_GB      = 0.006     (us/liter)
    theta_p_GB      = 0.001     (us/liter)
    rho0_GB         = 0         (1)
    : Misc
    synapseID       = 0
    verbose         = 0
    selected_for_report = 0
    conductance     = 0.0
    nc_type_param = 5
    minis_single_vesicle = 0   :// 0 -> no limit (old behavior)
    init_depleted = 0          :// 0 -> init full (old behavior)
    :sgid = -1
    :tgid = -1
}


VERBATIM
/**
 * This Verbatim block is needed to generate random numbers from a uniform
 * distribution U(0, 1).
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifndef NRN_VERSION_GTEQ_8_2_0
#include "nrnran123.h"
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);

#ifndef CORENEURON_BUILD
extern int ifarg(int iarg);

extern void* vector_arg(int iarg);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
#endif
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif


ENDVERBATIM


ASSIGNED {
    g_AMPA          (uS)    : AMPA Receptor
    g_NMDA          (uS)    : NMDA Receptor
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    rng                     : Random Number Generator
    usingR123               : TEMPORARY until mcellran4 completely deprecated
    ica_NMDA        (nA)    : NMDAR-mediated calcium current
    ica_VDCC        (nA)    : VDCC (R-type)
    : Long-term synaptic plasticity
    dep_GB          (1)
    pot_GB          (1)
    : Misc
    v               (mV)
    vsyn            (mV)
    i               (nA)

    : stuff for delayed connections
    delay_times
    delay_weights
    next_delay (ms)
}

STATE {
    : AMPA Receptor
    A_AMPA      (1)
    B_AMPA      (1)
    gmax_AMPA   (nS)
    : NMDA Receptor
    A_NMDA      (1)
    B_NMDA      (1)
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use_GB      (1)
    : VDCC (R-type)
    m_VDCC      (1)
    h_VDCC      (1)
    : Postsynaptic Ca2+ dynamics
    cai_CR      (mM)        <1e-6>
    : Long-term synaptic plasticity
    rho_GB      (1)
    effcai_GB   (us/liter)  <1e-3>
}

INITIAL{
    : AMPA Receptor
    A_AMPA      = 0
    B_AMPA      = 0
    gmax_AMPA   = gmax0_AMPA
    : NMDA Receptor
    A_NMDA      = 0
    B_NMDA      = 0
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use_GB      = Use
    : Postsynaptic Ca2+ dynamics
    cai_CR      = min_ca_CR
    : Long-term synaptic plasticity
    rho_GB      = rho0_GB
    effcai_GB   = 0
    dep_GB      = 0
    pot_GB      = 0
    : Delayed connection
    next_delay = -1

    : Initialize watchers
    net_send(0, 1)
}

PROCEDURE setup_delay_vecs() {
VERBATIM
#ifndef CORENEURON_BUILD
    void** vv_delay_times = (void**)(&_p_delay_times);
    void** vv_delay_weights = (void**)(&_p_delay_weights);
    *vv_delay_times = (void*)NULL;
    *vv_delay_weights = (void*)NULL;
    if (ifarg(1)) {
        *vv_delay_times = vector_arg(1);
    }
    if (ifarg(2)) {
        *vv_delay_weights = vector_arg(2);
    }
#endif
ENDVERBATIM
}


BREAKPOINT {
    LOCAL Eca_syn, mggate, i_AMPA, i_NMDA, Pf_NMDA, gca_bar_abs_VDCC, gca_VDCC
    SOLVE state METHOD euler
    : AMPA Receptor
    g_AMPA = (1e-3)*gmax_AMPA*(B_AMPA - A_AMPA)
    i_AMPA = g_AMPA*(v-E_AMPA)
    : NMDA Receptor
    mggate = 1 / (1 + exp(-slope_NMDA*v) * (mg/scale_NMDA))
    g_NMDA = (1e-3)*gmax_NMDA*mggate*(B_NMDA - A_NMDA)
    i_NMDA = g_NMDA*(v - E_NMDA)
    : NMDAR-mediated calcium current
    Pf_NMDA  = (4*cao_CR) / (4*cao_CR + (1/1.38) * 120 (mM)) * 0.6
    ica_NMDA = Pf_NMDA*g_NMDA*(v-40.0)
    : VDCC (R-type), assuming sphere for spine head
    gca_bar_abs_VDCC = gca_bar_VDCC * 4(um2)*PI*(3(1/um3)/4*volume_CR*1/PI)^(2/3)
    gca_VDCC = (1e-3) * gca_bar_abs_VDCC * m_VDCC * m_VDCC * h_VDCC
    Eca_syn = nernst(cai_CR, cao_CR, 2)
    ica_VDCC = gca_VDCC*(v-Eca_syn)
    : Update synaptic voltage (for recording convenience)
    vsyn = v
    : Update current
    i = i_AMPA + i_NMDA + ica_VDCC
}


DERIVATIVE state {
    LOCAL minf_VDCC, hinf_VDCC
    : AMPA Receptor
    A_AMPA'      = - A_AMPA/tau_r_AMPA
    B_AMPA'      = - B_AMPA/tau_d_AMPA
    gmax_AMPA'   = (gmax_d_AMPA + rho_GB*(gmax_p_AMPA - gmax_d_AMPA) - gmax_AMPA) / ((1e3)*tau_exp_GB)
    : NMDA Receptor
    A_NMDA'      = - A_NMDA/tau_r_NMDA
    B_NMDA'      = - B_NMDA/tau_d_NMDA
    : Stochastic Tsodyks-Markram Multi-Vesicular Release
    Use_GB'         = (Use_d + rho_GB*(Use_p - Use_d) - Use_GB) / ((1e3)*tau_exp_GB)
    : VDCC (R-type)
    minf_VDCC    = 1 / (1 + exp(((vhm_VDCC - ljp_VDCC) - v) / km_VDCC))
    hinf_VDCC    = 1 / (1 + exp(((vhh_VDCC - ljp_VDCC) - v) / kh_VDCC))
    m_VDCC'      = (minf_VDCC-m_VDCC)/mtau_VDCC
    h_VDCC'      = (hinf_VDCC-h_VDCC)/htau_VDCC
    : Long-term synaptic plasticity
    effcai_GB'   = - effcai_GB/tau_effca_GB + (cai_CR - min_ca_CR)
    rho_GB'      = ( - rho_GB*(1 - rho_GB)*(rho_star_GB - rho_GB)
                     + pot_GB*gamma_p_GB*(1 - rho_GB)
                     - dep_GB*gamma_d_GB*rho_GB ) / ((1e3)*tau_ind_GB)
    : Postsynaptic Ca2+ dynamics
    cai_CR'      = - (1e-9)*(ica_NMDA + ica_VDCC)*gamma_ca_CR/((1e-15)*volume_CR*2*FARADAY)
                   - (cai_CR - min_ca_CR)/tau_ca_CR
}


NET_RECEIVE (weight, u, tsyn (ms), recovered, unrecovered, nc_type) {
    : nc_type: 0=presynaptic netcon, 1=spontmini, 2=replay
    LOCAL p_rec, released, tp, factor, rec
    INITIAL {
        weight = 1
        u = 0
        tsyn = 0 (ms)
        if (init_depleted){
            recovered = 0
            unrecovered = Nrrp
        } else {
            recovered = Nrrp
            unrecovered = 0
        }
        if (nc_type == 0) {   : pre-synaptic netcon
    VERBATIM
            // setup self events for delayed connections to change weights
            IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
            IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
            if (vv_delay_times && vector_capacity(vv_delay_times)>=1) {
                double* deltm_el = vector_vec(vv_delay_times);
                int delay_times_idx;
                next_delay = 0;
                for(delay_times_idx = 0; delay_times_idx < vector_capacity(vv_delay_times); ++delay_times_idx) {
                    double next_delay_t = deltm_el[delay_times_idx];
    ENDVERBATIM
                    net_send(next_delay_t, 10)  : use flag 10 to avoid interfering with GluSynapse logic
    VERBATIM
                }
            }
    ENDVERBATIM
        }
    }


    if(verbose > 0){ UNITSOFF printf("Time = %g ms, incoming spike at synapse %g\n", t, synapseID) UNITSON }
    if(flag == 0) {
        if(weight <= 0){
            : Do not perform any calculations if the synapse (netcon) is deactivated.
            : This avoids drawing from the random stream
            : WARNING In this model *weight* is only used to activate/deactivate the
            :         synapse. The conductance is stored in gmax_AMPA and gmax_NMDA.
            if(verbose > 0){ printf("Inactive synapse, weight = %g\n", weight) }
        } else {
            : Flag 0: Regular spike
            if(verbose > 0){ printf("Flag 0, Regular spike\n") }
            : Update facilitation variable as Eq. 2 in Fuhrmann et al. 2002
            u = Use_GB + u*(1 - Use_GB)*exp(-(t - tsyn)/Fac)
            if ( verbose > 0 ) { printf("\tVesicle release probability = %g\n", u) }
            : Recovery
            p_rec = 1 - exp(-(t - tsyn)/Dep)
            if ( verbose > 0 ) { printf("\tVesicle recovery probability = %g\n", p_rec) }
            if ( verbose > 0 ) { printf("\tVesicle available before recovery = %g\n", recovered) }
            recovered = recovered + brand(unrecovered, p_rec)
            if ( verbose > 0 ) { printf("\tVesicles available after recovery = %g\n", recovered) }
            : Release
            rec = recovered  : Make a copy so we can change it for single vesicle minis w/o messing with recovered
            : Consider only a single recovered vesicle for minis (if minis_single_vesicle flag is set to 1)
            if (rec > 1 && minis_single_vesicle && nc_type == 1) { rec = 1 }
            released = brand(rec, u)
            if ( verbose > 0 ) { printf("\tReleased %g vesicles out of %g\n", released, recovered) }
            : Update vesicle pool
            recovered = recovered - released
            unrecovered = Nrrp - recovered
            if ( verbose > 0 ) { printf("\tFinal vesicle count, Recovered = %g, Unrecovered = %g, Nrrp = %g\n", recovered, unrecovered, Nrrp) }
            : Update AMPA variables
            tp = (tau_r_AMPA*tau_d_AMPA)/(tau_d_AMPA-tau_r_AMPA)*log(tau_d_AMPA/tau_r_AMPA)  : Time to peak
            factor = 1 / (-exp(-tp/tau_r_AMPA)+exp(-tp/tau_d_AMPA))  : Normalization factor
            A_AMPA = A_AMPA + released/Nrrp*factor
            B_AMPA = B_AMPA + released/Nrrp*factor
            : Update NMDA variables
            tp = (tau_r_NMDA*tau_d_NMDA)/(tau_d_NMDA-tau_r_NMDA)*log(tau_d_NMDA/tau_r_NMDA)  : Time to peak
            factor = 1 / (-exp(-tp/tau_r_NMDA)+exp(-tp/tau_d_NMDA))  : Normalization factor
            A_NMDA = A_NMDA + released/Nrrp*factor
            B_NMDA = B_NMDA + released/Nrrp*factor
            : Update tsyn
            : tsyn knows about all spikes, not only those that released
            : i.e. each spike can increase the u, regardless of recovered state
            :      and each spike trigger an evaluation of recovery
            tsyn = t
        }
    } else if(flag == 1) {
        : Flag 1, Initialize watchers
        if(verbose > 0){ printf("Flag 1, Initialize watchers\n") }
        WATCH (effcai_GB > theta_d_GB) 2
        WATCH (effcai_GB < theta_d_GB) 3
        WATCH (effcai_GB > theta_p_GB) 4
        WATCH (effcai_GB < theta_p_GB) 5
    } else if(flag == 2) {
        : Flag 2, Activate depression mechanisms
        if(verbose > 0){ printf("Flag 2, Activate depression mechanisms\n") }
        dep_GB = 1
    } else if(flag == 3) {
        : Flag 3, Deactivate depression mechanisms
        if(verbose > 0){ printf("Flag 3, Deactivate depression mechanisms\n") }
        dep_GB = 0
    } else if(flag == 4) {
        : Flag 4, Activate potentiation mechanisms
        if(verbose > 0){ printf("Flag 4, Activate potentiation mechanisms\n") }
        pot_GB = 1
    } else if(flag == 5) {
        : Flag 5, Deactivate potentiation mechanisms
        if(verbose > 0){ printf("Flag 5, Deactivate potentiation mechanisms\n") }
        pot_GB = 0
    } else if(flag == 10) {
        : Flag 10, Handle delayed connection weight changes
    VERBATIM
        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {
            double* weights_v = vector_vec(vv_delay_weights);
            double next_delay_weight = weights_v[(int)next_delay];
    ENDVERBATIM
            weight = conductance * next_delay_weight
            next_delay = next_delay + 1
    VERBATIM
        }
    ENDVERBATIM
    }
}

FUNCTION nernst(ci(mM), co(mM), z) (mV) {
    nernst = (1000) * R * (celsius + 273.15) / (z*FARADAY) * log(co/ci)
    if(verbose > 1) { UNITSOFF printf("nernst:%g R:%g temperature (c):%g \n", nernst, R, celsius) UNITSON }
}

PROCEDURE setRNG() {
    VERBATIM
    #ifndef CORENEURON_BUILD
    // For compatibility, allow for either MCellRan4 or Random123
    // Distinguish by the arg types
    // Object => MCellRan4, seeds (double) => Random123
    usingR123 = 0;
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        uint32_t a2 = 0;
        uint32_t a3 = 0;
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        if (ifarg(2)) {
            a2 = (uint32_t)*getarg(2);
        }
        if (ifarg(3)) {
            a3 = (uint32_t)*getarg(3);
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), a2, a3);
        usingR123 = 1;
    } else if( ifarg(1) ) {   // not a double, so assume hoc object type
        void** pv = (void**)(&_p_rng);
        *pv = nrn_random_arg(1);
    } else {  // no arg, so clear pointer
        void** pv = (void**)(&_p_rng);
        *pv = (void*)0;
    }
    #endif
    ENDVERBATIM
}


PROCEDURE clearRNG() {
VERBATIM
    #ifndef CORENEURON_BUILD
    if (usingR123) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
    } else {
        void** pv = (void**)(&_p_rng);
        if (*pv) {
            *pv = (void*)0;
        }
    }
    #endif
ENDVERBATIM
}


FUNCTION urand() {
    VERBATIM
    double value;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(RANDCAST _p_rng);
        #endif
    } else {
        value = 0.0;
    }
    _lurand = value;
    ENDVERBATIM
}

FUNCTION brand(n, p) {
    LOCAL result, count, success
    success = 0
    FROM count = 0 TO (n - 1) {
        result = urand()
        if(result <= p) {
            success = success + 1
        }
    }
    brand = success
}


FUNCTION bbsavestate() {
    bbsavestate = 0
    VERBATIM
    #ifndef CORENEURON_BUILD
        /* first arg is direction (0 save, 1 restore), second is array*/
        /* if first arg is -1, fill xdir with the size of the array */
        double *xdir, *xval;
        #ifndef NRN_VERSION_GTEQ_8_2_0
        double *hoc_pgetarg();
        long nrn_get_random_sequence(void* r);
        void nrn_set_random_sequence(void* r, int val);
        #endif
        xdir = hoc_pgetarg(1);
        xval = hoc_pgetarg(2);
        if (_p_rng) {
            // tell how many items need saving
            if (*xdir == -1) {  // count items
                if( usingR123 ) {
                    *xdir = 2.0;
                } else {
                    *xdir = 1.0;
                }
                return 0.0;
            } else if(*xdir ==0 ) {  // save
                if( usingR123 ) {
                    uint32_t seq;
                    char which;
                    nrnran123_getseq( (nrnran123_State*)_p_rng, &seq, &which );
                    xval[0] = (double) seq;
                    xval[1] = (double) which;
                } else {
                    xval[0] = (double)nrn_get_random_sequence(RANDCAST _p_rng);
                }
            } else {  // restore
                if( usingR123 ) {
                    nrnran123_setseq( (nrnran123_State*)_p_rng, (uint32_t)xval[0], (char)xval[1] );
                } else {
                    nrn_set_random_sequence(RANDCAST _p_rng, (long)(xval[0]));
                }
            }
        }
    #endif
    ENDVERBATIM
}


VERBATIM
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
    IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

    // make sure offset array non-null
    if (iArray) {
        // get handle to random123 instance
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        // get location for storing ids
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        // retrieve/store identifier seeds
        nrnran123_getids3(*pv, ia, ia+1, ia+2);
        // retrieve/store stream sequence
        char which;
        nrnran123_getseq(*pv, ia+3, &which);
        ia[4] = (int)which;
    }

    // increment integer offset (2 identifier), no double data
    *ioffset += 5;
    *doffset += 0;

    // serialize connection delay vectors
    if (vv_delay_times && vv_delay_weights &&
       (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {
        if (iArray) {
            uint32_t* di = ((uint32_t*)iArray) + *ioffset;
            // store vector sizes for deserialization
            di[0] = vector_capacity(vv_delay_times);
            di[1] = vector_capacity(vv_delay_weights);
        }
        if (dArray) {
            double* delay_times_el = vector_vec(vv_delay_times);
            double* delay_weights_el = vector_vec(vv_delay_weights);
            double* x_i = dArray + *doffset;
            int delay_vecs_idx;
            int x_idx = 0;
            for(delay_vecs_idx = 0; delay_vecs_idx < vector_capacity(vv_delay_times); ++delay_vecs_idx) {
                 x_i[x_idx++] = delay_times_el[delay_vecs_idx];
                 x_i[x_idx++] = delay_weights_el[delay_vecs_idx];
            }
        }
        // reserve space for connection delay data on serialization buffer
        *doffset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);
    } else {
        if (iArray) {
            uint32_t* di = ((uint32_t*)iArray) + *ioffset;
            di[0] = 0;
            di[1] = 0;
        }
    }
    // reserve space for delay vectors (may be 0)
    *ioffset += 2;
}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
    // make sure it's not previously set
    assert(!_p_rng);
    assert(!_p_delay_times && !_p_delay_weights);

    uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
    // make sure non-zero identifier seeds
    if (ia[0] != 0 || ia[1] != 0 || ia[2] != 0) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
        // get new stream
        *pv = nrnran123_newstream3(ia[0], ia[1], ia[2]);
        // restore sequence
        nrnran123_setseq(*pv, ia[3], (char)ia[4]);
    }
    // increment intger offset (2 identifiers), no double data
    *ioffset += 5;

    int delay_times_sz = iArray[5];
    int delay_weights_sz = iArray[6];
    *ioffset += 2;

    if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {
        double* x_i = dArray + *doffset;

        // allocate vectors
        _p_delay_times = (double*)vector_new1(delay_times_sz);
        _p_delay_weights = (double*)vector_new1(delay_weights_sz);

        double* delay_times_el = vector_vec((IvocVect*)_p_delay_times);
        double* delay_weights_el = vector_vec((IvocVect*)_p_delay_weights);

        // copy data
        int x_idx;
        int vec_idx = 0;
        for(x_idx = 0; x_idx < delay_times_sz + delay_weights_sz; x_idx += 2) {
            delay_times_el[vec_idx] = x_i[x_idx];
            delay_weights_el[vec_idx++] = x_i[x_idx+1];
        }
        *doffset += delay_times_sz + delay_weights_sz;
    }
}
ENDVERBATIM
