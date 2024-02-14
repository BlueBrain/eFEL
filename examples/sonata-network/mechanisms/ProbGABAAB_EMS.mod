COMMENT
/**
 * @file ProbGABAAB.mod
 * @brief
 * @author king, muller
 * @date 2011-08-17
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
ENDCOMMENT

TITLE GABAAB receptor with presynaptic short-term plasticity


COMMENT
GABAA receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al, 2002
Implemented by Srikanth Ramaswamy, Blue Brain Project, March 2009

_EMS (Eilif Michael Srikanth)
Modification of ProbGABAA: 2-State model by Eilif Muller, Michael Reimann, Srikanth Ramaswamy, Blue Brain Project, August 2011
This new model was motivated by the following constraints:

1) No consumption on failure.
2) No release just after release until recovery.
3) Same ensemble averaged trace as deterministic/canonical Tsodyks-Markram
   using same parameters determined from experiment.
4) Same quantal size as present production probabilistic model.

To satisfy these constaints, the synapse is implemented as a
uni-vesicular (generalization to multi-vesicular should be
straight-forward) 2-state Markov process.  The states are
{1=recovered, 0=unrecovered}.

For a pre-synaptic spike or external spontaneous release trigger
event, the synapse will only release if it is in the recovered state,
and with probability u (which follows facilitation dynamics).  If it
releases, it will transition to the unrecovered state.  Recovery is as
a Poisson process with rate 1/Dep.

This model satisys all of (1)-(4).

ENDCOMMENT


NEURON {
    THREADSAFE
    POINT_PROCESS ProbGABAAB_EMS

    RANGE tau_r_GABAA, tau_d_GABAA
    GLOBAL tau_r_GABAB, tau_d_GABAB
    RANGE Use, u, Dep, Fac, u0, tsyn
    RANGE unoccupied, occupied, Nrrp

    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB, GABAB_ratio
    RANGE A_GABAA_step, B_GABAA_step, A_GABAB_step, B_GABAB_step

    NONSPECIFIC_CURRENT i
    BBCOREPOINTER rng
    RANGE synapseID, selected_for_report, verboseLevel, conductance
    RANGE next_delay
    BBCOREPOINTER delay_times, delay_weights
    GLOBAL nc_type_param
    GLOBAL minis_single_vesicle
    GLOBAL init_depleted

    :RANGE sgid, tgid  : For debugging
}

PARAMETER {
    tau_r_GABAA = 0.2   (ms)  : dual-exponential conductance profile
    tau_d_GABAA = 8     (ms)  : IMPORTANT: tau_r < tau_d
    tau_r_GABAB = 3.5   (ms)  : dual-exponential conductance profile :Placeholder value from hippocampal recordings SR
    tau_d_GABAB = 260.9 (ms)  : IMPORTANT: tau_r < tau_d  :Placeholder value from hippocampal recordings
    Use = 1.0           (1)   : Utilization of synaptic efficacy (just initial values! Use, Dep and Fac are overwritten by BlueBuilder assigned values)
    Dep = 100           (ms)  : relaxation time constant from depression
    Fac = 10            (ms)  : relaxation time constant from facilitation
    e_GABAA = -80       (mV)  : GABAA reversal potential
    e_GABAB = -97       (mV)  : GABAB reversal potential
    gmax = .001         (uS)  : weight conversion factor (from nS to uS)
    u0   = 0                  : initial value of u, which is the running value of release probability
    Nrrp = 1            (1)   : Number of total release sites for given contact

    synapseID = 0
    verboseLevel = 0
    selected_for_report = 0
    GABAB_ratio = 0     (1)   : The ratio of GABAB to GABAA
    conductance = 0.0
    nc_type_param = 4
    minis_single_vesicle = 0   :// 0 - no limit (old behavior)
    init_depleted = 0          :// 0 - init full (old behavior)
}

COMMENT
The Verbatim block is needed to generate random nos. from a uniform distribution between 0 and 1
for comparison with Pr to decide whether to activate the synapse or not
ENDCOMMENT

VERBATIM

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#ifndef NRN_VERSION_GTEQ_8_2_0
#include "nrnran123.h"

#ifndef CORENEURON_BUILD
extern int ifarg(int iarg);

extern void* vector_arg(int iarg);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
#endif

double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif

ENDVERBATIM


ASSIGNED {
        v (mV)
        i (nA)
        i_GABAA (nA)
        i_GABAB (nA)
        g_GABAA (uS)
        g_GABAB (uS)
        g (uS)
        factor_GABAA
        factor_GABAB
        A_GABAA_step
        B_GABAA_step
        A_GABAB_step
        B_GABAB_step

        rng
        usingR123            : TEMPORARY until mcellran4 completely deprecated

        : MVR
        unoccupied (1) : no. of unoccupied sites following release event
        occupied   (1) : no. of occupied sites following one epoch of recovery
        tsyn (ms) : the time of the last spike
        u (1) : running release probability

        : stuff for delayed connections
        delay_times
        delay_weights
        next_delay (ms)
}

STATE {
        A_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAA
        B_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAA
        A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB
        B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB
}


INITIAL {
        LOCAL tp_GABAA, tp_GABAB

        tsyn = 0
        u=u0

        : MVR
        if ( init_depleted ) {
            unoccupied = Nrrp
            occupied = 0
        } else {
            unoccupied = 0
            occupied = Nrrp
        }

        A_GABAA = 0
        B_GABAA = 0

        A_GABAB = 0
        B_GABAB = 0

        tp_GABAA = (tau_r_GABAA*tau_d_GABAA)/(tau_d_GABAA-tau_r_GABAA)*log(tau_d_GABAA/tau_r_GABAA) :time to peak of the conductance
        tp_GABAB = (tau_r_GABAB*tau_d_GABAB)/(tau_d_GABAB-tau_r_GABAB)*log(tau_d_GABAB/tau_r_GABAB) :time to peak of the conductance

        factor_GABAA = -exp(-tp_GABAA/tau_r_GABAA)+exp(-tp_GABAA/tau_d_GABAA) :GABAA Normalization factor - so that when t = tp_GABAA, gsyn = gpeak
        factor_GABAA = 1/factor_GABAA

        factor_GABAB = -exp(-tp_GABAB/tau_r_GABAB)+exp(-tp_GABAB/tau_d_GABAB) :GABAB Normalization factor - so that when t = tp_GABAB, gsyn = gpeak
        factor_GABAB = 1/factor_GABAB

        A_GABAA_step = exp(dt*(( - 1.0 ) / tau_r_GABAA))
        B_GABAA_step = exp(dt*(( - 1.0 ) / tau_d_GABAA))
        A_GABAB_step = exp(dt*(( - 1.0 ) / tau_r_GABAB))
        B_GABAB_step = exp(dt*(( - 1.0 ) / tau_d_GABAB))

    VERBATIM
        if( usingR123 ) {
            nrnran123_setseq((nrnran123_State*)_p_rng, 0, 0);
        }
    ENDVERBATIM

        next_delay = -1

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
        SOLVE state

        g_GABAA = gmax*(B_GABAA-A_GABAA) :compute time varying conductance as the difference of state variables B_GABAA and A_GABAA
        g_GABAB = gmax*(B_GABAB-A_GABAB) :compute time varying conductance as the difference of state variables B_GABAB and A_GABAB
        g = g_GABAA + g_GABAB
        i_GABAA = g_GABAA*(v-e_GABAA) :compute the GABAA driving force based on the time varying conductance, membrane potential, and GABAA reversal
        i_GABAB = g_GABAB*(v-e_GABAB) :compute the GABAB driving force based on the time varying conductance, membrane potential, and GABAB reversal
        i = i_GABAA + i_GABAB
}

PROCEDURE state() {
        A_GABAA = A_GABAA*A_GABAA_step
        B_GABAA = B_GABAA*B_GABAA_step
        A_GABAB = A_GABAB*A_GABAB_step
        B_GABAB = B_GABAB*B_GABAB_step
}



NET_RECEIVE (weight, weight_GABAA, weight_GABAB, Psurv, nc_type) {
    : Psurv - survival probability of unrecovered state
    : nc_type:
    :   0 = presynaptic netcon
    :   1 = spontmini netcon
    :   2 = replay netcon

    LOCAL result, ves, occu
    weight_GABAA = weight
    weight_GABAB = weight * GABAB_ratio

    INITIAL {
        if (nc_type == 0) {  :// presynaptic netcon
    VERBATIM
            // setup self events for delayed connections to change weights
            IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
            IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

            if (vv_delay_times && vector_capacity(vv_delay_times)>=1) {
                double* deltm_el = vector_vec(vv_delay_times);
                int delay_times_idx;
                next_delay = 0;
                for (delay_times_idx = 0; delay_times_idx < vector_capacity(vv_delay_times); ++delay_times_idx) {
                    double next_delay_t = deltm_el[delay_times_idx];
    ENDVERBATIM
                    net_send(next_delay_t, 1)
    VERBATIM
                }
            }
    ENDVERBATIM
        }
    }

    if (flag == 1) {  :// self event to set next weight at
        :UNITSOFF
        :    printf( "gaba: self event at synapse %f time %g\n", synapseID, t)
        :UNITSON
    VERBATIM
        // setup self events for delayed connections to change weights
        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {
            double* weights_v = vector_vec(vv_delay_weights);
            double next_delay_weight = weights_v[(int)next_delay];
    ENDVERBATIM
            weight = conductance * next_delay_weight
            next_delay = next_delay + 1
    VERBATIM
        }
        return;
    ENDVERBATIM
    }

    : [flag == 0] Handle a spike which arrived
    :UNITSOFF
    :printf("[Syn %.0f] Received! (%f -> %f) with weight %g at time %g\n", synapseID, sgid, tgid, weight, t)
    :UNITSON

    : Do not perform any calculations if the synapse (netcon) is deactivated. This avoids drawing from
    : random number stream. Also, disable in case of t < 0 (in case of ForwardSkip) which causes numerical
    : instability if synapses are activated.
    if ( weight <= 0 || t < 0 ) {
    VERBATIM
        return;
    ENDVERBATIM
    }

    : calc u at event-
    if (Fac > 0) {
        u = u * exp(-(t - tsyn)/Fac)  :// update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    } else {
        u = Use
    }
    if(Fac > 0){
        u = u + Use*(1-u)  :// update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    }

    : recovery
    FROM counter = 0 TO (unoccupied - 1) {
        : Iterate over all unoccupied sites and compute how many recover
        Psurv = exp(-(t-tsyn)/Dep)
        result = urand()
        if (result>Psurv) {
            occupied = occupied + 1     :// recover a previously unoccupied site
            if ( verboseLevel > 0 ) {
                UNITSOFF
                printf("[Syn %.0f] Recovered! t = %g, Psurv = %g, urand = %g\n", synapseID, t, Psurv, result)
                UNITSON
            }
        }
    }

    ves = 0                  :// Initialize the number of released vesicles to 0
    occu = occupied          :// Make a copy, so we can update occupied in the loop
    if (occu > 1 && minis_single_vesicle && nc_type == 1) {    : // if nc_type is spont_mini consider single vesicle
        occu = 1
    }
    FROM counter = 0 TO (occu - 1) {
        : iterate over all occupied sites and compute how many release
        result = urand()
        if (result < u) {
            : release a single site!
            occupied = occupied - 1  :// decrease the number of occupied sites by 1
            ves = ves + 1            :// increase number of relesed vesicles by 1
        }
    }

    : Update number of unoccupied sites
    unoccupied = Nrrp - occupied

    : Update tsyn
    : tsyn knows about all spikes, not only those that released
    : i.e. each spike can increase the u, regardless of recovered state.
    :      and each spike trigger an evaluation of recovery
    tsyn = t

    if (ves > 0) { :no need to evaluate unless we have vesicle release
        A_GABAA = A_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA
        B_GABAA = B_GABAA + ves/Nrrp*weight_GABAA*factor_GABAA
        A_GABAB = A_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB
        B_GABAB = B_GABAB + ves/Nrrp*weight_GABAB*factor_GABAB

        if ( verboseLevel > 0 ) {
            UNITSOFF
            printf("[Syn %.0f] Release! t = %g, vals: %g %g %g %g\n",
                   synapseID, t, A_GABAA, weight_GABAA, factor_GABAA, weight)
            UNITSON
        }

    } else {
        : total release failure
        if ( verboseLevel > 0 ) {
            UNITSOFF
            printf("[Syn %.0f] Failure! t = %g: urand = %g\n", synapseID, t, result)
            UNITSON
        }
    }
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
    double value = 0.0;
    if ( usingR123 ) {
        value = nrnran123_dblpick((nrnran123_State*)_p_rng);
    } else if (_p_rng) {
        #ifndef CORENEURON_BUILD
        value = nrn_random_pick(RANDCAST _p_rng);
        #endif
    } else {
        // Note: prior versions used scop_random(1), but since we never use this model without configuring the rng.  Maybe should throw error?
        value = 0.0;
    }
    _lurand = value;
ENDVERBATIM
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

FUNCTION toggleVerbose() {
    verboseLevel = 1 - verboseLevel
}


VERBATIM
static void bbcore_write(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {
  IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
  IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

  if (d) {
    uint32_t* di = ((uint32_t*)d) + *d_offset;
    nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
    nrnran123_getids3(*pv, di, di+1, di+2);

    // write strem sequence
    char which;
    nrnran123_getseq(*pv, di+3, &which);
    di[4] = (int)which;
    //printf("ProbGABAAB_EMS bbcore_write %d %d %d\n", di[0], di[1], di[2]);

  }
  // reserve random123 parameters on serialization buffer
  *d_offset += 5;

  // serialize connection delay vectors
  if (vv_delay_times && vv_delay_weights &&
     (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {
    if (d) {
      uint32_t* di = ((uint32_t*)d) + *d_offset;
      // store vector sizes for deserialization
      di[0] = vector_capacity(vv_delay_times);
      di[1] = vector_capacity(vv_delay_weights);
    }
    if (x) {
      double* delay_times_el = vector_vec(vv_delay_times);
      double* delay_weights_el = vector_vec(vv_delay_weights);
      double* x_i = x + *x_offset;
      int delay_vecs_idx;
      int x_idx = 0;
      for(delay_vecs_idx = 0; delay_vecs_idx < vector_capacity(vv_delay_times); ++delay_vecs_idx) {
         x_i[x_idx++] = delay_times_el[delay_vecs_idx];
         x_i[x_idx++] = delay_weights_el[delay_vecs_idx];
      }
    }
    // reserve space for connection delay data on serialization buffer
    *x_offset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);
  } else {
    if (d) {
      uint32_t* di = ((uint32_t*)d) + *d_offset;
      di[0] = 0;
      di[1] = 0;
    }

  }
  // reserve space for delay vectors (may be 0)
  *d_offset += 2;

}

static void bbcore_read(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {
  // deserialize random123 data
  uint32_t* di = ((uint32_t*)d) + *d_offset;
  if (di[0] != 0 || di[1] != 0 || di[2] != 0) {
      nrnran123_State** pv = (nrnran123_State**)(&_p_rng);
#if !NRNBBCORE
      if(*pv) {
          nrnran123_deletestream(*pv);
      }
#endif
      *pv = nrnran123_newstream3(di[0], di[1], di[2]);
      // restore stream sequence
      char which = (char)di[4];
      nrnran123_setseq(*pv, di[3], which);
  }
  //printf("ProbAMPANMDA_EMS bbcore_read %d %d %d\n", di[0], di[1], di[2]);

  int delay_times_sz = di[5];
  int delay_weights_sz = di[6];
  *d_offset += 7;

  if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {
    double* x_i = x + *x_offset;

    // allocate vectors
    if (!_p_delay_times) {
      _p_delay_times = (double*)vector_new1(delay_times_sz);
    }
    assert(delay_times_sz == vector_capacity((IvocVect*)_p_delay_times));
    if (!_p_delay_weights) {
      _p_delay_weights = (double*)vector_new1(delay_weights_sz);
    }
    assert(delay_weights_sz == vector_capacity((IvocVect*)_p_delay_weights));

    double* delay_times_el = vector_vec((IvocVect*)_p_delay_times);
    double* delay_weights_el = vector_vec((IvocVect*)_p_delay_weights);

    // copy data
    int x_idx;
    int vec_idx = 0;
    for(x_idx = 0; x_idx < delay_times_sz + delay_weights_sz; x_idx += 2) {
      delay_times_el[vec_idx] = x_i[x_idx];
      delay_weights_el[vec_idx++] = x_i[x_idx+1];
    }
    *x_offset += delay_times_sz + delay_weights_sz;

  }
}
ENDVERBATIM
