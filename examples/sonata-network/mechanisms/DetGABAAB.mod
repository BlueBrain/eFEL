COMMENT
/**
 * @file DetGABAAB.mod
 * @brief Adapted from ProbGABAA_EMS.mod by Eilif, Michael and Srikanth
 * @author chindemi
 * @date 2014-05-25
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


TITLE GABAA and GABAB receptor with presynaptic short-term plasticity


COMMENT
GABAA and GABAB receptor conductance using a dual-exponential profile
presynaptic short-term plasticity based on Fuhrmann et al. 2002, deterministic
version.
ENDCOMMENT


NEURON {
    THREADSAFE

    POINT_PROCESS DetGABAAB
    RANGE tau_r_GABAA, tau_d_GABAA, tau_r_GABAB, tau_d_GABAB
    RANGE Use, u, Dep, Fac, u0, GABAB_ratio
    RANGE i, i_GABAA, i_GABAB, g_GABAA, g_GABAB, g, e_GABAA, e_GABAB
    NONSPECIFIC_CURRENT i
    RANGE synapseID, verboseLevel
    RANGE conductance
    RANGE next_delay
    BBCOREPOINTER delay_times, delay_weights
    GLOBAL nc_type_param
    : For debugging
    :RANGE sgid, tgid
}


PARAMETER {
    tau_r_GABAA  = 0.2   (ms) : dual-exponential conductance profile
    tau_d_GABAA  = 8     (ms) : IMPORTANT: tau_r < tau_d
    tau_r_GABAB  = 3.5   (ms) : dual-exponential conductance profile :Placeholder value from hippocampal recordings SR
    tau_d_GABAB  = 260.9 (ms) : IMPORTANT: tau_r < tau_d  :Placeholder value from hippocampal recordings
    Use          = 1.0   (1)  : Utilization of synaptic efficacy
    Dep          = 100   (ms) : relaxation time constant from depression
    Fac          = 10    (ms) :  relaxation time constant from facilitation
    e_GABAA      = -80   (mV) : GABAA reversal potential
    e_GABAB      = -97   (mV) : GABAB reversal potential
    gmax         = .001  (uS) : weight conversion factor (from nS to uS)
    u0           = 0          :initial value of u, which is the running value of release probability
    synapseID    = 0
    verboseLevel = 0
    GABAB_ratio  = 0     (1)  : The ratio of GABAB to GABAA
    conductance  = 0.0
    nc_type_param = 7
}


VERBATIM

#ifndef CORENEURON_BUILD
extern int ifarg(int iarg);
#ifndef NRN_VERSION_GTEQ_8_2_0
extern void* vector_arg(int iarg);
extern double* vector_vec(void* vv);
extern int vector_capacity(void* vv);
#endif
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

    : stuff for delayed connections
    delay_times
    delay_weights
    next_delay (ms)
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


STATE {
    A_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAA
    B_GABAA       : GABAA state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAA
    A_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_r_GABAB
    B_GABAB       : GABAB state variable to construct the dual-exponential profile - decays with conductance tau_d_GABAB
}


INITIAL{
    LOCAL tp_GABAA, tp_GABAB

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

    next_delay = -1

}


BREAKPOINT {
    SOLVE state METHOD cnexp
    g_GABAA = gmax*(B_GABAA-A_GABAA) :compute time varying conductance as the difference of state variables B_GABAA and A_GABAA
    g_GABAB = gmax*(B_GABAB-A_GABAB) :compute time varying conductance as the difference of state variables B_GABAB and A_GABAB
    g = g_GABAA + g_GABAB
    i_GABAA = g_GABAA*(v-e_GABAA) :compute the GABAA driving force based on the time varying conductance, membrane potential, and GABAA reversal
    i_GABAB = g_GABAB*(v-e_GABAB) :compute the GABAB driving force based on the time varying conductance, membrane potential, and GABAB reversal
    i = i_GABAA + i_GABAB
}


DERIVATIVE state{
    A_GABAA' = -A_GABAA/tau_r_GABAA
    B_GABAA' = -B_GABAA/tau_d_GABAA
    A_GABAB' = -A_GABAB/tau_r_GABAB
    B_GABAB' = -B_GABAB/tau_d_GABAB
}


NET_RECEIVE (weight, weight_GABAA, weight_GABAB, R, Pr, u, tsyn (ms), nc_type){
    LOCAL result
    weight_GABAA = weight
    weight_GABAB = weight * GABAB_ratio

    INITIAL{
        R=1
        u=u0
        tsyn=t

        if (nc_type == 0) {
            : nc_type {
            :   0 = presynaptic netcon
            :   1 = spontmini netcon
            :   2 = replay netcon
            : }
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
                net_send(next_delay_t, 1)
    VERBATIM
              }
            }
    ENDVERBATIM
        }
    }

    : Disable in case of t < 0 (in case of ForwardSkip) which causes numerical
    : instability if synapses are activated.
    if(t < 0 ) {
    VERBATIM
        return;
    ENDVERBATIM
    }

    if (flag == 1) {
        : self event to set next weight at delay
    VERBATIM
        // setup self events for delayed connections to change weights
        IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));
        if (vv_delay_weights && vector_capacity(vv_delay_weights)>=next_delay) {
          double* weights_v = vector_vec(vv_delay_weights);
          double next_delay_weight = weights_v[(int)next_delay];
    ENDVERBATIM
          weight = conductance*next_delay_weight
          next_delay = next_delay + 1
    VERBATIM
        }
        return;
    ENDVERBATIM
    }
    : flag == 0, i.e. a spike has arrived

    : calc u at event-
    if (Fac > 0) {
        u = u*exp(-(t - tsyn)/Fac) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    } else {
        u = Use
    }
    if(Fac > 0){
        u = u + Use*(1-u) :update facilitation variable if Fac>0 Eq. 2 in Fuhrmann et al.
    }

    R  = 1 - (1-R) * exp(-(t-tsyn)/Dep) :Probability R for a vesicle to be available for release, analogous to the pool of synaptic
                                        :resources available for release in the deterministic model. Eq. 3 in Fuhrmann et al.
    Pr  = u * R                         :Pr is calculated as R * u (running value of Use)
    R  = R - u * R                      :update R as per Eq. 3 in Fuhrmann et al.

    if( verboseLevel > 0 ) {
        printf("Synapse %f at time %g: R = %g Pr = %g erand = %g\n", synapseID, t, R, Pr, result )
    }

    tsyn = t

    A_GABAA = A_GABAA + Pr*weight_GABAA*factor_GABAA
    B_GABAA = B_GABAA + Pr*weight_GABAA*factor_GABAA
    A_GABAB = A_GABAB + Pr*weight_GABAB*factor_GABAB
    B_GABAB = B_GABAB + Pr*weight_GABAB*factor_GABAB

    if( verboseLevel > 0 ) {
        printf( " vals %g %g %g %g\n", A_GABAA, weight_GABAA, factor_GABAA, weight )
    }
}


FUNCTION toggleVerbose() {
    verboseLevel = 1-verboseLevel
}


VERBATIM
static void bbcore_write(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {

  IvocVect *vv_delay_times = *((IvocVect**)(&_p_delay_times));
  IvocVect *vv_delay_weights = *((IvocVect**)(&_p_delay_weights));

  // serialize connection delay vectors
  if (vv_delay_times && vv_delay_weights &&
     (vector_capacity(vv_delay_times) >= 1) && (vector_capacity(vv_delay_weights) >= 1)) {
    if (d && x) {
      int* d_i = d + *d_offset;
      // store vector sizes for deserialization
      d_i[0] = vector_capacity(vv_delay_times);
      d_i[1] = vector_capacity(vv_delay_weights);

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
  } else {
    if (d) {
      int* d_i = d + *d_offset;
      d_i[0] = 0;
      d_i[1] = 0;
    }
  }
  // reserve space for delay connection vector sizes on serialization buffer
  *d_offset += 2;

  // reserve space for connection delay data on serialization buffer
  if (vv_delay_times && vv_delay_weights) {
    *x_offset += vector_capacity(vv_delay_times) + vector_capacity(vv_delay_weights);
  }
}

static void bbcore_read(double* x, int* d, int* x_offset, int* d_offset, _threadargsproto_) {
  assert(!_p_delay_times && !_p_delay_weights);

  // first get delay vector sizes
  int* d_i = d + *d_offset;
  int delay_times_sz = d_i[0];
  int delay_weights_sz = d_i[1];
  *d_offset += 2;

  if ((delay_times_sz > 0) && (delay_weights_sz > 0)) {
    double* x_i = x + *x_offset;

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
    *x_offset += delay_times_sz + delay_weights_sz;
  }
}
ENDVERBATIM
