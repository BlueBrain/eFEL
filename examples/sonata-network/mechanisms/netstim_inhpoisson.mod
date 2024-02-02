COMMENT
/**
 * @file netstim_inhpoisson.mod
 * @brief Inhibitory poisson generator by the thinning method.
 * @author Eilif Muller
 * @date 2011-03-16
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
 *  Based on vecstim.mod and netstim2.mod shipped with PyNN. See
 *   Muller, Buesing, Schemmel, Meier (2007). "Spike-Frequency Adapting
 *   Neural Ensembles: Beyond Mean Adaptation and Renewal Theories",
 *   Neural Computation 19:11, 2958-3010. doi:10.1162/neco.2007.19.11.2958
 */
ENDCOMMENT

NEURON {
THREADSAFE
  ARTIFICIAL_CELL InhPoissonStim
  RANGE rmax
  RANGE duration
  BBCOREPOINTER uniform_rng, exp_rng, vecRate, vecTbins
  :THREADSAFE : only true if every instance has its own distinct Random
}

VERBATIM
#if defined(NRN_VERSION_GTEQ)
#if NRN_VERSION_GTEQ(9,0,0)
#define NRN_VERSION_GTEQ_9_0_0
#endif
#endif

#ifndef NRN_VERSION_GTEQ_8_2_0
extern int ifarg(int iarg);
#ifndef CORENEURON_BUILD
extern double* vector_vec(void* vv);
extern void* vector_new1(int _i);
extern int vector_capacity(void* vv);
extern void* vector_arg(int iarg);
double nrn_random_pick(void* r);
#endif
void* nrn_random_arg(int argpos);
#define RANDCAST 
#else
#define RANDCAST (Rand*)
#endif


#ifdef STIM_DEBUG
# define debug_printf(...) printf(__VA_ARGS__)
#else
# define debug_printf(...)
#endif

// constant used to indicate an event triggered after a restore to restart the main event loop
const int POST_RESTORE_RESTART_FLAG = -99;

ENDVERBATIM


PARAMETER {
  interval_min = 1.0  : average spike interval of surrogate Poisson process
  duration	= 1e6 (ms) <0,1e9>   : duration of firing (msec)
}

VERBATIM
#include "nrnran123.h"
ENDVERBATIM

ASSIGNED {
   vecRate
   vecTbins
   index
   curRate
   start (ms)
   event (ms)
   uniform_rng
   exp_rng
   usingR123
   rmax
   activeFlag
}

INITIAL {
   index = 0.
   activeFlag = 0.

   : determine start of spiking.
   VERBATIM
   IvocVect *vvTbins = *((IvocVect**)(&_p_vecTbins));
   double* px;

   if (vvTbins && vector_capacity(vvTbins)>=1) {
     px = vector_vec(vvTbins);
     start = px[0];
     if (start < 0.0) start=0.0;
   }
   else start = 0.0;

   /* first event is at the start
   TODO: This should draw from a more appropriate dist
   that has the surrogate process starting a t=-inf
   */
   event = start;

   /* set curRate */
   IvocVect *vvRate = *((IvocVect**)(&_p_vecRate));
   px = vector_vec(vvRate);

   /* set rmax */
   rmax = 0.0;
   int i;
   for (i=0;i<vector_capacity(vvRate);i++) {
      if (px[i]>rmax) rmax = px[i];
   }

   if (vvRate && vector_capacity(vvRate)>0) {
     curRate = px[0];
   }
   else {
      curRate = 1.0;
      rmax = 1.0;
   }

   /** after discussion with michael : rng streams should be set 0
     * in initial block. this is to make sure if initial block is
     * get called multiple times then the simulation should give the
     * same results. Otherwise this is an issue in coreneuron because
     * finitialized is get called twice in coreneuron (once from
     * neurodamus and then in coreneuron. But in general, initial state
     * should be callable multiple times.
     */
   if (_p_uniform_rng && usingR123) {
     nrnran123_setseq((nrnran123_State*)_p_uniform_rng, 0, 0);
   }
   if (_p_exp_rng && usingR123) {
     nrnran123_setseq((nrnran123_State*)_p_exp_rng, 0, 0);
   }

   ENDVERBATIM
   update_time()
   erand() : for some reason, the first erand() call seems
           : to give implausibly large values, so we discard it
   generate_next_event()
   : stop even producing surrogate events if we are past duration
   if (t+event < start+duration) {
VERBATIM
     debug_printf("InhPoisson: Initial event at t = %6.3f\n", t + event);
ENDVERBATIM
     net_send(event, activeFlag )
   }


}

: This procedure queues the next surrogate event in the
: poisson process (rate=ramx) to be thinned.
PROCEDURE generate_next_event() {
	event = 1000.0/rmax*erand()
	: but not earlier than 0
	if (event < 0) {
		event = 0
	}
}

: Supports multiple rng types: mcellran4, random123
: mcellran4:
: 1st arg: exp_rng
: 2nd arg: uniform_rng
: random123
: 3 exp seeds
: 3 uniform seeds
PROCEDURE setRNGs() {
VERBATIM
{
#ifndef CORENEURON_BUILD
    if( ifarg(1) && hoc_is_double_arg(1) ) {
        nrnran123_State** pv = (nrnran123_State**)(&_p_exp_rng);

        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(1), (uint32_t)*getarg(2), (uint32_t)*getarg(3));

        pv = (nrnran123_State**)(&_p_uniform_rng);
        if (*pv) {
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
        }
        *pv = nrnran123_newstream3((uint32_t)*getarg(4), (uint32_t)*getarg(5), (uint32_t)*getarg(6));

        usingR123 = 1;
    } else if( ifarg(1) ) {
        void** pv = (void**)(&_p_exp_rng);
        *pv = nrn_random_arg(1);

        pv = (void**)(&_p_uniform_rng);
        *pv = nrn_random_arg(2);
        usingR123 = 0;
    } else {
        if( usingR123 ) {
            nrnran123_State** pv = (nrnran123_State**)(&_p_exp_rng);
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
            pv = (nrnran123_State**)(&_p_uniform_rng);
            nrnran123_deletestream(*pv);
            *pv = (nrnran123_State*)0;
            //_p_exp_rng = (nrnran123_State*)0;
            //_p_uniform_rng = (nrnran123_State*)0;
        }
    }
#endif
}
ENDVERBATIM
}


FUNCTION urand() {
VERBATIM
	if (_p_uniform_rng) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.uniform(0,1)
		*/
            if( usingR123 ) {
		_lurand = nrnran123_dblpick((nrnran123_State*)_p_uniform_rng);
            } else {
#ifndef CORENEURON_BUILD
		_lurand = nrn_random_pick(RANDCAST _p_uniform_rng);
#endif
            }
	}else{
  	  hoc_execerror("multithread random in NetStim"," only via hoc Random");
	}
ENDVERBATIM
}

FUNCTION erand() {
VERBATIM
	if (_p_exp_rng) {
		/*
		:Supports separate independent but reproducible streams for
		: each instance. However, the corresponding hoc Random
		: distribution MUST be set to Random.negexp(1)
		*/
            if( usingR123 ) {
		_lerand = nrnran123_negexp((nrnran123_State*)_p_exp_rng);
            } else {
#ifndef CORENEURON_BUILD
		_lerand = nrn_random_pick(RANDCAST _p_exp_rng);
#endif
            }
	}else{
  	  hoc_execerror("multithread random in NetStim"," only via hoc Random");
	}
ENDVERBATIM
}





PROCEDURE setTbins() {
VERBATIM
  #ifndef CORENEURON_BUILD
  IvocVect** vv;
  vv = (IvocVect**)(&_p_vecTbins);
  *vv = (IvocVect*)0;

  if (ifarg(1)) {
    *vv = vector_arg(1);

    /*int size = vector_capacity(*vv);
    int i;
    double* px = vector_vec(*vv);
    for (i=0;i<size;i++) {
      printf("%f ", px[i]);
    }*/
  }
  #endif
ENDVERBATIM
}


PROCEDURE setRate() {
VERBATIM
  #ifndef CORENEURON_BUILD

  IvocVect** vv;
  vv = (IvocVect**)(&_p_vecRate);
  *vv = (IvocVect*)0;

  if (ifarg(1)) {
    *vv = vector_arg(1);

    int size = vector_capacity(*vv);
    int i;
    double max=0.0;
    double* px = vector_vec(*vv);
    for (i=0;i<size;i++) {
    	if (px[i]>max) max = px[i];
    }

    curRate = px[0];
    rmax = max;

    activeFlag = activeFlag + 1;
  }
  #endif
ENDVERBATIM
}

PROCEDURE update_time() {
VERBATIM
  IvocVect* vv; int i, i_prev, size; double* px;
  i = (int)index;
  i_prev = i;

  if (i >= 0) { // are we disabled?
    vv = *((IvocVect**)(&_p_vecTbins));
    if (vv) {
      size = vector_capacity(vv);
      px = vector_vec(vv);
      /* advance to current tbins without exceeding array bounds */
      while ((i+1 < size) && (t>=px[i+1])) {
	index += 1.;
	i += 1;
      }
      /* did the index change? */
      if (i!=i_prev) {
        /* advance curRate to next vecRate if possible */
        IvocVect *vvRate = *((IvocVect**)(&_p_vecRate));
        if (vvRate && vector_capacity(vvRate)>i) {
          px = vector_vec(vvRate);
          curRate = px[i];
        }
        else curRate = 1.0;
      }

      /* have we hit last bin? ... disable time advancing leaving curRate as it is*/
      if (i==size)
        index = -1.;

    } else { /* no vecTbins, use some defaults */
      rmax = 1.0;
      curRate = 1.0;
      index = -1.; /* no vecTbins ... disable time advancing & Poisson unit rate. */
    }
  }

ENDVERBATIM
}



COMMENT
/**
 * Upon a net_receive, we do up to two things.  The first is to determine the next time this artificial cell triggers
 * and sending a self event.  Second, we check to see if the synapse coupled to this artificial cell should be activated.
 * This second task is not done if we have just completed a state restore and only wish to restart the self event triggers.
 *
 * @param flag >= 0 for Typical activation, POST_RESTORE_RESTART_FLAG for only restarting the self event triggers 
 */
ENDCOMMENT
NET_RECEIVE (w) {
    : Note - if we have restored a sim from a saved state.  We need to restart the queue, but do not generate a spike now
    if ( flag == POST_RESTORE_RESTART_FLAG ) {
        if (t+event < start+duration) {
            net_send(event, activeFlag )
        }
    } else if( activeFlag == flag ) {
        update_time()
        generate_next_event()

        : stop even producing surrogate events if we are past duration
        if (t+event < start+duration) {
            net_send(event, activeFlag )
        }

        : check if we trigger event on coupled synapse
VERBATIM
        double u = (double)urand(_threadargs_);
        //printf("InhPoisson: spike time at time %g urand=%g curRate=%g, rmax=%g, curRate/rmax=%g \n",t, u, curRate, rmax, curRate/rmax);
        if (u<curRate/rmax) {
            debug_printf("\nInhPoisson: Spike time t = %g [urand=%g curRate=%g, rmax=%g]\n",
                         t, u, curRate, rmax);
ENDVERBATIM
            net_event(t)
VERBATIM
        }
ENDVERBATIM
    }
}


COMMENT
/**
 * Supply the POST_RESTORE_RESTART_FLAG.  For example, so a hoc program can call a NetCon.event with the proper event value
 *
 * @DEPRECATED Consider using restartEvent() for resuming the event loop
 * @return POST_RESTORE_RESTART_FLAG value for entities that wish to use its value
 */
ENDCOMMENT
FUNCTION getPostRestoreFlag() {
VERBATIM
    return POST_RESTORE_RESTART_FLAG;
ENDVERBATIM
}


COMMENT
/**
 * After a resume, send first event whose time is greater than the resume time.
 *
 * NOTE: Events generated right before the save time but scheduled for delivery afterwards
 *  will already be restored to the NetCon by the bbsavestate routines
 */
ENDCOMMENT
FUNCTION resumeEvent() {
    LOCAL elapsed_time
    : To be consistent with the previous run, it uses t=event as a starting point until it
    : reaches an elapsed_time >= resume_t.
    elapsed_time = event  : One event is always generated in the INITIAL block

    while( elapsed_time < t ) {
        update_time()
        generate_next_event()
        elapsed_time = elapsed_time + event
    }
    event = elapsed_time-t
    resumeEvent = elapsed_time
}

COMMENT
/**
 * Restart the event loop after a NEURON restore. It will discard events in the past
 */
ENDCOMMENT
PROCEDURE restartEvent() {
VERBATIM
#ifndef CORENEURON_BUILD
    double etime = resumeEvent(_threadargs_);
    if (etime < start+duration) {
        debug_printf("InhPoisson: First event after resume at t = %6.3f\n", etime);
        #if defined(NRN_VERSION_GTEQ_9_0_0)
        artcell_net_send(_tqitem, (double*)0, _ppvar[1].get<Point_process*>(), etime, activeFlag);
        #else
        artcell_net_send(_tqitem, (double*)0, (Point_process*)_ppvar[1]._pvoid, etime, activeFlag);
        #endif
    }
#endif
ENDVERBATIM
}


VERBATIM
static void bbcore_write(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        uint32_t dsize = 0;
        if (_p_vecRate)
        {
          dsize = (uint32_t)vector_capacity((IvocVect*)_p_vecRate);
        }
        if (iArray) {
                uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
                nrnran123_State** pv = (nrnran123_State**)(&_p_exp_rng);
                nrnran123_getids3(*pv, ia, ia+1, ia+2);

                // for stream sequence
                char which;

                nrnran123_getseq(*pv, ia+3, &which);
                ia[4] = (int)which;

                ia = ia + 5;
                pv = (nrnran123_State**)(&_p_uniform_rng);
                nrnran123_getids3( *pv, ia, ia+1, ia+2);

                nrnran123_getseq(*pv, ia+3, &which);
                ia[4] = (int)which;

                ia = ia + 5;
                IvocVect* vec = (IvocVect*)_p_vecRate;
                ia[0] = dsize;

                double *da = dArray + *doffset;
                double *dv;
                if(dsize)
                {
                  dv = vector_vec(vec);
                }
                int iInt;
                for (iInt = 0; iInt < dsize; ++iInt)
                {
                  da[iInt] = dv[iInt];
                }

                vec = (IvocVect*)_p_vecTbins;
                da = dArray + *doffset + dsize;
                if(dsize)
                {
                  dv = vector_vec(vec);
                }
                for (iInt = 0; iInt < dsize; ++iInt)
                {
                  da[iInt] = dv[iInt];
                }
        }
        *ioffset += 11;
        *doffset += 2*dsize;

}

static void bbcore_read(double* dArray, int* iArray, int* doffset, int* ioffset, _threadargsproto_) {
        uint32_t* ia = ((uint32_t*)iArray) + *ioffset;
        nrnran123_State** pv;
        if (ia[0] != 0 || ia[1] != 0)
        {
          pv = (nrnran123_State**)(&_p_exp_rng);
#if !NRNBBCORE
          if(*pv) {
              nrnran123_deletestream(*pv);
          }
#endif
          *pv = nrnran123_newstream3(ia[0], ia[1], ia[2] );
          nrnran123_setseq(*pv, ia[3], (char)ia[4]);
        }

        ia = ia + 5;
        if (ia[0] != 0 || ia[1] != 0)
        {
          pv = (nrnran123_State**)(&_p_uniform_rng);
#if !NRNBBCORE
          if(*pv) {
            nrnran123_deletestream(*pv);
          }
#endif
          *pv = nrnran123_newstream3(ia[0], ia[1], ia[2] );
          nrnran123_setseq(*pv, ia[2], (char)ia[3]);
        }

        ia = ia + 5;
        int dsize = ia[0];
        *ioffset += 11;

        double *da = dArray + *doffset;
        if(!_p_vecRate) {
          _p_vecRate = (double*)vector_new1(dsize);  /* works for dsize=0 */
        }
        assert(dsize == vector_capacity((IvocVect*)_p_vecRate));
        double *dv = vector_vec((IvocVect*)_p_vecRate);
        int iInt;
        for (iInt = 0; iInt < dsize; ++iInt)
        {
          dv[iInt] = da[iInt];
        }
        *doffset += dsize;

        da = dArray + *doffset;
        if(!_p_vecTbins) {
          _p_vecTbins = (double*)vector_new1(dsize);
        }
        assert(dsize == vector_capacity((IvocVect*)_p_vecTbins));
        dv = vector_vec((IvocVect*)_p_vecTbins);
        for (iInt = 0; iInt < dsize; ++iInt)
        {
          dv[iInt] = da[iInt];
        }
        *doffset += dsize;
}
ENDVERBATIM
