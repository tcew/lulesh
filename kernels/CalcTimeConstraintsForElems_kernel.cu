#include<kernels/kernelDefines.cu>

extern "C" __global__
#ifdef DOUBLE_PRECISION
__launch_bounds__(128,16)
#else
__launch_bounds__(128,16)
#endif
  void CalcTimeConstraintsForElems_kernel(occaKernelInfoArg,
    Index_t length,
    Real_t qqc2,
    Real_t dvovmax,
    Index_t *matElemlist,
    Real_t *ss,
    Real_t *vdov,
    Real_t *arealg,
    Real_t *dev_mindtcourant,
    Real_t *dev_mindthydro)
{
    int tid = threadIdx.x;
    int i=blockDim.x*blockIdx.x + tid;

    __shared__ volatile Real_t s_mindthydro[block_size];
    __shared__ volatile Real_t s_mindtcourant[block_size];

    Real_t mindthydro = Real_t(1.0e+20) ;
    Real_t mindtcourant = Real_t(1.0e+20) ;

    Real_t dthydro = mindthydro;
    Real_t dtcourant = mindtcourant;

    while (i<length) {

      Index_t indx = matElemlist[i] ;
      Real_t vdov_tmp = vdov[indx];

      // Computing dt_hydro
      if (vdov_tmp != Real_t(0.)) {
         Real_t dtdvov = dvovmax / (FABS(vdov_tmp)+Real_t(1.e-20)) ;
         if ( dthydro > dtdvov ) {
            dthydro = dtdvov ;
         }
      }
      if (dthydro < mindthydro)
        mindthydro = dthydro;

      // Computing dt_courant
      Real_t ss_tmp = ss[indx];
      Real_t area_tmp = arealg[indx];
      Real_t dtf = ss_tmp * ss_tmp ;

      dtf += ((vdov_tmp < 0.) ? qqc2*area_tmp*area_tmp*vdov_tmp*vdov_tmp : 0.);

      dtf = area_tmp / SQRT(dtf) ;

      /* determine minimum timestep with its corresponding elem */
      if (vdov_tmp != Real_t(0.) && dtf < dtcourant) {
        dtcourant = dtf ;
      }

      if (dtcourant< mindtcourant)
        mindtcourant= dtcourant;

      i += gridDim.x*blockDim.x;
    }

    s_mindthydro[tid]   = mindthydro;
    s_mindtcourant[tid] = mindtcourant;

    __syncthreads();

    // Do shared memory reduction
    if (block_size >= 1024) {
      if (tid < 512) {
        s_mindthydro[tid]   = min( s_mindthydro[tid]  , s_mindthydro[tid + 512]) ;
        s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid + 512]) ; }
      __syncthreads();  }

    if (block_size >=  512) {
      if (tid < 256) {
        s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid + 256]) ;
        s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid + 256]) ; }
      __syncthreads(); }

    if (block_size >=  256) {
      if (tid < 128) {
        s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid + 128]) ;
        s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid + 128]) ; }
      __syncthreads(); }

    if (block_size >=  128) {
      if (tid <  64) {
        s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +  64]) ;
        s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +  64]) ; }
      __syncthreads(); }

    if (tid <  32) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +  32]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +  32]) ;
    }

    if (tid <  16) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +  16]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +  16]) ;
    }
    if (tid <   8) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +   8]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +   8]) ;
    }
    if (tid <   4) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +   4]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +   4]) ;
    }
    if (tid <   2) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +   2]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +   2]) ;
    }
    if (tid <   1) {
      s_mindthydro[tid] = min( s_mindthydro[tid], s_mindthydro[tid +   1]) ;
      s_mindtcourant[tid] = min( s_mindtcourant[tid], s_mindtcourant[tid +   1]) ;
    }

    // Store in global memory
    if (tid==0) {
      dev_mindtcourant[blockIdx.x] = s_mindtcourant[0];
      dev_mindthydro[blockIdx.x] = s_mindthydro[0];
    }

}

/* occaKernel */
/* void CalcTimeConstraintsForElems_kernel(occaKernelInfoArg, */
/* 					const Index_t occaVariable length, */
/* 					const Real_t occaVariable qqc2, */
/* 					const Real_t occaVariable dvovmax, */
/* 					occaPointer Index_t *matElemlist, */
/* 					occaPointer Real_t *ss, */
/* 					occaPointer Real_t *vdov, */
/* 					occaPointer Real_t *arealg, */
/* 					occaPointer Real_t *dev_mindtcourant, */
/* 					occaPointer Real_t *dev_mindthydro){ */


/*   occaOuterFor0{ */

/*     // int tid = threadIdx.x; */
/*     // int i=blockDim.x*blockIdx.x + tid; */

/*     occaVolatile occaShared Real_t s_mindthydro[block_size] occaAligned; */
/*     occaVolatile occaShared Real_t s_mindtcourant[block_size] occaAligned; */

/*     Real_t mindthydro = Real_t(1.0e+20) ; */
/*     Real_t mindtcourant = Real_t(1.0e+20) ; */

/*     occaInnerFor0{ */
/*       Real_t dthydro = mindthydro; */
/*       Real_t dtcourant = mindtcourant; */

/*       int i = occaGlobalId0; */

/*       while (i<length) { */

/* 	Index_t indx = matElemlist[i] ; */
/* 	Real_t vdov_tmp = vdov[indx]; */

/* 	// Computing dt_hydro */
/* 	if (vdov_tmp != Real_t(0.)) { */
/* 	  Real_t dtdvov = dvovmax / (FABS(vdov_tmp)+Real_t(1.e-20)) ; */
/* 	  if ( dthydro > dtdvov ) { */
/*             dthydro = dtdvov ; */
/* 	  } */
/* 	} */
/* 	if (dthydro < mindthydro) */
/* 	  mindthydro = dthydro; */

/* 	// Computing dt_courant */
/* 	Real_t ss_tmp = ss[indx]; */
/* 	Real_t area_tmp = arealg[indx]; */
/* 	Real_t dtf = ss_tmp * ss_tmp ; */

/* 	dtf += ((vdov_tmp < 0.) ? qqc2*area_tmp*area_tmp*vdov_tmp*vdov_tmp : 0.); */

/* 	dtf = area_tmp / SQRT(dtf) ; */

/* 	/\* determine minimum timestep with its corresponding elem *\/ */
/* 	if (vdov_tmp != Real_t(0.) && dtf < dtcourant) { */
/* 	  dtcourant = dtf ; */
/* 	} */

/* 	if (dtcourant< mindtcourant) */
/* 	  mindtcourant= dtcourant; */

/* 	// i += gridDim.x*blockDim.x; */
/* 	i += occaGlobalDim0; */
/*       } */

/*       const int tid = occaInnerId0; */
/*       s_mindthydro[tid]   = mindthydro; */
/*       s_mindtcourant[tid] = mindtcourant; */
/*     } */


/*     occaBarrier(occaLocalMemFence); */

/* #if occaCPU==1 */
/*     Real_t mindtHydro   = s_mindthydro[0]; */
/*     Real_t mindtCourant = s_mindtcourant[0]; */

/*     for(int tid=1; tid<block_size; tid++){ */
/*       mindtHydro = MIN(mindtHydro, s_mindthydro[tid]); */
/*       mindtCourant = MIN(mindtCourant, s_mindtcourant[tid]); */
/*     } */

/*     dev_mindtcourant[occaOuterId0] = mindtCourant; */
/*     dev_mindthydro[occaOuterId0] = mindtHydro; */

/* #else */
/*     const int tid = occaInnerId0; */
/*     // Do shared memory reduction */
/*     if (block_size >= 1024) { */
/*       if (tid < 512) { */
/*         s_mindthydro[tid]   = MIN( s_mindthydro[tid]  , s_mindthydro[tid + 512]) ; */
/*         s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid + 512]) ; } */
/*       occaBarrier(occaLocalMemFence); */
/*     } */

/*     if (block_size >=  512) { */
/*       if (tid < 256) { */
/*         s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid + 256]) ; */
/*         s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid + 256]) ; } */
/*       occaBarrier(occaLocalMemFence); */
/*     } */

/*     if (block_size >=  256) { */
/*       if (tid < 128) { */
/*         s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid + 128]) ; */
/*         s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid + 128]) ; } */
/*       occaBarrier(occaLocalMemFence); */
/*       //__syncthreads(); */
/*     } */

/*     if (block_size >=  128) { */
/*       if (tid <  64) { */
/*         s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +  64]) ; */
/*         s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +  64]) ; } */
/*       occaBarrier(occaLocalMemFence); */
/*       //  __syncthreads(); */
/*     } */

/*     if (tid <  32) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +  32]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +  32]) ; */
/*       if(simdWidth < 32) occaBarrier(occaLocalMemFence); */
/*     } */

/*     if (tid <  16) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +  16]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +  16]) ; */
/*       if(simdWidth < 16) occaBarrier(occaLocalMemFence); */
/*     } */
/*     if (tid <   8) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +   8]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +   8]) ; */
/*       if(simdWidth < 8) occaBarrier(occaLocalMemFence); */
/*     } */
/*     if (tid <   4) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +   4]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +   4]) ; */
/*       if(simdWidth < 4) occaBarrier(occaLocalMemFence); */
/*     } */
/*     if (tid <   2) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +   2]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +   2]) ; */
/*       if(simdWidth < 2) occaBarrier(occaLocalMemFence); */
/*     } */
/*     if (tid <   1) { */
/*       s_mindthydro[tid] = MIN( s_mindthydro[tid], s_mindthydro[tid +   1]) ; */
/*       s_mindtcourant[tid] = MIN( s_mindtcourant[tid], s_mindtcourant[tid +   1]) ; */
/*     } */

/*     // Store in global memory */
/*     if (tid==0) { */
/*       dev_mindtcourant[occaOuterId0] = s_mindtcourant[0]; */
/*       dev_mindthydro[occaOuterId0] = s_mindthydro[0]; */
/*     } */
/* #endif */

/*   } */
/* } */
