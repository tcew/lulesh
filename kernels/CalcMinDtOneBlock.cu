#include<kernelDefines.cu>

extern "C" __global__
void CalcMinDtOneBlock(occaKernelInfoArg,
		       Real_t* dev_mindthydro,
		       Real_t* dev_mindtcourant,
		       Real_t* dtcourant,
		       Real_t* dthydro,
		       Index_t shared_array_size)
{

  volatile __shared__ Real_t s_data[block_size];
  int tid = threadIdx.x;

  if (blockIdx.x==0)
  {
    if (tid < shared_array_size)
      s_data[tid] = dev_mindtcourant[tid];
    else
      s_data[tid] = 1.0e20;

    __syncthreads();

    if (block_size >= 1024) { if (tid < 512) { s_data[tid] = min(s_data[tid],s_data[tid + 512]); } __syncthreads(); }
    if (block_size >=  512) { if (tid < 256) { s_data[tid] = min(s_data[tid],s_data[tid + 256]); } __syncthreads(); }
    if (block_size >=  256) { if (tid < 128) { s_data[tid] = min(s_data[tid],s_data[tid + 128]); } __syncthreads(); }
    if (block_size >=  128) { if (tid <  64) { s_data[tid] = min(s_data[tid],s_data[tid +  64]); } __syncthreads(); }
    if (tid <  32) { s_data[tid] = min(s_data[tid],s_data[tid +  32]); }
    if (tid <  16) { s_data[tid] = min(s_data[tid],s_data[tid +  16]); }
    if (tid <   8) { s_data[tid] = min(s_data[tid],s_data[tid +   8]); }
    if (tid <   4) { s_data[tid] = min(s_data[tid],s_data[tid +   4]); }
    if (tid <   2) { s_data[tid] = min(s_data[tid],s_data[tid +   2]); }
    if (tid <   1) { s_data[tid] = min(s_data[tid],s_data[tid +   1]); }

    if (tid<1)
    {
      *(dtcourant)= s_data[0];
    }
  }
  else if (blockIdx.x==1)
  {
    if (tid < shared_array_size)
      s_data[tid] = dev_mindthydro[tid];
    else
      s_data[tid] = 1.0e20;

    __syncthreads();

    if (block_size >= 1024) { if (tid < 512) { s_data[tid] = min(s_data[tid],s_data[tid + 512]); } __syncthreads(); }
    if (block_size >=  512) { if (tid < 256) { s_data[tid] = min(s_data[tid],s_data[tid + 256]); } __syncthreads(); }
    if (block_size >=  256) { if (tid < 128) { s_data[tid] = min(s_data[tid],s_data[tid + 128]); } __syncthreads(); }
    if (block_size >=  128) { if (tid <  64) { s_data[tid] = min(s_data[tid],s_data[tid +  64]); } __syncthreads(); }
    if (tid <  32) { s_data[tid] = min(s_data[tid],s_data[tid +  32]); }
    if (tid <  16) { s_data[tid] = min(s_data[tid],s_data[tid +  16]); }
    if (tid <   8) { s_data[tid] = min(s_data[tid],s_data[tid +   8]); }
    if (tid <   4) { s_data[tid] = min(s_data[tid],s_data[tid +   4]); }
    if (tid <   2) { s_data[tid] = min(s_data[tid],s_data[tid +   2]); }
    if (tid <   1) { s_data[tid] = min(s_data[tid],s_data[tid +   1]); }

    if (tid<1)
    {
      *(dthydro) = s_data[0];
    }
  }
}


/* #include<kernelDefines.occa> */

/* //template <int block_size> */
/* occaKernel void CalcMinDtOneBlock(occaKernelInfoArg, */
/* 				  occaPointer Real_t* dev_mindthydro, */
/* 				  occaPointer Real_t* dev_mindtcourant, */
/* 				  occaPointer Real_t* dtcourant, */
/* 				  occaPointer Real_t* dthydro, */
/* 				  const Index_t occaVariable shared_array_size){ */

/*   occaOuterFor0{ */

/*     occaVolatile occaShared Real_t s_data[block_size] occaAligned; */

/*     //    int tid = threadIdx.x; */

/*     if (occaOuterId0 == 0){ */

/*       occaInnerFor0{ */
/* 	const int tid = occaInnerId0; */

/* 	if (tid < shared_array_size) */
/* 	  s_data[tid] = dev_mindtcourant[tid]; */
/* 	else */
/* 	  s_data[tid] = 1.0e20; */
/*       } */

/*       occaBarrier(occaLocalMemFence); */

/* #if occaCPU==1 */
/*       (*dtcourant) = s_data[0]; */
/*       for(int tid=1; tid<block_size; tid++) */
/* 	(*dtcourant) = MIN((*dtcourant), s_data[tid]); */
/* #else */
/*       const int tid = occaInnerId0; */

/*       if (block_size >= 1024) { */
/* 	if (tid < 512) { s_data[tid] = MIN(s_data[tid],s_data[tid + 512]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  512) { */
/* 	if (tid < 256) { s_data[tid] = MIN(s_data[tid],s_data[tid + 256]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  256) { */
/* 	if (tid < 128) { s_data[tid] = MIN(s_data[tid],s_data[tid + 128]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  128) { */
/* 	if (tid <  64) { s_data[tid] = MIN(s_data[tid],s_data[tid +  64]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <  32) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +  32]); */
/* 	if(simdWidth < 32)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <  16) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +  16]); */
/* 	if(simdWidth < 32)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   8) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   8]); */
/* 	if(simdWidth < 8)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   4) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   4]); */
/* 	if(simdWidth < 4)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   2) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   2]); */
/* 	if(simdWidth < 2)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   1) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   1]); */
/* 	if(simdWidth < 1)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid<1){ */
/* 	*(dtcourant)= s_data[0]; */
/*       } */
/* #endif */
/*     } */
/*     else if (occaOuterId0==1){ */

/*       occaInnerFor0{ */
/* 	const int tid = occaInnerId0; */

/* 	if (tid < shared_array_size) */
/* 	  s_data[tid] = dev_mindthydro[tid]; */
/* 	else */
/* 	  s_data[tid] = 1.0e20; */
/*       } */

/*       occaBarrier(occaLocalMemFence); */
/* #if occaCPU==1 */
/*       (*dthydro) = s_data[0]; */
/*       for(int tid=1; tid<block_size; tid++) */
/* 	(*dthydro) = MIN((*dthydro), s_data[tid]); */

/* #else */
/*       const int tid = occaInnerId0; */
/*       if (block_size >= 1024) { */
/* 	if (tid < 512) { s_data[tid] = MIN(s_data[tid],s_data[tid + 512]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  512) { */
/* 	if (tid < 256) { s_data[tid] = MIN(s_data[tid],s_data[tid + 256]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  256) { */
/* 	if (tid < 128) { s_data[tid] = MIN(s_data[tid],s_data[tid + 128]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (block_size >=  128) { */
/* 	if (tid <  64) { s_data[tid] = MIN(s_data[tid],s_data[tid +  64]); } */
/* 	occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <  32) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +  32]); */
/* 	if(simdWidth < 32)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <  16) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +  16]); */
/* 	if(simdWidth < 32)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   8) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   8]); */
/* 	if(simdWidth < 8)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   4) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   4]); */
/* 	if(simdWidth < 4)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   2) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   2]); */
/* 	if(simdWidth < 2)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid <   1) { */
/* 	s_data[tid] = MIN(s_data[tid],s_data[tid +   1]); */
/* 	if(simdWidth < 1)  occaBarrier(occaLocalMemFence); */
/*       } */

/*       if (tid<1){ */
/* 	*(dthydro) = s_data[0]; */
/*       } */

/* #endif */
/*     } */
/*   } */
/* } */
