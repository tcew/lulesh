#include<kernels/kernelDefines.cu>

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
