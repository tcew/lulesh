extern "C" __global__
void CalcAccelerationForNodes_kernel(occaKernelInfoArg, int numNode,
                                     Real_t *xdd, Real_t *ydd, Real_t *zdd,
                                     Real_t *fx, Real_t *fy, Real_t *fz,
                                     Real_t *nodalMass)
{
  int tid=blockDim.x*blockIdx.x+threadIdx.x;
  if (tid < numNode)
  {
      Real_t one_over_nMass = Real_t(1.)/nodalMass[tid];
      xdd[tid]=fx[tid]*one_over_nMass;
      ydd[tid]=fy[tid]*one_over_nMass;
      zdd[tid]=fz[tid]*one_over_nMass;
  }
}

/* occaKernel */
/* void CalcAccelerationForNodes_kernel(occaKernelInfoArg, */
/* 				     const int occaVariable numNode, */
/*                                      occaPointer Real_t *xdd, */
/* 				     occaPointer Real_t *ydd, */
/* 				     occaPointer Real_t *zdd, */
/*                                      const occaPointer Real_t *fx, */
/* 				     const occaPointer Real_t *fy, */
/* 				     const occaPointer Real_t *fz, */
/*                                      const occaPointer Real_t *nodalMass){ */

/*   occaGlobalFor0{ */
/*     // int tid=blockDim.x*blockIdx.x+threadIdx.x; */
/*     const int tid = occaGlobalId0; */

/*     if (tid < numNode){ */
/*       const Real_t one_over_nMass = Real_t(1.)/nodalMass[tid]; */
/*       xdd[tid]=fx[tid]*one_over_nMass; */
/*       ydd[tid]=fy[tid]*one_over_nMass; */
/*       zdd[tid]=fz[tid]*one_over_nMass; */
/*     } */
/*   } */
/* } */
