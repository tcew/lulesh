extern "C" __global__
void ApplyAccelerationBoundaryConditionsForNodes_kernel(occaKernelInfoArg,
    int numNodeBC, Real_t *xyzdd,
    Index_t *symm)
{
    int i=blockDim.x*blockIdx.x+threadIdx.x;
    if (i < numNodeBC)
    {
      if (i<numNodeBC) {
          xyzdd[symm[i]] = Real_t(0.0) ;
      }
    }
}

/* occaKernel void */
/* ApplyAccelerationBoundaryConditionsForNodes_kernel(occaKernelInfoArg, */
/* 						   const int occaVariable numNodeBC, */
/* 						   occaPointer Real_t *xyzdd, */
/* 						   const occaPointer Index_t *symm){ */

/*   occaGlobalFor0{ */
/*     // int i=blockDim.x*blockIdx.x+threadIdx.x; */
/*     const int i = occaGlobalId0; */
/*     if (i < numNodeBC){ */
/*       if (i<numNodeBC) { */
/* 	xyzdd[symm[i]] = Real_t(0.0) ; */
/*       } */
/*     } */
/*   } */
/* } */
