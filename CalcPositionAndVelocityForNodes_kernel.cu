#include<kernelDefines.cu>

extern "C" __global__
void CalcPositionAndVelocityForNodes_kernel(occaKernelInfoArg, int numNode,
    const Real_t deltatime,
    const Real_t u_cut,
    Real_t* __restrict__ x,  Real_t* __restrict__ y,  Real_t* __restrict__ z,
    Real_t* __restrict__ xd, Real_t* __restrict__ yd, Real_t* __restrict__ zd,
    const Real_t* __restrict__ xdd, const Real_t* __restrict__ ydd, const Real_t* __restrict__ zdd)
{
    int i=blockDim.x*blockIdx.x+threadIdx.x;
    if (i < numNode)
    {
      Real_t xdtmp, ydtmp, zdtmp, dt;
      dt = deltatime;

      xdtmp = xd[i] + xdd[i] * dt ;
      ydtmp = yd[i] + ydd[i] * dt ;
      zdtmp = zd[i] + zdd[i] * dt ;

      if( FABS(xdtmp) < u_cut ) xdtmp = 0.0;
      if( FABS(ydtmp) < u_cut ) ydtmp = 0.0;
      if( FABS(zdtmp) < u_cut ) zdtmp = 0.0;

      x[i] += xdtmp * dt;
      y[i] += ydtmp * dt;
      z[i] += zdtmp * dt;

      xd[i] = xdtmp;
      yd[i] = ydtmp;
      zd[i] = zdtmp;
    }
}


/* #include<kernelDefines.occa> */

/* occaKernel void */
/* CalcPositionAndVelocityForNodes_kernel(occaKernelInfoArg, */
/* 				       const int occaVariable numNode, */
/* 				       const Real_t occaVariable deltatime, */
/* 				       const Real_t occaVariable u_cut, */
/* 				       occaPointer Real_t* occaRestrict x, */
/* 				       occaPointer Real_t* occaRestrict y, */
/* 				       occaPointer Real_t* occaRestrict z, */
/* 				       occaPointer Real_t* occaRestrict xd, */
/* 				       occaPointer Real_t* occaRestrict yd, */
/* 				       occaPointer Real_t* occaRestrict zd, */
/* 				       const occaPointer Real_t* occaRestrict xdd, */
/* 				       const occaPointer Real_t* occaRestrict ydd, */
/* 				       const occaPointer Real_t* occaRestrict zdd){ */

/*   occaGlobalFor0{ */
/*     // int i=blockDim.x*blockIdx.x+threadIdx.x; */
/*     const int i = occaGlobalId0; */
/*     if (i < numNode){ */

/*       Real_t xdtmp, ydtmp, zdtmp, dt; */
/*       dt = deltatime; */

/*       xdtmp = xd[i] + xdd[i] * dt ; */
/*       ydtmp = yd[i] + ydd[i] * dt ; */
/*       zdtmp = zd[i] + zdd[i] * dt ; */

/*       if( FABS(xdtmp) < u_cut ) xdtmp = 0.0; */
/*       if( FABS(ydtmp) < u_cut ) ydtmp = 0.0; */
/*       if( FABS(zdtmp) < u_cut ) zdtmp = 0.0; */

/*       x[i] += xdtmp * dt; */
/*       y[i] += ydtmp * dt; */
/*       z[i] += zdtmp * dt; */

/*       xd[i] = xdtmp; */
/*       yd[i] = ydtmp; */
/*       zd[i] = zdtmp; */
/*     } */
/*   } */
/* } */
