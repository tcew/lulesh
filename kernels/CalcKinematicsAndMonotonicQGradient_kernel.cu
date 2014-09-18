#include<kernels/kernelDefines.cu>

static __device__ __forceinline__
void CalcMonoGradient(Real_t *x, Real_t *y, Real_t *z,
                      Real_t *xv, Real_t *yv, Real_t *zv,
                      Real_t vol,
                      Real_t *delx_zeta,
                      Real_t *delv_zeta,
                      Real_t *delx_xi,
                      Real_t *delv_xi,
                      Real_t *delx_eta,
                      Real_t *delv_eta)
{

   #define SUM4(a,b,c,d) (a + b + c + d)
   const Real_t ptiny = Real_t(1.e-36) ;
   Real_t ax,ay,az ;
   Real_t dxv,dyv,dzv ;

   Real_t norm = Real_t(1.0) / ( vol + ptiny ) ;

   Real_t dxj = Real_t(-0.25)*(SUM4(x[0],x[1],x[5],x[4]) - SUM4(x[3],x[2],x[6],x[7])) ;
   Real_t dyj = Real_t(-0.25)*(SUM4(y[0],y[1],y[5],y[4]) - SUM4(y[3],y[2],y[6],y[7])) ;
   Real_t dzj = Real_t(-0.25)*(SUM4(z[0],z[1],z[5],z[4]) - SUM4(z[3],z[2],z[6],z[7])) ;

   Real_t dxi = Real_t( 0.25)*(SUM4(x[1],x[2],x[6],x[5]) - SUM4(x[0],x[3],x[7],x[4])) ;
   Real_t dyi = Real_t( 0.25)*(SUM4(y[1],y[2],y[6],y[5]) - SUM4(y[0],y[3],y[7],y[4])) ;
   Real_t dzi = Real_t( 0.25)*(SUM4(z[1],z[2],z[6],z[5]) - SUM4(z[0],z[3],z[7],z[4])) ;

   Real_t dxk = Real_t( 0.25)*(SUM4(x[4],x[5],x[6],x[7]) - SUM4(x[0],x[1],x[2],x[3])) ;
   Real_t dyk = Real_t( 0.25)*(SUM4(y[4],y[5],y[6],y[7]) - SUM4(y[0],y[1],y[2],y[3])) ;
   Real_t dzk = Real_t( 0.25)*(SUM4(z[4],z[5],z[6],z[7]) - SUM4(z[0],z[1],z[2],z[3])) ;

   /* find delvk and delxk ( i cross j ) */
   ax = dyi*dzj - dzi*dyj ;
   ay = dzi*dxj - dxi*dzj ;
   az = dxi*dyj - dyi*dxj ;

   *delx_zeta = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   ax *= norm ;
   ay *= norm ;
   az *= norm ;

   dxv = Real_t(0.25)*(SUM4(xv[4],xv[5],xv[6],xv[7]) - SUM4(xv[0],xv[1],xv[2],xv[3])) ;
   dyv = Real_t(0.25)*(SUM4(yv[4],yv[5],yv[6],yv[7]) - SUM4(yv[0],yv[1],yv[2],yv[3])) ;
   dzv = Real_t(0.25)*(SUM4(zv[4],zv[5],zv[6],zv[7]) - SUM4(zv[0],zv[1],zv[2],zv[3])) ;

   *delv_zeta = ax*dxv + ay*dyv + az*dzv ;

   /* find delxi and delvi ( j cross k ) */

   ax = dyj*dzk - dzj*dyk ;
   ay = dzj*dxk - dxj*dzk ;
   az = dxj*dyk - dyj*dxk ;

   *delx_xi = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   ax *= norm ;
   ay *= norm ;
   az *= norm ;

   dxv = Real_t(0.25)*(SUM4(xv[1],xv[2],xv[6],xv[5]) - SUM4(xv[0],xv[3],xv[7],xv[4])) ;
   dyv = Real_t(0.25)*(SUM4(yv[1],yv[2],yv[6],yv[5]) - SUM4(yv[0],yv[3],yv[7],yv[4])) ;
   dzv = Real_t(0.25)*(SUM4(zv[1],zv[2],zv[6],zv[5]) - SUM4(zv[0],zv[3],zv[7],zv[4])) ;

   *delv_xi = ax*dxv + ay*dyv + az*dzv ;

   /* find delxj and delvj ( k cross i ) */

   ax = dyk*dzi - dzk*dyi ;
   ay = dzk*dxi - dxk*dzi ;
   az = dxk*dyi - dyk*dxi ;

   *delx_eta = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

   ax *= norm ;
   ay *= norm ;
   az *= norm ;

   dxv = Real_t(-0.25)*(SUM4(xv[0],xv[1],xv[5],xv[4]) - SUM4(xv[3],xv[2],xv[6],xv[7])) ;
   dyv = Real_t(-0.25)*(SUM4(yv[0],yv[1],yv[5],yv[4]) - SUM4(yv[3],yv[2],yv[6],yv[7])) ;
   dzv = Real_t(-0.25)*(SUM4(zv[0],zv[1],zv[5],zv[4]) - SUM4(zv[3],zv[2],zv[6],zv[7])) ;

   *delv_eta = ax*dxv + ay*dyv + az*dzv ;
#undef SUM4
}


__device__
static
__forceinline__
void CalcElemVelocityGradient( const Real_t* const xvel,
                                const Real_t* const yvel,
                                const Real_t* const zvel,
                                const Real_t b[][8],
                                const Real_t detJ,
                                Real_t* const d )
{
  const Real_t inv_detJ = Real_t(1.0) / detJ ;
  Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
  const Real_t* const pfx = b[0];
  const Real_t* const pfy = b[1];
  const Real_t* const pfz = b[2];

  d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
                     + pfx[1] * (xvel[1]-xvel[7])
                     + pfx[2] * (xvel[2]-xvel[4])
                     + pfx[3] * (xvel[3]-xvel[5]) );

  d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
                     + pfy[1] * (yvel[1]-yvel[7])
                     + pfy[2] * (yvel[2]-yvel[4])
                     + pfy[3] * (yvel[3]-yvel[5]) );

  d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
                     + pfz[1] * (zvel[1]-zvel[7])
                     + pfz[2] * (zvel[2]-zvel[4])
                     + pfz[3] * (zvel[3]-zvel[5]) );

  dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
                      + pfx[1] * (yvel[1]-yvel[7])
                      + pfx[2] * (yvel[2]-yvel[4])
                      + pfx[3] * (yvel[3]-yvel[5]) );

  dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
                      + pfy[1] * (xvel[1]-xvel[7])
                      + pfy[2] * (xvel[2]-xvel[4])
                      + pfy[3] * (xvel[3]-xvel[5]) );

  dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
                      + pfx[1] * (zvel[1]-zvel[7])
                      + pfx[2] * (zvel[2]-zvel[4])
                      + pfx[3] * (zvel[3]-zvel[5]) );

  dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
                      + pfz[1] * (xvel[1]-xvel[7])
                      + pfz[2] * (xvel[2]-xvel[4])
                      + pfz[3] * (xvel[3]-xvel[5]) );

  dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
                      + pfy[1] * (zvel[1]-zvel[7])
                      + pfy[2] * (zvel[2]-zvel[4])
                      + pfy[3] * (zvel[3]-zvel[5]) );

  dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
                      + pfz[1] * (yvel[1]-yvel[7])
                      + pfz[2] * (yvel[2]-yvel[4])
                      + pfz[3] * (yvel[3]-yvel[5]) );
  d[5]  = Real_t( .5) * ( dxddy + dyddx );
  d[4]  = Real_t( .5) * ( dxddz + dzddx );
  d[3]  = Real_t( .5) * ( dzddy + dyddz );
}


__device__
static inline
Real_t AreaFace( const Real_t x0, const Real_t x1,
                 const Real_t x2, const Real_t x3,
                 const Real_t y0, const Real_t y1,
                 const Real_t y2, const Real_t y3,
                 const Real_t z0, const Real_t z1,
                 const Real_t z2, const Real_t z3)
{
   Real_t fx = (x2 - x0) - (x3 - x1);
   Real_t fy = (y2 - y0) - (y3 - y1);
   Real_t fz = (z2 - z0) - (z3 - z1);
   Real_t gx = (x2 - x0) + (x3 - x1);
   Real_t gy = (y2 - y0) + (y3 - y1);
   Real_t gz = (z2 - z0) + (z3 - z1);
   Real_t area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   return area ;
}


__device__
static inline
Real_t CalcElemCharacteristicLength( const Real_t x[8],
                                     const Real_t y[8],
                                     const Real_t z[8],
                                     const Real_t volume)
{
   Real_t a, charLength = Real_t(0.0);

   a = AreaFace(x[0],x[1],x[2],x[3],
                y[0],y[1],y[2],y[3],
                z[0],z[1],z[2],z[3]) ; // 38
   charLength = FMAX(a,charLength) ;

   a = AreaFace(x[4],x[5],x[6],x[7],
                y[4],y[5],y[6],y[7],
                z[4],z[5],z[6],z[7]) ;
   charLength = FMAX(a,charLength) ;

   a = AreaFace(x[0],x[1],x[5],x[4],
                y[0],y[1],y[5],y[4],
                z[0],z[1],z[5],z[4]) ;
   charLength = FMAX(a,charLength) ;

   a = AreaFace(x[1],x[2],x[6],x[5],
                y[1],y[2],y[6],y[5],
                z[1],z[2],z[6],z[5]) ;
   charLength = FMAX(a,charLength) ;

   a = AreaFace(x[2],x[3],x[7],x[6],
                y[2],y[3],y[7],y[6],
                z[2],z[3],z[7],z[6]) ;
   charLength = FMAX(a,charLength) ;

   a = AreaFace(x[3],x[0],x[4],x[7],
                y[3],y[0],y[4],y[7],
                z[3],z[0],z[4],z[7]) ;
   charLength = FMAX(a,charLength) ;

   charLength = Real_t(4.0) * volume / SQRT(charLength);

   return charLength;
}

extern "C" __global__
#ifdef DOUBLE_PRECISION
__launch_bounds__(64,8) // 64-bit
#else
__launch_bounds__(64,16) // 32-bit
#endif
  void CalcKinematicsAndMonotonicQGradient_kernel(occaKernelInfoArg,
    Index_t numElem, Index_t padded_numElem, const Real_t dt,
    const Index_t* __restrict__ nodelist, const Real_t* __restrict__ volo, const Real_t* __restrict__ v,
						  Real_t * __restrict__ x, //TextureObj<Real_t> x,
						  Real_t * __restrict__ y, //TextureObj<Real_t> y,
						  Real_t * __restrict__ z, //TextureObj<Real_t> z,
						  Real_t * __restrict__ xd,//TextureObj<Real_t> xd,
						  Real_t * __restrict__ yd,//TextureObj<Real_t> yd,
						  Real_t * __restrict__ zd,//TextureObj<Real_t> zd,
    Real_t* __restrict__ vnew,
    Real_t* __restrict__ delv,
    Real_t* __restrict__ arealg,
    Real_t* __restrict__ dxx,
    Real_t* __restrict__ dyy,
    Real_t* __restrict__ dzz,
    Real_t* __restrict__ vdov,
    Real_t* __restrict__ delx_zeta,
    Real_t* __restrict__ delv_zeta,
    Real_t* __restrict__ delx_xi,
    Real_t* __restrict__ delv_xi,
    Real_t* __restrict__ delx_eta,
    Real_t* __restrict__ delv_eta,
    Index_t* __restrict__ bad_vol,
    const Index_t num_threads
    )
{

  Real_t B[3][8] ; /** shape function derivatives */
  Index_t nodes[8] ;
  Real_t x_local[8] ;
  Real_t y_local[8] ;
  Real_t z_local[8] ;
  Real_t xd_local[8] ;
  Real_t yd_local[8] ;
  Real_t zd_local[8] ;
  Real_t D[6];

  int k=blockDim.x*blockIdx.x+threadIdx.x;

  if ( k < num_threads) {

    Real_t volume ;
    Real_t relativeVolume ;

    // get nodal coordinates from global arrays and copy into local arrays.
    #pragma unroll
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = nodelist[k+lnode*padded_numElem];
      nodes[lnode] = gnode;
      x_local[lnode] = x[gnode];
      y_local[lnode] = y[gnode];
      z_local[lnode] = z[gnode];
    }

    // volume calculations
    volume = CalcElemVolume(x_local, y_local, z_local );

    relativeVolume = volume / volo[k] ;
    vnew[k] = relativeVolume ;

    delv[k] = relativeVolume - v[k] ;
    // set characteristic length
    arealg[k] = CalcElemCharacteristicLength(x_local,y_local,z_local,volume);

    // get nodal velocities from global array and copy into local arrays.
    #pragma unroll
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = nodes[lnode];
      xd_local[lnode] = xd[gnode];
      yd_local[lnode] = yd[gnode];
      zd_local[lnode] = zd[gnode];
    }

    Real_t dt2 = Real_t(0.5) * dt;

    #pragma unroll
    for ( Index_t j=0 ; j<8 ; ++j )
    {
       x_local[j] -= dt2 * xd_local[j];
       y_local[j] -= dt2 * yd_local[j];
       z_local[j] -= dt2 * zd_local[j];
    }

    Real_t detJ;

    CalcElemShapeFunctionDerivatives(x_local,y_local,z_local,B,&detJ );

    CalcElemVelocityGradient(xd_local,yd_local,zd_local,B,detJ,D);

    // ------------------------
    // CALC LAGRANGE ELEM 2
    // ------------------------

    // calc strain rate and apply as constraint (only done in FB element)
    Real_t vdovNew = D[0] + D[1] + D[2];
    Real_t vdovthird = vdovNew/Real_t(3.0) ;

    // make the rate of deformation tensor deviatoric
    vdov[k] = vdovNew ;
    dxx[k] = D[0] - vdovthird ;
    dyy[k] = D[1] - vdovthird ;
    dzz[k] = D[2] - vdovthird ;

    // ------------------------
    // CALC MONOTONIC Q GRADIENT
    // ------------------------
    Real_t vol = volo[k]*vnew[k];

   // Undo x_local update
    #pragma unroll
    for ( Index_t j=0 ; j<8 ; ++j ) {
       x_local[j] += dt2 * xd_local[j];
       y_local[j] += dt2 * yd_local[j];
       z_local[j] += dt2 * zd_local[j];
    }

   CalcMonoGradient(x_local,y_local,z_local,xd_local,yd_local,zd_local,
                          vol,
                          &delx_zeta[k],&delv_zeta[k],&delx_xi[k],
                          &delv_xi[k], &delx_eta[k], &delv_eta[k]);

  //Check for bad volume
  if (relativeVolume < 0)
    *bad_vol = k;
  }
}
