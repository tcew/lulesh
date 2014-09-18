#include<kernels/kernelDefines.cu>

extern "C" __global__
#ifdef DOUBLE_PRECISION
__launch_bounds__(128,16)
#else
__launch_bounds__(128,16)
#endif
  void CalcMonotonicQRegionForElems_kernel(occaKernelInfoArg,
    Real_t qlc_monoq,
    Real_t qqc_monoq,
    Real_t monoq_limiter_mult,
    Real_t monoq_max_slope,
    Real_t ptiny,

    // the elementset length
    Index_t elength,

    Index_t *matElemlist,
    Index_t *elemBC,
    Index_t *lxim,
    Index_t *lxip,
    Index_t *letam,
    Index_t *letap,
    Index_t *lzetam,
    Index_t *lzetap,
    Real_t *delv_xi,
    Real_t *delv_eta,
    Real_t *delv_zeta,
    Real_t *delx_xi,
    Real_t *delx_eta,
    Real_t *delx_zeta,
    Real_t *vdov,Real_t *elemMass,Real_t *volo,Real_t *vnew,
    Real_t *qq, Real_t *ql,
    Real_t *q,
    Real_t qstop,
    Index_t* bad_q
    )
{
    int ielem=blockDim.x*blockIdx.x + threadIdx.x;

    if (ielem<elength) {
      Real_t qlin, qquad ;
      Real_t phixi, phieta, phizeta ;
      Index_t i = matElemlist[ielem];
      Int_t bcMask = elemBC[i] ;
      Real_t delvm, delvp ;

      /*  phixi     */
      Real_t norm = Real_t(1.) / ( delv_xi[i] + ptiny ) ;

      switch (bcMask & XI_M) {
         case 0:         delvm = delv_xi[lxim[i]] ; break ;
         case XI_M_SYMM: delvm = delv_xi[i] ;            break ;
         case XI_M_FREE: delvm = Real_t(0.0) ;                break ;
         default:        /* ERROR */ ;                        break ;
      }
      switch (bcMask & XI_P) {
         case 0:         delvp = delv_xi[lxip[i]] ; break ;
         case XI_P_SYMM: delvp = delv_xi[i] ;            break ;
         case XI_P_FREE: delvp = Real_t(0.0) ;                break ;
         default:        /* ERROR */ ;                        break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phixi = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm < phixi ) phixi = delvm ;
      if ( delvp < phixi ) phixi = delvp ;
      if ( phixi < Real_t(0.)) phixi = Real_t(0.) ;
      if ( phixi > monoq_max_slope) phixi = monoq_max_slope;


      /*  phieta     */
      norm = Real_t(1.) / ( delv_eta[i] + ptiny ) ;

      switch (bcMask & ETA_M) {
         case 0:          delvm = delv_eta[letam[i]] ; break ;
         case ETA_M_SYMM: delvm = delv_eta[i] ;             break ;
         case ETA_M_FREE: delvm = Real_t(0.0) ;                  break ;
         default:         /* ERROR */ ;                          break ;
      }
      switch (bcMask & ETA_P) {
         case 0:          delvp = delv_eta[letap[i]] ; break ;
         case ETA_P_SYMM: delvp = delv_eta[i] ;             break ;
         case ETA_P_FREE: delvp = Real_t(0.0) ;                  break ;
         default:         /* ERROR */ ;                          break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phieta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm  < phieta ) phieta = delvm ;
      if ( delvp  < phieta ) phieta = delvp ;
      if ( phieta < Real_t(0.)) phieta = Real_t(0.) ;
      if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

      /*  phizeta     */
      norm = Real_t(1.) / ( delv_zeta[i] + ptiny ) ;

      switch (bcMask & ZETA_M) {
         case 0:           delvm = delv_zeta[lzetam[i]] ; break ;
         case ZETA_M_SYMM: delvm = delv_zeta[i] ;              break ;
         case ZETA_M_FREE: delvm = Real_t(0.0) ;                    break ;
         default:          /* ERROR */ ;                            break ;
      }
      switch (bcMask & ZETA_P) {
         case 0:           delvp = delv_zeta[lzetap[i]] ; break ;
         case ZETA_P_SYMM: delvp = delv_zeta[i] ;              break ;
         case ZETA_P_FREE: delvp = Real_t(0.0) ;                    break ;
         default:          /* ERROR */ ;                            break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phizeta = Real_t(.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm   < phizeta ) phizeta = delvm ;
      if ( delvp   < phizeta ) phizeta = delvp ;
      if ( phizeta < Real_t(0.)) phizeta = Real_t(0.);
      if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

      /* Remove length scale */

      if ( vdov[i] > Real_t(0.) )  {
         qlin  = Real_t(0.) ;
         qquad = Real_t(0.) ;
      }
      else {
         Real_t delvxxi   = delv_xi[i]   * delx_xi[i]   ;
         Real_t delvxeta  = delv_eta[i]  * delx_eta[i]  ;
         Real_t delvxzeta = delv_zeta[i] * delx_zeta[i] ;

         if ( delvxxi   > Real_t(0.) ) delvxxi   = Real_t(0.) ;
         if ( delvxeta  > Real_t(0.) ) delvxeta  = Real_t(0.) ;
         if ( delvxzeta > Real_t(0.) ) delvxzeta = Real_t(0.) ;

         Real_t rho = elemMass[i] / (volo[i] * vnew[i]) ;

         qlin = -qlc_monoq * rho *
            (  delvxxi   * (Real_t(1.) - phixi) +
               delvxeta  * (Real_t(1.) - phieta) +
               delvxzeta * (Real_t(1.) - phizeta)  ) ;

         qquad = qqc_monoq * rho *
            (  delvxxi*delvxxi     * (Real_t(1.) - phixi*phixi) +
               delvxeta*delvxeta   * (Real_t(1.) - phieta*phieta) +
               delvxzeta*delvxzeta * (Real_t(1.) - phizeta*phizeta)  ) ;
      }

      qq[i] = qquad ;
      ql[i] = qlin  ;

      // Don't allow excessive artificial viscosity
      if (q[ielem] > qstop)
        *(bad_q) = ielem;

   }
}
