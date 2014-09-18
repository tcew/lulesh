#include<kernels/kernelDefines.cu>

static
__device__
__forceinline__
void ApplyMaterialPropertiesForElems_device(
    Real_t& eosvmin, Real_t& eosvmax,
    Real_t* vnew, Real_t *v,
    Real_t& vnewc, Index_t* bad_vol, Index_t i, Index_t zn)
{
  vnewc = vnew[zn] ;

  if (eosvmin != Real_t(0.)) {
      if (vnewc < eosvmin)
          vnewc = eosvmin ;
  }

  if (eosvmax != Real_t(0.)) {
      if (vnewc > eosvmax)
          vnewc = eosvmax ;
  }

  // Now check for valid volume
  Real_t vc = v[zn];
  if (eosvmin != Real_t(0.)) {
     if (vc < eosvmin)
        vc = eosvmin ;
  }
  if (eosvmax != Real_t(0.)) {
     if (vc > eosvmax)
        vc = eosvmax ;
  }
  if (vc <= 0.) {
     *bad_vol = i;
  }

}


static
__device__ __forceinline__
void CalcSoundSpeedForElems_device(Real_t& vnewc, Real_t rho0, Real_t &enewc,
                            Real_t &pnewc, Real_t &pbvc,
                            Real_t &bvc, Real_t ss4o3, Index_t nz,
                            Real_t *ss, Index_t i, Index_t iz)
{
  Real_t ssTmp = (pbvc * enewc + vnewc * vnewc *
             bvc * pnewc) / rho0;
  if (ssTmp <= Real_t(1.111111e-36)) {
     ssTmp = Real_t(1.111111e-36);
  }
  else {
    ssTmp = SQRT(ssTmp) ;
  }
  ss[iz] = ssTmp;
}

static
__device__ __forceinline__
void CalcPressureForElems_device(
                      Real_t& p_new, Real_t& bvc,
                      Real_t& pbvc, Real_t& e_old,
                      Real_t& compression, Real_t& vnewc,
                      Real_t pmin,
                      Real_t p_cut, Real_t eosvmax)
{

      Real_t c1s = Real_t(2.0)/Real_t(3.0);
      Real_t p_temp = p_new;

      bvc = c1s * (compression + Real_t(1.));
      pbvc = c1s;

      p_temp = bvc * e_old ;

      if ( FABS(p_temp) <  p_cut )
        p_temp = Real_t(0.0) ;

      if ( vnewc >= eosvmax ) /* impossible condition here? */
        p_temp = Real_t(0.0) ;

      if (p_temp < pmin)
        p_temp = pmin ;

      p_new = p_temp;

}

static
__device__
__forceinline__
void UpdateVolumesForElems_device(Index_t numElem, Real_t& v_cut,
                                  Real_t *vnew,
                                  Real_t *v,
                                  int i)
{
   Real_t tmpV ;
   tmpV = vnew[i] ;

   if ( FABS(tmpV - Real_t(1.0)) < v_cut )
      tmpV = Real_t(1.0) ;
   v[i] = tmpV ;
}

static
__device__
__forceinline__
void CalcEnergyForElems_device(Real_t& p_new, Real_t& e_new, Real_t& q_new,
                            Real_t& bvc, Real_t& pbvc,
                            Real_t& p_old, Real_t& e_old, Real_t& q_old,
                            Real_t& compression, Real_t& compHalfStep,
                            Real_t& vnewc, Real_t& work, Real_t& delvc, Real_t pmin,
                            Real_t p_cut, Real_t e_cut, Real_t q_cut, Real_t emin,
                            Real_t& qq, Real_t& ql,
                            Real_t& rho0,
                            Real_t& eosvmax,
                            Index_t length)
{
   const Real_t sixth = Real_t(1.0) / Real_t(6.0) ;
   Real_t pHalfStep;

   e_new = e_old - Real_t(0.5) * delvc * (p_old + q_old)
      + Real_t(0.5) * work;

   if (e_new  < emin ) {
      e_new = emin ;
   }

   CalcPressureForElems_device(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                   pmin, p_cut, eosvmax);

   Real_t vhalf = Real_t(1.) / (Real_t(1.) + compHalfStep) ;

   if ( delvc > Real_t(0.) ) {
      q_new /* = qq = ql */ = Real_t(0.) ;
   }
   else {
      Real_t ssc = ( pbvc * e_new
              + vhalf * vhalf * bvc * pHalfStep ) / rho0 ;

      if ( ssc <= Real_t(0.) ) {
         ssc =Real_t(.333333e-36) ;
      } else {
         ssc = SQRT(ssc) ;
      }

      q_new = (ssc*ql + qq) ;
   }

   e_new = e_new + Real_t(0.5) * delvc
      * (  Real_t(3.0)*(p_old     + q_old)
           - Real_t(4.0)*(pHalfStep + q_new)) ;

   e_new += Real_t(0.5) * work;

   if (FABS(e_new) < e_cut) {
      e_new = Real_t(0.)  ;
   }
   if (     e_new  < emin ) {
      e_new = emin ;
   }

   CalcPressureForElems_device(p_new, bvc, pbvc, e_new, compression, vnewc,
                   pmin, p_cut, eosvmax);

   Real_t q_tilde ;

   if (delvc > Real_t(0.)) {
      q_tilde = Real_t(0.) ;
   }
   else {
      Real_t ssc = ( pbvc * e_new
              + vnewc * vnewc * bvc * p_new ) / rho0 ;

      if ( ssc <= Real_t(0.) ) {
         ssc = Real_t(.333333e-36) ;
      } else {
         ssc = SQRT(ssc) ;
      }

      q_tilde = (ssc*ql + qq) ;
   }

   e_new = e_new - (  Real_t(7.0)*(p_old     + q_old)
                            - Real_t(8.0)*(pHalfStep + q_new)
                            + (p_new + q_tilde)) * delvc*sixth ;

   if (FABS(e_new) < e_cut) {
      e_new = Real_t(0.)  ;
   }
   if ( e_new  < emin ) {
      e_new = emin ;
   }

   CalcPressureForElems_device(p_new, bvc, pbvc, e_new, compression, vnewc,
                   pmin, p_cut, eosvmax);


   if ( delvc <= Real_t(0.) ) {
      Real_t ssc = ( pbvc * e_new
              + vnewc * vnewc * bvc * p_new ) / rho0 ;

      if ( ssc <= Real_t(0.) ) {
         ssc = Real_t(.333333e-36) ;
      } else {
         ssc = SQRT(ssc) ;
      }

      q_new = (ssc*ql + qq) ;

      if (FABS(q_new) < q_cut) q_new = Real_t(0.) ;
   }

   return ;
}


extern "C" __global__
void ApplyMaterialPropertiesAndUpdateVolume_kernel(occaKernelInfoArg,
        Index_t length,
        Real_t rho0,
        Real_t e_cut,
        Real_t emin,
        Real_t* ql,
        Real_t* qq,
        Real_t* vnew,
        Real_t* v,
        Real_t pmin,
        Real_t p_cut,
        Real_t q_cut,
        Real_t eosvmin,
        Real_t eosvmax,
        Index_t* matElemlist,
        Real_t* e,
        Real_t* delv,
        Real_t* p,
        Real_t* q,
        Real_t ss4o3,
        Real_t* ss,
        Real_t v_cut,
        Index_t* bad_vol
)
{

  Real_t e_old, delvc, p_old, q_old;
  Real_t compression, compHalfStep;
  Real_t qq_old, ql_old, work;
  Real_t p_new, e_new, q_new;
  Real_t bvc, pbvc, vnewc;

  Index_t i=blockDim.x*blockIdx.x + threadIdx.x;

  if (i<length) {

    Index_t zidx  = matElemlist[i] ;

    ApplyMaterialPropertiesForElems_device
      (eosvmin,eosvmax,vnew,v,vnewc,bad_vol,i,zidx);

    e_old = e[zidx];
    delvc = delv[zidx];
    p_old = p[zidx];
    q_old = q[zidx];

    Real_t vchalf ;
    compression = Real_t(1.) / vnewc - Real_t(1.);
    vchalf = vnewc - delvc * Real_t(.5);
    compHalfStep = Real_t(1.) / vchalf - Real_t(1.);

    if ( eosvmin != Real_t(0.) ) {
        if (vnewc <= eosvmin) { /* impossible due to calling func? */
            compHalfStep = compression ;
        }
    }
    if ( eosvmax != Real_t(0.) ) {
        if (vnewc >= eosvmax) { /* impossible due to calling func? */
            p_old        = Real_t(0.) ;
            compression  = Real_t(0.) ;
            compHalfStep = Real_t(0.) ;
        }
    }

    qq_old = qq[zidx] ;
    ql_old = ql[zidx] ;
    work = Real_t(0.) ;

    CalcEnergyForElems_device(p_new, e_new, q_new, bvc, pbvc,
                 p_old, e_old,  q_old, compression, compHalfStep,
                 vnewc, work,  delvc, pmin,
                 p_cut, e_cut, q_cut, emin,
                 qq_old, ql_old, rho0, eosvmax, length);

    p[zidx] = p_new ;
    e[zidx] = e_new ;
    q[zidx] = q_new ;

    CalcSoundSpeedForElems_device
       (vnewc,rho0,e_new,p_new,pbvc,bvc,ss4o3,length,ss,i,zidx);

    UpdateVolumesForElems_device(length,v_cut,vnew,v,i);

  }
}
