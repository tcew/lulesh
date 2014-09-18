#include <occa.hpp>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <assert.h>

#define LULESH_SHOW_PROGRESS 1
//#define DOUBLE_PRECISION
//#define SAMI

const int max_dimGrid = 256;

enum {
  VolumeError = -1,
  QStopError = -2,
  LFileError = -3
} ;

/****************************************************/
/* Allow flexibility for arithmetic representations */
/****************************************************/

/* Could also support fixed point and interval arithmetic types */
typedef float        real4 ;
typedef double       real8 ;

typedef int    Index_t ; /* array subscript and loop index */
typedef int    Int_t ;   /* integer representation */
#ifdef DOUBLE_PRECISION
typedef real8  Real_t ;  /* floating point representation */
#else
typedef real4  Real_t ;  /* floating point representation */
#endif

#define MAX(a, b) ( ((a) > (b)) ? (a) : (b))

/* Stuff needed for boundary conditions */
/* 2 BCs on each of 6 hexahedral faces (12 bits) */
#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002

#define XI_P        0x00038
#define XI_P_SYMM   0x00008
#define XI_P_FREE   0x00010

#define ETA_M       0x001c0
#define ETA_M_SYMM  0x00040
#define ETA_M_FREE  0x00080

#define ETA_P       0x00e00
#define ETA_P_SYMM  0x00200
#define ETA_P_FREE  0x00400

#define ZETA_M      0x07000
#define ZETA_M_SYMM 0x01000
#define ZETA_M_FREE 0x02000

#define ZETA_P      0x38000
#define ZETA_P_SYMM 0x08000
#define ZETA_P_FREE 0x10000

/* Given a number of bytes, nbytes, and a byte alignment, align, (e.g., 2,
 * 4, 8, or 16), return the smallest integer that is larger than nbytes and
 * a multiple of align.
 */
#define PAD_DIV(nbytes, align)  (((nbytes) + (align) - 1) / (align))
#define PAD(nbytes, align)  (PAD_DIV((nbytes),(align)) * (align))

template <typename T>
void occaCheck(occa::memory &a){

#ifndef NDEBUG
  const int n = a.bytes()/sizeof(T);

  if(n){
    std::vector<T> testA(n);

    a.copyTo(&(testA[0]));

    for(int i=0; i<n; i++){

      if(testA[i] != testA[i])
	std::cout<<"a["<<i<<"] = "<< testA[i] << std::endl;

      assert(testA[i] == testA[i]);
    }
  }
#endif
}


template <typename T>
void occaCompare(occa::memory &o_a, occa::memory &o_b){

#ifndef NDEBUG
  int len = o_a.bytes()/sizeof(T);

  assert(len == o_b.bytes()/sizeof(T));

  std::vector<T> a(len);
  std::vector<T> b(len);

  o_a.copyTo(&(a[0]));
  o_b.copyTo(&(b[0]));


  for(int i=0; i<len; i++)
    assert(a[i] == b[i]);
#endif

}

Real_t CalcElemVolumeTemp( const Real_t x0, const Real_t x1,
			   const Real_t x2, const Real_t x3,
			   const Real_t x4, const Real_t x5,
			   const Real_t x6, const Real_t x7,
			   const Real_t y0, const Real_t y1,
			   const Real_t y2, const Real_t y3,
			   const Real_t y4, const Real_t y5,
			   const Real_t y6, const Real_t y7,
			   const Real_t z0, const Real_t z1,
			   const Real_t z2, const Real_t z3,
			   const Real_t z4, const Real_t z5,
			   const Real_t z6, const Real_t z7){

  Real_t twelveth = ((Real_t) 1.0)/((Real_t) 12.0);

  Real_t dx61 = x6 - x1;
  Real_t dy61 = y6 - y1;
  Real_t dz61 = z6 - z1;

  Real_t dx70 = x7 - x0;
  Real_t dy70 = y7 - y0;
  Real_t dz70 = z7 - z0;

  Real_t dx63 = x6 - x3;
  Real_t dy63 = y6 - y3;
  Real_t dz63 = z6 - z3;

  Real_t dx20 = x2 - x0;
  Real_t dy20 = y2 - y0;
  Real_t dz20 = z2 - z0;

  Real_t dx50 = x5 - x0;
  Real_t dy50 = y5 - y0;
  Real_t dz50 = z5 - z0;

  Real_t dx64 = x6 - x4;
  Real_t dy64 = y6 - y4;
  Real_t dz64 = z6 - z4;

  Real_t dx31 = x3 - x1;
  Real_t dy31 = y3 - y1;
  Real_t dz31 = z3 - z1;

  Real_t dx72 = x7 - x2;
  Real_t dy72 = y7 - y2;
  Real_t dz72 = z7 - z2;

  Real_t dx43 = x4 - x3;
  Real_t dy43 = y4 - y3;
  Real_t dz43 = z4 - z3;

  Real_t dx57 = x5 - x7;
  Real_t dy57 = y5 - y7;
  Real_t dz57 = z5 - z7;

  Real_t dx14 = x1 - x4;
  Real_t dy14 = y1 - y4;
  Real_t dz14 = z1 - z4;

  Real_t dx25 = x2 - x5;
  Real_t dy25 = y2 - y5;
  Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3)		\
  ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  // 11 + 3*14
  Real_t volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
		   dy31 + dy72, dy63, dy20,
		   dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
		   dy43 + dy57, dy64, dy70,
		   dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
		   dy14 + dy25, dy61, dy50,
		   dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

  volume *= twelveth;

  return volume;
}


Real_t CalcElemVolume( const Real_t x[8],
		       const Real_t y[8],
		       const Real_t z[8]){

  Real_t volume =
    CalcElemVolumeTemp( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
			y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
			z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);

  return volume;

}




#if defined(_WIN64) || defined(__LP64__)
// 64-bit pointer operand constraint for inlined asm
#define _ASM_PTR_ "l"
#else
// 32-bit pointer operand constraint for inlined asm
#define _ASM_PTR_ "r"
#endif


occa::device occaHandle;
occa::kernel CalcVolumeForceForElems_kernel;
occa::kernel AddNodeForcesFromElems_kernel;
occa::kernel CalcAccelerationForNodes_kernel;
occa::kernel ApplyAccelerationBoundaryConditionsForNodes_kernel;
occa::kernel CalcPositionAndVelocityForNodes_kernel;
occa::kernel CalcKinematicsAndMonotonicQGradient_kernel;
occa::kernel CalcMonotonicQRegionForElems_kernel;
occa::kernel ApplyMaterialPropertiesAndUpdateVolume_kernel;
occa::kernel CalcTimeConstraintsForElems_kernel;
occa::kernel CalcMinDtOneBlock;
occa::kernel FillKernel;
occa::kernel ZeroKernel;

class Domain
{

public:

  // Index_t max_streams;
  // std::vector<cudaStream_t> streams;

  /* Elem-centered */

  // Vector_d<Index_t> matElemlist ; /* material indexset */
  // Vector_d<Index_t> nodelist ;    /* elemToNode connectivity */

  occa::memory matElemlist, nodelist;

  // Vector_d<Index_t> lxim ;        /* element connectivity through face */
  // Vector_d<Index_t> lxip ;
  // Vector_d<Index_t> letam ;
  // Vector_d<Index_t> letap ;
  // Vector_d<Index_t> lzetam ;
  // Vector_d<Index_t> lzetap ;

  occa::memory lxim, lxip, letam, letap, lzetam, lzetap;

  //  Vector_d<Int_t> elemBC ;        /* elem face symm/free-surf flag */

  occa::memory elemBC;

  // Vector_d<Real_t> e ;            /* energy */

  // Vector_d<Real_t> p ;            /* pressure */

  // Vector_d<Real_t> q ;            /* q */
  // Vector_d<Real_t> ql ;           /* linear term for q */
  // Vector_d<Real_t> qq ;           /* quadratic term for q */

  // Vector_d<Real_t> v ;            /* relative volume */

  // Vector_d<Real_t> volo ;         /* reference volume */
  // Vector_d<Real_t> delv ;         /* m_vnew - m_v */
  // Vector_d<Real_t> vdov ;         /* volume derivative over volume */

  // Vector_d<Real_t> arealg ;       /* char length of an element */

  // Vector_d<Real_t> ss ;           /* "sound speed" */

  // Vector_d<Real_t> elemMass ;     /* mass */

  occa::memory e, p, q, ql, qq, v, volo, delv, vdov, arealg, ss, elemMass;


  // Vector_d<Real_t>* vnew ;         /* new relative volume -- temporary */

  // Vector_d<Real_t>* delv_xi ;      /* velocity gradient -- temporary */
  // Vector_d<Real_t>* delv_eta ;
  // Vector_d<Real_t>* delv_zeta ;

  // Vector_d<Real_t>* delx_xi ;      /* coordinate gradient -- temporary */
  // Vector_d<Real_t>* delx_eta ;
  // Vector_d<Real_t>* delx_zeta ;

  // Vector_d<Real_t>* dxx ;          /* principal strains -- temporary */
  // Vector_d<Real_t>* dyy ;
  // Vector_d<Real_t>* dzz ;

  occa::memory vnew, delv_xi, delv_eta, delv_zeta, delx_xi, delx_eta, delx_zeta, dxx, dyy, dzz;

  /* Node-centered */

  // Vector_d<Real_t> x ;            /* coordinates */
  // Vector_d<Real_t> y ;
  // Vector_d<Real_t> z ;

  occa::memory x, y, z;

  // TextureObj<Real_t> tex_x;
  // TextureObj<Real_t> tex_y;
  // TextureObj<Real_t> tex_z;

  // TODO: use textures instead of global
  //  occa::textMem tex_x, tex_y, tex_z;
  occa::memory tex_x, tex_y, tex_z;


  // Vector_d<Real_t> xd ;           /* velocities */
  // Vector_d<Real_t> yd ;
  // Vector_d<Real_t> zd ;
  occa::memory xd, yd, zd;

  // TextureObj<Real_t> tex_xd;
  // TextureObj<Real_t> tex_yd;
  // TextureObj<Real_t> tex_zd;

  // TODO: use textures instead of global
  // occa::textMem tex_xd, tex_yd, tex_zd;
  occa::memory tex_xd, tex_yd, tex_zd;


  // Vector_d<Real_t> xdd ;          /* accelerations */
  // Vector_d<Real_t> ydd ;
  // Vector_d<Real_t> zdd ;
  occa::memory xdd, ydd, zdd;

  // Vector_d<Real_t> fx ;           /* forces */
  // Vector_d<Real_t> fy ;
  // Vector_d<Real_t> fz ;
  occa::memory fx, fy, fz;

  occa::memory fx_elem, fy_elem, fz_elem;

  // Vector_d<Real_t> nodalMass ;    /* mass */
  occa::memory nodalMass;

  /* Boundary nodesets */

  // Vector_d<Index_t> symmX ;       /* symmetry plane nodesets */
  // Vector_d<Index_t> symmY ;
  // Vector_d<Index_t> symmZ ;
  occa::memory symmX, symmY, symmZ;

  // Vector_d<Int_t> nodeElemCount ;
  // Vector_d<Int_t> nodeElemStart;
  // Vector_d<Index_t> nodeElemCornerList ;
  occa::memory nodeElemCount, nodeElemStart, nodeElemCornerList;

  /* Parameters */

  Real_t dtfixed ;               /* fixed time increment */
  Real_t deltatimemultlb ;
  Real_t deltatimemultub ;
  Real_t stoptime ;              /* end time for simulation */
  Real_t dtmax ;                 /* maximum allowable time increment */
  Int_t cycle ;                  /* iteration count for simulation */

  Real_t* dthydro_h;             /* hydro time constraint */
  Real_t* dtcourant_h;           /* courant time constraint */
  Index_t* bad_q_h;              /* flag to indicate Q error */
  Index_t* bad_vol_h;            /* flag to indicate volume error */

  /* cuda Events to indicate completion of certain kernels */
  //  cudaEvent_t time_constraint_computed;

  Real_t time_h ;               /* current time */
  Real_t deltatime_h ;          /* variable time increment */

  Real_t u_cut ;                /* velocity tolerance */
  Real_t hgcoef ;               /* hourglass control */
  Real_t qstop ;                /* excessive q indicator */
  Real_t monoq_max_slope ;
  Real_t monoq_limiter_mult ;
  Real_t e_cut ;                /* energy tolerance */
  Real_t p_cut ;                /* pressure tolerance */
  Real_t ss4o3 ;
  Real_t q_cut ;                /* q tolerance */
  Real_t v_cut ;                /* relative volume tolerance */
  Real_t qlc_monoq ;            /* linear term coef for q */
  Real_t qqc_monoq ;            /* quadratic term coef for q */
  Real_t qqc ;
  Real_t eosvmax ;
  Real_t eosvmin ;
  Real_t pmin ;                 /* pressure floor */
  Real_t emin ;                 /* energy floor */
  Real_t dvovmax ;              /* maximum allowable volume change */
  Real_t refdens ;              /* reference density */

  Index_t sizeX ;
  Index_t sizeY ;
  Index_t sizeZ ;
  Index_t maxPlaneSize ;

  Index_t numElem ;
  Index_t padded_numElem ;

  Index_t numNode;
  Index_t padded_numNode ;

  Index_t numSymmX ;
  Index_t numSymmY ;
  Index_t numSymmZ ;

  Index_t octantCorner;

} ;


void occaCheckDomain(Domain *domain){

#ifndef NDEBUG
  occaCheck<Index_t>(domain->matElemlist);
  occaCheck<Index_t>(domain->nodelist);
  occaCheck<Index_t>(domain->lxim);
  occaCheck<Index_t>(domain->lxip);
  occaCheck<Index_t>(domain->letam);
  occaCheck<Index_t>(domain->letap);
  occaCheck<Index_t>(domain->lzetam);
  occaCheck<Index_t>(domain->lzetap);

  occaCheck<Int_t>(domain->elemBC);

  occaCheck<Real_t>(domain->e);
  occaCheck<Real_t>(domain->p);
  occaCheck<Real_t>(domain->q);
  occaCheck<Real_t>(domain->ql);
  occaCheck<Real_t>(domain->qq);
  occaCheck<Real_t>(domain->v);
  occaCheck<Real_t>(domain->volo);
  occaCheck<Real_t>(domain->delv);
  occaCheck<Real_t>(domain->vdov);
  occaCheck<Real_t>(domain->arealg);
  occaCheck<Real_t>(domain->ss);
  occaCheck<Real_t>(domain->elemMass);

  occaCheck<Real_t>(domain->vnew);
  occaCheck<Real_t>(domain->delv_xi);
  occaCheck<Real_t>(domain->delv_eta);
  occaCheck<Real_t>(domain->delv_zeta);
  occaCheck<Real_t>(domain->delx_xi);
  occaCheck<Real_t>(domain->delx_eta);
  occaCheck<Real_t>(domain->delx_zeta);
  occaCheck<Real_t>(domain->dxx);
  occaCheck<Real_t>(domain->dyy);

  occaCheck<Real_t>(domain->x);
  occaCheck<Real_t>(domain->y);
  occaCheck<Real_t>(domain->z);


  occaCheck<Real_t>(domain->tex_x);
  occaCheck<Real_t>(domain->tex_y);
  occaCheck<Real_t>(domain->tex_z);

  occaCheck<Real_t>(domain->xd);
  occaCheck<Real_t>(domain->yd);
  occaCheck<Real_t>(domain->zd);


  occaCheck<Real_t>(domain->tex_xd);
  occaCheck<Real_t>(domain->tex_yd);
  occaCheck<Real_t>(domain->tex_zd);

  occaCheck<Real_t>(domain->xdd);
  occaCheck<Real_t>(domain->ydd);
  occaCheck<Real_t>(domain->zdd);

  occaCheck<Real_t>(domain->fx);
  occaCheck<Real_t>(domain->fy);
  occaCheck<Real_t>(domain->fz);

  occaCheck<Real_t>(domain->nodalMass);

  occaCheck<Index_t>(domain->symmX);
  occaCheck<Index_t>(domain->symmY);
  occaCheck<Index_t>(domain->symmZ);

  occaCheck<Int_t>(domain->nodeElemCount);
  occaCheck<Int_t>(domain->nodeElemStart);
  occaCheck<Index_t>(domain->nodeElemCornerList);
#endif
}


// void cuda_init()
// {
//     Int_t deviceCount, dev;
//     cudaDeviceProp cuda_deviceProp;

//     cudaSafeCall( cudaGetDeviceCount(&deviceCount) );
//     if (deviceCount == 0) {
//         fprintf(stderr, "cuda_init(): no devices supporting CUDA.\n");
//         exit(1);
//     }

//     dev = 0;

//     if ((dev < 0) || (dev > deviceCount-1)) {
//         fprintf(stderr, "cuda_init(): requested device (%d) out of range [%d,%d]\n",
//                 dev, 0, deviceCount-1);
//         exit(1);
//     }

//     printf("Setting CUDA device %d\n",dev);
//     cudaSafeCall( cudaSetDevice(dev) );

//     cudaSafeCall( cudaGetDeviceProperties(&cuda_deviceProp, dev) );
//     if (cuda_deviceProp.major < 3) {
//         fprintf(stderr, "cuda_init(): This implementation of Lulesh requires device SM 3.0+.\n", dev);
//         exit(1);
//     }

// #if CUDART_VERSION < 5000
//    fprintf(stderr,"cuda_init(): This implementation of Lulesh uses texture objects, which is requires Cuda 5.0+.\n");
//    exit(1);
// #endif

// }

#if 0
static void buildLuleshKernels(){

  occa::kernelInfo defs;

#ifdef DOUBLE_PRECISION
  defs.addDefine("Real_t", "double");
  defs.addDefine("DOUBLE_PRECISION", 1);
#else
  defs.addDefine("Real_t", "float");
  defs.addDefine("DOUBLE_PRECISION", 0);
#endif
  defs.addDefine("Index_t", "int");
  defs.addDefine("Int_t", "int");

  AddNodeForcesFromElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/AddNodeForcesFromElems_kernel.cu",
     "AddNodeForcesFromElems_kernel", defs);

  CalcAccelerationForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcAccelerationForNodes_kernel.cu",
     "CalcAccelerationForNodes_kernel", defs);

  ApplyAccelerationBoundaryConditionsForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/ApplyAccelerationBoundaryConditionsForNodes_kernel.cu",
     "ApplyAccelerationBoundaryConditionsForNodes_kernel",
     defs);

  CalcPositionAndVelocityForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcPositionAndVelocityForNodes_kernel.cu",
     "CalcPositionAndVelocityForNodes_kernel", defs);

  CalcKinematicsAndMonotonicQGradient_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcKinematicsAndMonotonicQGradient_kernel.cu",
     "CalcKinematicsAndMonotonicQGradient_kernel", defs);

  CalcMonotonicQRegionForElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcMonotonicQRegionForElems_kernel.cu",
     "CalcMonotonicQRegionForElems_kernel", defs);

  ApplyMaterialPropertiesAndUpdateVolume_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/ApplyMaterialPropertiesAndUpdateVolume_kernel.cu",
     "ApplyMaterialPropertiesAndUpdateVolume_kernel", defs);

  const int simdWidth = 32;
  {
    const int block_size = 128;

    occa::kernelInfo defs1 = defs;

    defs1.addDefine("block_size", block_size);
    defs1.addDefine("simdWidth", simdWidth);
    CalcTimeConstraintsForElems_kernel =
      occaHandle.buildKernelFromSource
      ("kernels/CalcTimeConstraintsForElems_kernel.cu",
       "CalcTimeConstraintsForElems_kernel", defs1);
  }

  {

    //    const int max_dimGrid = 1024;
    occa::kernelInfo defs1 = defs;

    defs1.addDefine("block_size", max_dimGrid);
    defs1.addDefine("simdWidth", simdWidth);

    CalcMinDtOneBlock =
      occaHandle.buildKernelFromSource
      ("kernels/CalcMinDtOneBlock.cu",
       "CalcMinDtOneBlock", defs1);
  }


  FillKernel =
    occaHandle.buildKernelFromSource
    ("kernels/Fill_kernel.cu", "Fill_kernel", defs);

  ZeroKernel =
    occaHandle.buildKernelFromSource
    ("kernels/Zero_kernel.cu", "Zero_kernel");

  CalcVolumeForceForElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcVolumeForceForElems_kernel.cu",
     "CalcVolumeForceForElems_kernel", defs);

}

#else
static void buildLuleshKernels(){

  occa::kernelInfo defs;

#ifdef DOUBLE_PRECISION
  defs.addDefine("Real_t", "double");
  defs.addDefine("DOUBLE_PRECISION", 1);
#else
  defs.addDefine("Real_t", "float");
  defs.addDefine("DOUBLE_PRECISION", 0);
#endif
  defs.addDefine("Index_t", "int");
  defs.addDefine("Int_t", "int");

  AddNodeForcesFromElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/AddNodeForcesFromElems_kernel.occa",
     "AddNodeForcesFromElems_kernel", defs);

  CalcAccelerationForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcAccelerationForNodes_kernel.occa",
     "CalcAccelerationForNodes_kernel", defs);

  ApplyAccelerationBoundaryConditionsForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/ApplyAccelerationBoundaryConditionsForNodes_kernel.occa",
     "ApplyAccelerationBoundaryConditionsForNodes_kernel",
     defs);

  CalcPositionAndVelocityForNodes_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcPositionAndVelocityForNodes_kernel.occa",
     "CalcPositionAndVelocityForNodes_kernel", defs);

  CalcKinematicsAndMonotonicQGradient_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcKinematicsAndMonotonicQGradient_kernel.occa",
     "CalcKinematicsAndMonotonicQGradient_kernel", defs);

  CalcMonotonicQRegionForElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcMonotonicQRegionForElems_kernel.occa",
     "CalcMonotonicQRegionForElems_kernel", defs);

  ApplyMaterialPropertiesAndUpdateVolume_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/ApplyMaterialPropertiesAndUpdateVolume_kernel.occa",
     "ApplyMaterialPropertiesAndUpdateVolume_kernel", defs);

  const int simdWidth = 32;
  {
    const int block_size = 128;

    occa::kernelInfo defs1 = defs;

    defs1.addDefine("block_size", block_size);
    defs1.addDefine("simdWidth", simdWidth);
    CalcTimeConstraintsForElems_kernel =
      occaHandle.buildKernelFromSource
      ("kernels/CalcTimeConstraintsForElems_kernel.occa",
       "CalcTimeConstraintsForElems_kernel", defs1);
  }

  {

    //    const int max_dimGrid = 1024;
    occa::kernelInfo defs1 = defs;

    defs1.addDefine("block_size", max_dimGrid);
    defs1.addDefine("simdWidth", simdWidth);

    CalcMinDtOneBlock =
      occaHandle.buildKernelFromSource
      ("kernels/CalcMinDtOneBlock.occa",
       "CalcMinDtOneBlock", defs1);
  }


  FillKernel =
    occaHandle.buildKernelFromSource
    ("kernels/Fill_kernel.occa", "Fill_kernel", defs);

  ZeroKernel =
    occaHandle.buildKernelFromSource
    ("kernels/Zero_kernel.occa", "Zero_kernel");


  CalcVolumeForceForElems_kernel =
    occaHandle.buildKernelFromSource
    ("kernels/CalcVolumeForceForElems_kernel.occa",
     "CalcVolumeForceForElems_kernel", defs);

}

#endif

static void occa_init(char model[], const int plat, const int dev){

  occa::availableDevices<occa::OpenCL>();

  occaHandle.setup(model, plat, dev);

  buildLuleshKernels();

}

#define MAX_THREAD_BLOCKS 65535

occa::memory occaMalloc(const size_t nbytes,
			void *host_ptr = NULL){

  occa::memory *a = new occa::memory;

  if(host_ptr != NULL){

    *a = occaHandle.malloc(nbytes, host_ptr);

  }
  else{

    *a = occaHandle.malloc(nbytes);

    size_t dim = 1;
    occa::dim inner(max_dimGrid);

    int nblocks = (nbytes + max_dimGrid-1)/max_dimGrid;

    if(nblocks > MAX_THREAD_BLOCKS)
      nblocks = MAX_THREAD_BLOCKS;

    occa::dim outer(nblocks);

    ZeroKernel.setWorkingDims(dim, inner, outer);

    ZeroKernel(nbytes, *a);
  }

  return *a;

}

void fill(occa::memory &a, Real_t val){

  const int block_size = max_dimGrid;
  const int n = a.bytes()/sizeof(Real_t);

  if(n){
    int numBlocks = (n+block_size-1)/block_size;

    if(numBlocks > MAX_THREAD_BLOCKS)
      numBlocks = MAX_THREAD_BLOCKS;

    size_t dims = 1;
    occa::dim inner(block_size);
    occa::dim outer(numBlocks);

    FillKernel.setWorkingDims(dims, inner, outer);

    FillKernel(n, val, a);
  }
}



void AllocateNodalPersistent(Domain* domain,
			     size_t domNodes){
  // domain->x.resize(domNodes) ;  /* coordinates */
  // domain->y.resize(domNodes) ;
  // domain->z.resize(domNodes) ;
  domain->x = occaMalloc(domNodes*sizeof(Real_t));
  domain->y = occaMalloc(domNodes*sizeof(Real_t));
  domain->z = occaMalloc(domNodes*sizeof(Real_t));

  // TODO:
  // domain->tex_x.initialize(domain->x.raw(),domNodes);
  // domain->tex_y.initialize(domain->y.raw(),domNodes);
  // domain->tex_z.initialize(domain->z.raw(),domNodes);

  // domain->tex_x = domain->x;
  // domain->tex_y = domain->y;
  // domain->tex_z = domain->z;

  // std::vector<float> tempX(domNodes);
  // std::vector<float> tempY(domNodes);
  // std::vector<float> tempZ(domNodes);

  // domain->x.copyTo(&(tempX[0]));
  // domain->y.copyTo(&(tempY[0]));
  // domain->z.copyTo(&(tempZ[0]));

  // domain->tex_x = occaHandle.talloc(1, occa::dim(domNodes), &(tempX[0]), occa::floatFormat, occa::readWrite);
  // domain->tex_y = occaHandle.talloc(1, occa::dim(domNodes), &(tempY[0]), occa::floatFormat, occa::readWrite);
  // domain->tex_z = occaHandle.talloc(1, occa::dim(domNodes), &(tempZ[0]), occa::floatFormat, occa::readWrite);

  domain->tex_x = occaMalloc(domNodes*sizeof(Real_t));
  domain->tex_y = occaMalloc(domNodes*sizeof(Real_t));
  domain->tex_z = occaMalloc(domNodes*sizeof(Real_t));

  // occaCompare<Real_t>(domain->x, domain->tex_x);
  // occaCompare<Real_t>(domain->y, domain->tex_y);
  // occaCompare<Real_t>(domain->z, domain->tex_z);



  // domain->xd.resize(domNodes) ; /* velocities */
  // domain->yd.resize(domNodes) ;
  // domain->zd.resize(domNodes) ;
  domain->xd = occaMalloc(domNodes*sizeof(Real_t));
  domain->yd = occaMalloc(domNodes*sizeof(Real_t));
  domain->zd = occaMalloc(domNodes*sizeof(Real_t));

  // TODO: use textures instead of global
  // domain->tex_xd.initialize(domain->xd.raw(),domNodes);
  // domain->tex_yd.initialize(domain->yd.raw(),domNodes);
  // domain->tex_zd.initialize(domain->zd.raw(),domNodes);

  // domain->tex_xd = domain->xd;
  // domain->tex_yd = domain->yd;
  // domain->tex_zd = domain->zd;
  // #ifndef DOUBLE_PRECISION
  // domain->xd.copyTo(&(tempX[0]));
  // domain->yd.copyTo(&(tempY[0]));
  // domain->zd.copyTo(&(tempZ[0]));

  // domain->tex_xd = occaHandle.talloc(1, occa::dim(domNodes), &(tempX[0]), occa::floatFormat, occa::readWrite);
  // domain->tex_yd = occaHandle.talloc(1, occa::dim(domNodes), &(tempY[0]), occa::floatFormat, occa::readWrite);
  // domain->tex_zd = occaHandle.talloc(1, occa::dim(domNodes), &(tempZ[0]), occa::floatFormat, occa::readWrite);

  domain->tex_xd = occaMalloc(domNodes*sizeof(Real_t));
  domain->tex_yd = occaMalloc(domNodes*sizeof(Real_t));
  domain->tex_zd = occaMalloc(domNodes*sizeof(Real_t));

  // occaCompare<Real_t>(domain->xd, domain->tex_xd);
  // occaCompare<Real_t>(domain->yd, domain->tex_yd);
  // occaCompare<Real_t>(domain->zd, domain->tex_zd);
  //#endif
  // domain->xdd.resize(domNodes) ; /* accelerations */
  // domain->ydd.resize(domNodes) ;
  // domain->zdd.resize(domNodes) ;
  domain->xdd = occaMalloc(domNodes*sizeof(Real_t));
  domain->ydd = occaMalloc(domNodes*sizeof(Real_t));
  domain->zdd = occaMalloc(domNodes*sizeof(Real_t));

  // domain->fx.resize(domNodes) ;  /* forces */
  // domain->fy.resize(domNodes) ;
  // domain->fz.resize(domNodes) ;
  domain->fx = occaMalloc(domNodes*sizeof(Real_t));
  domain->fy = occaMalloc(domNodes*sizeof(Real_t));
  domain->fz = occaMalloc(domNodes*sizeof(Real_t));

  //  domain->nodalMass.resize(domNodes) ;  /* mass */
  domain->nodalMass = occaMalloc(domNodes*sizeof(Real_t));
}

void AllocateElemPersistent(Domain* domain, size_t domElems, size_t padded_domElems){
  // domain->matElemlist.resize(domElems) ;  /* material indexset */
  // domain->nodelist.resize(8*padded_domElems) ;   /* elemToNode connectivity */
  domain->matElemlist = occaMalloc(domElems*sizeof(Index_t));
  domain->nodelist = occaMalloc(8*padded_domElems*sizeof(Index_t));

  // domain->lxim.resize(domElems) ; /* elem connectivity through face */
  // domain->lxip.resize(domElems) ;
  // domain->letam.resize(domElems) ;
  // domain->letap.resize(domElems) ;
  // domain->lzetam.resize(domElems) ;
  // domain->lzetap.resize(domElems) ;
  domain->lxim = occaMalloc(domElems*sizeof(Index_t));
  domain->lxip = occaMalloc(domElems*sizeof(Index_t));
  domain->letam = occaMalloc(domElems*sizeof(Index_t));
  domain->letap = occaMalloc(domElems*sizeof(Index_t));
  domain->lzetam = occaMalloc(domElems*sizeof(Index_t));
  domain->lzetap = occaMalloc(domElems*sizeof(Index_t));

  // domain->elemBC.resize(domElems) ;  /* elem face symm/free-surf flag */
  domain->elemBC = occaMalloc(domElems*sizeof(Int_t));

  // domain->e.resize(domElems) ;   /* energy */
  // domain->p.resize(domElems) ;   /* pressure */

  // domain->q.resize(domElems) ;   /* q */
  // domain->ql.resize(domElems) ;  /* linear term for q */
  // domain->qq.resize(domElems) ;  /* quadratic term for q */

  // domain->v.resize(domElems) ;     /* relative volume */

  // domain->volo.resize(domElems) ;  /* reference volume */
  // domain->delv.resize(domElems) ;  /* m_vnew - m_v */
  // domain->vdov.resize(domElems) ;  /* volume derivative over volume */

  // domain->arealg.resize(domElems) ;  /* elem characteristic length */

  // domain->ss.resize(domElems) ;      /* "sound speed" */

  // domain->elemMass.resize(domElems) ;  /* mass */

  domain->e = occaMalloc(domElems*sizeof(Real_t));
  domain->p = occaMalloc(domElems*sizeof(Real_t));
  domain->q = occaMalloc(domElems*sizeof(Real_t));
  domain->ql = occaMalloc(domElems*sizeof(Real_t));
  domain->qq = occaMalloc(domElems*sizeof(Real_t));
  domain->v = occaMalloc(domElems*sizeof(Real_t));
  domain->volo = occaMalloc(domElems*sizeof(Real_t));
  domain->delv = occaMalloc(domElems*sizeof(Real_t));
  domain->vdov = occaMalloc(domElems*sizeof(Real_t));
  domain->arealg = occaMalloc(domElems*sizeof(Real_t));
  domain->ss = occaMalloc(domElems*sizeof(Real_t));
  domain->elemMass = occaMalloc(domElems*sizeof(Real_t));


}

void AllocateSymmX(Domain* domain, size_t size){
  // domain->symmX.resize(size) ;
  domain->symmX = occaMalloc(size*sizeof(Index_t));
}

void AllocateSymmY(Domain* domain, size_t size){
  // domain->symmY.resize(size) ;
  domain->symmY = occaMalloc(size*sizeof(Index_t));
}

void AllocateSymmZ(Domain* domain, size_t size){
  // domain->symmZ.resize(size) ;
  domain->symmZ = occaMalloc(size*sizeof(Index_t));
}

void InitializeFields(Domain* domain){
  /* Basic Field Initialization */
  // thrust::fill(domain->ss.begin(),domain->ss.end(),0.);
  // thrust::fill(domain->e.begin(),domain->e.end(),0.);
  // thrust::fill(domain->p.begin(),domain->p.end(),0.);
  // thrust::fill(domain->q.begin(),domain->q.end(),0.);
  // thrust::fill(domain->v.begin(),domain->v.end(),1.);
  fill(domain->ss, 0.);
  fill(domain->e, 0.);
  fill(domain->p, 0.);
  fill(domain->q, 0.);
  fill(domain->v, 1.);

  // thrust::fill(domain->xd.begin(),domain->xd.end(),0.);
  // thrust::fill(domain->yd.begin(),domain->yd.end(),0.);
  // thrust::fill(domain->zd.begin(),domain->zd.end(),0.);
  fill(domain->xd, 0.);
  fill(domain->yd, 0.);
  fill(domain->zd, 0.);

  // thrust::fill(domain->xdd.begin(),domain->xdd.end(),0.);
  // thrust::fill(domain->ydd.begin(),domain->ydd.end(),0.);
  // thrust::fill(domain->zdd.begin(),domain->zdd.end(),0.);
  fill(domain->xdd, 0.);
  fill(domain->ydd, 0.);
  fill(domain->zdd, 0.);

  // thrust::fill(domain->nodalMass.begin(),domain->nodalMass.end(),0.);
  fill(domain->nodalMass, 0.);






  // fill(domain->ql, 0.);
  // fill(domain->qq, 0.);
  // fill(domain->volo, 0.);
  // fill(domain->delv, 0.);
  // fill(domain->vdov, 0.);
  // fill(domain->arealg, 0.);
  // fill(domain->elemMass, 0.);
  // fill(domain->vnew, 0.);
  // fill(domain->delv_xi, 0.);
  // fill(domain->delv_eta, 0.);
  // fill(domain->delv_zeta, 0.);
  // fill(domain->dxx, 0.);
  // fill(domain->dyy, 0.);

  // fill(domain->x, 0.);
  // fill(domain->y, 0.);
  // fill(domain->z, 0.);

  // fill(domain->tex_x, 0.);
  // fill(domain->tex_y, 0.);
  // fill(domain->tex_z, 0.);


  // fill(domain->tex_xd, 0.);
  // fill(domain->tex_yd, 0.);
  // fill(domain->tex_zd, 0.);

  // fill(domain->fx, 0.);
  // fill(domain->fy, 0.);
  // fill(domain->fz, 0.);
}

Domain *NewDomain(char* argv[], Index_t nx, bool structured)
{

  Domain *domain = new Domain ;

  // domain->max_streams = 32;
  // domain->streams.resize(domain->max_streams);

  // for (Int_t i=0;i<domain->max_streams;i++)
  //   cudaStreamCreate(&(domain->streams[i]));

  // cudaEventCreateWithFlags(&domain->time_constraint_computed,cudaEventDisableTiming);

  Index_t domElems;
  Index_t domNodes;
  Index_t padded_domElems;

  // Vector_h<Index_t> nodelist_h;
  // Vector_h<Real_t> x_h;
  // Vector_h<Real_t> y_h;
  // Vector_h<Real_t> z_h;

  std::vector<Index_t> nodelist_h;
  std::vector<Real_t> x_h;
  std::vector<Real_t> y_h;
  std::vector<Real_t> z_h;

  if (structured)
    {
      Real_t tx, ty, tz ;
      Index_t nidx, zidx;

      Index_t edgeElems = nx ;
      Index_t edgeNodes = edgeElems+1 ;

      domain->sizeX = edgeElems ;
      domain->sizeY = edgeElems ;
      domain->sizeZ = edgeElems ;
      domain->numElem = domain->sizeX*domain->sizeY*domain->sizeZ ;
      domain->padded_numElem = PAD(domain->numElem,32);

      domain->numNode = (domain->sizeX+1)*(domain->sizeY+1)*(domain->sizeZ+1) ;
      domain->padded_numNode = PAD(domain->numNode,32);

      domElems = domain->numElem ;
      domNodes = domain->numNode ;
      padded_domElems = domain->padded_numElem ;

      AllocateElemPersistent(domain,domElems,padded_domElems);
      AllocateNodalPersistent(domain,domNodes);

      InitializeFields(domain);

      /* initialize nodal coordinates */

      x_h.resize(domNodes);
      y_h.resize(domNodes);
      z_h.resize(domNodes);

      nidx = 0 ;
      tz = ((Real_t) 0.0) ;
      for (Index_t plane=0; plane<edgeNodes; ++plane) {
	ty = ((Real_t) 0.0) ;
	for (Index_t row=0; row<edgeNodes; ++row) {
          tx = ((Real_t) 0.0) ;
          for (Index_t col=0; col<edgeNodes; ++col) {
	    x_h[nidx] = tx ;
	    y_h[nidx] = ty ;
	    z_h[nidx] = tz ;
	    ++nidx ;
	    tx = ((Real_t) 1.125)*((Real_t) col+1)/((Real_t) nx) ;
          }
          ty = ((Real_t) 1.125)*((Real_t) row+1)/((Real_t) nx) ;
	}
	tz = ((Real_t) 1.125)*((Real_t) plane+1)/((Real_t) nx) ;
      }

      // domain->x = x_h;
      // domain->y = y_h;
      // domain->z = z_h;
      domain->x.copyFrom(&(x_h[0]));
      domain->y.copyFrom(&(y_h[0]));
      domain->z.copyFrom(&(z_h[0]));

      domain->tex_x.copyFrom(&(x_h[0]));
      domain->tex_y.copyFrom(&(y_h[0]));
      domain->tex_z.copyFrom(&(z_h[0]));

      /* embed hexehedral elements in nodal point lattice */

      nodelist_h.resize(padded_domElems*8);
      nidx = 0 ;
      zidx = 0 ;
      for (Index_t plane=0; plane<edgeElems; ++plane) {
	for (Index_t row=0; row<edgeElems; ++row) {
          for (Index_t col=0; col<edgeElems; ++col) {
	    nodelist_h[0*padded_domElems+zidx] = nidx                                       ;
	    nodelist_h[1*padded_domElems+zidx] = nidx                                   + 1 ;
	    nodelist_h[2*padded_domElems+zidx] = nidx                       + edgeNodes + 1 ;
	    nodelist_h[3*padded_domElems+zidx] = nidx                       + edgeNodes     ;
	    nodelist_h[4*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes                 ;
	    nodelist_h[5*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes             + 1 ;
	    nodelist_h[6*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
	    nodelist_h[7*padded_domElems+zidx] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
	    ++zidx ;
	    ++nidx ;
          }
          ++nidx ;
	}
	nidx += edgeNodes ;
      }
      //    domain->nodelist = nodelist_h;
      domain->nodelist.copyFrom(&(nodelist_h[0]));

      domain->numSymmX = (edgeElems+1)*(edgeElems+1) ;
      domain->numSymmY = (edgeElems+1)*(edgeElems+1) ;
      domain->numSymmZ = (edgeElems+1)*(edgeElems+1) ;

      AllocateSymmX(domain,edgeNodes*edgeNodes);
      AllocateSymmY(domain,edgeNodes*edgeNodes);
      AllocateSymmZ(domain,edgeNodes*edgeNodes);

      /* set up symmetry nodesets */

      // Vector_h<Index_t> symmX_h(domain->symmX.size());
      // Vector_h<Index_t> symmY_h(domain->symmY.size());
      // Vector_h<Index_t> symmZ_h(domain->symmZ.size());

      std::vector<Index_t> symmX_h(domain->symmX.bytes()/sizeof(Index_t));
      std::vector<Index_t> symmY_h(domain->symmY.bytes()/sizeof(Index_t));
      std::vector<Index_t> symmZ_h(domain->symmZ.bytes()/sizeof(Index_t));

      nidx = 0 ;
      for (Index_t i=0; i<edgeNodes; ++i) {
	Index_t planeInc = i*edgeNodes*edgeNodes ;
	Index_t rowInc   = i*edgeNodes ;
	for (Index_t j=0; j<edgeNodes; ++j) {
          symmX_h[nidx] = planeInc + j*edgeNodes ;
          symmY_h[nidx] = planeInc + j ;
          symmZ_h[nidx] = rowInc   + j ;
          ++nidx ;
	}
      }

      // domain->symmX = symmX_h;
      // domain->symmY = symmY_h;
      // domain->symmZ = symmZ_h;
      domain->symmX.copyFrom(&(symmX_h[0]));
      domain->symmY.copyFrom(&(symmY_h[0]));
      domain->symmZ.copyFrom(&(symmZ_h[0]));


      // Vector_h<Index_t> lxim_h(domElems);
      // Vector_h<Index_t> lxip_h(domElems);
      // Vector_h<Index_t> letam_h(domElems);
      // Vector_h<Index_t> letap_h(domElems);
      // Vector_h<Index_t> lzetam_h(domElems);
      // Vector_h<Index_t> lzetap_h(domElems);

      std::vector<Index_t> lxim_h(domElems);
      std::vector<Index_t> lxip_h(domElems);
      std::vector<Index_t> letam_h(domElems);
      std::vector<Index_t> letap_h(domElems);
      std::vector<Index_t> lzetam_h(domElems);
      std::vector<Index_t> lzetap_h(domElems);

      /* set up elemement connectivity information */
      lxim_h[0] = 0 ;
      for (Index_t i=1; i<domElems; ++i) {
	lxim_h[i]   = i-1 ;
	lxip_h[i-1] = i ;
      }
      lxip_h[domElems-1] = domElems-1 ;

      for (Index_t i=0; i<edgeElems; ++i) {
	letam_h[i] = i ;
	letap_h[domElems-edgeElems+i] = domElems-edgeElems+i ;
      }
      for (Index_t i=edgeElems; i<domElems; ++i) {
	letam_h[i] = i-edgeElems ;
	letap_h[i-edgeElems] = i ;
      }

      for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
	lzetam_h[i] = i ;
	lzetap_h[domElems-edgeElems*edgeElems+i] = domElems-edgeElems*edgeElems+i ;
      }
      for (Index_t i=edgeElems*edgeElems; i<domElems; ++i) {
	lzetam_h[i] = i - edgeElems*edgeElems ;
	lzetap_h[i-edgeElems*edgeElems] = i ;
      }

      /* set up boundary condition information */
      // Vector_h<Index_t> elemBC_h(domElems);
      std::vector<Index_t> elemBC_h(domElems);
      for (Index_t i=0; i<domElems; ++i) {
	elemBC_h[i] = 0 ;  /* clear BCs by default */
      }

      /* symmetry plane or free surface BCs */
      for (Index_t i=0; i<edgeElems; ++i) {
	Index_t planeInc = i*edgeElems*edgeElems ;
	Index_t rowInc   = i*edgeElems ;
	for (Index_t j=0; j<edgeElems; ++j) {

	  elemBC_h[planeInc+j*edgeElems] |= XI_M_SYMM ;
	  elemBC_h[planeInc+j*edgeElems+edgeElems-1] |= XI_P_FREE ;
	  elemBC_h[planeInc+j] |= ETA_M_SYMM ;
	  elemBC_h[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE ;
	  elemBC_h[rowInc+j] |= ZETA_M_SYMM ;
	  elemBC_h[rowInc+j+domElems-edgeElems*edgeElems] |= ZETA_P_FREE ;
	}
      }

      // domain->lxim = lxim_h;
      // domain->lxip = lxip_h;
      // domain->letam = letam_h;
      // domain->letap = letap_h;
      // domain->lzetam = lzetam_h;
      // domain->lzetap = lzetap_h;
      // domain->elemBC = elemBC_h;
      domain->lxim.copyFrom(&(lxim_h[0]));
      domain->lxip.copyFrom(&(lxip_h[0]));
      domain->letam.copyFrom(&(letam_h[0]));
      domain->letap.copyFrom(&(letap_h[0]));
      domain->lzetam.copyFrom(&(lzetam_h[0]));
      domain->lzetap.copyFrom(&(lzetap_h[0]));
      domain->elemBC.copyFrom(&(elemBC_h[0]));

      /* deposit energy */
      domain->octantCorner = 0;
      // domain->e[domain->octantCorner] = ((Real_t) 3.948746e+7) ;

      Real_t val = ((Real_t) 3.948746e+7);

      domain->e.copyFrom(&val, sizeof(Real_t), 0);
    }
  else
    {
      FILE *fp;
      int ee, en;

      if ((fp = fopen(argv[2], "r")) == 0) {
	printf("could not open file %s\n", argv[2]) ;
	exit( LFileError ) ;
      }

      bool fsuccess;
      fsuccess = fscanf(fp, "%d %d", &ee, &en) ;
      domain->numElem = Index_t(ee);
      domain->padded_numElem = PAD(domain->numElem,32);

      domain->numNode = Index_t(en);
      domain->padded_numNode = PAD(domain->numNode,32);

      domElems = domain->numElem ;
      domNodes = domain->numNode ;
      padded_domElems = domain->padded_numElem ;

      AllocateElemPersistent(domain,domElems,padded_domElems);
      AllocateNodalPersistent(domain,domNodes);

      InitializeFields(domain);

      /* initialize nodal coordinates */
      x_h.resize(domNodes);
      y_h.resize(domNodes);
      z_h.resize(domNodes);

      for (Index_t i=0; i<domNodes; ++i) {
	double px, py, pz ;
	fsuccess = fscanf(fp, "%lf %lf %lf", &px, &py, &pz) ;
	x_h[i] = ((Real_t) px) ;
	y_h[i] = ((Real_t) py) ;
	z_h[i] = ((Real_t) pz) ;
      }
      // domain->x = x_h;
      // domain->y = y_h;
      // domain->z = z_h;
      domain->x.copyFrom(&(x_h[0]));
      domain->y.copyFrom(&(y_h[0]));
      domain->z.copyFrom(&(z_h[0]));

      domain->tex_x.copyFrom(&(x_h[0]));
      domain->tex_y.copyFrom(&(y_h[0]));
      domain->tex_z.copyFrom(&(z_h[0]));

      /* embed hexehedral elements in nodal point lattice */
      nodelist_h.resize(padded_domElems*8);
      for (Index_t zidx=0; zidx<domElems; ++zidx) {
	for (Index_t ni=0; ni<Index_t(8); ++ni) {
          int n ;
          fsuccess = fscanf(fp, "%d", &n) ;
          nodelist_h[ni*padded_domElems+zidx] = Index_t(n);
	}
      }
      // domain->nodelist = nodelist_h;
      domain->nodelist.copyFrom(&(nodelist_h[0]));

      /* set up face-based element neighbors */
      // Vector_h<Index_t> lxim_h(domElems);
      // Vector_h<Index_t> lxip_h(domElems);
      // Vector_h<Index_t> letam_h(domElems);
      // Vector_h<Index_t> letap_h(domElems);
      // Vector_h<Index_t> lzetam_h(domElems);
      // Vector_h<Index_t> lzetap_h(domElems);

      std::vector<Index_t> lxim_h(domElems);
      std::vector<Index_t> lxip_h(domElems);
      std::vector<Index_t> letam_h(domElems);
      std::vector<Index_t> letap_h(domElems);
      std::vector<Index_t> lzetam_h(domElems);
      std::vector<Index_t> lzetap_h(domElems);

      for (Index_t i=0; i<domElems; ++i) {
	int xi_m, xi_p, eta_m, eta_p, zeta_m, zeta_p ;
	fsuccess = fscanf(fp, "%d %d %d %d %d %d",
			  &xi_m, &xi_p, &eta_m, &eta_p, &zeta_m, &zeta_p) ;

	lxim_h[i]   = Index_t(xi_m) ;
	lxip_h[i]   = Index_t(xi_p) ;
	letam_h[i]  = Index_t(eta_m) ;
	letap_h[i]  = Index_t(eta_p) ;
	lzetam_h[i] = Index_t(zeta_m) ;
	lzetap_h[i] = Index_t(zeta_p) ;
      }

      // domain->lxim = lxim_h;
      // domain->lxip = lxip_h;
      // domain->letam = letam_h;
      // domain->letap = letap_h;
      // domain->lzetam = lzetam_h;
      // domain->lzetap = lzetap_h;
      domain->lxim.copyFrom(&(lxim_h[0]));
      domain->lxip.copyFrom(&(lxip_h[0]));
      domain->letam.copyFrom(&(letam_h[0]));
      domain->letap.copyFrom(&(letap_h[0]));
      domain->lzetam.copyFrom(&(lzetam_h[0]));
      domain->lzetap.copyFrom(&(lzetap_h[0]));

      /* set up X symmetry nodeset */

      fsuccess = fscanf(fp, "%d", &domain->numSymmX) ;
      // Vector_h<Index_t> symmX_h(domain->numSymmX);
      std::vector<Index_t> symmX_h(domain->numSymmX);
      for (Index_t i=0; i<domain->numSymmX; ++i) {
	int n ;
	fsuccess = fscanf(fp, "%d", &n) ;
	symmX_h[i] = Index_t(n) ;
      }
      // domain->symmX = symmX_h;
      domain->symmX.copyFrom(&(symmX_h[0]));

      fsuccess = fscanf(fp, "%d", &domain->numSymmY) ;
      // Vector_h<Index_t> symmY_h(domain->numSymmY);
      std::vector<Index_t> symmY_h(domain->numSymmY);
      for (Index_t i=0; i<domain->numSymmY; ++i) {
	int n ;
	fsuccess = fscanf(fp, "%d", &n) ;
	symmY_h[i] = Index_t(n) ;
      }
      // domain->symmY = symmY_h;
      domain->symmY.copyFrom(&(symmY_h[0]));

      fsuccess = fscanf(fp, "%d", &domain->numSymmZ) ;
      // Vector_h<Index_t> symmZ_h(domain->numSymmZ);
      std::vector<Index_t> symmZ_h(domain->numSymmZ);
      for (Index_t i=0; i<domain->numSymmZ; ++i) {
	int n ;
	fsuccess = fscanf(fp, "%d", &n) ;
	symmZ_h[i] = Index_t(n) ;
      }
      // domain->symmZ = symmZ_h;
      domain->symmZ.copyFrom(&(symmZ_h[0]));
      /* set up free surface nodeset */
      Index_t numFreeSurf;
      fsuccess = fscanf(fp, "%d", &numFreeSurf) ;
      // Vector_h<Index_t> freeSurf_h(numFreeSurf);
      std::vector<Index_t> freeSurf_h(numFreeSurf);
      for (Index_t i=0; i<numFreeSurf; ++i) {
	int n ;
	fsuccess = fscanf(fp, "%d", &n) ;
	freeSurf_h[i] = Index_t(n) ;
      }

      fclose(fp);

      /* set up boundary condition information */
      // Vector_h<Index_t> elemBC_h(domElems);
      // Vector_h<Index_t> surfaceNode_h(domNodes);

      std::vector<Index_t> elemBC_h(domElems);
      std::vector<Index_t> surfaceNode_h(domNodes);

      for (Index_t i=0; i<domain->numElem; ++i) {
	elemBC_h[i] = 0 ;
      }

      for (Index_t i=0; i<domain->numNode; ++i) {
	surfaceNode_h[i] = 0 ;
      }

      for (Index_t i=0; i<domain->numSymmX; ++i) {
	surfaceNode_h[symmX_h[i]] = 1 ;
      }

      for (Index_t i=0; i<domain->numSymmY; ++i) {
	surfaceNode_h[symmY_h[i]] = 1 ;
      }

      for (Index_t i=0; i<domain->numSymmZ; ++i) {
	surfaceNode_h[symmZ_h[i]] = 1 ;
      }

      for (Index_t zidx=0; zidx<domain->numElem; ++zidx) {
	Int_t mask = 0 ;

	for (Index_t ni=0; ni<8; ++ni) {
          mask |= (surfaceNode_h[nodelist_h[ni*domain->padded_numElem+zidx]] << ni) ;
	}

	if ((mask & 0x0f) == 0x0f) elemBC_h[zidx] |= ZETA_M_SYMM ;
	if ((mask & 0xf0) == 0xf0) elemBC_h[zidx] |= ZETA_P_SYMM ;
	if ((mask & 0x33) == 0x33) elemBC_h[zidx] |= ETA_M_SYMM ;
	if ((mask & 0xcc) == 0xcc) elemBC_h[zidx] |= ETA_P_SYMM ;
	if ((mask & 0x99) == 0x99) elemBC_h[zidx] |= XI_M_SYMM ;
	if ((mask & 0x66) == 0x66) elemBC_h[zidx] |= XI_P_SYMM ;
      }

      for (Index_t zidx=0; zidx<domain->numElem; ++zidx) {
	if (elemBC_h[zidx] == (XI_M_SYMM | ETA_M_SYMM | ZETA_M_SYMM)) {
          domain->octantCorner = zidx ;
          break ;
	}
      }

      for (Index_t i=0; i<domain->numNode; ++i) {
	surfaceNode_h[i] = 0 ;
      }

      for (Index_t i=0; i<numFreeSurf; ++i) {
	surfaceNode_h[freeSurf_h[i]] = 1 ;
      }

      for (Index_t zidx=0; zidx<domain->numElem; ++zidx) {
	Int_t mask = 0 ;

	for (Index_t ni=0; ni<8; ++ni) {
          mask |= (surfaceNode_h[nodelist_h[ni*domain->padded_numElem+zidx]] << ni) ;
	}

	if ((mask & 0x0f) == 0x0f) elemBC_h[zidx] |= ZETA_M_SYMM ;
	if ((mask & 0xf0) == 0xf0) elemBC_h[zidx] |= ZETA_P_SYMM ;
	if ((mask & 0x33) == 0x33) elemBC_h[zidx] |= ETA_M_SYMM ;
	if ((mask & 0xcc) == 0xcc) elemBC_h[zidx] |= ETA_P_SYMM ;
	if ((mask & 0x99) == 0x99) elemBC_h[zidx] |= XI_M_SYMM ;
	if ((mask & 0x66) == 0x66) elemBC_h[zidx] |= XI_P_SYMM ;
      }

      // domain->elemBC = elemBC_h;
      domain->elemBC.copyFrom(&(elemBC_h[0]));
      /* deposit energy */
      //    domain->e[domain->octantCorner] = ((Real_t) 3.948746e+7) ;
      Real_t val = ((Real_t) 3.948746e+7) ;

      domain->e.copyFrom(&val, sizeof(Real_t), 0);
    }

  /* set up node-centered indexing of elements */
  // Vector_h<Index_t> nodeElemCount_h(domNodes);
  std::vector<Index_t> nodeElemCount_h(domNodes);

  for (Index_t i=0; i<domNodes; ++i) {
    nodeElemCount_h[i] = 0 ;
  }

  for (Index_t i=0; i<domElems; ++i) {
    for (Index_t j=0; j < 8; ++j) {
      ++(nodeElemCount_h[nodelist_h[j*padded_domElems+i]]);
    }
  }

  // Vector_h<Index_t> nodeElemStart_h(domNodes);
  std::vector<Index_t> nodeElemStart_h(domNodes);

  nodeElemStart_h[0] = 0;
  for (Index_t i=1; i < domNodes; ++i) {
    nodeElemStart_h[i] =
      nodeElemStart_h[i-1] + nodeElemCount_h[i-1] ;
  }

  // Vector_h<Index_t> nodeElemCornerList_h(nodeElemStart_h[domNodes-1] +
  //                nodeElemCount_h[domNodes-1] );

  std::vector<Index_t> nodeElemCornerList_h(nodeElemStart_h[domNodes-1] +
					    nodeElemCount_h[domNodes-1] );

  for (Index_t i=0; i < domNodes; ++i) {
    nodeElemCount_h[i] = 0;
  }

  for (Index_t j=0; j < 8; ++j) {
    for (Index_t i=0; i < domElems; ++i) {
      Index_t m = nodelist_h[padded_domElems*j+i];
      Index_t k = padded_domElems*j + i ;
      Index_t offset = nodeElemStart_h[m] +
	nodeElemCount_h[m] ;
      nodeElemCornerList_h[offset] = k;
      ++(nodeElemCount_h[m]) ;
    }
  }

  Index_t clSize = nodeElemStart_h[domNodes-1] +
    nodeElemCount_h[domNodes-1] ;
  for (Index_t i=0; i < clSize; ++i) {
    Index_t clv = nodeElemCornerList_h[i] ;
    if ((clv < 0) || (clv > padded_domElems*8)) {
      fprintf(stderr,
	      "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
      exit(1);
    }
  }

  // domain->nodeElemStart = nodeElemStart_h;
  // domain->nodeElemCount = nodeElemCount_h;
  // domain->nodeElemCornerList = nodeElemCornerList_h;

  domain->nodeElemStart = occaMalloc(domNodes*sizeof(Index_t), &(nodeElemStart_h[0]));
  domain->nodeElemCount = occaMalloc(domNodes*sizeof(Index_t), &(nodeElemCount_h[0]));
  domain->nodeElemCornerList = occaMalloc(nodeElemCornerList_h.size()*sizeof(Index_t), &(nodeElemCornerList_h[0]));


  /* Create a material IndexSet (entire domain same material for now) */
  // Vector_h<Index_t> matElemlist_h(domElems);
  std::vector<Index_t> matElemlist_h(domElems);
  for (Index_t i=0; i<domElems; ++i) {
    matElemlist_h[i] = i ;
  }
  //  domain->matElemlist = matElemlist_h;
  domain->matElemlist = occaMalloc(domElems*sizeof(Index_t), &(matElemlist_h[0]));

  // cudaMallocHost(&domain->dtcourant_h,sizeof(Real_t),0);
  // cudaMallocHost(&domain->dthydro_h,sizeof(Real_t),0);
  // cudaMallocHost(&domain->bad_vol_h,sizeof(Index_t),0);
  // cudaMallocHost(&domain->bad_q_h,sizeof(Index_t),0);

  domain->dtcourant_h = (Real_t*) malloc(sizeof(Real_t));
  domain->dthydro_h = (Real_t*) malloc(sizeof(Real_t));
  domain->bad_vol_h = (Index_t*) malloc(sizeof(Index_t));
  domain->bad_q_h = (Index_t*) malloc(sizeof(Index_t));

  *(domain->bad_vol_h)=-1;
  *(domain->bad_q_h)=-1;
  *(domain->dthydro_h)=1e20;
  *(domain->dtcourant_h)=1e20;

  /* initialize material parameters */
  domain->deltatime_h = ((Real_t) 1.0e-7) ;
  domain->time_h      = ((Real_t) 0.) ;
  domain->dtfixed = ((Real_t) -1.0e-7) ;
  domain->deltatimemultlb = ((Real_t) 1.1) ;
  domain->deltatimemultub = ((Real_t) 1.2) ;
  domain->stoptime  = ((Real_t) 1.0e-2) ;
  domain->dtmax     = ((Real_t) 1.0e-2) ;
  domain->cycle   = 0 ;

  domain->e_cut = ((Real_t) 1.0e-7) ;
  domain->p_cut = ((Real_t) 1.0e-7) ;
  domain->q_cut = ((Real_t) 1.0e-7) ;
  domain->u_cut = ((Real_t) 1.0e-7) ;
  domain->v_cut = ((Real_t) 1.0e-10) ;

  domain->hgcoef      = ((Real_t) 3.0) ;
  domain->ss4o3       = ((Real_t) 4.0)/((Real_t) 3.0) ;

  domain->qstop              =  ((Real_t) 1.0e+12) ;
  domain->monoq_max_slope    =  ((Real_t) 1.0) ;
  domain->monoq_limiter_mult =  ((Real_t) 2.0) ;
  domain->qlc_monoq          = ((Real_t) 0.5) ;
  domain->qqc_monoq          = ((Real_t) 2.0)/((Real_t) 3.0) ;
  domain->qqc                = ((Real_t) 2.0) ;

  domain->pmin =  ((Real_t) 0.) ;
  domain->emin = ((Real_t) -1.0e+15) ;

  domain->dvovmax =  ((Real_t) 0.1) ;

  domain->eosvmax =  ((Real_t) 1.0e+9) ;
  domain->eosvmin =  ((Real_t) 1.0e-9) ;

  domain->refdens =  ((Real_t) 1.0) ;

  /* initialize field data */
  // Vector_h<Real_t> nodalMass_h(domNodes);
  // Vector_h<Real_t> volo_h(domElems);
  // Vector_h<Real_t> elemMass_h(domElems);

  std::vector<Real_t> nodalMass_h(domNodes);
  std::vector<Real_t> volo_h(domElems);
  std::vector<Real_t> elemMass_h(domElems);

  for (Index_t i=0; i<domElems; ++i) {
    Real_t x_local[8], y_local[8], z_local[8] ;
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
      {
	Index_t gnode = nodelist_h[lnode*padded_domElems+i];
	x_local[lnode] = x_h[gnode];
	y_local[lnode] = y_h[gnode];
	z_local[lnode] = z_h[gnode];
      }

    // volume calculations
    Real_t volume = CalcElemVolume(x_local, y_local, z_local );
    volo_h[i] = volume ;
    elemMass_h[i] = volume ;
    for (Index_t j=0; j<8; ++j) {
      Index_t gnode = nodelist_h[j*padded_domElems+i];
      nodalMass_h[gnode] += volume / ((Real_t) 8.0) ;

      assert(nodalMass_h[gnode] > 0.);
    }
  }

  // domain->nodalMass = nodalMass_h;
  // domain->volo = volo_h;
  // domain->elemMass= elemMass_h;
  domain->nodalMass.copyFrom(&(nodalMass_h[0]));
  domain->volo.copyFrom(&(volo_h[0]));
  domain->elemMass.copyFrom(&(elemMass_h[0]));


  Index_t padded_numElem = domain->padded_numElem;

  domain->fx_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));
  domain->fy_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));
  domain->fz_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));


  int allElem = domain->numElem +  /* local elem */
    2*domain->sizeX*domain->sizeY ; /* plane ghosts */

  domain->vnew = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->dxx = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->dyy = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->dzz = occaMalloc(domain->numElem*sizeof(Real_t));


  domain->delx_xi    = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->delx_eta   = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->delx_zeta  = occaMalloc(domain->numElem*sizeof(Real_t));
  domain->delv_xi    = occaMalloc(allElem*sizeof(Real_t));
  domain->delv_eta   = occaMalloc(allElem*sizeof(Real_t));
  domain->delv_zeta  = occaMalloc(allElem*sizeof(Real_t));


  return domain ;
}

static inline
void TimeIncrement(Domain* domain){
  // To make sure dtcourant and dthydro have been updated on host
  // cudaEventSynchronize(domain->time_constraint_computed);
  Real_t targetdt = domain->stoptime - domain->time_h;

  if ((domain->dtfixed <= ((Real_t) 0.0)) && (domain->cycle != Int_t(0))) {

    Real_t ratio ;

    /* This will require a reduction in parallel */
    Real_t newdt = ((Real_t) 1.0e+20) ;

    if ( *(domain->dtcourant_h) < newdt) {
      newdt = *(domain->dtcourant_h) / ((Real_t) 2.0) ;
    }
    if ( *(domain->dthydro_h) < newdt) {
      newdt = *(domain->dthydro_h) * ((Real_t) 2.0) / ((Real_t) 3.0) ;
    }

    Real_t olddt = domain->deltatime_h;
    ratio = newdt / olddt ;
    if (ratio >= ((Real_t) 1.0)) {
      if (ratio < domain->deltatimemultlb) {
	newdt = olddt ;
      }
      else if (ratio > domain->deltatimemultub) {
	newdt = olddt*domain->deltatimemultub ;
      }
    }

    if (newdt > domain->dtmax) {
      newdt = domain->dtmax ;
    }
    domain->deltatime_h = newdt ;
  }

  /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
  if ((targetdt > domain->deltatime_h) &&
      (targetdt < (((Real_t) 4.0) * domain->deltatime_h / ((Real_t) 3.0))) ) {
    targetdt = ((Real_t) 2.0) * domain->deltatime_h / ((Real_t) 3.0) ;
  }

  if (targetdt < domain->deltatime_h) {
    domain->deltatime_h = targetdt ;
  }

  domain->time_h += domain->deltatime_h ;

  ++domain->cycle ;
}



static inline
void CalcVolumeForceForElems(const Real_t hgcoef,Domain *domain)
{
  Index_t numElem = domain->numElem ;
  Index_t padded_numElem = domain->padded_numElem;

#ifdef DOUBLE_PRECISION
  // Vector_d<Real_t>* fx_elem = Allocator< Vector_d<Real_t> >::allocate(padded_numElem*8);
  // Vector_d<Real_t>* fy_elem = Allocator< Vector_d<Real_t> >::allocate(padded_numElem*8);
  // Vector_d<Real_t>* fz_elem = Allocator< Vector_d<Real_t> >::allocate(padded_numElem*8);

  // occa::memory* fx_elem = Allocator<occa::memeory >::allocate(padded_numElem*8);
  // occa::memory* fy_elem = Allocator<occa::memeory >::allocate(padded_numElem*8);
  // occa::memory* fz_elem = Allocator<occa::memeory >::allocate(padded_numElem*8);

  // occa::memory fx_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));
  // occa::memory fy_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));
  // occa::memory fz_elem = occaMalloc(padded_numElem*8*sizeof(Real_t));

#else
  // thrust::fill(domain->fx.begin(),domain->fx.end(),0.);
  // thrust::fill(domain->fy.begin(),domain->fy.end(),0.);
  // thrust::fill(domain->fz.begin(),domain->fz.end(),0.);

  occa::tic("fill");
  fill(domain->fx, 0.);
  fill(domain->fy, 0.);
  fill(domain->fz, 0.);
  occa::toc("fill");
#endif


  occa::tic("volume force kernel");
  int num_threads = numElem ;
  const int block_size = 64;
  int dimGrid = PAD_DIV(num_threads,block_size);

  size_t dims = 1;
  occa::dim inner(block_size);
  occa::dim outer(dimGrid);

  bool hourg_gt_zero = hgcoef > ((Real_t) 0.0);

  if (hourg_gt_zero)
    {
      //       CalcVolumeForceForElems_kernel<true> <<<dimGrid,block_size>>>
      //       ( domain->volo.raw(),
      //         domain->v.raw(),
      //         domain->p.raw(),
      //         domain->q.raw(),
      // 	      hgcoef, numElem, padded_numElem,
      //         domain->nodelist.raw(),
      //         domain->ss.raw(),
      //         domain->elemMass.raw(),
      //         domain->tex_x, domain->tex_y, domain->tex_z, domain->tex_xd, domain->tex_yd, domain->tex_zd,
      // #ifdef DOUBLE_PRECISION
      //         fx_elem->raw(),
      //         fy_elem->raw(),
      //         fz_elem->raw() ,
      // #else
      //         domain->fx.raw(),
      //         domain->fy.raw(),
      //         domain->fz.raw(),
      // #endif
      //         domain->bad_vol_h,
      //         num_threads
      //       );


      occaCheckDomain(domain);

      occa::memory bad_vol = occaMalloc(sizeof(Index_t), domain->bad_vol_h);

      CalcVolumeForceForElems_kernel.setWorkingDims(dims, inner, outer);

      CalcVolumeForceForElems_kernel
	( domain->volo,
	  domain->v,
	  domain->p,
	  domain->q,
	  hgcoef, numElem, padded_numElem,
	  domain->nodelist,
	  domain->ss,
	  domain->elemMass,
	  domain->tex_x, domain->tex_y, domain->tex_z, domain->tex_xd, domain->tex_yd, domain->tex_zd,
#ifdef DOUBLE_PRECISION
	  domain->fx_elem,
	  domain->fy_elem,
	  domain->fz_elem,
#else
	  domain->fx,
	  domain->fy,
	  domain->fz,
#endif
	  bad_vol,
	  num_threads, true);

      bad_vol.copyTo(domain->bad_vol_h);

      bad_vol.free();

      occaCheckDomain(domain);

    }
  else
    {
      //       CalcVolumeForceForElems_kernel<false> <<<dimGrid,block_size>>>
      //       ( domain->volo.raw(),
      //         domain->v.raw(),
      //         domain->p.raw(),
      //         domain->q.raw(),
      // 	      hgcoef, numElem, padded_numElem,
      //         domain->nodelist.raw(),
      //         domain->ss.raw(),
      //         domain->elemMass.raw(),
      //         domain->tex_x, domain->tex_y, domain->tex_z, domain->tex_xd, domain->tex_yd, domain->tex_zd,
      // #ifdef DOUBLE_PRECISION
      //         fx_elem->raw(),
      //         fy_elem->raw(),
      //         fz_elem->raw() ,
      // #else
      //         domain->fx.raw(),
      //         domain->fy.raw(),
      //         domain->fz.raw(),
      // #endif
      //         domain->bad_vol_h,
      //         num_threads
      //       );


      occaCheck<Real_t>(domain->fx_elem);
      occaCheck<Real_t>(domain->fy_elem);
      occaCheck<Real_t>(domain->fz_elem);

      occaCheckDomain(domain);

      occa::memory bad_vol = occaMalloc(sizeof(Index_t), domain->bad_vol_h);

      CalcVolumeForceForElems_kernel.setWorkingDims(dims, inner, outer);
      CalcVolumeForceForElems_kernel
	( domain->volo,
	  domain->v,
	  domain->p,
	  domain->q,
	  hgcoef, numElem, padded_numElem,
	  domain->nodelist,
	  domain->ss,
	  domain->elemMass,
	  domain->tex_x, domain->tex_y, domain->tex_z, domain->tex_xd, domain->tex_yd, domain->tex_zd,
#ifdef DOUBLE_PRECISION
	  domain->fx_elem,
	  domain->fy_elem,
	  domain->fz_elem,
#else
	  domain->fx,
	  domain->fy,
	  domain->fz,
#endif
	  bad_vol,
	  num_threads, false);

      bad_vol.copyTo(domain->bad_vol_h);

      bad_vol.free();

      occaCheckDomain(domain);
    }

  occa::toc("volume force kernel", CalcVolumeForceForElems_kernel);

#ifdef DOUBLE_PRECISION
  num_threads = domain->numNode;

  // Launch boundary nodes first
  dimGrid= PAD_DIV(num_threads,block_size);

  // AddNodeForcesFromElems_kernel<<<dimGrid,block_size>>>
  // ( domain->numNode,
  //   domain->padded_numNode,
  //   domain->nodeElemCount.raw(),
  //   domain->nodeElemStart.raw(),
  //   domain->nodeElemCornerList.raw(),
  //   fx_elem->raw(),
  //   fy_elem->raw(),
  //   fz_elem->raw(),
  //   domain->fx.raw(),
  //   domain->fy.raw(),
  //   domain->fz.raw(),
  //   num_threads
  // );

  occa::tic("add node forces kernel");
  occaCheck<Real_t>(domain->fx_elem);
  occaCheck<Real_t>(domain->fy_elem);
  occaCheck<Real_t>(domain->fz_elem);


  dims = 1;
  inner.x = block_size;
  outer.x = dimGrid;

  AddNodeForcesFromElems_kernel.setWorkingDims(dims, inner, outer);

  AddNodeForcesFromElems_kernel
    ( domain->numNode,
      domain->padded_numNode,
      domain->nodeElemCount,
      domain->nodeElemStart,
      domain->nodeElemCornerList,
      domain->fx_elem,
      domain->fy_elem,
      domain->fz_elem,
      domain->fx,
      domain->fy,
      domain->fz,
      num_threads);


  occaCheck<Real_t>(domain->fx);
  occaCheck<Real_t>(domain->fy);

  occa::toc("add node forces kernel", AddNodeForcesFromElems_kernel);
  //cudaDeviceSynchronize();
  //cudaCheckError();

  // Allocator<Vector_d<Real_t> >::free(fx_elem,padded_numElem*8);
  // Allocator<Vector_d<Real_t> >::free(fy_elem,padded_numElem*8);
  // Allocator<Vector_d<Real_t> >::free(fz_elem,padded_numElem*8);

  // Allocator<occa::memory >::free(fx_elem, padded_numElem*8);
  // Allocator<occa::memory >::free(fy_elem, padded_numElem*8);
  // Allocator<occa::memory >::free(fz_elem, padded_numElem*8);

  // fx_elem.free();
  // fy_elem.free();
  // fz_elem.free();
#endif // ifdef DOUBLE_PRECISION
  return ;
}


static inline
void CalcVolumeForceForElems(Domain* domain)
{
  const Real_t hgcoef = domain->hgcoef ;

  CalcVolumeForceForElems(hgcoef,domain);
}

static inline void checkErrors(Domain* domain,int its)
{
  if (*(domain->bad_vol_h) != -1)
    {
      printf("Volume Error in cell %d at iteration %d\n",*(domain->bad_vol_h),its);
      exit(VolumeError);
    }

  if (*(domain->bad_q_h) != -1)
    {
      printf("Q Error in cell %d at iteration %d\n",*(domain->bad_q_h),its);
      exit(QStopError);
    }
}

static inline void CalcForceForNodes(Domain *domain)
{
  CalcVolumeForceForElems(domain);
}


static inline
void CalcAccelerationForNodes(Domain *domain)
{
  Index_t dimBlock = 128;
  Index_t dimGrid = PAD_DIV(domain->numNode,dimBlock);


  occa::tic("node forces kernel");
  size_t dims = 1;
  occa::dim inner(dimBlock);
  occa::dim outer(dimGrid);

  CalcAccelerationForNodes_kernel.setWorkingDims(dims, inner, outer);

  // CalcAccelerationForNodes_kernel<<<dimGrid, dimBlock>>>
  //     (domain->numNode,
  //      domain->xdd.raw(),domain->ydd.raw(),domain->zdd.raw(),
  //      domain->fx.raw(),domain->fy.raw(),domain->fz.raw(),
  //      domain->nodalMass.raw());

  occaCheck<Real_t>(domain->fx);
  occaCheck<Real_t>(domain->fy);
  occaCheck<Real_t>(domain->fz);
  occaCheck<Real_t>(domain->nodalMass);

  CalcAccelerationForNodes_kernel
    (domain->numNode,
     domain->xdd,
     domain->ydd,
     domain->zdd,
     domain->fx,
     domain->fy,
     domain->fz,
     domain->nodalMass);

  occa::toc("node forces kernel", CalcAccelerationForNodes_kernel);
  //cudaDeviceSynchronize();
  //cudaCheckError();
}

static inline
void ApplyAccelerationBoundaryConditionsForNodes(Domain *domain)
{

  Index_t dimBlock = 128;

  Index_t dimGrid = PAD_DIV(domain->numSymmX,dimBlock);

  occa::tic("acceleration bc kernel");
  size_t dims = 1;
  occa::dim inner(dimBlock);
  occa::dim outer(dimGrid);


  // ApplyAccelerationBoundaryConditionsForNodes_kernel<<<dimGrid, dimBlock>>>
  //     (domain->numSymmX,
  //      domain->xdd.raw(),
  //      domain->symmX.raw());

  ApplyAccelerationBoundaryConditionsForNodes_kernel.setWorkingDims(dims, inner, outer);

  ApplyAccelerationBoundaryConditionsForNodes_kernel
    (domain->numSymmX,
     domain->xdd,
     domain->symmX);

  occa::toc("acceleration bc kernel", ApplyAccelerationBoundaryConditionsForNodes_kernel);

  occa::tic("acceleration bc kernel");
  dimGrid = PAD_DIV(domain->numSymmY,dimBlock);
  outer.x = dimGrid;

  // ApplyAccelerationBoundaryConditionsForNodes_kernel<<<dimGrid, dimBlock>>>
  //     (domain->numSymmY,
  //      domain->ydd.raw(),
  //      domain->symmY.raw());

  ApplyAccelerationBoundaryConditionsForNodes_kernel.setWorkingDims(dims, inner, outer);
  ApplyAccelerationBoundaryConditionsForNodes_kernel
    (domain->numSymmY,
     domain->ydd,
     domain->symmY);

  occa::toc("acceleration bc kernel", ApplyAccelerationBoundaryConditionsForNodes_kernel);

  occa::tic("acceleration bc kernel");
  dimGrid = PAD_DIV(domain->numSymmZ,dimBlock);
  outer.x = dimGrid;
  // ApplyAccelerationBoundaryConditionsForNodes_kernel<<<dimGrid, dimBlock>>>
  //     (domain->numSymmZ,
  //      domain->zdd.raw(),
  //      domain->symmZ.raw());

  ApplyAccelerationBoundaryConditionsForNodes_kernel.setWorkingDims(dims, inner, outer);
  ApplyAccelerationBoundaryConditionsForNodes_kernel
    (domain->numSymmZ,
     domain->zdd,
     domain->symmZ);

  occa::toc("acceleration bc kernel", ApplyAccelerationBoundaryConditionsForNodes_kernel);
}



static inline
void CalcPositionAndVelocityForNodes(const Real_t u_cut, Domain* domain)
{
  Index_t dimBlock = 128;
  Index_t dimGrid = PAD_DIV(domain->numNode,dimBlock);

  occa::tic("pos n vel kernel");
  size_t dims = 1;
  occa::dim inner(dimBlock);
  occa::dim outer(dimGrid);

  // CalcPositionAndVelocityForNodes_kernel<<<dimGrid, dimBlock>>>
  //     (domain->numNode,domain->deltatime_h,u_cut,
  //      domain->x.raw(),domain->y.raw(),domain->z.raw(),
  //      domain->xd.raw(),domain->yd.raw(),domain->zd.raw(),
  //      domain->xdd.raw(),domain->ydd.raw(),domain->zdd.raw());

  // occaCompare<Real_t>(domain->x, domain->tex_x);
  // occaCompare<Real_t>(domain->y, domain->tex_y);
  // occaCompare<Real_t>(domain->z, domain->tex_z);

  // occaCompare<Real_t>(domain->xd, domain->tex_xd);
  // occaCompare<Real_t>(domain->yd, domain->tex_yd);
  // occaCompare<Real_t>(domain->zd, domain->tex_zd);

  CalcPositionAndVelocityForNodes_kernel.setWorkingDims(dims, inner, outer);
  CalcPositionAndVelocityForNodes_kernel
    (domain->numNode,
     domain->deltatime_h,
     u_cut,
     domain->tex_x,
     domain->tex_y,
     domain->tex_z,
     domain->tex_xd,
     domain->tex_yd,
     domain->tex_zd,
     domain->xdd,
     domain->ydd,
     domain->zdd);

  occa::toc("pos n vel kernel", CalcPositionAndVelocityForNodes_kernel);
  //cudaDeviceSynchronize();
  //cudaCheckError();
}

static inline
void LagrangeNodal(Domain *domain)
{

  Real_t u_cut = domain->u_cut ;

  occa::tic("calc force nodes");
  occaCheckDomain(domain);

  CalcForceForNodes(domain);
  occa::toc("calc force nodes");

  occa::tic("calc time increment");
  occaCheckDomain(domain);

  TimeIncrement(domain);
  occa::toc("calc time increment");

  occa::tic("calc acceleration");
  occaCheckDomain(domain);

  CalcAccelerationForNodes(domain);
  occa::toc("calc acceleration");

  occa::tic("acceleration bcs");
  occaCheckDomain(domain);

  ApplyAccelerationBoundaryConditionsForNodes(domain);
  occa::toc("acceleration bcs");

  occa::tic("position and velocity");
  occaCheckDomain(domain);

  CalcPositionAndVelocityForNodes(u_cut, domain);

  occaCheckDomain(domain);
  occa::toc("position and velocity");

  return;
}


static inline
void CalcKinematicsAndMonotonicQGradient(Domain *domain)
{
  Index_t numElem = domain->numElem ;
  Index_t padded_numElem = domain->padded_numElem;

  int num_threads = numElem;

  const int block_size = 64;
  int dimGrid = PAD_DIV(num_threads,block_size);

  size_t dims = 1;
  occa::dim inner(block_size);
  occa::dim outer(dimGrid);

  // CalcKinematicsAndMonotonicQGradient_kernel<<<dimGrid,block_size>>>
  // (  numElem,padded_numElem, domain->deltatime_h,
  //    domain->nodelist.raw(),
  //    domain->volo.raw(),
  //    domain->v.raw(),
  //    domain->tex_x,domain->tex_y,domain->tex_z,domain->tex_xd,domain->tex_yd,domain->tex_zd,
  //    domain->vnew->raw(),
  //    domain->delv.raw(),
  //    domain->arealg.raw(),
  //    domain->dxx->raw(),
  //    domain->dyy->raw(),
  //    domain->dzz->raw(),
  //    domain->vdov.raw(),
  //    domain->delx_zeta->raw(),
  //    domain->delv_zeta->raw(),
  //    domain->delx_xi->raw(),
  //    domain->delv_xi->raw(),
  //    domain->delx_eta->raw(),
  //    domain->delv_eta->raw(),
  //    domain->bad_vol_h,
  //    num_threads
  // );


  occaCheckDomain(domain);

  occa::memory bad_vol = occaMalloc(sizeof(Index_t),
				    domain->bad_vol_h);

  CalcKinematicsAndMonotonicQGradient_kernel.setWorkingDims(dims, inner, outer);

  CalcKinematicsAndMonotonicQGradient_kernel
    (  numElem,padded_numElem, domain->deltatime_h,
       domain->nodelist,
       domain->volo,
       domain->v,
       domain->tex_x,domain->tex_y,domain->tex_z,domain->tex_xd,domain->tex_yd,domain->tex_zd,
       domain->vnew,
       domain->delv,
       domain->arealg,
       domain->dxx,
       domain->dyy,
       domain->dzz,
       domain->vdov,
       domain->delx_zeta,
       domain->delv_zeta,
       domain->delx_xi,
       domain->delv_xi,
       domain->delx_eta,
       domain->delv_eta,
       bad_vol,
       num_threads);

  bad_vol.copyTo(domain->bad_vol_h);

  occaCheck<Real_t>(domain->vnew);
  occaCheck<Real_t>(domain->delv);
  occaCheck<Real_t>(domain->arealg);
  occaCheck<Real_t>(domain->dxx);
  occaCheck<Real_t>(domain->dyy);
  occaCheck<Real_t>(domain->dzz);
  occaCheck<Real_t>(domain->vdov);
  occaCheck<Real_t>(domain->delx_zeta);
  occaCheck<Real_t>(domain->delv_zeta);
  occaCheck<Real_t>(domain->delx_xi);
  occaCheck<Real_t>(domain->delv_xi);
  occaCheck<Real_t>(domain->delx_eta);
  occaCheck<Real_t>(domain->delv_eta);

  //cudaDeviceSynchronize();
  //cudaCheckError();
}

static inline
void CalcMonotonicQRegionForElems(Domain *domain)
{

  const Real_t ptiny        = ((Real_t) 1.e-36) ;
  Real_t monoq_max_slope    = domain->monoq_max_slope ;
  Real_t monoq_limiter_mult = domain->monoq_limiter_mult ;

  Real_t qlc_monoq = domain->qlc_monoq;
  Real_t qqc_monoq = domain->qqc_monoq;
  Index_t elength = domain->numElem;

  Index_t dimBlock= 128;
  Index_t dimGrid = PAD_DIV(elength,dimBlock);

  size_t dims = 1;
  occa::dim inner(dimBlock);
  occa::dim outer(dimGrid);

  // CalcMonotonicQRegionForElems_kernel<<<dimGrid,dimBlock>>>
  // ( qlc_monoq,qqc_monoq,monoq_limiter_mult,monoq_max_slope,ptiny,elength,
  //   domain->matElemlist.raw(),domain->elemBC.raw(),
  //   domain->lxim.raw(),domain->lxip.raw(),
  //   domain->letam.raw(),domain->letap.raw(),
  //   domain->lzetam.raw(),domain->lzetap.raw(),
  //   domain->delv_xi->raw(),domain->delv_eta->raw(),domain->delv_zeta->raw(),
  //   domain->delx_xi->raw(),domain->delx_eta->raw(),domain->delx_zeta->raw(),
  //   domain->vdov.raw(),domain->elemMass.raw(),domain->volo.raw(),domain->vnew->raw(),
  //   domain->qq.raw(),domain->ql.raw(),
  //   domain->q.raw(),
  //   domain->qstop,
  //   domain->bad_q_h
  // );

  occa::memory bad_q = occaMalloc(sizeof(Index_t), domain->bad_q_h);

  CalcMonotonicQRegionForElems_kernel.setWorkingDims(dims, inner, outer);

  CalcMonotonicQRegionForElems_kernel
    ( qlc_monoq,qqc_monoq,
      monoq_limiter_mult,
      monoq_max_slope,
      ptiny,
      elength,
      domain->matElemlist,
      domain->elemBC,
      domain->lxim,
      domain->lxip,
      domain->letam,
      domain->letap,
      domain->lzetam,
      domain->lzetap,
      domain->delv_xi,
      domain->delv_eta,
      domain->delv_zeta,
      domain->delx_xi,
      domain->delx_eta,
      domain->delx_zeta,
      domain->vdov,
      domain->elemMass,
      domain->volo,
      domain->vnew,
      domain->qq,
      domain->ql,
      domain->q,
      domain->qstop,
      bad_q);

  bad_q.copyTo(domain->bad_q_h);

  bad_q.free();

  //cudaDeviceSynchronize();
  //cudaCheckError();
}


static inline
void ApplyMaterialPropertiesAndUpdateVolume(Domain *domain)
{
  Index_t length = domain->numElem ;

  if (length != 0) {

    Index_t dimBlock = 128;
    Index_t dimGrid = PAD_DIV(length,dimBlock);

    size_t dims = 1;
    occa::dim inner(dimBlock);
    occa::dim outer(dimGrid);

    // ApplyMaterialPropertiesAndUpdateVolume_kernel<<<dimGrid,dimBlock>>>
    //     (length,
    //      domain->refdens,
    //      domain->e_cut,
    //      domain->emin,
    //      domain->ql.raw(),
    //      domain->qq.raw(),
    //      domain->vnew->raw(),
    //      domain->v.raw(),
    //      domain->pmin,
    //      domain->p_cut,
    //      domain->q_cut,
    //      domain->eosvmin,
    //      domain->eosvmax,
    //      domain->matElemlist.raw(),
    //      domain->e.raw(),
    //      domain->delv.raw(),
    //      domain->p.raw(),
    //      domain->q.raw(),
    //      domain->ss4o3,
    //      domain->ss.raw(),
    //      domain->v_cut,
    //      domain->bad_vol_h
    //      );

    occa::memory bad_vol = occaMalloc(sizeof(Index_t), domain->bad_vol_h);

    ApplyMaterialPropertiesAndUpdateVolume_kernel.setWorkingDims(dims, inner, outer);

    ApplyMaterialPropertiesAndUpdateVolume_kernel
      (length,
       domain->refdens,
       domain->e_cut,
       domain->emin,
       domain->ql,
       domain->qq,
       domain->vnew,
       domain->v,
       domain->pmin,
       domain->p_cut,
       domain->q_cut,
       domain->eosvmin,
       domain->eosvmax,
       domain->matElemlist,
       domain->e,
       domain->delv,
       domain->p,
       domain->q,
       domain->ss4o3,
       domain->ss,
       domain->v_cut,
       bad_vol);

    bad_vol.copyTo(domain->bad_vol_h);
    bad_vol.free();

    //    occaCheck<Real_t>(domain->e);

    //cudaDeviceSynchronize();
    //cudaCheckError();
  }
}



static inline
void LagrangeElements(Domain *domain)
{

  int allElem = domain->numElem +  /* local elem */
    2*domain->sizeX*domain->sizeY ; /* plane ghosts */

  // domain->vnew = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->dxx  = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->dyy  = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->dzz  = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);

  // domain->vnew = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->dxx = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->dyy = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->dzz = occaMalloc(domain->numElem*sizeof(Real_t));

  // domain->delx_xi    = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->delx_eta   = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->delx_zeta  = Allocator< Vector_d<Real_t> >::allocate(domain->numElem);
  // domain->delv_xi    = Allocator< Vector_d<Real_t> >::allocate(allElem);
  // domain->delv_eta   = Allocator< Vector_d<Real_t> >::allocate(allElem);
  // domain->delv_zeta  = Allocator< Vector_d<Real_t> >::allocate(allElem);


  // domain->delx_xi    = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->delx_eta   = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->delx_zeta  = occaMalloc(domain->numElem*sizeof(Real_t));
  // domain->delv_xi    = occaMalloc(allElem*sizeof(Real_t));
  // domain->delv_eta   = occaMalloc(allElem*sizeof(Real_t));
  // domain->delv_zeta  = occaMalloc(allElem*sizeof(Real_t));

  /*********************************************/
  /*  Calc Kinematics and Monotic Q Gradient   */
  /*********************************************/
  occa::tic("kinematics n monotonic Q grad");
  occaCheckDomain(domain);
  CalcKinematicsAndMonotonicQGradient(domain);
  occaCheckDomain(domain);
  occa::toc("kinematics n monotonic Q grad");

  // Allocator<Vector_d<Real_t> >::free(domain->dxx,domain->numElem);
  // Allocator<Vector_d<Real_t> >::free(domain->dyy,domain->numElem);
  // Allocator<Vector_d<Real_t> >::free(domain->dzz,domain->numElem);

  // Allocator< occa::memory >::free(domain->dxx,domain->numElem);
  // Allocator< occa::memory >::free(domain->dyy,domain->numElem);
  // Allocator< occa::memory >::free(domain->dzz,domain->numElem);
  // domain->dxx.free();
  // domain->dyy.free();
  // domain->dzz.free();


  /***********************************************************/
  /* Transfer veloctiy gradients in the first order elements */
  /* problem->commElements->Transfer(CommElements::monoQ) ;  */
  /***********************************************************/

  occa::tic("monotonic Q region");
  /**********************************
   *    Calc Monotic Q Region
   **********************************/
  occaCheckDomain(domain);
  CalcMonotonicQRegionForElems(domain);
  occaCheckDomain(domain);
  occa::toc("monotonic Q region");

  // Allocator<Vector_d<Real_t> >::free(domain->delx_xi,domain->numElem);
  // Allocator<Vector_d<Real_t> >::free(domain->delx_eta,domain->numElem);
  // Allocator<Vector_d<Real_t> >::free(domain->delx_zeta,domain->numElem);
  // Allocator<Vector_d<Real_t> >::free(domain->delv_xi,allElem);
  // Allocator<Vector_d<Real_t> >::free(domain->delv_eta,allElem);
  // Allocator<Vector_d<Real_t> >::free(domain->delv_zeta,allElem);

  // Allocator<occa::memory >::free(domain->delx_xi,domain->numElem);
  // Allocator<occa::memory >::free(domain->delx_eta,domain->numElem);
  // Allocator<occa::memory >::free(domain->delx_zeta,domain->numElem);
  // Allocator<occa::memory >::free(domain->delv_xi,allElem);
  // Allocator<occa::memory >::free(domain->delv_eta,allElem);
  // Allocator<occa::memory >::free(domain->delv_zeta,allElem);

  // domain->delx_xi.free();
  // domain->delx_eta.free();
  // domain->delx_zeta.free();
  // domain->delv_xi.free();
  // domain->delv_eta.free();
  // domain->delv_zeta.free();

  occa::tic("update volume");
  occaCheckDomain(domain);
  ApplyMaterialPropertiesAndUpdateVolume(domain) ;
  occaCheckDomain(domain);
  occa::toc("update volume");
  //  Allocator<Vector_d<Real_t> >::free(domain->vnew,domain->numElem);
  // Allocator<occa::memory >::free(domain->vnew,domain->numElem);
  // domain->vnew.free();
}



static inline
void CalcTimeConstraintsForElems(Domain* domain){

  Real_t qqc = domain->qqc;
  Real_t qqc2 = ((Real_t) 64.0) * qqc * qqc ;
  Real_t dvovmax = domain->dvovmax ;

  const Index_t length = domain->numElem;

  //  const int max_dimGrid = 1024;
  const int dimBlock = 128;
  int dimGrid=std::min(max_dimGrid,PAD_DIV(length,dimBlock));

  size_t dims = 1;
  occa::dim inner(dimBlock);
  occa::dim outer(dimGrid);

  // TODO: whats this?
  // cudaFuncSetCacheConfig(CalcTimeConstraintsForElems_kernel<dimBlock>, cudaFuncCachePreferShared);

  // Vector_d<Real_t>* dev_mindtcourant= Allocator< Vector_d<Real_t> >::allocate(dimGrid);
  // Vector_d<Real_t>* dev_mindthydro  = Allocator< Vector_d<Real_t> >::allocate(dimGrid);

  occa::tic("time constraints kernel");

  // occa::memory* dev_mindtcourant= Allocator< occa::memory >::allocate(dimGrid);
  // occa::memory* dev_mindthydro  = Allocator< occa::memory >::allocate(dimGrid);
  occa::memory dev_mindtcourant = occaMalloc(dimGrid*sizeof(Real_t));
  occa::memory dev_mindthydro = occaMalloc(dimGrid*sizeof(Real_t));

  // CalcTimeConstraintsForElems_kernel<dimBlock> <<<dimGrid,dimBlock>>>
  //     (length,qqc2,dvovmax,
  //      domain->matElemlist.raw(),domain->ss.raw(),domain->vdov.raw(),domain->arealg.raw(),
  //      dev_mindtcourant->raw(),dev_mindthydro->raw());


  // TODO: careful with dimBlock template variable
  CalcTimeConstraintsForElems_kernel.setWorkingDims(dims, inner, outer);

  occaCheckDomain(domain);

  CalcTimeConstraintsForElems_kernel
    (length,qqc2,dvovmax,
     domain->matElemlist,
     domain->ss,
     domain->vdov,
     domain->arealg,
     dev_mindtcourant,
     dev_mindthydro);

  occaCheckDomain(domain);

  occa::toc("time constraints kernel", CalcTimeConstraintsForElems_kernel);

  // TODO: if dimGrid < 1024, should launch less threads
  // CalcMinDtOneBlock<max_dimGrid> <<<2,max_dimGrid, max_dimGrid*sizeof(Real_t), domain->streams[1]>>>(dev_mindthydro->raw(),dev_mindtcourant->raw(),domain->dtcourant_h,domain->dthydro_h, dimGrid);


  occa::tic("min Dt kernel");
  dims = 1;
  inner.x = max_dimGrid;
  outer.x = 2;
  CalcMinDtOneBlock.setWorkingDims(dims, inner, outer);

  occa::memory dtcourant_d = occaMalloc(sizeof(Real_t), domain->dtcourant_h);
  occa::memory dthydro_d = occaMalloc(sizeof(Real_t), domain->dthydro_h);

  occaCheckDomain(domain);

  // TODO: add occa get stream and set stream
  CalcMinDtOneBlock(dev_mindthydro,
		    dev_mindtcourant,
		    dtcourant_d,
		    dthydro_d, dimGrid);


  dtcourant_d.copyTo(domain->dtcourant_h);
  dthydro_d.copyTo(domain->dthydro_h);

  occaCheckDomain(domain);
  // cudaEventRecord(domain->time_constraint_computed,domain->streams[1]);

  // Allocator<Vector_d<Real_t> >::free(dev_mindtcourant,dimGrid);
  // Allocator<Vector_d<Real_t> >::free(dev_mindthydro,dimGrid);

  // Allocator<occa::memory >::free(dev_mindtcourant,dimGrid);
  // Allocator<occa::memory >::free(dev_mindthydro,dimGrid);

  dev_mindtcourant.free();
  dev_mindthydro.free();

  dtcourant_d.free();
  dthydro_d.free();

  occa::toc("min Dt kernel", CalcMinDtOneBlock);
}


static inline
void LagrangeLeapFrog(Domain* domain)
{

  /* calculate nodal forces, accelerations, velocities, positions, with
   * applied boundary conditions and slide surface considerations */
  occa::tic("lagrange nodal");
  LagrangeNodal(domain);
  occa::toc("lagrange nodal");

  /* calculate element quantities (i.e. velocity gradient & q), and update
   * material states */
  occa::tic("lagrange elements");
  LagrangeElements(domain);
  occa::toc("lagrange elements");

  occa::tic("time constraints");
  CalcTimeConstraintsForElems(domain);
  occa::toc("time constraints");

}

void printUsage(char* argv[])
{
  printf("Usage: \n");
  printf("Unstructured grid:  %s -u <file.lmesh> ThreadModel Platform Device \n", argv[0]) ;
  printf("Structured grid:    %s -s numEdgeElems ThreadModel Platform Device \n", argv[0]) ;
  printf("\nExamples:\n") ;
  printf("%s -s 45 OpenCL 0 1 \n", argv[0]) ;
  printf("%s -u sedov15oct.lmesh CUDA 1 1 \n", argv[0]) ;
}


#ifdef SAMI

#ifdef __cplusplus
extern "C" {
#endif
#include "silo.h"
#ifdef __cplusplus
}
#endif

#define MAX_LEN_SAMI_HEADER  10

#define SAMI_HDR_NUMBRICK     0
#define SAMI_HDR_NUMNODES     3
#define SAMI_HDR_NUMMATERIAL  4
#define SAMI_HDR_INDEX_START  6
#define SAMI_HDR_MESHDIM      7

#define MAX_ADJACENCY  14  /* must be 14 or greater */

void DumpSAMI(Domain *domain, char *name)
{
  DBfile *fp ;
  int headerLen = MAX_LEN_SAMI_HEADER ;
  int headerInfo[MAX_LEN_SAMI_HEADER];
  char varName[] = "brick_nd0";
  char coordName[] = "x";
  int version = 121 ;
  int numElem = int(domain->numElem) ;
  int numNode = int(domain->numNode) ;
  int count ;

  int *materialID ;
  int *nodeConnect ;
  double *nodeCoord ;

  if ((fp = DBCreate(name, DB_CLOBBER, DB_LOCAL,
		     NULL, DB_PDB)) == NULL)
    {
      printf("Couldn't create file %s\n", name) ;
      exit(1);
    }

  for (int i=0; i<MAX_LEN_SAMI_HEADER; ++i) {
    headerInfo[i] = 0 ;
  }
  headerInfo[SAMI_HDR_NUMBRICK]    = numElem ;
  headerInfo[SAMI_HDR_NUMNODES]    = numNode ;
  headerInfo[SAMI_HDR_NUMMATERIAL] = 1 ;
  headerInfo[SAMI_HDR_INDEX_START] = 1 ;
  headerInfo[SAMI_HDR_MESHDIM]     = 3 ;

  DBWrite(fp, "mesh_data", headerInfo, &headerLen, 1, DB_INT) ;

  count = 1 ;
  DBWrite(fp, "version", &version, &count, 1, DB_INT) ;

  nodeConnect = new int[numElem] ;

  // Vector_h<Index_t> nodelist_h = domain->nodelist;
  std::vector<Index_t> nodelist_h = domain->nodelist;

  for (Index_t i=0; i<8; ++i)
    {
      for (Index_t j=0; j<numElem; ++j) {
	nodeConnect[j] = int(nodelist_h[i*domain->padded_numElem + j]) + 1 ;
      }
      varName[8] = '0' + i;
      DBWrite(fp, varName, nodeConnect, &numElem, 1, DB_INT) ;
    }

  delete [] nodeConnect ;

  nodeCoord = new double[numNode] ;

  // Vector_h<Real_t> x_h = domain->x;
  // Vector_h<Real_t> y_h = domain->y;
  // Vector_h<Real_t> z_h = domain->z;

  std::vector<Real_t> x_h = domain->x;
  std::vector<Real_t> y_h = domain->y;
  std::vector<Real_t> z_h = domain->z;

  for (Index_t i=0; i<3; ++i)
    {
      for (Index_t j=0; j<numNode; ++j) {
	Real_t coordVal ;
	switch(i) {
	case 0: coordVal = double(x_h[j]) ; break ;
	case 1: coordVal = double(y_h[j]) ; break ;
	case 2: coordVal = double(z_h[j]) ; break ;
	}
	nodeCoord[j] = coordVal ;
      }
      coordName[0] = 'x' + i ;
      DBWrite(fp, coordName, nodeCoord, &numNode, 1, DB_DOUBLE) ;
    }

  delete [] nodeCoord ;

  materialID = new int[numElem] ;

  for (Index_t i=0; i<numElem; ++i)
    materialID[i] = 1 ;

  DBWrite(fp, "brick_material", materialID, &numElem, 1, DB_INT) ;

  delete [] materialID ;

  DBClose(fp);
}
#endif

#ifdef SAMI
void DumpDomain(Domain *domain)
{
  char meshName[64] ;
  printf("Dumping SAMI file\n");
  sprintf(meshName, "sedov_%d.sami", int(domain->cycle)) ;

  DumpSAMI(domain, meshName) ;

}
#endif

void write_solution(Domain* locDom)
{
  // Vector_h<Real_t> x_h = locDom->x;
  // Vector_h<Real_t> y_h = locDom->y;
  // Vector_h<Real_t> z_h = locDom->z;

  size_t domNodes = locDom->x.bytes()/sizeof(Real_t);

  std::vector<Real_t> x_h(domNodes);
  std::vector<Real_t> y_h(domNodes);
  std::vector<Real_t> z_h(domNodes);

  locDom->x.copyTo(&(x_h[0]));
  locDom->y.copyTo(&(y_h[0]));
  locDom->z.copyTo(&(z_h[0]));

  printf("Writing solution to file xyz.asc\n");
  std::stringstream filename;
  filename << "xyz.asc";

  FILE *fout = fopen(filename.str().c_str(),"wb");

  for (Index_t i=0; i<locDom->numNode; i++) {
    fprintf(fout,"%10d\n",i);
    fprintf(fout,"%.10f\n",x_h[i]);
    fprintf(fout,"%.10f\n",y_h[i]);
    fprintf(fout,"%.10f\n",z_h[i]);
  }
  fclose(fout);
}

int main(int argc, char *argv[])
{
  if (argc != 6) {
    printUsage(argv);
    exit( LFileError );
  }

  if (  strcmp(argv[1],"-u") != 0 && strcmp(argv[1],"-s") != 0 )
    {
      printUsage(argv);
      exit( LFileError ) ;
    }

  bool structured = ( strcmp(argv[1],"-s") == 0 );

  occa::tic("initialization");
  // cuda_init();
  occa_init(argv[3], atoi(argv[4]), atoi(argv[5]));

  occa::initTimer(occaHandle);

  printf("OCCA device and kernels are initialized \n");

  /* assume cube subdomain geometry for now */
  Index_t nx = atoi(argv[2]);

  Domain *locDom ;
  locDom = NewDomain(argv,nx,structured) ;

  // TODO: whats this?
  //  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);

  /* timestep to solution */
  int its=0;

  if (structured)
    printf("Running until t=%f, Problem size=%dx%dx%d\n",locDom->stoptime,nx,nx,nx);
  else
    printf("Running until t=%f, Problem size=%d \n",locDom->stoptime,locDom->numElem);

  occa::toc("initialization");

  // cudaEvent_t timer_start, timer_stop;
  // cudaEventCreate(&timer_start);
  // cudaEventCreate(&timer_stop);
  // cudaEventRecord( timer_start );
  double timer_start = occa::currentTime();

  occa::tic("time stepping");

  while(locDom->time_h < locDom->stoptime)
    {
      // Time increment has been moved after computation of volume forces to hide launch latencies
      //TimeIncrement(locDom) ;

      occa::tic("lagrange leap frog");
      LagrangeLeapFrog(locDom) ;
      occa::toc("lagrange leap frog");

      checkErrors(locDom,its);

#if LULESH_SHOW_PROGRESS
      printf("time = %e, dt=%e\n", double(locDom->time_h), double(locDom->deltatime_h) ) ;
#endif
      its++;
    }

  double timer_stop = occa::currentTime();

  occa::toc("time stepping");
  // float elapsed_time;
  // // cudaEventRecord( timer_stop );
  // // cudaEventSynchronize( timer_stop);
  // // cudaEventElapsedTime( &elapsed_time, timer_start, timer_stop );
  // elapsed_time*=1.e-3f;

  double elapsed_time = timer_stop - timer_start;


  occa::printTimer();


  printf("Run completed:  \n");
  printf("   Elapsed Time        =  %8.4e seconds\n",elapsed_time);
  if (structured)
    printf("   Problem size        =  %ix%ix%i \n",    nx,nx,nx);
  else
    printf("   Problem size        =  %i \n",    locDom->numElem);
  printf("   Iteration count     =  %i \n",    its);

  Real_t e_zero;

  // TODO:
  // Real_t* d_ezero_ptr = locDom->e.raw() + locDom->octantCorner;
  // cudaMemcpy(&e_zero, d_ezero_ptr, sizeof(Real_t),cudaMemcpyDeviceToHost) ;

  locDom->e.copyTo(&e_zero, sizeof(Real_t), locDom->octantCorner *sizeof(Real_t));

  printf("   Final Origin Energy =  %16.10e \n", e_zero);

  size_t free_mem, total_mem, used_mem;
  // TODO:
  // cudaMemGetInfo(&free_mem, &total_mem);
  used_mem= total_mem - free_mem;
  printf("   Used Memory         =  %8.4f Mb\n", used_mem / (1024.*1024.) );

  bool write_solution_flag=true;
  if (write_solution_flag) {
    write_solution(locDom);
  }

#ifdef SAMI
  DumpDomain(locDom) ;
#endif

  return 0 ;
}
