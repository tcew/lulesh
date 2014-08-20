
extern occa::device occaHandle;
extern occa::kernel CalcVolumeForceForElems_kernel_true;
extern occa::kernel CalcVolumeForceForElems_kernel_false;
extern occa::kernel AddNodeForcesFromElems_kernel;
extern occa::kernel CalcAccelerationForNodes_kernel;
extern occa::kernel ApplyAccelerationBoundaryConditionsForNodes_kernel;
extern occa::kernel CalcPositionAndVelocityForNodes_kernel;
extern occa::kernel CalcKinematicsAndMonotonicQGradient_kernel;
extern occa::kernel CalcMonotonicQRegionForElems_kernel;
extern occa::kernel ApplyMaterialPropertiesAndUpdateVolume_kernel;
extern occa::kernel CalcTimeConstraintsForElems_kernel;
extern occa::kernel CalcMinDtOneBlock;
extern occa::kernel FillKernel;

void buildLuleshKernels();

