occaKernel void Fill_kernel(occaKernelInfoArg,
			    const Index_t occaVariable n,
			    const Real_t occaVariable alpha,
			    occaPointer Real_t *x){

  occaGlobalFor0{

    int i = occaGlobalId0;

    while(i < n){
      x[i] = alpha;

      i += occaGlobalDim0;

    }
  }
}
