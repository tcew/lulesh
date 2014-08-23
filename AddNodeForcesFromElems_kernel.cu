extern "C" __global__
void AddNodeForcesFromElems_kernel( occaKernelInfoArg, Index_t numNode,
                                    Index_t padded_numNode,
                                    const Int_t* nodeElemCount,
                                    const Int_t* nodeElemStart,
                                    const Index_t* nodeElemCornerList,
                                    const Real_t* fx_elem,
                                    const Real_t* fy_elem,
                                    const Real_t* fz_elem,
                                    Real_t* fx_node,
                                    Real_t* fy_node,
                                    Real_t* fz_node,
                                    const Int_t num_threads)
{
    int tid=blockDim.x*blockIdx.x+threadIdx.x;
    if (tid < num_threads)
    {
      Index_t g_i = tid;
      Int_t count=nodeElemCount[g_i];
      Int_t start=nodeElemStart[g_i];
      Real_t fx,fy,fz;
      fx=fy=fz=Real_t(0.0);

      for (int j=0;j<count;j++)
      {
          Index_t pos=nodeElemCornerList[start+j]; // Uncoalesced access here
          fx += fx_elem[pos];
          fy += fy_elem[pos];
          fz += fz_elem[pos];
      }


      fx_node[g_i]=fx;
      fy_node[g_i]=fy;
      fz_node[g_i]=fz;
    }
}

/* occaKernel void */
/* AddNodeForcesFromElems_kernel(occaKernelInfoArg, */
/* 			      const Index_t occaVariable numNode, */
/* 			      const Index_t occaVariable padded_numNode, */
/* 			      const occaPointer Int_t* nodeElemCount, */
/* 			      const occaPointer Int_t* nodeElemStart, */
/* 			      const occaPointer Index_t* nodeElemCornerList, */
/* 			      const occaPointer Real_t* fx_elem, */
/* 			      const occaPointer Real_t* fy_elem, */
/* 			      const occaPointer Real_t* fz_elem, */
/* 			      occaPointer Real_t* fx_node, */
/* 			      occaPointer Real_t* fy_node, */
/* 			      occaPointer Real_t* fz_node, */
/* 			      const Int_t occaVariable num_threads){ */

/*   occaGlobalFor0{ */

/*     // int tid=blockDim.x*blockIdx.x+threadIdx.x; */
/*     const int tid = occaGlobalId0; */

/*     if (tid < num_threads){ */
/*       const Index_t g_i = tid; */
/*       const Int_t count=nodeElemCount[g_i]; */
/*       const Int_t start=nodeElemStart[g_i]; */
/*       Real_t fx,fy,fz; */
/*       fx=fy=fz=Real_t(0.0); */

/*       for (int j=0;j<count;j++){ */
/* 	const Index_t pos=nodeElemCornerList[start+j]; // Uncoalesced access here */
/* 	fx += fx_elem[pos]; */
/* 	fy += fy_elem[pos]; */
/* 	fz += fz_elem[pos]; */
/*       } */


/*       fx_node[g_i]=fx; */
/*       fy_node[g_i]=fy; */
/*       fz_node[g_i]=fz; */
/*     } */
/*   } */
/* } */
