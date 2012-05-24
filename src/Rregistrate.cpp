#include "gvector_gmatrix.h"
#include "flowPeaks.h"
#include <R_ext/Rdynload.h>
static const R_CMethodDef cMethods[]={
    {"Rpack_get_flowpeaks",(DL_FUNC)&Rpack_get_flowpeaks,12},
    {"Rpack_get_flowpeaks2",(DL_FUNC)&Rpack_get_flowpeaks2,10},
    {"Rpack_raster_image",(DL_FUNC)&Rpack_raster_image,7},
    {"Rpack_kmeans",(DL_FUNC)&Rpack_kmeans,11},
    {"Rpack_assign_kmeans",(DL_FUNC)&Rpack_assign_kmeans,6},
    {"Rpack_relevel",(DL_FUNC)&Rpack_relevel,6},
    {"Rpack_voronoi",(DL_FUNC)&Rpack_voronoi,8},
    {"Rpack_summarize_cluster",(DL_FUNC)&Rpack_summarize_cluster,9},
    {"Rpack_kmeans_center",(DL_FUNC)&Rpack_kmeans_center,10},
    {"Rpack_seedplusplus",(DL_FUNC)&Rpack_seedplusplus,6},
    {NULL,NULL,0}
};
extern "C" void R_init_flowPeaks(DllInfo *info)
{
    R_registerRoutines(info,cMethods,NULL,NULL,NULL);
}    
