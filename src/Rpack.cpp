/* 
* The author of this software is Yongchao Ge.
* Permission to use, copy, modify, and distribute this software for any
* purpose without fee is hereby granted, provided that this entire notice
* is included in all copies of any software which is or includes a copy
* or modification of this software and in all copies of the supporting
* documentation for such software.
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, THE AUTHOR DOES NOT MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
* OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*/


//We need to redefine the printing and error handling differently
#include "gvector_gmatrix.h"
#include "flowPeaks.h"
#include "func_collect.h"
#include "kd_tree.h"
#include "VoronoiDiagramGenerator.h"
#include <map>
#include <R.h>
class C_Call_error{};
void Rpack_stream_handler(const char * label, const char * file,
			  int line, const char * reason)
{
    if(line==-1){
	Rprintf("%s\n",reason); //no need to print extra info when its -1
    }else{
	Rprintf("%s:%d: %s: %s\n", file, line, label, reason);
    }
}

void Rpack_error_handler(const char *reason, const char* file, 
		   int line, int gsl_errno)
{
    Rpack_stream_handler("Error",file,line,reason);
    throw C_Call_error();
}
class Rpack{
public:
    gsl_error_handler_t* old_error_handler;
    gsl_stream_handler_t* old_stream_handler;
    Rpack(){
	old_error_handler=gsl_set_error_handler(&Rpack_error_handler);
	old_stream_handler=gsl_set_stream_handler(&Rpack_stream_handler);
	//setting to the same seed for every times it calls this library
	g_rng.set_seed(0);
    }
    ~Rpack(){
	gsl_set_error_handler(old_error_handler);
	gsl_set_stream_handler(old_stream_handler);
    }
};
void Rpack_get_flowpeaks(double *Cw,double* Cm, double* CS,int *pK, int *pp, 
			 int *Nb,double *Cpeaks, double *Cpeaks_f, double * Cpeaks_df, 
			 int* cfound,int *cid, double* ptol)
{
    Rpack rpack;

    try{
	get_flowpeaks(Cw,Cm,CS,pK,pp,Nb,Cpeaks,Cpeaks_f,Cpeaks_df,cfound,cid, ptol);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}
void Rpack_get_flowpeaks2(double *Cw,double* Cm, double* CS, int *pK,int *pp, 
			  int *Nb,double* Cm_f,
			   double *Cpeaks, int *cid, double *ptol)
{
    Rpack rpack;
    try{
	get_flowpeaks2(Cw,Cm,CS,pK,pp,Nb,Cm_f,Cpeaks,cid, ptol);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}

void Rpack_kmeans(double *A, int *pn, int*pp, int*pK, 
		  int*IC, double *C,
		  int *nc,double *S, int *Nb, double*twss,double* stime)
{
    Rpack rpack;
    try{
	get_kmeans(A,pn,pp,pK,IC,C,nc,S,Nb,twss,stime);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}

void Rpack_kmeans_center(double *A, int*pn, int *pp, int*pK,
			 double *initC,double *C,double *err, int* iter_max,
			 double*twss,double* stime)
{
    Rpack rpack;
    try{
	get_kmeans_center(A,pn,pp,pK,initC,C,err,iter_max,twss,true,stime);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}
void Rpack_raster_image(double *raw, int *rawid,
		      int *pn, int* pres, double*grid, 
		      int *grid_id, int*pngrid)
{
    Rpack rpack;
    try{
	raster_image(raw, rawid,
		     pn, pres, grid, 
		     grid_id, pngrid);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}

void Rpack_voronoi(int *n, double *xv, double *yv, int* autobox, 
		   double*box,int *nedge, double *coordv, int *sitev)
{
    Rpack rpack;
    try{
	voronoi(n,xv,yv,autobox,box,nedge,coordv,sitev);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }

}

void Rpack_relevel(int *A, int*B, int*pn, int *L_old, int *L_new, int *pk)
{
    int n=*pn;
    int k=*pk;
    map<int,int> levels;
    for(int i=0;i<k;i++){
	levels[L_old[i]]=L_new[i];
	if(levels.size()<(i+1)){
	    error("There are duplicated values in your old level settings\n");
	}
    }
    map<int,int>::iterator it;
    for(int i=0;i<n;i++){
	it=levels.find(A[i]);
	if(it==levels.end()){
	    error("The data does not belong to the old levels \n");
	}
	B[i]=(*it).second;
    }
}

void Rpack_summarize_cluster(double *A, int *pn, int*pp, int *pK, 
			     int *IC, int*ret_nc, double* ret_C, 
			     double* ret_S,double *ret_twss)
{
    Rpack rpack;
    try{
	get_summarize(A, pn, pp, pK, IC,
		      ret_nc, ret_C, ret_S, ret_twss,false);
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}
void Rpack_assign_kmeans(double *A, int *pn, int *pp, double *C,
			 int *pK, int *ret_IC)
{
    Rpack rpack;
    try{
	int n=*pn;
	int K=*pK;
	int p=*pp;
	for(int i=0;i<n;i++){
	    ret_IC[i]=get_IC(A+i*p,p,K,C);
	}
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}

//allow outliers detections.
void Rpack_assign_flowPeaks(double *A, int *pn, int *pp, double *Cw,double*Cm,
			    double *CS,int *pK, int* cid, double* ftol,
			    double *ffc,
			    int *ret_IC)
//if ret_IC=-1, means NA
{
    Rpack rpack;
    try{
	assign_flowPeaks(A,pn,pp,Cw,Cm,
			 CS,pK, cid, ftol,ffc,
			 ret_IC);	
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}

void Rpack_seedplusplus(double *A, int *pn, int *pp, int *pK,
			double *C, double* ptwss)
{
    Rpack rpack;
    try{
	double time_begin=getRunningTime(true);
	*ptwss=SeedPlusPlus(A,*pn,*pp,*pK,C);
	Rprintf("The tot.twss with the initial seeds is %.5E after %.2f sec\n",
		*ptwss,getRunningTime());
    }
    catch (...){
	error("The underlying C++ function has produced errors\n");
    }
}
