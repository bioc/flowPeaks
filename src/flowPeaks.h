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

extern "C"{
    void Rpack_get_flowpeaks(double *Cw,double* Cm, double* CS,int *pK, int *pp, 
		       int *Nb,double *Cpeaks, double *Cpeaks_f, double * Cpeaks_df, 
		       int* cfound,int *cid, double* ptol);
    void Rpack_get_flowpeaks2(double *Cw,double* Cm, double* CS, int *pK,int *pp, 
			       int *Nb, double* Cm_f,
			       double *Cpeaks, int *cid, double *ptol);
    void Rpack_raster_image(double *raw, int *rawid,
			    int *pn, int* pres, double*grid, 
			    int *grid_id, int*pngrid);
    void Rpack_kmeans(double *A, int *pn, int*pp, int*pK,
		      int*IC, double *C,
		      int *nc,double *S, int *Nb,
		      double*twss,double* stime);
    void Rpack_assign_kmeans(double *A, int *pn, int *pp, double *C,
			     int *pK, int *ret_IC);
    void Rpack_assign_flowPeaks(double *A, int *pn, int *pp, double *Cw,double*Cm,
		      double *CS,int *pK, int* cid, double* ftol,double *ffc,
		      int *ret_IC);
    void Rpack_relevel(int *A, int*B, int*pn, int* L_old, int*L_new, int *pk);
    void Rpack_voronoi(int *n, double *xv, double *yv, int* autobox, 
		       double*box,int *nedge, double *coordv, int *sitev);
    //assuming the values of A belonging to L_old[0..(k-1)]
    void Rpack_summarize_cluster(double *A, int *pn, int*pp, int *pK, 
				 int *IC, int*ret_nc, double* ret_C,
				 double* ret_S,double *ret_twss);
    void Rpack_kmeans_center(double *A, int*pn, int *pp, int*pK,
			 double *initC,double *C,double *err, int* iter_max,
			     double*twss, double*stime);
    void Rpack_seedplusplus(double *A, int *pn, int *pp, int *pK,
			    double *C, double* ptwss);
}
void get_flowpeaks(double *Cw,double* Cm, double* CS,int *pK, int *pp, int *Nb,
		   double *Cpeaks, double *Cpeaks_f, double * Cpeaks_df, 
		   int* cfound,int *cid, double* ptol);
void get_flowpeaks2(double *Cw,double* Cm, double* CS, int *pK,int *pp, 
		    int *Nb, double* Cm_f,
		     double *Cpeaks, int *cid, double *ptol);
void raster_image(double *raw, int *rawid,
		  int *pn, int* pres, double*grid, 
		  int *grid_id, int*pngrid);
void get_kmeans(double *A, int *pn, int*pp, int*pK, 
		int*ret_IC,double *ret_C,int*ret_nc,double *ret_S, 
		int *ret_Nb,double *re_twss,double* stime);
void get_kmeans_center(double *A, int*pn, int *pp, int*pK,
			 double *initC,double *C,double *err, int* iter_max,
		       double*twss,bool transpose=true,double *stime=NULL);
void voronoi(int *n, double *xv, double *yv, int* autobox, 
	     double*box,int *nedge, double *coordv, int *sitev);
void get_Var_local(const double *A, const int n, const int p,
		   const int* id, const int nid, double *sV);
void compute_Nb(const int *IC1, const int *IC2, const int n,
		const int K, int* Nb);
void assign_flowPeaks(double *A, int *pn, int *pp, double *Cw,double*Cm,
		      double *CS,int *pK, int* cid, double* ftol,double *ffc,
		      int *ret_IC);
