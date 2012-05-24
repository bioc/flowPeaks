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

#ifndef FUNC_COLLECT__
#define FUNC_COLLECT__
void get_summarize(double *A, int *pn, int*pp, int *pK, int *IC,
		   int*ret_nc, double* ret_C, double* ret_S, double *ret_twss,bool computeSonly);
double summarize_IC1(const double *A, const int n, const int p, const int K,
		     const int* IC1, 
		     int *nc, double*C, double *D);
double get_IC1_IC2(const double *a, const int p,
		 const int K, const double *C,
		 int* pIC1, int*pIC2);
int get_IC(const double *a, const int p, const int K, const double *C, double* pd=NULL);
void doublecopy2lower(double *A, const int p);

//for kmeans
double SeedPlusPlus(const double *A, const int n, const int p, const int K, double *C);
double KMeans_HW_plain(double *A, int n, int p, int K, double* init_C,
		       int *ret_IC, double*ret_C, int *ret_nc,
		       double eps,int iter_max, int &iter_real, int*ret_IC2=NULL);
void  Kmeans_HW_once(const double *A, const int n, const int p, const int K,
		     double *C, int *IC1, int *IC2, double* D, int *nc, 
		     double& twss,const double eps, const int iter_max, 
		     int* piter_real);
void sample_nK(int n, int K, int* p);
void read_file(const char *filename, int* plength, Vchar& buffv);
void headline(const char *filename, const char sep, Vstring& vline,int skip=0);
void scanfast(const char *filename, const int skiprows,const int skipcols,
	      const char sep,
	      vector<double>& xv, int* pnrows) ;
double getRunningTime(bool init=false,double time0_=0.0);
void doubletranspose(const double *A, const int n, const int p,
		     double *B);
void adjustS(double *S, double *x,
	     double *w,double h, double h0, 
	     int n, int p,int K);
//if init=false, time0_ is useless
#endif
