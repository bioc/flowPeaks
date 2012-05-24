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

#include "gvector_gmatrix.h"
#include "quick_op.h"
#include "func_collect.h"
#include "flowPeaks.h"

void get_id_smin(const double *A, const int n, const int p,const double *c0,
		 const double smin,int* id, int&nid)
{//id needs to be preallocated with space of n
    nid=0;
    for(int i=0;i<n;i++){
	if(L2dist(c0,A+i*p,p)<=smin){
	    id[nid]=i;
	    nid++;
	}
    }
}

void get_Var_local(const double *A, const int n, const int p,
		   const int* id, const int nid, double *sV)
{
    //compute the mean and covariance matrix
    gvector sm(p);
    sm.set_zero();
    for(int i=0;i<nid;i++){
	const double *a=A+id[i]*p;
	doubleadd(sm.data,a,p);
    }
    doublemul(sm.data,1/(nid+0.0),p);
    
    doublecopy(sV,0,p*p);
    for(int i=0;i<nid;i++){
	const double *a=A+id[i]*p;
	for(int j=0;j<p;j++){
	    for(int l=j;l<p;l++){
		sV[j*p+l]+=(a[j]-sm.data[j])*(a[l]-sm.data[l]);
	    }
	}
    }
    if(nid>1){
	doublemul(sV,1/(nid-1.0),p*p);
    }
    //copy back to the lower triangle.
    doublecopy2lower(sV,p);
}

/*get the actual neighboring clusters, which is more restrictive than 
  the voronoi diagrams*/
void compute_Nb(const int *IC1, const int *IC2, const int n,
		const int K, int* Nb)
{     
    fill(Nb,Nb+K*K,0);
    for(int i=0;i<n;i++){
	int ic1=IC1[i];
	int ic2=IC2[i];
	Nb[ic1*K+ic2]++;
    }
}

