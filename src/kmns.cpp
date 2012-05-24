/* 
* The author of this software is Yongchao Ge.
* Permission to use, copy, modify, and distribute this software for any
* purpose without fee is hereby granted, provided that this entire notice
* is included in all copies of any software which is or includes a copy
* or modification of this software and in all copies of the supporting
* documentation for such software.
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, THE AUTHOR  DOES NOT MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
* OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*/

#include "gvector_gmatrix.h"
#include "quick_op.h"
#include "func_collect.h"
//#include "kd_tree.h"
//#include "flowPeaks.h"
using namespace std;
/*A rewrite of the Hartigan-Wong's Kmeans algorithm in C++ so that I can have 
  the flexibility to change the behavior---Yongchao Ge 11/05/2010*/
void Optim_Transfer(const double *A, const int n, const int p, const int K,
		    double *C, int *IC1, int *IC2, double* D, int *nc, double& twss,
		    double* C1, double* C2,int* isCqtran, int* C_last, int* live, int &Optim_good);
void Quick_Transfer(const double *A, const int n, const int p, const int K,
		    double *C, int *IC1, int *IC2, double* D, int *nc, double& twss,
		    double* C1, double* C2, int* isCqtran, int* C_last, int &Optim_good);

void Kmeans_HW_init(const double *A, const int n, const int p, const int K,
		    double *C, int *IC1, int *IC2, double *D,int *nc,
		      double &twss)
{
    //determine the memebership
    for(int i=0;i<n;i++){
	get_IC1_IC2(A+i*p,p,K,C,&IC1[i],&IC2[i]);
    }
    //recompute the center
    twss=summarize_IC1(A,n,p,K,IC1,nc,C,D);
}

//Input, A,m,n,K, and C, IC1, IC2, D, nc,twss, these can be computed from another aloorithm based on k-d tree
//Output, C, IC1, IC2, D, nc, twss are updated
//options, eps, and iter_max             
void  Kmeans_HW_once(const double *A, const int n, const int p, const int K,
		  double *C, int *IC1, int *IC2, double* D, int *nc, double& twss,
		  const double eps, const int iter_max, int* piter_real)
{
    VAdouble C1(K),C2(K);
    for(int k=0;k<K;k++){
	updateC12(nc[k],C1[k],C2[k]);
    }
    VAint isCqtran(1,K);//is it updated in quick transfer
    VAint  C_last(0,K);/*the last time it is changed, for Quick_Transfer, it is n+step,
		       for Optim_Transfer, it is step, where at step, the 
		       cluster center has been updated */
    
 
    //starting on the iterations
    VAint live(K);
    int Optim_good=0;
    int iter;
    for(iter=0;iter<iter_max;iter++){
	//cerr<<iter<<":"<<twss<<"\t";
	double oldtwss=twss;
	Optim_Transfer(A,n,p,K,
		       C,IC1,IC2,D,nc,twss,
		       &C1[0],&C2[0],&isCqtran[0],&C_last[0],&live[0],Optim_good);
	if(Optim_good==n) break;
	Quick_Transfer(A,n,p,K,
		       C,IC1,IC2,D,nc,twss,
		       &C1[0],&C2[0],&isCqtran[0],&C_last[0],Optim_good);
	if(K==2) break;
	if(oldtwss-twss<eps*twss) break;
	
	C_last=-1; //reset so that we can check whether it is changed
    }
    if(iter==iter_max){
	*piter_real=iter_max;
    }else{
	*piter_real=iter+1;
    }
}

void Optim_Transfer(const double *A, const int n, const int p, const int K,
		    double *C, int *IC1, int *IC2, double* D, int *nc, double& twss,
		    double* C1, double* C2, int* isCqtran, int* C_last, int* live, int &Optim_good)
{
    //initialize the liveset
    for(int k=0;k<K;k++){
	if (isCqtran[k]==1) {
	    live[k]=n;
	}
    }

    for(int i=0;i<n;i++){
	Optim_good++;
	int k1=IC1[i];
	int k2=IC2[i];
	if(nc[k1]==1) continue; //no use to swap a data from a singleton
	
	if(C_last[k1]!=-1){
	    //update D[i] only when the cluster has been updated
	    D[i]=L2dist(A+i*p,C+k1*p,p)*C1[k1];
	}
	//find the cluster with minimum R2
	double R2=L2dist(A+i*p,C+k2*p,p)*C2[k2];
	for(int k=0;k<K;k++){
	    if(k!=k1 && k!=k2 && (i<live[k1] || i<live[k])){
		double R2k=C2[k]*L2dist(A+i*p,C+k*p,p);
		if(R2k<R2){
		    R2=R2k;
		    k2=k;
		}
	    }
	}
	
	if(R2>=D[i]){
	    IC2[i]=k2; //no swap of the data, but needs to update IC2
	}else{//swap i from k1 to k2
	    //update the center, etc.
	    Optim_good=0;
	    live[k1]=n+i;
	    live[k2]=n+i;
	    C_last[k1]=i;
	    C_last[k2]=i;
	    doubleweightsum(C+k1*p,A+i*p,p,-1/(nc[k1]-1.0));
	    doubleweightsum(C+k2*p,A+i*p,p,1/(nc[k2]+1.0));
	    twss+=R2-D[i];
	    nc[k1]--;
	    nc[k2]++;
	    updateC12(nc[k1],C1[k1],C2[k1]);
	    updateC12(nc[k2],C1[k2],C2[k2]);
	    IC1[i]=k2;
	    IC2[i]=k1;
	}
    }
    
    if(Optim_good==n)return;
    live-=n; //to be used for the next call of Optim_Transfer
}
void Quick_Transfer(const double *A, const int n, const int p, const int K,
		    double *C, int *IC1, int *IC2, double* D, int *nc, double& twss,
		    double* C1, double* C2,int* isCqtran, int* C_last, int &Optim_good)
{

    for(int k=0;k<K;k++){
	isCqtran[k]=0;
    }
    int istep=0;
    int Quick_good=0;
    while(Quick_good<n+1){
	for(int i=0;i<n && Quick_good<n+1; i++,istep++,Quick_good++){
	    int k1=IC1[i];
	    int k2=IC2[i];
	    if(nc[k1]==1) continue;
	    if(istep<=C_last[k1]){//needs to recompute D[i] for at most n steps before
		D[i]=L2dist(A+i*p,C+k1*p,p)*C1[k1];
	    }
	    if(istep>=C_last[k1] && istep>=C_last[k2]) continue;
	    double R2=C2[k2]*L2dist(A+i*p,C+k2*p,p);
	    if(R2<D[i]){
		Quick_good=0;
		Optim_good=0;
		isCqtran[k1]=1;
		isCqtran[k2]=1;
		C_last[k1]=n+istep; 
		C_last[k2]=n+istep;
		doubleweightsum(C+k1*p,A+i*p,p,-1/(nc[k1]-1.0));
		doubleweightsum(C+k2*p,A+i*p,p,1/(nc[k2]+1.0));
		twss+=R2-D[i];
		nc[k1]--;
		nc[k2]++;
		updateC12(nc[k1],C1[k1],C2[k1]);
		updateC12(nc[k2],C1[k2],C2[k2]);
		IC1[i]=k2;
		IC2[i]=k1;
	    }
	}
    }
    //cerr<<"\t"<<istep/(n+0.0)<<"\n";
}

double KMeans_HW_plain(double *A, int n, int p, int K, double* init_C,
		       int *ret_IC, double*ret_C, int *ret_nc,
		       double eps,int iter_max, int &iter_real, int *ret_IC2)
{
    gmatrix C(K,p);
    VAint nc(K);
    VAint IC1(n), IC2(n);
    VAdouble D(n);
    double s;
    if(init_C){
	doublecopy(init_C,C.data,K*p);
    }else{
	/*VAint permu(K);
	sample_nK(n,K,permu);
	for(int j=0;j<K;j++){
	    doublecopy(A+permu[j]*p,C.data + j*p,p);
	    }*/
	SeedPlusPlus(A,n,p,K,C.data);
    }

    Kmeans_HW_init(A,n,p,K,C.data,&IC1[0],&IC2[0],&D[0],&nc[0],s);
    Kmeans_HW_once(A, n, p, K, 
		   C.data, &IC1[0], &IC2[0],&D[0], &nc[0], s,
		   eps,iter_max, &iter_real);
    if(ret_IC){
	copy(&IC1[0],&IC1[0]+n,ret_IC);
    }
    if(ret_IC2){
	copy(&IC2[0],&IC2[0]+n,ret_IC2);
    }
    if(ret_C){
	doublecopy(C.data,ret_C,K*p);
    }
    if(ret_nc){
	copy(&nc[0],&nc[0]+K,ret_nc);
    }
    return s;
}


