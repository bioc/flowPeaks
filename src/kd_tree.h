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

#ifndef KD_TREE_H__
#define KD_TREE_H__
#include "gvector_gmatrix.h"

struct Node {
    int n;
    int begin_index;
    Node* left;
    Node* right;  //children nodes, using index of the node_data
    double wss; //with sum of square of deviations to the center mu
    double* mu,*middle,*rad; 
    mutable int ic1; //the best cluster assignment
};


class KD_Tree{
public:
    KD_Tree(const int n, const int p, const double* A);
    ~KD_Tree();
    void RunKMeans_EM(const int K, double *C, double* Cnew, int *nc, 
		      double &s, double eps, int iter_max, int*piter_eral=NULL);
    double RunKMeans_GE(const int K, double eps, int iter_max, double *ret_C, 
			int* ret_IC, int *ret_nc);
    double compute_newCenter(const int K, double *C, double*Cnew,int *nc) const; 
    void  compute_IC2(const int K, const double*C, const int *nc, 
		      int*IC2) const;
    void summarize_IC1(int *IC1) const;
    double summarize_twss(const double *C) const;
//data
    int n_, p_; //root is always at index zero
    const double* A_; //a copy of A with possible rearrange orders
    Node* noderoot_;// the two are exactly the same address.

    char* nodedata_; //store the data for mu, middle and radius and the Node
    int* index_;
    double* tmpV1_,*tmpV2_;

//datastorage
    //VAdouble A_v;
    VAchar nodedata_v;
    VAint index_v;
    VAdouble tmpV1_v;
    VAdouble tmpV2_v;

private:
    Node* BuildNodes(const double *A, const int begin_index, 
		     const int end_index, 
		     char* & nodenext);
    double compute_twss(const Node* node, const double* center) const;
    double compute_newCenter(const Node* node, const int* cand_C, const int K,
			     double *C, double *Cnew, int *nc) const;
    void compute_IC2(const Node* node,const int* cand_C, const int K, 
		     const double *C, const int *nc,int *IC2) const;
    bool ShouldBePruned(const double* box_middle, const double *box_radius,
			const double *C,
			const int best_index, const int test_index) const;
    void summarize_IC1(const Node* node, int *IC1 ) const;
    double summarize_twss(const Node* node, const double *C) const;
    void quick_transfer(const int K, double *C, int *nc, int*IC1, int*IC2, double*D, 
			double& twss, int iter_max) const;
};
double KMeans_EM(double *A, int n, int p, int K,
	       int trials=1,bool ppseed=true,
	       int *ret_IC=NULL, double*ret_C=NULL, int* ret_nc=NULL, 
	       double eps=0.0001,int iter_max=100);
#endif
