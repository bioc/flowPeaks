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

#include "kd_tree.h"
#include "quick_op.h"
#include "flowPeaks.h"
#include "func_collect.h"
#include <sstream>
#include <ctime>
#include <cstring>
using namespace std;
const double EPS=1e-6;
KD_Tree::KD_Tree(const int n, const int p, const double *A): n_(n), p_(p) {
    try{
	int N=2*n-1;
	//allocate all data in the same place to avoid page-in and out
	nodedata_v.resize(N*(sizeof(Node)+3*p*sizeof(double)));
	nodedata_=&nodedata_v[0];
	//takes about 6+10/p times more memory than the data itself
	// 5/p/ is for the extra info for the node
	//also need to consider 2 times of the data in R use
	// The extra for R can be avoided by using a command line tool
	index_v.resize(n);
	index_=&index_v[0];
	for (int i = 0; i < n; i++){
	    index_[i] = i;
	}
	tmpV1_v.resize(p);
	tmpV1_=&tmpV1_v[0];
	
	tmpV2_v.resize(p);
	tmpV2_=&tmpV2_v[0];
	
	noderoot_=(Node*) nodedata_;
	char* nodenow=nodedata_;
	BuildNodes(A,0,n,nodenow);
	//rearrange the data according to index_
	A_=A;
	//A_v.resize(n*p);
	//A_=&A_v[0];
	
	/*doublecopy(A,A_,n*p);
	for(int i=0;i<n;i++){
	    doublecopy(A+index_[i]*p,A_+i*p,p);
	    }*/
    }
    catch (bad_alloc& ba){
	gsl_error("failed to allocate space in constructing kd tree",__FILE__,__LINE__,GSL_ENOMEM);
    }
} 

KD_Tree::~KD_Tree(){
}
Node* KD_Tree::BuildNodes(const double *A, const int begin_index, const int end_index,
			char* & nodenext) 
{
    //nodenext is updated whenever this function is executed
    Node* node=(Node*)nodenext;
    node->mu=(double*)(nodenext+sizeof(Node));
    node->middle=(double*)(nodenext+sizeof(Node)+p_*sizeof(double));
    node->rad=(double*)(nodenext+sizeof(Node)+2*p_*sizeof(double));
    nodenext+=sizeof(Node)+3*p_*sizeof(double);

    node->n = (end_index - begin_index);
    node->begin_index = begin_index;
    int  sizeL = 0;
    { 
	// Calculate the bounding box
	//need to hide these data for recursive algorithm,
	const double *first=A+index_[begin_index]*p_;
	double* dmin=tmpV1_;
	double* dmax=tmpV2_;
	doublecopy(first,dmin,p_);
	doublecopy(first,dmax,p_);
	for (int i = begin_index+1; i < end_index; i++){
	    for (int j = 0; j < p_; j++) {
		double c = A[index_[i]*p_ + j];
		if (c<dmin[j]) dmin[j] = c;
		if (c>dmax[j]) dmax[j] = c;
	    }
	}
	double max_rad = -1;
	int split_p = -1;
	for (int j = 0; j < p_; j++) {
	    node->middle[j] = (dmin[j] + dmax[j]) / 2;
	    node->rad[j] = (dmax[j] - dmin[j]) / 2;
	    if (node->rad[j] > max_rad) {
		max_rad = node->rad[j];
		split_p = j;
	    }
	}
	
	// If the max spread is 0, make this a leaf node
	if (max_rad <EPS) {
	    node->left = node->right = NULL;
	    doublecopy(first,node->mu,p_);
	    node->wss = 0;
	    return node;
	}
	
	double split_pos = node->middle[split_p];
	int iL = begin_index, iR = end_index-1;
	while (iL <= iR) {
	    bool is_iL_good = (A[index_[iL]*p_ + split_p] < split_pos);
	    bool is_iR_good = (A[index_[iR]*p_ + split_p] >= split_pos);
	    if (!is_iL_good && !is_iR_good) {
		swap(index_[iL],index_[iR]);
		is_iL_good = is_iR_good = true;
	    }
	    if (is_iL_good) {
		iL++;
		sizeL++;
	    }
	    if (is_iR_good) {
		iR--;
	    }
	}
    }
    node->left = BuildNodes(A,begin_index, begin_index + sizeL, nodenext);
    node->right = BuildNodes(A,begin_index + sizeL, end_index, nodenext);

    /*it is probably twice as fast with this updating strategy than 
      a naive recomputations*/
    double w=(node->left)->n/(node->n+0.0);
    for(int j=0;j<p_;j++){
	node->mu[j]=node->left->mu[j]*w+node->right->mu[j]*(1-w);
    }
    node->wss=compute_twss(node->left,node->mu)+
	compute_twss(node->right,node->mu);
    return node;
}

//Note that IC1, nc, twss are always based in the old Center, 
void KD_Tree::summarize_IC1(int *IC1) const{
    summarize_IC1(noderoot_,IC1);
}
void KD_Tree::summarize_IC1(const Node* node, int *IC1 ) const 
{
    if (node->ic1!=-1) {
	for(int i=node->begin_index;i<node->begin_index+node->n;i++){
	    IC1[i]=node->ic1;
	}
	return;
    }
    if(node->left!=NULL){
	summarize_IC1(node->left,IC1);
	summarize_IC1(node->right,IC1);
    }
}


double KD_Tree::summarize_twss(const double *C) const{
    return summarize_twss(noderoot_,C);
}
double KD_Tree::summarize_twss(const Node* node, const double *C) const 
{
    if (node->ic1!=-1) {
	return compute_twss(node,C+node->ic1*p_);
    }
    if(node->left!=NULL){
	return summarize_twss(node->left,C)+
	    summarize_twss(node->right,C);
    }
    return -1;//please check this should be impossible
}



double KD_Tree::compute_newCenter(const int K, double *C, double *Cnew,
				  int *nc) const
{
    VAint cand_C(K);
    for(int i=0;i<K;i++){
	cand_C[i]=i;
	nc[i]=0;
    }
    for(int i=0;i<K*p_;i++){
	Cnew[i]=0;
    }
    double twss=compute_newCenter(noderoot_,&cand_C[0],K,C,Cnew,nc);
    for(int i=0;i<K;i++){
	if(nc[i]==0){
	    /*randomly assign a new center,
	      A better choice would be using Seed*/
	    const double* Aloc=A_+gsl_rng_uniform_int(g_rng.r,n_)*p_;
	    doublecopy(Aloc,Cnew+i*p_,p_);
	    gsl_stream_printf("Warning",__FILE__, __LINE__,
			      "Empty clusters, you need to check with your data");
	}
    }
    return twss;
}

double KD_Tree::compute_twss(const Node* node, const double* center) const {
    double s=L2dist(node->mu,center,p_);
    return node->wss + node->n*s;
}

void KD_Tree::quick_transfer(const int K, double *C, int *nc, int*IC1, int*IC2, double*D, 
			double& twss, int iter_max) const //good choice of itermax maybe 10
{
    const double *A=A_;
    int n=n_;
    int p=p_;

    //update IC1, IC2 from the currentCenter, not the newCenter
    //twss,ic1, nc are corrected assigned by the current center
    summarize_IC1(noderoot_,IC1);
    compute_IC2(K,C,nc,IC2);

    VAdouble C1(K);
    VAdouble C2(K);
    VAint C_last(K);
    for(int k=0;k<K;k++){
	C_last[k]=n-1;//to make sure D[i] is computed during the first time
	updateC12(nc[k],C1[k],C2[k]);
    };
    
    int istep=0;
    int quick_good=0;
    for(int iter=0;iter<iter_max && quick_good<n+1;iter++){
	for(int i=0;i<n&& quick_good<n+1;i++,quick_good++,istep++){
	    int k1=IC1[i];
	    int k2=IC2[i];
	    if(nc[k1]==1) continue;
	    if(C_last[k1]>=istep){ //update D[i] even when exactly n steps ago as we didn't update D[i] at that time
		D[i]=L2dist(A+i*p,C+k1*p,p)*C1[k1];
	    }
	    if(C_last[k1]<=istep && C_last[k2]<=istep) continue;
	    double R2=C2[k2]*L2dist(A+i*p,C+k2*p,p);
	    if(R2>=D[i]) continue;
	    quick_good=0; //reset the count of good data
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
    //cerr<<"    Qt_iter="<<istep/(n+0.0);
}

void KD_Tree::compute_IC2(const int K, const double*C, const int *nc, 
		      int*IC2) const
{
    Vint cand_C(K);
    for(int k=0;k<K;k++){
	cand_C[k]=k;
    }
    compute_IC2(noderoot_,&cand_C[0],K,C,nc,IC2);
}

void KD_Tree::compute_IC2(const Node* node,const int* cand_C, const int K, 
			  const double *C, const int *nc,int *IC2) const 
{
    if(node->ic1==-1){
	compute_IC2(node->left,cand_C, K,C,nc,IC2);
	compute_IC2(node->right,cand_C,K,C,nc,IC2);
	return;
    }

    int ic1=node->ic1;
    int ic2=cand_C[0];
    if(ic1==ic2) ic2=cand_C[1];

    //find the minimu, note that K>2
    if(K>2){
	double s2min=L2dist(node->middle,C+ic2*p_, p_);
	for(int i=1;i<K;i++){
	    if(cand_C[i]==ic1 || cand_C[i]==ic2) continue;
	    double s=L2dist(node->middle,C+cand_C[i]*p_, p_);
	    if (s < s2min) {
		s2min = s;
		ic2 = cand_C[i];
	    }
	}
    }
    if(node->left==NULL || K==2){
	for (int i = node->begin_index; i < node->begin_index + node->n; i++){
	    IC2[i] = ic2;
	}
	return;
    }
    int* cand_Cnew=mymalloc(K,int);//(int*)malloc(sizeof(int)*K);
    cand_Cnew[0]=ic2;
    int Knew=1;
    for (int i = 0; i < K; i++){
	if(ic1==cand_C[i] || ic2==cand_C[i]) continue;
	if (!ShouldBePruned(node->middle, node->rad, C, ic2, cand_C[i])){
	    cand_Cnew[Knew]=cand_C[i];
	    Knew++;
	}
    }
    if (Knew > 1) {
	//add back ic1 for the candidate
	cand_Cnew[Knew]=ic1;
	Knew++;
	//propaget the IC1 to the children
	node->left->ic1=ic1;
	compute_IC2(node->left,cand_Cnew,Knew,C,nc,IC2);
	node->right->ic1=ic1;
	compute_IC2(node->right,cand_Cnew,Knew,C,nc,IC2);
    }else{
	for (int i = node->begin_index; i < node->begin_index + node->n; i++){
	    IC2[i] = ic2;
	}
    }
    myfree(cand_Cnew);
    return;

}

//not using Vint for Cand_C to save the amount of saved local
//variables for recursive aglorithm
double KD_Tree::compute_newCenter(const Node* node, const int* cand_C, const int K,
				  double *C, double *Cnew, int *nc) const 
{
    double smin =L2dist(node->middle,C+cand_C[0]*p_, p_);
    int ic1=cand_C[0];
    for (int i = 1; i < K; i++) {
	double s=L2dist(node->middle,C+cand_C[i]*p_, p_);
	if (s < smin) {
	    smin = s;
	    ic1 = cand_C[i];
	}
    }
    // If this is a non-leaf node, recurse if necessary
    if (node->left != NULL) {
	int *cand_Cnew=mymalloc(K,int);//(int*)malloc(K*sizeof(int));
	cand_Cnew[0]=ic1;
	int Knew=1;

	for (int i = 0; i < K; i++){
	    if(ic1==cand_C[i]) continue;
	    if (!ShouldBePruned(node->middle, node->rad, C, ic1, cand_C[i])){
		cand_Cnew[Knew]=cand_C[i];
		Knew++;
	    }
	}
	if (Knew > 1) {
	    node->ic1=-1;
	    double result=compute_newCenter(node->left,cand_Cnew,Knew,C,Cnew,nc) +
		compute_newCenter(node->right,cand_Cnew,Knew,C,Cnew,nc);
	    myfree(cand_Cnew);
	    return result;
	}else{
	    myfree(cand_Cnew);
	}
    }
    
    //update Cnew
    node->ic1=ic1;
    double w=node->n/(nc[ic1]+node->n+0.0);
    doubleweightsum(Cnew+ic1*p_,node->mu,p_,w);
    /*for(int j=0;j<p_;j++){
	Cnew[ic1*p_+j]+=(node->mu[j]-Cnew[ic1*p_+j])*w;
	}*/
    nc[ic1] += node->n;
    return compute_twss(node, C + ic1*p_);
}
bool KD_Tree::ShouldBePruned(const double* box_middle, 
			     const double *box_radius, const double *C,
			     const int best_index, const int test_index) const 
{
    const double *best=C+best_index*p_;
    const double *test =C+test_index*p_;
    double s=0;
    for (int i = 0; i < p_; i++) {
	double dif=best[i]-test[i];
	double sum=best[i]+test[i];
	int sgn=dif<0?1:-1;
	s+=dif*(2*(box_middle[i]+sgn*box_radius[i])-sum);
	/*compute 2*(x-y) \dot ( z-(x+y)/2)*/
    }
    return (s>=0);
}

double KD_Tree::RunKMeans_GE(const int K,double eps,int iter_max, double *ret_C, int *ret_IC, int *ret_nc)
{
    int p=p_;
    int n=n_;
    gmatrix C(K,p);
    gmatrix Cnew(K,p);
    VAint nc(K);
    VAint IC1(n);
    VAint IC2(n);
    VAdouble D(n);

    double s;
    s=SeedPlusPlus(A_,n_,p_,K,C.data);
    RunKMeans_EM(K,C.data,Cnew.data,&nc[0],s,eps,100);
    summarize_IC1(&IC1[0]);
    compute_IC2(K,C.data,&nc[0],&IC2[0]);
    VAdouble C1(K);
    VAdouble C2(K);
    VAint C_last(K);
    for(int k=0;k<K;k++){
	C_last[k]=n-1;//to make sure D[i] is computed during the first time
	updateC12(nc[k],C1[k],C2[k]);
    };
    int iter_real;
    //cerr<<"\n";
    Kmeans_HW_once(A_,n,p,K,
		   C.data,&IC1[0],&IC2[0],&D[0],&nc[0],s,eps,3,&iter_real);
    
    if(ret_C){
	doublecopy(C.data,ret_C,K*p);
    }
    if(ret_nc){
	copy(&nc[0],&nc[0]+K,ret_nc);
    }
    if(ret_IC){
	summarize_IC1(ret_IC);
    }
    return s;
}


void KD_Tree::RunKMeans_EM(const int K, double *C, double* Cnew, int *nc, 
			   double& s,double eps,int iter_max, int* piter_real) 
{
    double olds=s;
    int iter_real=0;
    for(int i=0;i<iter_max; i++){
	s=compute_newCenter(K,C,Cnew,nc);
	if((i>0 && olds-s<eps*s) || (i==iter_max-1))//the first time this s is the same as the olds 
	{
	    iter_real=i+1;//to get the correct number
	    doublecopy(Cnew,C,K*p_);
	    break;
	}
	doublecopy(Cnew,C,K*p_);//update C
	olds=s;
    }
    if(piter_real!=NULL){
	*piter_real=iter_real;
    }
    //update s, since ic1, and Cnew have been updated.
    //s=summarize_twss(C); this seems useless as ic1 needs to be based Cnew, not //oldC
}

double KMeans_EM(double *A, int n, int p, int K,
		 int trials, bool ppseed, int *ret_IC1, double*ret_C, 
		 int* ret_nc, double eps,int iter_max)
{
    gmatrix C(K,p);
    gmatrix Cnew(K,p);
    VAint nc(K);
    KD_Tree tree(n,p,A);
    double smin=GSL_POSINF,s=GSL_POSINF;
	    
    for(int i=0;i<trials;i++){
	//determine the seed
	if(ppseed){
	    s=SeedPlusPlus(A,n,p,K,C.data);
	    ostringstream os;
	    os<<"        step 0, set the intial seeds, tot.wss="<<s;
	    string str=os.str();
	    gsl_stream_printf("","",-1,str.c_str());
	}else{
	    VAint permuv(K);
	    int *permu=&permuv[0];
	    sample_nK(n,K,permu);
	    for(int j=0;j<K;j++){
		doublecopy(A+permu[j]*p,C.data + j*p,p);
	    }
	}
	tree.RunKMeans_EM(K,C.data,Cnew.data,&nc[0],s,eps,iter_max);
	//updates
	if(s<smin){
	    smin=s;
	    if(ret_IC1){
		tree.summarize_IC1(ret_IC1);
	    }
	    if(ret_C){
		copy(C.data,C.data+K*p,ret_C); //note that s, IC1 are for the current C
	    }
	    if(ret_nc){
		copy(&nc[0],&nc[0]+K,ret_nc);
	    }
	}
	//cerr<<"\n Finished random trial number "<<i+1<<" with smin="<<smin;
	//cerr<<" with the current s="<<s<<"\n";
    }
    return s;
}

void collect_ic1_ic2(const int n, const int *ic1, const int *ic2, const int k, 
		     int *id, int &nid)
{
    nid=0;
    for(int i=0;i<n;i++){
	if(ic1[i]==k||ic2[i]==k){
	    id[nid]=i;
	    nid++;
	}
    }
}

void transposeA(double *A, int n, int p)
{
 
    VAdouble Acv(n*p);
    double *Ac=&Acv[0];
    memcpy(Ac,A,n*p*sizeof(double));
    for(int i=0;i<n;i++){
	for(int j=0;j<p;j++){
	    A[i*p+j]=Ac[j*n+i];
	}
    }
}
void get_kmeans_center(double *A, int*pn, int *pp, int*pK,
		       double *initC,double *C,double *eps, int* iter_max,
		       double*twss,bool transpose, double *stime)
{
    //for this version, A is tranposed in C function, and .C with dup=true
    //no danger to modofi A
    int n=*pn;
    int p=*pp;
    int K=*pK;
    if(transpose) transposeA(A,n,p);//A has been modified
    if(stime==NULL){
	getRunningTime(true);
    }else{
	getRunningTime(true,*stime);
    }
    VAint nc(K);
    int iter_real;
    double s;
    {
	KD_Tree tree(n,p,A);
	tree.RunKMeans_EM(K,initC,C,&nc[0],s,eps[0],iter_max[0],&iter_real);
    }
    char msgstr[1000];
    sprintf(msgstr,"Finished kd-tree at %d iterations with tot.wss=%.5e at %5.2f seconds\n",iter_real,s,getRunningTime());
    gsl_stream_printf("","",-1,msgstr);
    //note that initC and C are identical
    //now works with Hartigan and Wong
    s=KMeans_HW_plain(A,n,p,K,C,
		      NULL,C,NULL,
		      eps[1],iter_max[1],iter_real,NULL);
    sprintf(msgstr,"Finished Hartigan_wong at %d iterations with tot.wss=%.5e at %5.2f seconds\n",iter_real,s,getRunningTime());
    gsl_stream_printf("","",-1,msgstr);

    *twss=s;
    if(stime!=NULL){
	*stime=getRunningTime();
    }
}

void get_kmeans(double *A, int *pn, int*pp, int*pK, 
		int* ret_IC, double *ret_C, int *ret_nc,
		double *ret_S, int *ret_Nb,double *ret_twss, double*stime)
{
    int n=*pn;
    int p=*pp;
    int K=*pK;
    if(stime==NULL){
	getRunningTime(true);
    }else{
	getRunningTime(true,*stime);
    }
    double eps=min(0.01/K,0.0001);
    *ret_twss=KMeans_EM(A,n,p,K,1,true,NULL,ret_C,NULL,eps,100);
    ostringstream os;
    os<<"        step 1, do the rough EM, tot.wss="<<*ret_twss<<
	" at "<<getRunningTime()<<" sec";
    string str=os.str();
    gsl_stream_printf("","",-1,str.c_str());
    int iter_real=10;
    VAint IC2(n);
    eps=eps*10; //get a quick convergences
    *ret_twss=KMeans_HW_plain(A,n,p,K,ret_C,
			      ret_IC,ret_C,ret_nc,
			      eps,10,iter_real,&IC2[0]);
    compute_Nb(ret_IC,&IC2[0],n,K,ret_Nb);
    os.str("");
    os<<"        step 2, do the fine transfer of Hartigan-Wong Algorithm\n"<<
       "                 tot.wss="<<*ret_twss<<" at "<<getRunningTime()<<" sec";
    str=os.str();
    gsl_stream_printf("","",-1,str.c_str());
    //work on the computing the covariance matrix
    if(ret_S!=NULL){
	get_summarize(A,pn,pp,pK,ret_IC,ret_nc,ret_C,ret_S,NULL,true);
	//delegate the adjustment them into R
	/*VAint id(n);
	int nid;
	for(int k=0;k<K;k++){
	    //compute the local variances for the data that has k 
	    //as the first or second best choice
	    gmatrix sV(p,p);
	    collect_ic1_ic2(n,ret_IC,&IC2[0],k,&id[0],nid);
	    get_Var_local(A,n,p,&id[0],nid,sV.data);
	    double wV=K/(nid+K+0.0);
	    doubleweightsum(sV.data,S0,p*p,wV);
	    double* s=ret_S+k*p*p;
	    
	    double ws=K/(ret_nc[k]+K+0.0);
	    doublemul(s,1.5,p*p);
	    doubleweightsum(s,S0,p*p,ws); //only globally
	}*/
    }
    if(stime!=NULL){
	*stime=getRunningTime();
    }
}
