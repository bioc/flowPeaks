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
#include "flowPeaks.h"
#include "func_collect.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <algorithm>
#include <set>
//#include <fstream>, it caused conflict for Mac computers.
using namespace std;

typedef struct
{
  int iter;
  double step;
  double max_step;
  double tol;
  gsl_vector *x1;
  gsl_vector *dx1;
  gsl_vector *x2;
  double pnorm;
  gsl_vector *p;
  double g0norm;
  gsl_vector *g0;
}
conjugate_fr_state_t;

typedef struct
{
  double step;
  double max_step;
  double tol;
  gsl_vector *x1;
  gsl_vector *g1;
}
steepest_descent_state_t;


class GMM{ //Gaussian mixutre models
public:
    vector<gvector> Cm;
    gmatrix Cm_mat; //the same as Cm, but in matrix format.
    vector<gmatrix> CS; 
    VAdouble  Csigma_min;
    gvector Cw; //the mean, Variance matrix, and the weight of each group
    size_t getK() const{return Cw.size;}; //the number of components
    size_t getp()const {return Cm[0].size;}//the dimension of the data
    double getCCmin() const{return CCmin;}
    const gvector& getCC() const{ return CC;};
    double my_f(const gsl_vector *v);
    double pseudo_maxf();
    int ck(const gsl_vector& v) const; //the best cluster
    void my_df(const gsl_vector *v,gsl_vector *df);
    double my_fdf(const gsl_vector *v,gsl_vector *df);
    void linedens_change_xy(const gvector &x, const gvector &y, const int n,
			    VAdouble& A,VAdouble& B, VAdouble &C);
    GMM(double *Cw,double* Cm, double* CS, int *pK, int *pp){
	Init(Cw,Cm,CS,pK,pp);
    }
    void Init(double *Cw,double* Cm, double* CS, int *pK, int *pp);
    GMM(){};
	
protected:
    vector<gmatrix> CSinv;
    vector<gmatrix> CShalfinv;
    gvector CC;  /*for the convience of computing the 
		    density and density derive functions*/
    double CCmin;
    gvector vc;
    gvector vd;
    gvector ve;/*temporary objects in computing the density and 
		   density derive functions*/
    void resize(int K, int p){
	Cm.resize(K);
	Cm_mat.resize(K,p);
	CS.resize(K);
	Csigma_min.resize(K);
	Cw.resize(K);
	CSinv.resize(K);
	CShalfinv.resize(K);
	CC.resize(K);
	vc.resize(p);
	vd.resize(p);
	ve.resize(p);
	for(int i=0;i<K;i++){
	    Cm[i].resize(p);
	    CS[i].resize(p,p);
	    CSinv[i].resize(p,p);
	    CShalfinv[i].resize(p,p);
	}
    }
};
void GMM::Init(double *pCw, double* pCm, double* pCS, int *pK, int *pp)
{
    int K=*pK;
    int p=*pp;
    resize(K,p);
    copy(pCm,pCm+K*p,Cm_mat.data);
    for(int i=0;i<K;i++){
	Cm[i].memcpy(gvector_view(pCm+i*p,p));
	CS[i].memcpy(gmatrix_view(pCS+i*(p*p),p,p));
	gmatrix &Smat=CS[i];
	double smin=sqrt(Smat(0,0));
	for(int j=1;j<p;j++){
	    double s=sqrt(Smat(j,j));
	    if(s<smin){
		s=smin;
	    }
	}
	Csigma_min[i]=smin;
	Cw[i]=pCw[i];
	CSinv[i].inverse(CS[i]);
	CShalfinv[i].half(CSinv[i]);
	CC[i]=CShalfinv[i].log_det()-p*log(sqrt(2*M_PI))+
	    log(Cw[i]);
    }
    CCmin=CC.min();
    CC-=CCmin;
}
int GMM::ck(const gsl_vector& v) const{
    return get_IC(v.data,getp(),getK(),Cm_mat.data);
}

double GMM::my_f(const gsl_vector *v){
    double s=0;
    int K=getK();
    for(int i=0;i<K;i++){
	vc.memcpy_fast(v->data);
	vc.sub_fast(Cm[i].data);
	vd.mv_fast(CShalfinv[i].data,vc.data,vc.size);
	s+=exp(CC[i]-vd.dot_fast(vd.data)/2.0);
    }
    return -s;//the negative of the log density function.
}
double GMM::pseudo_maxf()
{
    int K=getK();
    double s=GSL_NEGINF;
    for(int k=0;k<K;k++){
	double sk=-my_f(&Cm[k]);
	if(sk>s) s=sk;
    }
    return s;
}
void GMM::my_df(const gsl_vector *v,gsl_vector *df)
{
    int K=getK();
    doublecopy(df->data,0.0,df->size);
    for(int i=0;i<K;i++){
	vc.memcpy_fast(v->data);
	vc.sub_fast(Cm[i].data);
	vd.mv_fast(CShalfinv[i].data,vc.data,vc.size);
	double f=exp(CC[i]-vd.dot_fast(vd.data)/2.0);
	vd.mv_fast(CSinv[i].data,vc.data,vc.size);
	vd.mul_fast(f);
	doubleadd(df->data,vd.data,vd.size);
    }
}


double GMM::my_fdf(const gsl_vector *v,gsl_vector *df)
{
    int K=getK();
    //gsl_vector_set_zero(df);
    doublecopy(df->data,0.0,df->size);
    double f=0;
    for(int i=0;i<K;i++){
	/*vc.memcpy(*v);
	vc-=Cm[i];
	vd.mv(CShalfinv[i],vc);
	double fi=exp(CC[i]-vd.dot(vd)/2.0);
	f+=fi;
	vd.mv(CSinv[i],vc);
	vd*=fi;
	gsl_vector_add(df,&vd);*/
	vc.memcpy_fast(v->data);
	vc.sub_fast(Cm[i].data);
	vd.mv_fast(CShalfinv[i].data,vc.data,vc.size);
	double fi=exp(CC[i]-vd.dot_fast(vd.data)/2.0);
	f+=fi;
	vd.mv_fast(CSinv[i].data,vc.data,vc.size);
	vd.mul_fast(fi);
	doubleadd(df->data,vd.data,vd.size);
    }
    return -f;
}



double my_f(const gsl_vector *v, void *params)
{
    GMM *gmm = (GMM *)params;
    return gmm->my_f(v);
}
/* The gradient of f*/
void my_df(const gsl_vector *v, void *params,
	    gsl_vector *df)
{
    GMM *gmm = (GMM *)params;
    gmm->my_df(v,df);
}
/* Compute both f and df together. could be improved
   for the moment, lazy function calls*/
void my_fdf (const gsl_vector *x, void *params,
	     double *f, gsl_vector *df)
{
    //*f = my_f(x, params);
    //my_df(x, params, df); save 50% of time
    GMM *gmm = (GMM *)params;
    *f=gmm->my_fdf(x,df);
}

//speed up the computations for linedev
class LineDens{
public:
    LineDens(const int K){
	A.resize(K);
	B.resize(K);
	C.resize(K);
    };
    ~LineDens(){};
    void change_xy(const gvector & x, const gvector& y, GMM& gmm, const int n){
	gmm.linedens_change_xy(x,y,n,A,B,C);}
    double operator()(const int i, GMM& gmm);
    
private:
    VAdouble A;
    VAdouble B;
    VAdouble C;
};
void GMM::linedens_change_xy(const gvector &x, const gvector &y, const int n, VAdouble& A,
		  VAdouble& B, VAdouble &C)
{
    int K=getK();
    gvector& step=ve;
    step.memcpy(y);
    step-=x;
    step*=1.0/n;//step=(y-x)/n
    for(int i=0;i<K;i++){
	vc.memcpy(x);
	vc-=Cm[i];
	vd.mv(CShalfinv[i],vc);
	A[i]=vd.dot(vd);
	vc.mv(CShalfinv[i],step);
	B[i]=2*vc.dot(vd);
	C[i]=vc.dot(vc);
    }
}
    
double LineDens::operator()(const int x, GMM& gmm)
{
    double s=0;
    int K=gmm.getK();
    for(int i=0;i<K;i++){
	s+=exp(gmm.getCC()[i]-(A[i]+B[i]*x+C[i]*x*x)/2);
    }
    return s;//positive sign
}
double linedev(const gvector& x, const gvector& y, GMM& gmm, LineDens& ld, 
	       int N, int dir=0,int cri=0,int n=100)
//dir=0, absolute, dir=-1, change sign.
//cri=0, the relative difference to the fitted line
//cri=1, the relative for the peaks and valleys.
{
    {
	static int ilinedev=1;
	//cerr<<"ilinedev="<<ilinedev<<'\n';
	ilinedev++;
    }
    ld.change_xy(x,y,gmm,n);
    int Nx=ceil(gmm.Cw[gmm.ck(x)]*N);
    int Ny=ceil(gmm.Cw[gmm.ck(y)]*N);
    //double sf=1/pow(Nx+Ny,1.0/2)*3.0; //penalize to mix the 
    //peaks with sparse density and think density
    //relatively to the average size of the clusters
    double sf=pow(2.0*N/(gmm.getK()+0.0)/(Nx+Ny),1.0/2);
    double px= ld(0,gmm);
    double py= ld(n,gmm);
    double dev=0,devi,pz;

    if(cri==0){
	for(int i=1;i<n;i++){
	    pz=ld(i,gmm);
	    double f_fit=px+(i+0.0)/n*(py-px);
	    devi=(f_fit-pz)/f_fit;
	    if(dir==0){
		devi=fabs(devi);
	    }else if(dir==-1){
		devi=-devi;
	    }
	    if(devi>dev){
		dev=devi;
	    }
	}
	dev/=sf;
	//just use the counts,it didn't work for the neighboring kmeans clusters
	/*gvector z(x.size);
	z.memcpy(x);
	int Nmin=min(Nx,Ny);
	int Nmax=max(Nx,Ny);
	for(int i=1;i<n;i++){
	    doubleweightsum(z.data,y.data,z.size,1.0/(n-i+1.0));
	    int Nz=ceil(gmm.Cw[gmm.ck(z)]*N);
	    if(Nz>Nmax){
		Nmax=Nz;
	    }else if(Nz<Nmin){
		Nmin=Nz;
	    }
	}
	dev=(Nmax-Nmin)/((sqrt(Nmax)+sqrt(Nmin))*3.0/2.0);*/
    }else{
	//more intelligently searching 
	Vdouble Vz(n+1);
	Vz[0]=px;
	Vz[n]=py;
	for(int i=1;i<n;i++){
	    Vz[i]=ld(i,gmm);
	}
	//find the valley;
	int vi=min_element(&Vz[1],&Vz[n])-&Vz[0];
	//find the two peaks in the two sides
	double peak1=*max_element(&Vz[0],&Vz[vi]);
	double peak2=*max_element(&Vz[vi+1],&Vz[n+1]);
	double valley=Vz[vi];
	double peakmin=min(peak1,peak2);
	devi=(peakmin-valley)/peakmin;
	if(devi>dev){
		dev=devi;
	    }
    }
    return dev;
}
void computeSmatTol(const gmatrix& Cpeaks, const GMM& gmm, gmatrix& S)
{
    int K=S.nrows();
    int p=Cpeaks.ncols();
    VAint ck2(K);
    VAdouble s2(K);
    for(int i=0;i<K;i++){
	const gvector &Ci=Cpeaks[i];
	int ic1, ic2;
	get_IC1_IC2(Ci.data,p,gmm.getK(),gmm.Cm_mat.data,&ic1,&ic2);
	s2[i]=L2dist(Ci.data,gmm.Cm[ic2].data,p);
	ck2[i]=ic2;
    }
    for(int i=0;i<K-1;i++){
	for(int j=(i+1);j<K;j++){
	    double s=sqrt(s2[i])+sqrt(s2[j]);
	    S(i,j)=s*s*4; //note that is has been squared
	    S(j,i)=S(i,j);
	}
    }
}

void unique_peaks( const gmatrix & SD, const double tol,vector<Vint>& groups)
{
    //first divide the peaks into groups with the almost identical data points
    Vint vtmp(1);
    int loc=0;
    int K=SD.nrows();
    while(loc<K){
	unsigned int a;
	for(a=0;a<groups.size();a++){
	    if(SD(loc,groups[a][0])<tol){
		groups[a].push_back(loc);
		loc++;
		break;
	    }
	}
	if(a==groups.size()){
	    //build a new group
	    vtmp[0]=loc;
	    groups.push_back(vtmp);
	    loc++;
	}
    }
}


void MatDevLine(GMM& gmm,const gmatrix &M, gmatrix &FD, int *Nb, 
		int dir=0,int cri=0)
{
    int K=gmm.getK();
    int p=gmm.getp();
    LineDens ld(K);

    //M and FD are having data of uK rows
    int uK=M.nrows();
    gvector_view x,y;
    
    //compute n.
    double n=0;
    for(int i=0;i<K*K;i++){
	n+=Nb[i];
    }
    //In the future, we should clean this up with just n rather than with Nb

    //compute the nearest centers
    VAint Ck(uK);
    VAdouble smin(uK);
    for(int i=0;i<uK;i++){
	x.change_view(M[i]);
	Ck[i]=gmm.ck(x);
	int ic1,ic2;
	get_IC1_IC2(gmm.Cm[Ck[i]].data,p,K,gmm.Cm_mat.data,&ic1,&ic2);
	smin[i]=L2dist(gmm.Cm[ic1].data,gmm.Cm[ic2].data,p);
    }

    //compute the normalizatiion factors
    Vdouble Mtolv(K*K);
    double *Mtol=&Mtolv[0];
    for(int i=0;i<K-1;i++){
	for(int j=(i+1);j<K;j++){
	    int Nx=ceil(gmm.Cw[i]*n);
	    int Ny=ceil(gmm.Cw[j]*n);
	    double px=-gmm.my_f(&gmm.Cm[i]);
	    double py=-gmm.my_f(&gmm.Cm[j]);
	    Mtol[i*K+j]=fabs(px-py)/(px+py)/2.0/
		(1/pow(Nx+Ny,1.0/2)*3.0);
	    Mtol[j*K+i]=Mtol[i*K+j];
	}
	Mtol[i*K+i]=GSL_POSINF;
    }
    //collect the smallest Mtol for each i
    Vdouble Vtolv(K);
    double* Vtol=&Vtolv[0];
    for(int i=0;i<K;i++){
	Vtol[i]=*min_element(Mtol+i*K,Mtol+(i+1)*K);
    }
    //find the median for Vtol
    int halfK=K/2;
    nth_element(Vtol,Vtol+K,Vtol+halfK);
    //double tolsf=*(Vtol+halfK);

    for(int i=0;i<uK-1;i++){
	FD(i,i)=0;
	x.change_view(M[i]);
	for(int j=i+1;j<uK;j++){
	    y.change_view(M[j]);
	    double d=linedev(x,y,gmm,ld,n,dir,cri); ///tolsf;
	    /*double sij=L2dist(x.data,y.data,p);
	     if(sij<=min(smin[i],smin[j])){
		FD(i,j)=d/2.0;
		FD(j,i)=d/2.0;
		continue;	
		}else{*/
	    
	    FD(i,j)=d;
	    FD(j,i)=d;
	    //}
	}
    }
}

void Norm2(const gmatrix &M, gmatrix &D)
{
    int K=D.nrows();
    gvector_view x,y;
    for(int i=0;i<K-1;i++){
	x.change_view(M[i]);
	D(i,i)=0;
	for(int j=(i+1);j<K;j++){
	    y.change_view(M[j]);
	    D(i,j)=L2dist(x,y);
	    D(j,i)=D(i,j);
	}
    }
}

bool minpair(const gmatrix& SD,const gmatrix& SmatTol, const gmatrix& FD,
	     const double Ftol,const gmatrix& checked,
	     const int K,
	     int& a, int &b)
{
//note that SD and FD's dimension may be greater than K*K as I am
//only using the submatrix
//The SD has higher priority
    a=0;b=1;
    bool found=false;
    double d=1e10;//setting it to be a very big number
    for(int i=0;i<K-1;i++){
	for(int j=i+1;j<K;j++){
	    if(checked(i,j)>0 || SD(i,j)>SmatTol(i,j)||FD(i,j)>Ftol) continue;
	    if(SD(i,j)<d){
		d=SD(i,j);
		a=i;
		b=j;
		found=true;
	    }
	}
    }
    return found;
}

void merge_matrix(gmatrix& D, int gsize, int a, int b)
{//assuming a<b<gsize, D is a size K matrix (k>gsize)
 //We are only concerned for the submatrix
    for(int k=0;k<gsize;k++){
	if(k!=a && k!=b){
	    D(a,k)=min(D(a,k),D(b,k));
	    D(k,a)=D(a,k);
	}
    }
    D(a,a)=min(D(a,a),D(b,b));//mostly unnecessary
    //shift data after b-row and b-col
    for(int k=b;k<gsize-1;k++){
	for(int h=0;h<=k;h++){
	    if(h<b){
		D(k,h)=D(k+1,h);
	    }else{
		D(k,h)=D(k+1,h+1);
	    }
	    D(h,k)=D(k,h);
	}
    }
}


void peaks_merge(double *Sd,double *Fd, double*pFtol, double* SmatTold,
		 int *pK, int *cid)
{
    double Ftol=*pFtol; 
    int K=*pK;

    gmatrix_view SD_org, FD_org,SmatTol_org;
    SD_org.change_view(Sd,K,K);
    FD_org.change_view(Fd,K,K);
    SmatTol_org.change_view(SmatTold,K,K);


    gmatrix SD(K,K);
    gmatrix FD(K,K);
    gmatrix SmatTol(K,K);
    SD.memcpy(SD_org);
    FD.memcpy(FD_org);
    SmatTol.memcpy(SmatTol_org);


    vector<Vint> peaks;
    peaks.resize(K);
    for(int i=0;i<K;i++){
	peaks[i].push_back(i);
    }
    gmatrix checked(K,K);
    checked.set_all(-1);
    while(true){
	int gsize=peaks.size();
	if(gsize==1) break;//since everything is in one peaks

	int a,b;
	if(!minpair(SD,SmatTol,FD, Ftol,checked,gsize,a,b)) break; //can not find suitable pairs
	Vint &ga=peaks[a];
	Vint &gb=peaks[b];
	//need to check if all of them are correct
	double verified=false;
	for(unsigned int i=0;i<ga.size();i++){
	    for(unsigned int j=0;j<gb.size();j++){
		int k1=ga[i];
		int k2=gb[j];
		if(SD_org(k1,k2)<=SmatTol_org(k1,k2) && FD_org(k1,k2)<=Ftol){
		    verified=true;
		    break;
		}
	    }
	    if(verified) break;
	} 
	if(!verified){
	    checked(a,b)=1;
	    continue; //find the next possible pair
	}
	//peaks a and b together
	ga.insert(ga.end(),gb.begin(),gb.end());
	peaks.erase(peaks.begin()+b); //erase the peaks b
	//recompute the distance matrix
	merge_matrix(SD,gsize,a,b);
	merge_matrix(FD,gsize,a,b);
	merge_matrix(SmatTol,gsize,a,b);
	//reset the checked
	checked.set_all(-1);
    }
    //output the peaks into cid
    for(unsigned int i=0;i<peaks.size();i++){
	for(unsigned int j=0;j<peaks[i].size();j++){
	    cid[peaks[i][j]]=i;
	}
    }
}

//useful for the get the unique peaks
//and peaks finding
double get_maxstepsize(GMM& gmm,int n=10)
{
    int p=gmm.getp();
    int K=gmm.getK();

    VAdouble h(p);
    for(int i=0;i<p;i++){
	double xmin= 1e10;
	double xmax= -1e10;
	for(int j=0;j<K;j++){
	    double m=gmm.Cm[j][i];
	    double threesd=3*sqrt(gmm.CS[j](i,i));
	    if(m+threesd>xmax) xmax=m+threesd;
	    if(m-threesd<xmin) xmin=m-threesd;
	}
	h[i]=(xmax-xmin)/n;
    }

    double* hloc=&h[0];
    int halfp=p/2;
    nth_element(hloc,hloc+p,hloc+halfp);
    return hloc[halfp];
}

//it is too convervative
double medianbinh(GMM& gmm)
{
    int p=gmm.getp();
    int K=gmm.getK();

    gmatrix M(p,K);
    
    for(int i=0;i<p;i++){
	for(int j=0;j<K;j++){
	    M(i,j)=gmm.CS[j](i,i);
	}
    }

    VAdouble h(p);
    for(int i=0;i<p;i++){
	double *loc=&M(i,0);
	int halfK=K/2;
	nth_element(loc,loc+halfK,loc+K);
	h[i]=*(loc+halfK);
    }
    double* hloc=&h[0];
    int halfp=p/2;
    nth_element(hloc,hloc+p,hloc+halfp);
    double binh=*(hloc+halfp);
    return sqrt(binh)/K/2;
}

double get_min(const gvector& x, GMM & gmm,
	       gvector & y, gvector& yd, bool & found,
	       const double maxstepsize=1e8)
{
    int n=x.size;
    int ck=gmm.ck(x);
    double kstepsize=gmm.Csigma_min[ck]/3;
    gvector grad1(n);

    gsl_multimin_function_fdf my_func;
    my_func.n = n;
    my_func.f = &my_f;
    my_func.df = &my_df;
    my_func.fdf = &my_fdf;
    my_func.params = (void *)&gmm;

    const gsl_multimin_fdfminimizer_type *T;
    T= gsl_multimin_fdfminimizer_steepest_descent;
    //T=gsl_multimin_fdfminimizer_conjugate_fr;

    gsl_multimin_fdfminimizer *s;
    s = gsl_multimin_fdfminimizer_alloc (T,n );
    gsl_multimin_fdfminimizer_set (s, &my_func, &x, kstepsize/10, 0.25);
    /*when it succeeds, move twice as much,
      when it fails, shrink to one quarter*/
    found=false;
    size_t iter = 0;
    int status;
    do
    {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate (s);
	if (status)
	    break;
	double dxsize=sqrt(L2dist(s->dx->data,n));
	if(dxsize<kstepsize*2e-4){
	    /*equivalently about 6 times shrinking from kstepsize,
	  since 1/4^6=2.4e-4 */
	    found=true; //converges, not wanting to waste time on it;
	    break;
	}
	steepest_descent_state_t *state = (steepest_descent_state_t *) s->state;

	int newck=gmm.ck(*(s->x)); 
        //This computation is much less than computing one f and df, so 
	//it can safely ignored to handle performance details
	
	if(L2dist(gmm.Cm[newck].data,s->x->data,n)<
	   L2dist(gmm.Cm[ck].data,s->x->data,n)){
	    //update ck and kstepsize
	    ck=newck;
	    kstepsize=gmm.Csigma_min[ck]/3;
	    state->step=kstepsize/10;
	    //moved to a new cluster, starting to work 
	    //with the new cluster center
	    double f1=gmm.my_fdf(&gmm.Cm[newck],&grad1);
	    if(f1 < s->f){//improves, move to the center
		s->f=f1;
		doublecopy(gmm.Cm[ck].data,s->x->data,n);
		doublecopy(grad1.data,s->gradient->data,n);
	    }//else stay where they are
	}else{
	    if(state->step > kstepsize){
		state->step=kstepsize; //constrain the movements of the points
	    }
	}
    }while (iter < 10000);
    
    y.memcpy(*s->x);
    yd.memcpy(*s->gradient);
    double fmin=s->f;
    gsl_multimin_fdfminimizer_free (s);
    return fmin;
}

void get_flowpeaks(double *Cw,double* Cm, double* CS, int *pK, int *pp, 
	      int *Nb,double *Cpeaks, double *Cpeaks_f, double * Cpeaks_df, 
		   int* cfound,int *cid,  double* ptol)
{
    //step 0, collect information.
    GMM gmm(Cw,Cm,CS,pK,pp);
    int p=*pp;
    int K=*pK;

    //step 1, obtain the peaks
    gvector_view y,yd;
    double maxstepsize=get_maxstepsize(gmm);//medianbinh(gmm);
    //set<int> obrit;
    for(int i=0;i<K;i++){
	yd.change_view(Cpeaks_df+i*p,p);
	y.change_view(Cpeaks+i*p,p);
	bool found;
	char msg[40];
	sprintf(msg,"i=%d\t",i);
	//gsl_stream_printf("Information",__FILE__,__LINE__,msg);
	Cpeaks_f[i]=get_min(gmm.Cm[i],gmm,y,yd,found,maxstepsize);
	cfound[i]=found;
    }
    //step 2, peaks merging

    //step 2.1 Find the unique peaks
    gmatrix_view Cpeaks_M(Cpeaks,K,p);
    gmatrix SD(K,K);
    Norm2(Cpeaks_M,SD);
    double Stol=maxstepsize*maxstepsize/100;
    vector<Vint> peak_groups;
    unique_peaks(SD,Stol,peak_groups);
    //step 2.2 Update the uCpeaks,and the original output with the optimum peaks,
    int uK=peak_groups.size();
    gmatrix uCpeaks(uK,p);
    for(int i=0;i<uK;i++){
	Vint& gi=peak_groups[i];
	int loc=gi[0];
	double fi=Cpeaks_f[loc];
		
	for(unsigned int j=1;j<gi.size();j++){
	    double fj=Cpeaks_f[gi[j]];
	    if(fj<fi){//a higher peak because of the negative sign
		fi=fj;
		loc=gi[j];
	    }
	}
	uCpeaks[i].memcpy(Cpeaks_M[loc]);
	if(gi.size()==1) continue;
	//update the original output from the peaks finding
	for(unsigned int j=0;j<gi.size();j++){
	    int gloc=gi[j];
	    Cpeaks_M[gloc].memcpy(Cpeaks_M[loc]);
	    Cpeaks_f[gloc]=Cpeaks_f[loc];
	    cfound[gloc]=cfound[loc];
	    copy(Cpeaks_df+loc*p,Cpeaks_df+(loc+1)*p,Cpeaks_df+gloc*p);
	}
    }

    //step 3 recompute the uSD, and uFD, and merge
    gmatrix uSD(uK,uK);
    Norm2(uCpeaks,uSD);//could use a submatrix to speed up, the computation time maybe too little to be taken care of
    gmatrix uFD(uK,uK); //the number of unique peaks.
    MatDevLine(gmm,uCpeaks,uFD,Nb,0,0);
    Vint ucid(uK);
    gmatrix SmatTol(uK,uK);
    computeSmatTol(uCpeaks,gmm,SmatTol);
    //check what's wrong with merging of 3
    /*fstream ofile1("Smat.txt",fstream::out);
    ofile1<<SmatTol;
    fstream ofile2("Ftol.txt",fstream::out);
    ofile2<<uFD;
    fstream ofile3("Cpeaks.txt",fstream::out);
    ofile3<<"\t"<<uCpeaks;
    fstream ofile4("uSD.txt",fstream::out);
    ofile4<<"\t"<<uSD;*/

    //find out possible breaks
    /*{
	fstream ofile("tol.txt",fstream::out);
	double tol1=uFD.data[1];
	double tol2=tol1;
	for(int i=0;i<uK-1;i++){
	    for(int j=(i+1);j<uK;j++){
		double tmp=uFD.data[i*uK+j];
		if(tmp>tol2){
		    tol2=tmp;
		}else if (tmp<tol1){
		    tol1=tmp;
		}
	    }
	}
	tol1*=0.9;
	tol2*=1.1;
	for(int i=0;i<=200;i++){
	    double tol=tol1*exp((log(tol2)-log(tol1))*i/200);
	    peaks_merge(uSD.data,uFD.data,&tol,SmatTol.data,&uK,&ucid[0]);
	    int nK=*max_element(&ucid[0],&ucid[0]+uK)+1;
	    ofile<<tol<<'\t'<<nK<<'\n';
	}
	ofile.close();
	}*/



    peaks_merge(uSD.data,uFD.data,ptol,SmatTol.data,&uK,&ucid[0]);
    //step 4 gets back to the original kmeans
    for(int i=0;i<uK;i++){
	Vint& gi=peak_groups[i];
	for(unsigned int j=0;j<gi.size();j++){
	    cid[gi[j]]=ucid[i];
	}
    }
}

void raster_image(double *raw, int *rawid,
		  int *pn, int* pres, double*grid, 
		  int *grid_id, int*pngrid) 
{
    int n=*pn;
    int nres=*pres;//should be between 400 to 2000
    if(nres<400 || nres>2000){
	gsl_error("The resolution is too high or too low",
		  __FILE__,__LINE__,GSL_EDOM);
	
    }
    gmatrix_view rawM;
    rawM.change_view(raw,n,2);
    gvector_view xcol,ycol;
    xcol=rawM(0);
    double xmin=xcol.min();
    double xmax=xcol.max();
    ycol=rawM(1);
    double ymin=ycol.min();
    double ymax=ycol.max();
    double xstep=(xmax-xmin)/(nres);
    double ystep=(ymax-ymin)/(nres);
    
    gmatrix gridM(nres+1,nres+1);
    for(int i=0;i<n;i++){
	int loci=round((rawM(i,0)-xmin)/xstep);
	int locj=round((rawM(i,1)-ymin)/ystep);
	gridM(loci,locj)=rawid[i];
    }
 
    gmatrix_view gridMV(grid,n,2);
    int ngrid=0;
    for(int i=0;i<=nres;i++){
	for(int j=0;j<=nres;j++){
	    if(gridM(i,j)>0){
		gridMV(ngrid,0)=i*xstep+xmin;
		gridMV(ngrid,1)=j*ystep+ymin;
		grid_id[ngrid]=round(gridM(i,j));
		ngrid++;
	    }
	}
    }
    *pngrid=ngrid;
}

//A second version, which only looks at the neibors
void get_flowpeaks2(double *Cw,double* Cm, double* CS, int *pK,int *pp, 
		     int *Nb, double* Cm_f, double* Cpeaks,
		     int *cid, double *ptol)
{
    int K=*pK;
    int p=*pp;
    
    GMM gmm(Cw,Cm,CS,pK,pp);
    VAint Nbk(K);
    double n=0;
    for(int k=0;k<K;k++){
	for(int j=0;j<K;j++){
	    n+=Nb[k*K+j];
	}
	Nbk[k]=max_element(&Nb[k*K],&Nb[(k+1)*K])-&Nb[k*K];
    }
    n/=2.0; //it should be the same as given.

    VAint cnext(K); //the better neighbors
    gvector midpt(p);
    //use the density to find neighbors
     LineDens ld(K);
    for(int k=0;k<K;k++){
	double f0=gmm.my_f(&gmm.Cm[k]);
	Cm_f[k]=f0;
	double df=0; 
	int ck=-1;
	double nmin=1.0/(K+0.0)*Cw[k]*n; //a little bit conservative
	//use the closest neighbors
	double smin=L2dist(gmm.Cm[k].data,gmm.Cm[Nbk[k]].data,p)*4;
	for(int l=0;l<K;l++){
	    if(Nb[k*K+l]>=nmin){//to make sure they share
		//enough of data points for these two clusters 
		double f1=gmm.my_f(&gmm.Cm[l]);
		if(f1>=f0) continue;
		double sl=L2dist(gmm.Cm[k].data,
				 gmm.Cm[l].data,p);
		if(sl>smin) continue;
		double tol=linedev(gmm.Cm[k],gmm.Cm[l],gmm,ld,n,1,1);
		if(tol>ptol[0]/2.0) continue; //ignore the data with deep valley
		double df1=(f1-f0)/sqrt(sl);
		if(df1<df){ 
		    df=df1;//with a better derivative
		    ck=l;
		}
	    }
	}
	if(ck==-1){
	    ck=k;
	}
	cnext[k]=ck;
    }

    VAint peaks_id(-1,K);
    Vint peaks_set;
    Vint gint;
    for(int k=0;k<K;k++){
	int l=k;
	int peak_old=-1;
	gint.clear();
	while(true){
	    gint.push_back(l);
	    if(peaks_id[l]>0){//found existing peak
		peak_old=peaks_id[l];
		break;
	    }
	    int lnext=cnext[l];
	    if(l==lnext) break;//found new peak
	    l=cnext[l];
	}
	if(peak_old==-1){
	    peak_old=l;
	    peaks_set.push_back(l);
	}
	for(unsigned int j=0;j<gint.size();j++){
	    peaks_id[gint[j]]=peak_old;
	}
    }
    //copy Cpeaks
    for(int k=0;k<K;k++){
	double *loc=gmm.Cm[peaks_id[k]].data;
	doublecopy(loc,Cpeaks+k*p,p);
    }

    int uK=peaks_set.size();
    gmatrix uCpeaks(uK,p);
    for(int uk=0;uk<uK;uk++){
	uCpeaks[uk].memcpy(gmm.Cm[peaks_set[uk]]);
    }
    gmatrix uSD(uK,uK);
    Norm2(uCpeaks,uSD);
    gmatrix uFD(uK,uK); //the number of unique peaks.
    MatDevLine(gmm,uCpeaks,uFD,Nb,0,0);

    gmatrix SmatTol(uK,uK);
    //fill the SmatTol to be very big number so that this parameter is useless
    double smax=*max_element(uSD.data,uSD.data+uK*uK);
    SmatTol.set_all(smax*1.5);
    computeSmatTol(uCpeaks,gmm,SmatTol);
    VAint ucid(uK);
    peaks_merge(uSD.data,uFD.data,ptol,SmatTol.data,&uK,&ucid[0]);

    VAint peaks2ucid(K);
    for(int uk=0;uk<uK;uk++){
	peaks2ucid[peaks_set[uk]]=ucid[uk];
    }

    for(int k=0;k<K;k++){
	cid[k]=peaks2ucid[peaks_id[k]];
    }
}

void assign_flowPeaks(double *A, int *pn, int *pp, double *Cw,double*Cm,
		      double *CS,int *pK, int* cid, double* ftol,double*ffc,
		      int *ret_IC)
{
    int n=*pn;
    int K=*pK;
    int p=*pp;
    //find out the number of clusters
    int nC=*max_element(cid,cid+K)+1;
    double tol=*ftol;
    double fc=*ffc;

    vector<GMM> peakinfo(nC);
    vector<double> peakdens(nC);
    {
	double* tmpCw=mymalloc(K, double);// (double*)malloc(K*sizeof(double));
	double* tmpCm=mymalloc(K*p, double);//(double*)malloc(K*p*sizeof(double));
	double* tmpCS=mymalloc(K*p*p,double); //(double*)malloc(K*p*p*sizeof(double));
	for(int i=0;i<nC;i++){
	    int tmpK=0;
	    for(int j=0;j<K;j++){
		if(cid[j]==i){
		    tmpCw[tmpK]=Cw[j];
		    doublecopy(Cm+j*p,tmpCm+tmpK*p,p);
		    doublecopy(CS+j*p*p,tmpCS+tmpK*p*p,p*p);
		    tmpK++;
		}
	    }
	    //renormlize tmpCw
	    double s=doublesum(tmpCw,tmpK);
	    doublemul(tmpCw,1/s,tmpK);
	    peakinfo[i].Init(tmpCw,tmpCm,tmpCS,&tmpK,pp);
	    peakdens[i]=peakinfo[i].pseudo_maxf();
	}
	myfree(tmpCw);
	myfree(tmpCm);
	myfree(tmpCS);
    }
    gvector prob(nC);
    gvector CC(nC);
    for(int j=0;j<nC;j++){
	CC[j]=peakinfo[j].getCCmin();
    }
    CC-=CC.min();;
    for(int i=0;i<n;i++){
	int k=cid[get_IC(A+i*p,p,K,Cm)];
	gvector_view a(A+i*p,p);
	double f=-peakinfo[k].my_f(&a);
	if(f<peakdens[k]*tol){
	    ret_IC[i]=-(nC+1);
	}else{
	    if(nC==1){
		ret_IC[i]=1;
		continue;
	    }
	    //check for the posterier distributions
	    for(int j=0;j<nC;j++){
		prob[j]=(-peakinfo[j].my_f(&a))*exp(CC[j]);
	    }
	    prob/=prob.sum();
	    /*int k1=0;//the best
	    int k2=1; //the second best
	    if (prob[1]>prob[0]){
		swap(k1,k2);
	    }
	    for(int j=2;j<nC;j++){
		if(prob[j]>prob[k1]){
		    k2=k1;
		    k1=j;
		}else if(prob[j]>prob[k2]){
		    k2=j;
		}
	    }
	    if(k1!=k){
		ret_IC[i]=-(k1+1);
		}else */
	    if(prob[k]<fc){
		ret_IC[i]=-(nC+2);
	    }else{
		ret_IC[i]=k+1;
	    }
	}
    }
}
