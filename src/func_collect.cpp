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
#include <ctime>
//#include <fstream>
#include <cstring>
#include <cctype>
#include <stdexcept>
#define EPS (1e-6)
double summarize_IC1(const double *A, const int n, const int p, const int K,
		     const int* IC1, 
		     int *nc, double*C, double *D)
{
    for(int k=0;k<K;k++){
	for(int j=0;j<p;j++){
	    C[k*p+j]=0;
	}
	nc[k]=0;
    }
    for(int i=0;i<n;i++){
	int k=IC1[i];
	nc[k]++;
	for(int j=0;j<p;j++){
	    C[k*p+j]+=A[i*p+j];
	}
    }
    for(int k=0;k<K;k++){
	if(nc[k]==0){
	    //reset the center
	    gsl_stream_printf("Warning",__FILE__, __LINE__,"Empty clusters");
	    const double* Aloc=A+gsl_rng_uniform_int(g_rng.r,n)*p;
	    doublecopy(Aloc,C+k*p,p);
	}else{
	    doublemul(C+k*p,1.0/(nc[k]),p);
	}
    }
    if(D==NULL) return GSL_POSINF; //no point to compute twss either
    //compute the intial twss
    double twss=0;
    for(int i=0;i<n;i++){
	int k=IC1[i];
	D[i]=L2dist(A+i*p,C+k*p,p);
	twss+=D[i];
	if(nc[k]!=1){
	    D[i]*=nc[k]/(nc[k]-1.0);
	}else{
	    D[i]=0;
	}
    }
    return twss;
}

int get_IC(const double *a, const int p, const int K, const double *C, double* pd)
{
    int IC=0;
    double d=L2dist(a,C,p);
    for(int k=1;k<K;k++){
	double dk=L2dist(a,C+k*p,p);
	if(dk<d){
	    IC=k;
	    d=dk;
	}
    }
    if(pd!=NULL){
	*pd=d;
    };
    return IC;
}
double get_IC1_IC2(const double *a, const int p,
		 const int K, const double *C,
		 int* pIC1, int*pIC2)
{
    int IC1=0;
    int IC2=1;
    double d1=L2dist(a,C,p);
    double d2=L2dist(a,C+p,p);
    if(d2<d1){
	swap(IC1,IC2);
	swap(d1,d2);
    }
    for(int k=2;k<K;k++){
	double dk=L2dist(a,C+k*p,p);
	if(dk<d1){
	    IC2=IC1;
	    d2=d1;
	    IC1=k;
	    d1=dk;
	}else if(dk<d2){
	    IC2=k;
	    d2=dk;
	}
    }
    *pIC1=IC1;
    *pIC2=IC2;
    return d1;
}
void get_summarize(double *A, int *pn, int*pp, int *pK, int *IC, 
		   int*ret_nc, double* ret_C, double* ret_S, double *ret_twss,bool computeSonly)
{
    int n=*pn;
    int p=*pp;
    int K=*pK;

    //first compute ret_C, ret_nc
    if(!computeSonly){
	VAdouble D(*pn);
	*ret_twss=summarize_IC1(A,*pn,*pp,*pK,IC,ret_nc,ret_C,&D[0]);
    }

    //compute S
    if(ret_S==NULL){
	return; //no computing of S
    }

    //initiaize to be zero
    double *S;
    S=ret_S;
    doublecopy(S,0,K*p*p);
    //compute S
    for(int i=0;i<n;i++){
	int k=IC[i];
	double *a=A+i*p;
	double *c=ret_C+k*p;
	double *s=S+k*p*p;
	for(int j=0;j<p;j++){
	    for(int l=j;l<p;l++){
		s[j*p+l]+=(a[j]-c[j])*(a[l]-c[l]);
	    }
	}
    }
    //clean up
    for(int k=0;k<K;k++){
	double *s=S+k*p*p;
	double w0=0;
	if(ret_nc[k]>1){
	    w0=1/(ret_nc[k]-1.0);
	}
	doublemul(s,w0,p*p);
	doublecopy2lower(s,p);
    }
}
void doublecopy2lower(double *A, const int p)
{
    for(int i=0;i<p;i++){
	for(int j=i+1;j<p;j++){
	    A[j*p+i]=A[i*p+j];
	}
    }
}

double computeD(const double *A, const double *C, const int n, 
		const int p, const int k, double* D)
{
    double s=0;
    for(int i=0;i<n;i++){
	double d=L2dist(A+i*p,C+k*p,p);
	if(d<D[i]){
	    D[i]=d;
	}
	s+=D[i];
    }
    return s;
}
double SeedPlusPlus(const double *A, int n, int p, int K, double *C)
{
    int loc=gsl_rng_uniform_int(g_rng.r,n);
    doublecopy(A+loc*p,C,p);
    
    gvector DV(n);
    double* D=DV.data;
    double s = 0;
    for (int i = 0; i < n; i++) {
	D[i] = L2dist(A+i*p, C, p);
	s+=D[i];
    }

    for (int k = 1; k < K; k++) {
	//I restrict that D[i] > s/(n_*5.0) to avoid too close centers
	double dmin=s/(n*5.0);
	int found=-1;
	while(found<0){
	    double cutoff = s*gsl_rng_uniform(g_rng.r)*(1-EPS);
	    double scur = 0;
	    for (int i = 0; i < n; i++) {
		scur += D[i];
		if (scur >= cutoff && D[i]>dmin){
		    found=i;
		    break;
		}
	    }
	    //when not found, we have to search another seed
	}
	doublecopy(A + found*p,C+k*p,p);
	s = computeD(A,C,n,p,k,D);
    }
    return s;
}
void sample_nK(int n, int K, int *p)
{
    VAint remainVv(n);
    int *remainV=&remainVv[0];
    for(int i=0;i<n;i++){
        remainV[i]=i;
    }
    int remainN = n;
    for (int j = 0; j < K; j++) {
        int h = gsl_rng_uniform_int(g_rng.r,remainN);
        remainN--;
        p[j]=remainV[h];
        remainV[h] = remainV[remainN];
    }
}
double getRunningTime(bool init, double time0_)
{
    static double time_at=0;//with a higher precision, but easily overflows
    static time_t sec_at=0;//good say seconds is greater than 10.
    static double time0=0;
    if(init) {
	time0=time0_;
	time_at=double(clock()) / CLOCKS_PER_SEC;
	sec_at=time(NULL);
	return time_at+time0;
    }
    double timelap=difftime(time(NULL),sec_at);
    if(timelap<100){
	return double(clock()) / CLOCKS_PER_SEC - time_at+time0;
    }else{
	return timelap+time0;
    }
}
void read_file(const char *filename, int* plength,Vchar& buffv){
    FILE* fp=fopen(filename,"rb");
    char msgstr[902];

    if(fp==NULL){
	 *plength=0;
	 snprintf(msgstr,900,"Can not opne file %s\n",filename);
	 throw(domain_error(msgstr));
     }
    fseek(fp,0,SEEK_END);
    int length=ftell(fp);
    buffv.resize(length+2);
    char *buff=&buffv[0];
    rewind(fp);
    int count=fread(buff,1,length,fp);
    fclose(fp);
    if(count!=length){
	snprintf(msgstr,900,"Can not read enough data in %s\n",filename);
	throw(domain_error(msgstr));
    }

    if(buff[length-1]!='\n'){//append an '\n' at the end of the file 
	buff[length++]='\n';
    }
    buff[length+1]=0;
    *plength=length;
}
/*void headline(const char *filename, const char sep, Vstring& vline, int skip)
{
    ifstream ifs (filename, ifstream::in|ifstream::binary );
    string line;
    for(int i=0;i<skip;i++){
	getline(ifs,line);
    }
    getline(ifs,line);
    split(line,vline,sep);
    ifs.close();
    }*/
void scanfast(const char *filename, const int skiprows,const int skipcols,
	      const char sep,
	      vector<double>&xv, int* pnrows) 
{
    char msgstr[902];
    xv.reserve(20000);//concentrate on large data files
    //skip lines
    int length;
    Vchar buffv;
    read_file(filename,&length,buffv);
    char* buff=&buffv[0];
    //skip lines
    char* p=buff;
    char* plim=buff+length;
    for(int i=0;i<skiprows;i++,p++){
	p=(char*)memchr(p,'\n',plim-p);
	if(p>=plim){
	    snprintf(msgstr,900,"there are only %d lines, can not skip %d rows\n",
		    i,skiprows);
	    throw(domain_error(msgstr));
	}
    }
    //read the data
    int nlines=skiprows;
    int ncolhead=0;
    while(p<plim){
	char* pline=(char*)memchr(p,'\n',plim-p);
	if(pline>=plim){
	    throw(domain_error("Something is wrong in the las line of the data"));
	};
	//one line at a time
	nlines++;
	
	//first skip columns
	for(int j=0;j<skipcols;j++){
	    p=(char*)memchr(p,sep,pline-p);
	    if(p==NULL||p==pline){
		snprintf(msgstr,900,"We can not skip %d columns at line %d\n",
			skipcols,nlines);
		throw(domain_error(msgstr));
	    }
	}

	//read the remaining of the dataset
	int ncol=0;
	for(;p<pline;p++){

	    char *pnext=(char*)memchr(p,sep,pline-p);
	    if(pnext==NULL){
		pnext=pline;
	    }

	    char *oldp=p;
	    double d=strtod(oldp,&p);
	   
	    if(p==oldp||p>pnext){
		snprintf(msgstr,900,"The data is incorrect number format at line %d\n",
			nlines);
		throw(domain_error(msgstr));
	    }
	    //read the numbers succesfully
	    ncol++;
	    xv.push_back(d);

	    //clean up the last columns
	    if(pnext==pline){
		//skip white spaces
		for(;p<pline;p++){
		    if(!isspace(*p)){
			snprintf(msgstr,900,"The data format is wrong for the last column at line %d\n",
				nlines);
			throw(domain_error(msgstr));
		    }
		}
	    }
	    if(p>=pline) 
	    {
		if(nlines==skiprows+1){
		    ncolhead=ncol;
		}else{
		    if(ncolhead!=ncol){
			snprintf(msgstr,900,"The data is has %d numbers at line %d\n",
				ncol,nlines);
			throw(domain_error(msgstr));	
		    }
		}
	    }
	}
    }
    *pnrows=nlines-skiprows;
}

void doubletranspose(const double *A, const int n, const int p,
		     double *B)
//if A is NULL, the transpose is copied to itself in B 
{
    if(A==B ||B==NULL){
	throw(domain_error("You can not set B to be the same as A or B to be NULL"));
    }
    Vdouble Acv;
    double* Ac;
    if(A==NULL){
	Acv.resize(n*p);
	Ac=&Acv[0];
	memcpy(Ac,B,n*p*sizeof(double));
    }else{
	Ac=const_cast<double*>(A);
    }
    for(int i=0;i<n;i++){
	for(int j=0;j<p;j++){
	    B[i*p+j]=Ac[j*n+i];
	}
    }
}
void adjustS(double *S, double *x,
	     double *w,double h0, double h, 
	     int n, int p,int K)
{
    //first compute S0;
    Vdouble S0v(p*p);
    double *S0=&S0v[0];
    doublecopy(S0,0.0,p*p);
    for(int j=0;j<p;j++){
	double zmin=x[j];
	double zmax=zmin;
	for(int i=1;i<n;i++){
	    double xnow=x[i*p+j];
	    if(xnow<zmin){
		zmin=xnow;
	    }else if(xnow>zmax){
		zmax=xnow;
	    }
	}
	S0[j*p+j]=pow((zmax-zmin)/pow(K+0.0,1.0/p),2.0)/3.0;
    }
    for(int k=0;k<K;k++){
	double* Sk=S+k*p*p;
	double lambda=K/(w[k]*n+K);
	for(int i=0;i<p*p;i++){
	    Sk[i]=Sk[i]*h+S0[i]*h0*lambda;
	}
    }
}
