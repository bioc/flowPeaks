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
#include <sstream>
#include <algorithm>
#include <stdexcept>
RNG g_rng;
void error_msg(int lineN);
CBLAS_TRANSPOSE_t get_transposeid(bool transpose){
    if(transpose) return CblasTrans;
    return CblasNoTrans;
}
int gvector::init(const size_t n, bool clean){
    //rewrite the gsl_vector calloc and alloc functions
    if (n==0)
    {
	gsl_error("vector dimension n must be positive integer",
		  __FILE__,__LINE__,GSL_EDOM);
    }
    block=gsl_block_alloc(n);
    if(block==0){
	gsl_error("failed to allocate space for block",__FILE__,__LINE__,GSL_ENOMEM);
    }
    data=block->data;
    size=n;
    stride=1;
    owner=1;
    if(clean){
	set_zero();
    }
    return 0;
}
int gvector::resize(size_t n, bool clean){
    if (!owner){
	gsl_error("You can't resize a vector view", __FILE__,__LINE__,
		  GSL_EINVAL);
    }
    if(n==size) {
	if(clean){
	    set_zero();
	    return 0;
	}
    }
    if(size>0 && owner) gsl_block_free(block);
    return init(n,clean);
}

/*gvector& gvector::operator=(const gsl_vector& other){
    if(this==&other) return *this;
    if (!owner || stride!=1){
	cerr<<"You can't copyassign into a gvector_view, please use memcpy or changeview function\n";
	exit(0);
    }
    if(size!=other.size){
	cerr<<"When you copy assign vectors, please be sure that \n";
	cerr<<"the two vectors have the same sizes\n";
	exit(0);
    }

    //resize(other.size,false); 
    memcpy(other);
    return *this;
    }*/
gvector& gvector::operator-(){
    for(unsigned int i=0;i<size;i++){
	(*this)[i]=-(*this)[i];
    }
    return *this;
}
void gvector_view::assign(const gsl_vector & other){
    block=other.block;
    data=other.data;
    size=other.size;
    stride=other.stride;
    owner=0;
} 



int gmatrix::init(size_t n1,size_t n2, bool clean){
    if (n1 == 0)
    {
	gsl_error("matrix dimension n1 must be positive integer",
		  __FILE__,__LINE__,GSL_EINVAL);
    }
    else if (n2 == 0)
    {
	gsl_error("matrix dimension n2 must be positive integer",
		  __FILE__,__LINE__,GSL_EINVAL);
    }
    
    block=gsl_block_alloc(n1*n2);
    if(block==0){
	gsl_error("failed to allocate space for block",
		  __FILE__,__LINE__,GSL_ENOMEM);
    }
    data=block->data;
    size1=n1;
    size2=n2;
    tda=n2;
    owner=1;
    if(clean){
	set_zero();
    }
    return 0;
}

int gmatrix::resize(const size_t n1, const size_t n2,bool clean){
    if (!owner){
	gsl_error("You can't resize a matrix view",
		  __FILE__,__LINE__,GSL_EINVAL);
    }
    if(n1==size1 && n2==size2) {
	if(clean){
	    set_zero();
	}
	return 0;
    }

    if(size1>0 && size2>0 &&owner==1) gsl_block_free(block);
    return init(n1,n2,clean);
}

int gmatrix::transpose(){
    //first store the data elsewhere and then copy it back
    gmatrix tmp(size1,size2);
    doublecopy(data,tmp.data,size1*size2);
    std::swap(size1,size2);
    for(size_t i=0;i<size1;i++){
	for(size_t j=0;j<size2;j++){
	    data[i*size2+j]=tmp.data[j*size1+i];
	}
    }
    return 1;
}

/*
gmatrix& gmatrix::operator=(const gsl_matrix& other){
    if(this==&other) return *this;
    if (!owner){
	cerr<<"You can't copy assign into a gmatrix_view, please use memcpy or changeview function\n";
	exit(0);
    }
    if(size1!=other.size1 ||size2!=other.size2){
	cerr<<"When you copy assign matrices, please be sure that \n";
	cerr<<"the two matrices have the same dimensions\n";
    }
    //resize(other.size1,other.size2,false);
    memcpy(other);
    return *this;
    }*/
gmatrix& gmatrix::operator-(){
    for(unsigned int i=0;i<size1;i++){
	for(unsigned int j=0;j<size2;j++){
	    (*this)(i,j)=-(*this)(i,j);
	}
    }
    return *this;
}
void gmatrix_view::assign(const gsl_matrix & other){
    block=other.block;
    data=other.data;
    size1=other.size1;
    size2=other.size2;
    tda=other.tda;
    owner=0;
} 

int gvector::solve(const gsl_matrix &other, const gsl_vector &b){
    gmatrix A(other.size1,other.size2);
    A.memcpy(other);
    gpermutation p(A.nrows());

    int signnum;
    gsl_linalg_LU_decomp(&A,&p,&signnum);
    return gsl_linalg_LU_solve(&A,&p,&b,this); //*this=A^{-1} b
}

 int gmatrix::inverse(const gsl_matrix &other){
    gmatrix A(other.size1,other.size2);
    A.memcpy(other);
    gpermutation p(A.nrows());

    int signnum;
    gsl_linalg_LU_decomp(&A,&p,&signnum);
    return gsl_linalg_LU_invert(&A,&p,this); //*this=A^{-1}
}

int gmatrix::half(const gsl_matrix &other){
    const gmatrix_view otherv(other);
    if(!otherv.is_square() || !otherv.is_symmetric()){
	gsl_stream_printf("ERROR",__FILE__,__LINE__,
			  "No square root matrix can be computed");
	gsl_error("The matrix is not square or symmetric",__FILE__,__LINE__,
		  GSL_EDOM);
    }
    int size=size1;
    gmatrix U(size,size);
    gmatrix V(size,size);
    gvector s(size);
    gmatrix S(size,size);
    otherv.svd(U,s,V);
    S.set_zero();
    for(int i=0;i<size;i++){
	S(i,i)=sqrt(s[i]); //set the diagnonal matrix;
    }
    //multiply the matrix to get the final product
    gmatrix A1(size,size); //avoid the aliasing
    A1.mm(U,S);
    mm(A1,V,false,true);
    return 0;
}
double gmatrix::log_det() const{
    gmatrix A(size1,size2);
    A.memcpy(*this);
    gpermutation p(A.nrows());

    int signnum;
    gsl_linalg_LU_decomp(&A,&p,&signnum);
    return gsl_linalg_LU_lndet(&A); //log(det(A))
}

int gmatrix::svd(gsl_matrix& U, gsl_vector& S,gsl_matrix& V) const{
    gvector work(size1);
    gsl_matrix_memcpy(&U,this);
    return gsl_linalg_SV_decomp(&U,&V,&S,&work);
}
bool gmatrix::is_symmetric(double tol) const{
    if(!is_square()) return false;
    for(unsigned int i=0;i<size1-1;i++){
	for(unsigned int j=i+1;j<size2;j++){
	    double dif=(*this)(i,j)-(*this)(j,i);
	    if(fabs(dif>tol)) return false;
	}
    }
    return true;
}

int gpermutation::init(const size_t n,bool clean)
{
    if (n == 0)
    {
	gsl_error("permutation length n must be positive integer",__FILE__,
		  __LINE__,GSL_EDOM);
    }
    data = mymalloc(n,size_t);//(size_t *) malloc (n * sizeof (size_t));
    if(data == 0){
        gsl_error("failed to allocate space for permutation data",__FILE__,
		  __LINE__,GSL_EDOM);
    }
    size = n;
    if(clean){
	set_identity();
    }
    return 0;
}

int gpermutation::resize(const size_t n,bool clean)
{
    if(n==size){
	if(clean){
	    set_identity();
	}
	return 0;
    }

    if(size>0) myfree(data);
    return init(n,clean);
}
/*gpermutation& gpermutation::operator=(const gsl_permutation& other)
{
    if(this==&other) return *this;
    if(size!=other.size){
	cerr<<"When you copy assign permutations, please be sure that \n";
	cerr<<"the two permutations have the same sizes\n";
	//resize(other.size);
    }
    memcpy(other);
    return *this;
    }*/

void split(const string& s, Vstring& sv,char sep)
{//it is possible to be much fater with C's implmenetation strtok.
//It's now about at least three times slower than C's results
//We should maintain the code as it is to be user friendly and sacrifice the speed.
    int oldl,send,sbegin,l,r,b;
    sv.clear();
    sbegin=0;
    send=s.size();
    oldl=sbegin;
    for(b=sbegin;b<send;b++){
	if(s[b]==sep){
	    l=oldl;
	    r=b-1;
	    oldl=b+1;
	}else{
	    if(b+1!=send) {
		continue;
	    }else{
		l=oldl;
		r=b;
	    }
	}
	//remove the white space
	while(l<=r){
	    if(isspace(s[l])){
		l++;
	    }else break;
	}
	while(l<=r){
	    if(isspace(s[r])){
		r--;
	    }else break;
	}

	if(l>r){
	    sv.push_back("");
	}else{
	    sv.push_back(s.substr(l-sbegin,r-l+1));
	}

	if(b+1==send && s[b]==sep){ //how to deal with that b is the last one, we need to compute both
	    sv.push_back("");
	}
    }
    if(send==0) sv.push_back("");
    //if s is empty, sv should be "" of length 1.    
}

bool string2double(string& s, double &d)
{
    istringstream is(s);
    is>>d;
    if(!is) return false;
    if(!is.eof()) return false; //there is some extra leftoever;
    return true;
}

bool readrow(Vstring& v, Vdouble & D, bool & label, int nlab,bool clear)
//nlab=-1, searching, otherwise the first nlab elements
//are used for labels.
{
    if(clear) D.clear();
    double d;
    bool issearch=(nlab==-1);
    if(nlab== -1){
	if (string2double(v[0],d)){
	    label=false;
	    D.push_back(d);
	}else{
	    label=true;
	}
	nlab=1;
    }
    for(unsigned int i=nlab;i<v.size();i++){
	if(string2double(v[i],d)){
	    D.push_back(d);
	}else{
	    if(!issearch){
		char msgstr[902];
		snprintf(msgstr,900,"Error in reading field %d as %s is \
not a number.\n",i,v[i].c_str());
		throw(domain_error(msgstr));
	    }
	    return false;
	}
    }
    return true;
}
void error_msg(int lineN)
{
    char msgstr[902];
    snprintf(msgstr,900,"The data for line %d is in incorrect format.\n",lineN);
    throw(domain_error(msgstr));
}

void gmatrix_frame::transpose()
{
    //std::swap(rownames,colnames);//not work for vararry
    VAstring tmp(rownames);
    rownames.resize(colnames.size());
    rownames=colnames;
    colnames.resize(tmp.size());
    colnames=tmp;
    gmatrix::transpose();
}

void gmatrix_frame::set_split_char(char ch){
    split_char=ch;
}

void gmatrix_frame::cleanformat(Vdouble& TD, Vstring& Trownames, Vstring& Tcolnames)
{
    //check if we need to remove the empty first entry
    if(!Trownames.empty()){
	if(Trownames[0]=="" && TD.size()==0){
	    Trownames.erase(Trownames.begin());
	    size1--;
	    gsl_stream_printf("Warning!",__FILE__,__LINE__,"the empty first entry is removed when the data have no column data");
	}
	rownames.resize(Trownames.size());
	copy(Trownames.begin(),Trownames.end(),&rownames[0]);
    }
    
    if(!Tcolnames.empty()){
	if(Tcolnames[0]=="" && TD.size()==0){
	    Tcolnames.erase(Tcolnames.begin());
	    size2--;
	    gsl_stream_printf("Warning!",__FILE__,__LINE__,"the empty first entry is removed when the data have no row data");
	}
	colnames.resize(Tcolnames.size());
	copy(Tcolnames.begin(),Tcolnames.end(),&colnames[0]);
    }
    //resize assume the data dimension was allocated corretly.
    //however, our setting of size1 and size2 have not allocated enough space
    //we should use init, but init is private.
    double size1_new=size1;
    double size2_new=size2;
    size1=0;
    size2=0;
    gmatrix::resize(size1_new,size2_new);
    copy(TD.begin(),TD.end(),&data[0]);
}

istream& operator>>(istream& input, gmatrix_frame& T)
{
    //first find out if there are rows and columns
    //and compute n and p.
    //read the firt two lines when necessary
    //assuming the data has at least one row
    Vdouble D;
    Vstring rownames;
    Vstring colnames;
    int lineN=0;
    string line;
    Vstring lineV1;
    getline(input,line); lineN++;
    split(line,lineV1,T.split_char);
    line=""; //reset
    Vstring lineV2;
    getline(input,line); lineN++;
    split(line,lineV2,T.split_char);
    bool labeltmp;
    if(line==""){
	if(input){
	    error_msg(2);
	}else{
	    if(readrow(lineV1,D,labeltmp)){
		if(labeltmp){
		    rownames.push_back(lineV1[0]);
		}
		T.size1=1;
		T.size2=D.size();
		T.cleanformat(D,rownames,colnames);
		return input;
	    }else{//all are column names
		T.size1=0;
		T.size2=lineV1.size();
		T.cleanformat(D,rownames,lineV1);
		return input;
	    }
	}
    }

    T.size1=0;
    bool rowlabel=false;

    if(lineV2.size()==lineV1.size()+1){
	rowlabel=true;
	colnames=lineV1;
    }else if(lineV2.size()!=lineV1.size()){
	gsl_error("The number of the fields are  unequal among the first two lines.",__FILE__,__LINE__,GSL_EDOM);
    }else{//the two rows has the same number of fields.
	if(readrow(lineV1,D,rowlabel)){//success means no column labels
	    if(rowlabel){
		rownames.push_back(lineV1[0]);
	    }
	    T.size1++; 
	}else {//failure means that already column labels
	    if(readrow(lineV2,D,rowlabel)){
		D.clear(); //prepare for the repeat visit of line 2.
		colnames.assign(lineV1.begin()+rowlabel,lineV1.end());
		if(lineV1[0]!="" && rowlabel){
		    string msg="Warning! the nonempty first entry "+lineV1[0]+
			" for the data with rows and columns is iggored";
		    gsl_stream_printf("Warning",__FILE__,__LINE__,msg.c_str());
		    
		}
	    }else{
		error_msg(2);
	    }
	}
    }

    T.size2=lineV2.size()-rowlabel;
    while(true){
	if(!readrow(lineV2,D,labeltmp,rowlabel,false)){
	    error_msg(lineN);
	}
	if(rowlabel){
	    rownames.push_back(lineV2[0]);
	}
	T.size1++;
	
	if(!getline(input,line)) break;
	lineN++;
	split(line,lineV2,T.split_char);
	if(lineV2.size()!=T.size2+rowlabel){
	    error_msg(lineN);
	}
	if(line==""){
	    if(input) error_msg(lineN);//check if it is the last empty line
	    break;
	}
    }
    T.cleanformat(D,rownames,colnames);
    return input;
}


ostream& operator<<(ostream& output, const gmatrix_frame& T)
{
    output<<T.size1<<" row x "<<T.size2<<" column matrix\n";
    bool print_rownames = (T.rownames.size()!=0);
    bool print_colnames = (T.colnames.size()!=0);
    if(print_colnames){
	if(print_rownames) output<<"\t";
	output<<T.colnames[0];
	for(unsigned int j=1;j<T.size2;j++){
	    output<<"\t"<<T.colnames[j];
	}
	output<<"\n";
    }
    for(unsigned int i=0;i<T.size1;i++){
	if(print_rownames) output<<T.rownames[i]<<"\t";
	if(T.size2>0) output<<T(i,0);
	for(unsigned int j=1;j<T.size2;j++){
	    output<<"\t"<<T(i,j);
	};
	output<<"\n";
    };
    return output;
}


void gvector::key_sort(VAint & V)
{

    gpermutation p(size);
    sort_index(p);
    permute_by(p);
    //permute v by p;
    VAint Vc(V);
    for(unsigned int i=0;i<size;i++){
	V[i]=Vc[p[i]];
    }
}
istream& operator>>(istream& input, gpermutation& p)
{
    //assuming know the size
    if(p.size==0){
	gsl_stream_printf("Warning",__FILE__, __LINE__,
			  "The size of gpermutation is zero, nothing to be read");
	return input;
    } 
    for(unsigned int i=0;i<p.size;i++){
	size_t temp;
	input>>temp;
	p[i]=temp;
    }
    return input;
}

ostream& operator<<(ostream& output, gpermutation& p)
{
    //assuming know the size
    if(p.size==0){
	gsl_stream_printf("Warning",__FILE__,__LINE__,
	    "The size of gpermutation is zero, nothing to be written");
	return output;
    } 
    output<<p[0];
    for(unsigned int i=1;i<p.size;i++){
	output<<'\t'<<p[i];
    }
    return output;
}

istream& operator>>(istream& input, gvector& v)
{
    //assuming know the size
    if(v.size==0){
	gsl_stream_printf("Warning",__FILE__,__LINE__,
			  "The vector is zero, nothing to be read");
	return input;
    } 
    for(unsigned int i=0;i<v.size;i++){
	double temp;
	input>>temp;
	v[i]=temp;
    }
    return input;
}

ostream& operator<<(ostream& output, gvector& v)
{
    //assuming know the size
    if(v.size==0){
	gsl_stream_printf("Warning",__FILE__,__LINE__,
			  "The vector size is zero, nothing to be written");
	return output;
    } 
    output<<v[0];
    for(unsigned int i=1;i<v.size;i++){
	output<<'\t'<<v[i];
    }
    return output;
}

istream& operator>>(istream& input, gmatrix& M)
{
    //assuming know the size
    if(M.isempty()){
	gsl_stream_printf("Warning",__FILE__,__LINE__,
			  "The matrix is empty, nothing to be read");
	return input;
    } 
    for(unsigned int i=0;i<M.size1;i++){
	for(unsigned int j=0;j<M.size2;j++){
	    input>>M(i,j);
	}
    }
    return input;
}

ostream& operator<<(ostream& output, const gmatrix& M)
{
    //assuming we know the size
    if(M.isempty()){
	gsl_stream_printf("Warning",__FILE__,__LINE__,
			  "The matrix size is zero, nothing to be written");
	return output;
    } 
    for(unsigned int i=0;i<M.size1;i++){
	output<<M(i,0);
	for(unsigned int j=1;j<M.size2;j++){
	    output<<'\t'<<M(i,j);
	}
	output<<'\n';
    }
    return output;
}

double L2dist(const gvector& v1,const gvector & v2){
    double s=0;
    double dif;
    for(unsigned int i=0;i<v1.size;i++){
	dif=v1[i]-v2[i];
	s+=dif*dif;
    }
    return s;
}


/*int main()
{
    gmatrix_frame M;
    cin>>M;
    cout<<M;
    gvector v(10);
    for(int i=0;i<10;i++){
	v[i]=i;
    }
    cout<<v.mean()<<"\t"<<v.sd()<<endl;;
}

*/
