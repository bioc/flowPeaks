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

#ifndef _GVECTOR_GMATRIX_H
#define _GVECTOR_GMATRIX_H

/* Please note that for copy assignment or other i/o, fread,fwritef, input<<, the user takes the responbility to make sure the sizes are correct. The only exception are the views, which are not deep copied, or the inout<< for gmatrix_frame, which figures out the dimnensions by reading the whole file into memory and then convert it to a gmatrix_frame */

#include <iostream>
#include <string>
#include <vector>
#include <valarray>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include "quick_op.h"
#include "memory_handle.h"
using namespace std;

typedef vector<int> Vint; 
typedef vector<char> Vchar; 
typedef vector<string> Vstring; 
typedef vector<double> Vdouble; 
typedef vector<bool> Vbool; 

typedef valarray<int> VAint; 
typedef valarray<char> VAchar; 
typedef valarray<string> VAstring; 
typedef valarray<double> VAdouble; 
typedef valarray<bool> VAbool; 

CBLAS_TRANSPOSE_t get_transposeid(bool transpose);
/*the definition of gsl_vector
typedef struct 
{
  size_t size;
  size_t stride;
  double *data;
  gsl_block *block;
  int owner;
} 
gsl_vector;

typedef struct
{
  gsl_vector vector;
} _gsl_vector_view;
*/
  /*
typedef struct 
{
  size_t size1;
  size_t size2;
  size_t tda;
  double * data;
  gsl_block * block;
  int owner;
} gsl_matrix;

typedef struct
{
  gsl_matrix matrix;
} _gsl_matrix_view;
*/

/*struct gsl_permutation_struct
{
  size_t size;
  size_t *data;
  };*/

class gpermutation: public gsl_permutation{
public:
    gpermutation(const size_t n,bool clean=true){
	init(n,clean);}
    ~gpermutation(){
	if(size>0) free(this->data);};

    int next(){
	return gsl_permutation_next(this);};

    int prev(){
	return gsl_permutation_prev(this);};

    void reverse(){
	gsl_permutation_reverse(this);};

    int inverse(gsl_permutation & inv){
	return gsl_permutation_inverse(&inv,this);}

    bool isempty(){ return size>0;}

    int memcpy(const gsl_permutation &other){
	return gsl_permutation_memcpy(this,&other);}
    void ran_shuffle(gsl_rng* r){
	set_identity();
	gsl_ran_shuffle(r,data,size,sizeof(size_t));}

    int resize(const size_t n,bool clean=true);
    bool isvalid(){
	return  gsl_permutation_valid(this);};

    void set_identity(){
	gsl_permutation_init(this);}

    int swap(size_t i, size_t j){
	return gsl_permutation_swap(this,i,j);}

    size_t& operator[](size_t i){
	return data[i];}

    const size_t& operator[](size_t i) const{
	return data[i];}


    friend istream& operator>>(istream& input, gpermutation& p);
    friend ostream& operator<<(ostream& output, const gpermutation& p);
 private:
    int init(const size_t n,bool clean=true);
    gpermutation& operator=(const gsl_permutation& other);
};



class gvector : public gsl_vector{
public:
    //section 8.3.1
    gvector(){size=0; stride=1;block=0;data=0;owner=1;}
    gvector(const size_t n, const bool clean=true){
	init(n,clean);}
    ~gvector(){if (size>0 && owner==1) gsl_block_free(block);}

    int resize(const size_t n, bool clean=true);

    //section 8.3.2 elements access
    double&  operator[](const size_t i){ 
	return *gsl_vector_ptr(this,i);}

    const double &operator[](const size_t i) const { 
	return *gsl_vector_const_ptr(this,i);}

    double& operator()(const size_t i){ 
	return *gsl_vector_ptr(this,i);}

    const double& operator()(const size_t i) const {
	return *gsl_vector_const_ptr(this,i);}

    //section 8.3.3
    void set_all(double x){
	gsl_vector_set_all (this,x);};

    void set_zero(){
	gsl_vector_set_zero (this);}

    int set_basis (const size_t i) {
	return gsl_vector_set_basis (this,i);}

    //section 8.3.4.i/o
    int fwrite (FILE * stream) const {
	return gsl_vector_fwrite (stream, this);}

    int fread (FILE * stream) {
	return gsl_vector_fread (stream, this);}

    int fprintf (FILE * stream, const char * format) const {
	return gsl_vector_fprintf (stream, this,format) ;}

    int fscanf (FILE * stream)  {
	return gsl_vector_fscanf (stream, this); }

    //section 8.3.6
    int memcpy(const gsl_vector &other){
	return gsl_vector_memcpy(this,&other);} 

    int swap(gsl_vector & other){
	return gsl_vector_swap(this,&other);}
	
    //section 8.3.7
    int swap(const size_t i,const size_t j){
	return gsl_vector_swap_elements(this,i,j);}
    int reverse(){
	return gsl_vector_reverse(this);}
	
    //section 8.3.8. operations
    int operator+=(double x) {
	return gsl_vector_add_constant (this,x);}
    int operator+=(const gsl_vector &other){
	return gsl_vector_add(this,&other);}

    int operator-=(double x) {
	return gsl_vector_add_constant (this,-x);}
    int operator-=(const gsl_vector &other){
	return gsl_vector_sub(this,&other);}
    

    int operator*=(double x) {
	return gsl_vector_scale (this, x);}
    int operator*=(const gsl_vector &other) {
	return gsl_vector_mul(this,&other);}

    int operator/=(const gsl_vector &other) {
	return gsl_vector_div (this, &other);}
    int operator/=(double x) {
	return gsl_vector_scale (this, 1/x);};
    
    // some fast compuations, assuming the slide is one
    // by default, all of these functions are inlined
    void memcpy_fast(const double *b){
	doublecopy(b,data,size);
    }
    void memcpy_fast(const double b){
	doublecopy(data,b,size);
    }
    void sub_fast(const double *b){
	doublesub(data,b,size);
    }
    void mul_fast(const double b){
	doublemul(data,b,size);
    }
    void add_fast(const double *b){
	doubleadd(data,b,size);
    }
    double dot_fast(const double *b) const{
	return doubledot(data,b,size);
    }
    void mv_fast(const double* A, const double* b, const int p){
        doublemv(A,b,size,p,data);
    }

    double max() const{
	return gsl_vector_max (this);}
    double min() const{
	return gsl_vector_min (this);}
    size_t max_index(){return gsl_vector_max_index (this);}
    size_t min_index(){return gsl_vector_min_index (this);}

    
    bool isempty() const {return size==0;} //isnull does not make sense
    bool ispos() const {return gsl_vector_ispos (this);}
    bool isneg() const {return gsl_vector_isneg (this);}
    bool isnoneg() const {return gsl_vector_isnonneg (this);}

    //more operators for convenience
    gvector& operator-();

    //convienet functions
    int permute_by(const gsl_permutation& other){
	return gsl_permute_vector(&other,this);}
    int permute_inverse_by(const gsl_permutation& other){
	return gsl_permute_vector_inverse(&other,this);}
    double dot(const gsl_vector& other) const{
	double d;
	gsl_blas_ddot(this,&other,&d);
	return d;}
    double L2norm() const{
	return gsl_blas_dnrm2(this);};
    double mean() const{
	return gsl_stats_mean(data,stride,size);}
    double sd() const{
	return gsl_stats_sd(data,stride,size);};
    double var() const{
	return gsl_stats_variance(data,stride,size);};
    double cov(gsl_vector& other) const{
	return gsl_stats_covariance (data, stride,other.data,other.stride,size);}
    void sort(){
	gsl_sort_vector(this);}
    int sort_index(gsl_permutation&other) const{
	return gsl_sort_vector_index(&other,this);}
    void key_sort(VAint & V);

    double quantile_from_sorted(double f){
	return gsl_stats_quantile_from_sorted_data(data,stride,size,f);}

    double sum(){
	double s=0;
	for(unsigned int i=0;i<size;i++){
	    s+=data[i];
	}
	return s;}

    double abssum(){
	return gsl_blas_dasum(this);}
    int axpy(double a,const gsl_vector& other){
	return gsl_blas_daxpy(a,&other,this); //*this+=a*other
    }
    int mv(const gsl_matrix& A, const gsl_vector& other, 
	   const bool TransA=false,double a=1,double b=0){
	return gsl_blas_dgemv(get_transposeid(TransA),a,&A,&other,b,this);}// *this=a*op(A)*other+b*(*this);

    int solve(const gsl_matrix &other,const gsl_vector &b);

    friend ostream& operator<<(ostream& output, const gvector& v);
    friend istream& operator>>(istream& input, gvector& v);
private:
    int init(const size_t n, bool clean=true);
    //forbid the copy assignment, use the memcpy instead
    gvector& operator=(const gsl_vector& other);

};

//for this class, owner is always 0.
class gvector_view : public gvector{
public:
    gvector_view(){owner=0;};
    //section 8.3.5. 
    /*subvector(size_t offset,size_t n,size_t stride=1);*/
    /*also works when other if the type of other is gvector and gvector_view*/
    /*base also works with the &base[0], where base is a valarray*/
    gvector_view(const gsl_vector& other, size_t n=0,size_t offset=0,size_t stride=1){
	if(n==0) n=other.size;
	change_view(other,offset,n,stride);};
    gvector_view(const double*base,size_t n, size_t stride=1){
	change_view(base,n,stride);}
    gvector_view(const gsl_vector_view&other){
	change_view(other.vector);
    }
    gvector_view(const gsl_vector_const_view&other){
	change_view(other.vector);
    }

    void change_view(const gsl_vector& other, size_t n=0,size_t offset=0,size_t stride=1){
	if(n==0) n=other.size;
	assign(gsl_vector_const_subvector_with_stride(&other,offset,stride,n).vector);}
    void change_view(const double*base,size_t n=0, size_t stride=1){
	assign(gsl_vector_const_view_array_with_stride(base,stride,n).vector);}
    
   //for view, no copy of the data, but just the information
    gvector_view& operator=(const gsl_vector& other){
	if(this==&other) return *this;
	if(owner==1){
	    cerr<<"Can not assign gvector unless it's a a gvector_view\n";
	    exit(0);
	}
	assign(other);
	return *this;} //no deep copy
//implmentation details
private:
    void assign(const gsl_vector & other); //no copy, set owner=0;
};


//define matrix and matrix view and some blas functions.
class gmatrix : public gsl_matrix{
public:
    gmatrix(){size1=0;size2=0;tda=0;data=0;block=0;owner=1;}
    gmatrix(size_t n1,size_t n2, bool clean=true){
	init(n1,n2,clean);}
    ~gmatrix(){if(size1>0 && size2>0 &&owner==1) gsl_block_free(block);}

    size_t nrows() const {return size1;}
    size_t ncols() const {return size2;}
    int resize(const size_t n1, const size_t n2, bool clean=true);

    //section 8.4.2
    const double &operator()( size_t row, size_t col) const {
	return *gsl_matrix_const_ptr(this,row, col) ;}
    double &operator()( size_t row, size_t col ) {
	return *gsl_matrix_ptr(this, row, col) ;}

    //section 8.4.3
    void set_all (double  x) {
	gsl_matrix_set_all(this,x);}
    void set_zero() {
	gsl_matrix_set_zero(this);}
    void set_identity(){
	gsl_matrix_set_identity(this);}

    //section 8.4.4
    int fwrite (FILE * stream) const {
	return gsl_matrix_fwrite (stream, this);}
    int fread (FILE * stream) {
	return gsl_matrix_fread (stream, this);}
    int fprintf (FILE * stream, const char * format) const {
	return gsl_matrix_fprintf (stream, this,format) ;}
    int fscanf (FILE * stream)  {
	return gsl_matrix_fscanf (stream, this); }

    friend ostream& operator<<(ostream& output, const gmatrix& M);
    friend istream& operator>>(istream& input, gmatrix & M);
    
    //section 8.4.6.
    //[] for row, () for column
    gvector_view operator()(size_t j){
	return gsl_matrix_column(this,j);}
    gvector_view operator[](size_t i){
	return gsl_matrix_row(this,i);}
    gvector_view diagonal(){
	return gsl_matrix_diagonal(this);}
    gvector_view subdiagonal(size_t k){
	return gsl_matrix_subdiagonal(this,k);}
    gvector_view supdiagonal(size_t k){
	return gsl_matrix_superdiagonal(this,k);}

    const gvector_view operator()(size_t j) const{
	return gsl_matrix_const_column(this,j);}
    const gvector_view operator[](size_t i) const{
	return gsl_matrix_const_row(this,i);}
    const gvector_view diagonal() const{
	return gsl_matrix_const_diagonal(this);}
    const gvector_view subdiagonal(size_t k)const{
	return gsl_matrix_const_subdiagonal(this,k);}
    const gvector_view supdiagonal(size_t k) const{
	return gsl_matrix_const_superdiagonal(this,k);}


    //section 8.4.7
    int memcpy(const gsl_matrix& other){
	return gsl_matrix_memcpy(this,&other);};
    int swap(gsl_matrix& other){
	return gsl_matrix_swap(this,&other);};
    int get_row(gsl_vector& v, size_t i){
	return gsl_matrix_get_row(&v,this,i);};
    int get_col(gsl_vector& v, size_t j){
	return gsl_matrix_get_col(&v,this,j);};
    int set_row(gsl_vector &v, size_t i){
	return gsl_matrix_set_row(this,i,&v);}
    int set_col(gsl_vector &v, size_t j){
	return gsl_matrix_set_col(this,j,&v);}
    

    //section 8.4.9
    int swap_rows(size_t i, size_t j){
	return gsl_matrix_swap_rows(this,i,j);}

    int swap_cols(size_t i, size_t j){
	return gsl_matrix_swap_columns(this,i,j);}

    int swap_rowcol(size_t i, size_t j){
	return gsl_matrix_swap_rowcol(this,i,j);}

    int transpose(); /*{
    return gsl_matrix_transpose(this);}, the standard do required both dimension are the same*/

    int transpose_memcpy(const gsl_matrix &other){
	return gsl_matrix_transpose_memcpy(this,&other);}

    //section 8.4.10
    int operator+=(double x){
	return gsl_matrix_add_constant(this,x);}
    int operator+=(const gsl_matrix &other){
	return gsl_matrix_add(this,&other);}

    int operator-=(double x){
	return gsl_matrix_add_constant(this,-x);}
    int operator-=(const gsl_matrix &other){
	return gsl_matrix_sub(this,&other);}

    int operator*=(double x){
	return gsl_matrix_scale(this,x);}
    int operator*=(const gsl_matrix &other){
	return gsl_matrix_mul_elements(this,&other);}

    int operator/=(double x){
	return gsl_matrix_scale(this,1/x);}
    int operator/=(const gsl_matrix &other){
	return gsl_matrix_div_elements(this,&other);}

    double max(){
	return gsl_matrix_max(this);}
    double min(){
	return gsl_matrix_min(this);}

    void max_index(size_t& i, size_t&j){
	return gsl_matrix_max_index(this,&i,&j);}

    void min_index(size_t& i, size_t&j){
	return gsl_matrix_min_index(this,&i,&j);}

    //section 8.4.12
    bool isempty() const{
	return size1==0||size2==0;}
    bool ispos() const {
	return gsl_matrix_ispos(this);}
    bool isneg() const {
	return gsl_matrix_isneg(this);}
    bool isnonneg() const {
	return gsl_matrix_isnonneg(this);}

//more operators for convenience
    gmatrix& operator-(); //change the sign

    // (*this)=aAB+b(*this)
    int mm(const gsl_matrix & A, const gsl_matrix &B, bool TransA=false, bool TransB=false,double a=1, double b=0){
	return gsl_blas_dgemm(get_transposeid(TransA),get_transposeid(TransB),
			      a,&A,&B,b,this);}
    
    int inverse(const gsl_matrix& A); /*this=A^{-1}*/
    int half(const gsl_matrix &A); //*this=A^{/2}*/
    int svd(gsl_matrix& U, gsl_vector& S,gsl_matrix& V) const;//singular valude dcomposation;
    bool is_square() const { return (size1==size2);};
    bool is_symmetric(double tol=1e-8) const ;
    double log_det() const;
private:
    int init(const size_t n1, const size_t n2, bool clean=true);
    //forbidden the deep copy, use memcpy instead
    gmatrix& operator=(const gsl_matrix& other); 
};

//for this class, owner is always zero
class gmatrix_view: public gmatrix{
public:
    gmatrix_view(){owner=0;}
    gmatrix_view(const gsl_matrix& other){
	change_view(other,0,0,other.size1,other.size2);}
    gmatrix_view(const gsl_matrix& other,size_t k1, size_t k2, size_t n1, size_t n2){
	change_view(other,k1,k2,n1,n2);}
    gmatrix_view(const double*base,size_t n1, size_t n2,size_t tda=0){
	if(tda==0) tda=n2;
	change_view(base,n1,n2,tda);}
    void change_view(const gsl_matrix& other,size_t k1, size_t k2, size_t n1, size_t n2){
	assign(gsl_matrix_const_submatrix(&other,k1,k2,n1,n2).matrix);}
    void change_view(const double*base,size_t n1, size_t n2,size_t tda=0){
	if(tda==0) tda=n2;
	assign(gsl_matrix_const_view_array_with_tda(base,n1,n2,tda).matrix);}

    //no deep copy
    gmatrix_view& operator=(const gsl_matrix& other){
	if(this==&other) return *this;
	if(owner==1){
	    cerr<<"Can not assign gmatrix unless it's a a gmatrix_view\n";
	    exit(0);
	}
	assign(other);
	return *this;}
//implmentation details
private:
    void assign(const gsl_matrix& other);
};

//data with rownames and colnames
class gmatrix_frame: public gmatrix{
public:
    VAstring rownames;
    VAstring colnames;
    void set_split_char(char ch);
    gmatrix_frame(){split_char='\t';};
    ~gmatrix_frame(){};
    
    void transpose();

    //intelligent read, automatically figured out the number of rows and columns
    //and the rownames, colnames.
    friend ostream& operator<<(ostream& output, const gmatrix_frame& T);

    friend istream& operator>>(istream& input, gmatrix_frame& T);

private:
    void cleanformat(Vdouble& TD, Vstring& Trownames, Vstring& Tcolnames);
    char split_char;
};

//for random number generators
class RNG{
public:
    gsl_rng* r;
    RNG(const gsl_rng_type* T=gsl_rng_ranlxs1){
	r=gsl_rng_alloc(T);
	//set a fixed seed for reproducible data
	gsl_rng_set(r,0);
    };
    void set_seed(unsigned long int s){
	gsl_rng_set(r,s);
    }
    ~RNG(){if(r!=0) gsl_rng_free(r);}
};
extern RNG g_rng;

//some useful helper functions
void split(const string& s, Vstring& sv,char sep='\t');
bool string2double(string &s, double &d);
bool readrow(Vstring& v, Vdouble & D, bool & label, int nlab=-1,bool clear=true);
double L2dist(const gvector& v1,const gvector & v2);
inline void updateC12(int nc,double &C1, double &C2)
{
    C2=nc/(nc+1.0);
    C1=GSL_POSINF;
    if(nc>1){
	C1=nc/(nc-1.0);
    }
}
#endif
