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


#ifndef QUICK_OP__
#define QUICK_OP__
/*This file provides the quick operations, the goal is for speed
  not for safety*/

inline double L2dist(const double*a, const double *b, const int n)
{
    double s=0;
    for(int i=0;i<n;i++){
	double dif=a[i]-b[i];
	s+=dif*dif;
    }
    return s;
}
inline double L2dist(const double*a, const int n)
{
    double s=0;
    for(int i=0;i<n;i++){
	s+=a[i]*a[i];
    }
    return s;
}

inline double doubledot(const double*a, const double *b, const int n)
{
    double s=0;
    for(int i=0;i<n;i++){
	s+=a[i]*b[i];
    }
    return s;
}
inline void doublecopy(const double *a, double * b, const int n)
{
    for(int i=0;i<n;i++){
	b[i]=a[i];
    }
}
inline void doublecopy(double *a, const double b, const int n)
{
    for(int i=0;i<n;i++){
	a[i]=b;
    }
}


inline void doubleweightsum(double *a, const double*b, const int n, double w)
{
    //a=(1-w)*a+w*b or a[i]+=(b-a)*w
    for(int i=0;i<n;i++){
	a[i]+=(b[i]-a[i])*w;
    }
}

inline void doubleadd(double *a, const double b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]+=b;
    }
}
inline void doubleadd(double *a, const double* b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]+=b[i];
    }
}

inline void doublemul(double *a, const double b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]*=b;
    }
}


inline void doublemul(double *a, const double* b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]*=b[i];
    }
}
inline void doublediv(double *a, const double* b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]/=b[i];
    }
}
inline void doublesub(double *a, const double* b, const int n)
{

    for(int i=0;i<n;i++){
	a[i]-=b[i];
    }
}
inline double doublesum(const double *a, const int n)
{
    double s=0;
    for(int i=0;i<n;i++){
	s+=a[i];
    }
    return s;
}

inline void doublemv(const double *A, const double *b, const int n, const int p, double *x)
{
    for(int i=0;i<n;i++){
	x[i]=doubledot(A+i*p,b,p);
    }
}

inline void doubleapply(double*a, const int n, double func(double))
{
    for(int i=0;i<n;i++){
	a[i]=func(a[i]);
    }
}
#endif
