#ifndef _MYLIB_H
#define _MYLIB_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex>
#include <lapacke.h>
#include <cblas.h>
using namespace std;
typedef  complex<double> Complex;
#define tol1 0.0000005
int find(double *srcArray,int len,double ref,int *dstArray);
int findi(int *srcArray,int len,int ref,int *dstArray);
complex<double> conj(complex<double> a);
complex<double>* conj_(int n,complex<double> *a);
void cat4matrix(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result);
void cat4matrix_Ac(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result,double h);
int* sort(int *a,int n);
double *sortD (double *a,int n);
double norm(int n,double *x);
void MyDgemv(int Trans,int A_row,int A_col,double alpha,double *A,int lda,double *b,int inc_b,double beta,double *c,int inc_c);
void MyDgesv(int Major,int A_row,int nhis,double *A,int lda,int *ipiv,double *b,int ldb);
void MyDgemm(int Major,int m,int k,int n,double alpha,double *A,int lda,double *B,int ldb,double beta,double *C,int ldc);
void MyDgeev(int Major,int A_row,double *A,int lda,double *wr,double *wi,double *vl,int n1,double *vr,int n2);
void MyZgemm(int Major,int m,int k,int n,Complex alpha,Complex *A,int lda,Complex *B,int ldb,Complex beta,Complex *C,int ldc);
void MyZgemv(int Trans,int A_row,int A_col,Complex alpha,Complex *A,int lda,Complex *b,int inc_b,Complex beta,Complex *c,int inc_c);
#endif
