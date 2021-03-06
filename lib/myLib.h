#ifndef _MYLIB_H
#define _MYLIB_H
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <complex>
using namespace std;
#define tol1 0.0000005
int find(double *srcArray,int len,double ref,int *dstArray);
int findi(int *srcArray,int len,int ref,int *dstArray);
complex<double> conj(complex<double> a);
complex<double>* conj_(int n,complex<double> *a);
void cat4matrix(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result);
void cat4matrix_Ac(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result,double h);
int* sort(int *a,int n);
double *sortD (double *a,int n);
#endif
