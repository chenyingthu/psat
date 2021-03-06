#include "myLib.h"
#include<lapacke.h>
int find(double *srcArray,int len,double ref,int *dstArray)
{
	int num=0;
	for (int i=0;i<len;++i)
		if(srcArray[i]>ref+tol1)
			dstArray[num++]=i;
	return num;
}
int findi(int *srcArray,int len,int ref,int *dstArray)
{
	int num=0;
	for (int i=0;i<len;++i)
		if(srcArray[i]>ref)
			dstArray[num++]=i;
	return num;
}
complex<double> conj(complex<double>a)
{
	a=complex<double>(a.real(),-a.imag());
	return a;
}
complex<double> *conj_(int n,complex<double>*a)
{
	complex<double> *b=new complex<double> [n];
	for (int i=0;i<n;++i)
		b[i]=complex<double>(a[i].real(),-a[i].imag());
	return b;
}
void cat4matrix(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result)
{
	int n=na+nc;
	int m=ma+mb;
	for (int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			if(i<na&&j<ma){
				Result[j+i*m]=A[i+j*na];
			}
			else if(i<nb&&j>=ma){
				Result[j+i*m]=B[i+(j-ma)*nb];
			}
			else if(i>=na&&j<mc){
				Result[j+i*m]=C[i-na+j*nc];
			}
			else if(i>=ma&&j>=mc){
				Result[j+i*m]=D[i-ma+(j-mc)*nd];
			}
		}
	}
}
void cat4matrix_Ac(double *A,int na,int ma,double *B,int nb,int mb,double *C,int nc,int mc,double *D,int nd,int md,double *Result,double h)
{
	int n=na+nc;
	int m=ma+mb;
	for (int i=0;i<n;++i){
		for(int j=0;j<m;++j){
			if(i<na&&j<ma){
				Result[j+i*m]=-h*A[i+j*na];
				if(i==j)
				  Result[j+i*m]=1+Result[j+i*m];
				else
				  Result[j+i*m]=Result[j+i*m];
			}
			else if(i<nb&&j>=ma){
				Result[j+i*m]=B[i+(j-ma)*nb];
			}
			else if(i>=na&&j<mc){
				Result[j+i*m]=-h*C[i-na+j*nc];
			}
			else if(i>=ma&&j>=mc){
				Result[j+i*m]=D[i-ma+(j-mc)*nd];
			}
		}
	}
}
int cmp(const void *a, const void *b)
{
     return(*(int *)a-*(int *)b);
}
int* sort(int *s,int n)
{
	int *a=new int [n];
	for(int i=0;i<n;++i)
		a[i]=s[i];
	qsort(a,n,sizeof(a[0]),cmp);
	return a;
}
int cmpD(const void * a, const void * b)
{
	return((*(double*)a-*(double*)b>0)?1:-1);
}
double* sortD(double *s,int n)
{
	double *a=new double [n];
	for(int i=0;i<n;++i)
		a[i]=s[i];
	qsort(a,n,sizeof(a[0]),cmpD);
	return a;
}
