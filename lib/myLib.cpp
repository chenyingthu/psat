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
void main1(){
//	culaInitialize();
	double A[4]={1,2,3,4};
	double B[4]={5,6,7,8};
	double C[4]={11,12,13,14};
	double D[4]={17,18,19,20};
	double A1[4]={1,17,18,19};
	double B1[6]={2,3,4,20,21,22};
	double C1[6]={5,6,7,23,24,25};
	double D1[9]={8,9,10,11,12,13,14,15,16};
	double wr[5];
	double wi[5];
	double vl[25];
	double vr[25];
	double Result[25];
	cat4matrix(A1,2,2,C1,2,3,B1,3,2,D1,3,3,Result);
	for(int i=0;i<5;++i){
		for(int j=0;j<5;++j){
			printf("%lf\t",Result[j+i*5]);
		}
		printf("\n");
	}
	LAPACKE_dgeev(LAPACK_COL_MAJOR,'N','N',5,Result,5,wr,wi,vl,5,vr,5);
	for(int i=0;i<5;++i)
		printf("%lf\t%lf\n",wr[i],wi[i]);
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
