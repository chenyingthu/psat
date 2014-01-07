#include "myLib.h"
//#include<lapacke.h>
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
double norm(int n,double *x){
	double temp=0;
	for ( int i = 0; i < n; i += 1 ) {
		temp+=x[i]*x[i];
	}
	return sqrt(temp);
}
void MyDgemv(int Trans,int A_row,int A_col,double alpha,double *A,int lda,double *b,int inc_b,double beta,double *c,int inc_c){
	if(Trans==0)//No Trans
		culaDgemv('N',A_row,A_col,alpha,A,lda,b,inc_b,beta,c,inc_c);
	else if (Trans==1) //Trans
		culaDgemv('T',A_row,A_col,alpha,A,lda,b,inc_b,beta,c,inc_c);
	else
		printf("Wrong in MyDgemv\n");
}
void MyDgesv(int Major,int A_row,int nhis,double *A,int lda,int *ipiv,double *b,int ldb){
	if(Major==0)//Col Major
		culaDgesv(A_row,nhis,A,lda,ipiv,b,ldb);
	else if(Major==1)//Row Major
		culaDgesv(A_row,nhis,A,lda,ipiv,b,ldb);
	else
		printf("Wrong in MyDgesv\n");
}
void MyDgemm(int Major,int m,int k,int n,double alpha,double *A,int lda,double *B,int ldb,double beta,double *C,int ldc){
	if(Major==0)//Col Major
		culaDgemm('N','N',m,k,n,alpha,A,lda,B,ldb,beta,C,ldc);
	else if(Major==1)
		culaDgemm('N','N',m,k,n,alpha,A,lda,B,ldb,beta,C,ldc);
}
void MyDgeev(int Major,int A_row,double *A,int lda,double *wr,double *wi,double *vl,int n1,double *vr,int n2){
	if(Major==0)//Col Major
		culaDgeev('N','N',A_row,A,lda,wr,wi,vl,n1,vr,n2);
	else if(Major==1)//Row Major
		culaDgeev('N','N',A_row,A,lda,wr,wi,vl,n1,vr,n2);
	else
		printf("Worng in MyZgeev\n");
}
void MyZgemm(int Major,int m,int k,int n,Complex alpha,Complex *A,int lda,Complex *B,int ldb,Complex beta,Complex *C,int ldc){
	if(Major==0)//Col Major
		culaZgemm('N','N',m,k,n,alpha,A,lda,B,ldb,beta,C,ldc);
	else if(Major==1)
		culaZgemm('N','N',m,k,n,alpha,A,lda,B,ldb,beta,C,ldc);
}
void MyZgemv(int Trans,int A_row,int A_col,Complex alpha,Complex *A,int lda,Complex *b,int inc_b,Complex beta,Complex *c,int inc_c){
	if(Trans==0)//No Trans
		culaZgemv('N',A_row,A_col,alpha,A,lda,b,inc_b,beta,c,inc_c);
	else if (Trans==1) //Trans
		culaZgemv('T',A_row,A_col,alpha,A,lda,b,inc_b,beta,c,inc_c);
	else
		printf("Wrong in MyZgemv\n");
}


void MyDeviceDgemv(int Trans,int A_row,int A_col,double alpha,double *A,int lda,double *b,int inc_b,double beta,double *c,int inc_c){
	double *A_dev;
	cudaMalloc((void**)&A_dev, A_row*A_col * sizeof(double));
	double *b_dev;
	cudaMalloc((void**)&b_dev, A_col * sizeof(double));
	double *c_dev;
	cudaMalloc((void**)&c_dev, A_row * sizeof(double));
	cudaMemcpy(A_dev, A, A_row*A_col * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(b_dev, b, A_col * sizeof(double), cudaMemcpyHostToDevice);
	if(Trans==0)//No Trans
		culaDeviceDgemv('N',A_row,A_col,alpha,A_dev,lda,b_dev,inc_b,beta,c_dev,inc_c);
	else if (Trans==1) //Trans
		culaDeviceDgemv('T',A_row,A_col,alpha,A_dev,lda,b_dev,inc_b,beta,c_dev,inc_c);
	else
		printf("Wrong in MyDeviceDgemv\n");
	cudaMemcpy(c, c_dev, A_row * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(A_dev);
	cudaFree(b_dev);
	cudaFree(c_dev);
}
void MyDeviceDgesv(int Major,int A_row,int nhis,double *A,int lda,int *ipiv,double *b,int ldb){
	double *A_dev;
	cudaMalloc((void**)&A_dev, A_row*A_row * sizeof(double));
	double *b_dev;
	cudaMalloc((void**)&b_dev, A_row * sizeof(double));
	int *ipiv_dev;
	cudaMalloc((void**)&ipiv_dev, A_row * sizeof(double));
	cudaMemcpy(A_dev, A, A_row*A_row * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(b_dev, b, A_row * sizeof(double), cudaMemcpyHostToDevice);
	if(Major==0)//Col Major
		culaDeviceDgesv(A_row,nhis,A_dev,lda,ipiv_dev,b_dev,ldb);
	else if(Major==1)//Row Major
		culaDeviceDgesv(A_row,nhis,A_dev,lda,ipiv_dev,b_dev,ldb);
	else
		printf("Wrong in MyDeviceDgesv\n");
	cudaMemcpy(b, b_dev, A_row * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(A_dev);
	cudaFree(b_dev);
	cudaFree(ipiv_dev);
}
void MyDeviceDgemm(int Major,int m,int k,int n,double alpha,double *A,int lda,double *B,int ldb,double beta,double *C,int ldc){
	double *A_dev;
	cudaMalloc((void**)&A_dev, m*n * sizeof(double));
	double *B_dev;
	cudaMalloc((void**)&B_dev, n*k * sizeof(double));
	double *C_dev;
	cudaMalloc((void**)&C_dev, m*k * sizeof(double));
	cudaMemcpy(A_dev, A, m*n * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(B_dev, B,n*k * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(C_dev, C,n*k * sizeof(double), cudaMemcpyHostToDevice);
	if(Major==0)//Col Major
		culaDeviceDgemm('N','N',m,k,n,alpha,A_dev,lda,B_dev,ldb,beta,C_dev,ldc);
	else if(Major==1)
		culaDeviceDgemm('N','N',m,k,n,alpha,A_dev,lda,B_dev,ldb,beta,C_dev,ldc);
	cudaMemcpy(C, C_dev, m*k * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(A_dev);
	cudaFree(B_dev);
	cudaFree(C_dev);
}
void MyDeviceDgeev(int Major,int A_row,double *A,int lda,double *wr,double *wi,double *vl,int n1,double *vr,int n2){
	if(Major==0)//Col Major
		culaDeviceDgeev('N','N',A_row,A,lda,wr,wi,vl,n1,vr,n2);
	else if(Major==1)//Row Major
		culaDeviceDgeev('N','N',A_row,A,lda,wr,wi,vl,n1,vr,n2);
	else
		printf("Worng in MyDeviceZgeev\n");
}
void MyDeviceZgemm(int Major,int m,int k,int n,Complex alpha,Complex *A,int lda,Complex *B,int ldb,Complex beta,Complex *C,int ldc){
	Complex *A_dev;
	cudaMalloc((void**)&A_dev, m*n * sizeof(Complex));
	Complex *B_dev;
	cudaMalloc((void**)&B_dev, n*k * sizeof(Complex));
	Complex *C_dev;
	cudaMalloc((void**)&C_dev, m*k * sizeof(Complex));
	cudaMemcpy(A_dev, A, m*n * sizeof(Complex), cudaMemcpyHostToDevice);
	cudaMemcpy(B_dev, B,n*k * sizeof(Complex), cudaMemcpyHostToDevice);
	cudaMemcpy(C_dev, C,n*k * sizeof(Complex), cudaMemcpyHostToDevice);
	if(Major==0)//Col Major
		culaDeviceZgemm('N','N',m,k,n,alpha,A_dev,lda,B_dev,ldb,beta,C_dev,ldc);
	else if(Major==1)
		culaDeviceZgemm('N','N',m,k,n,alpha,A_dev,lda,B_dev,ldb,beta,C_dev,ldc);
	cudaMemcpy(C, C_dev, m*k * sizeof(Complex), cudaMemcpyDeviceToHost);
	cudaFree(A_dev);
	cudaFree(B_dev);
	cudaFree(C_dev);
}
void MyDeviceZgemv(int Trans,int A_row,int A_col,Complex alpha,Complex *A,int lda,Complex *b,int inc_b,Complex beta,Complex *c,int inc_c){
	Complex *A_dev;
	Complex *b_dev;
	Complex *c_dev;
	cudaMalloc((void**)&A_dev, A_row*A_col * sizeof(Complex));
	cudaMalloc((void**)&b_dev, A_col * sizeof(Complex));
	cudaMalloc((void**)&c_dev, A_row * sizeof(Complex));
	cudaMemcpy(A_dev, A, A_row*A_col * sizeof(Complex), cudaMemcpyHostToDevice);
	cudaMemcpy(b_dev, b, A_col * sizeof(Complex), cudaMemcpyHostToDevice);
	if(Trans==0)//No Trans
		culaDeviceZgemv('N',A_row,A_col,alpha,A_dev,lda,b_dev,inc_b,beta,c_dev,inc_c);
	else if (Trans==1) //Trans
		culaDeviceZgemv('T',A_row,A_col,alpha,A_dev,lda,b_dev,inc_b,beta,c_dev,inc_c);
	else
		printf("Wrong in MyDeviceZgemv\n");
	cudaMemcpy(c, c_dev, A_row * sizeof(Complex), cudaMemcpyDeviceToHost);
	cudaFree(A_dev);
	cudaFree(b_dev);
	cudaFree(c_dev);
}
