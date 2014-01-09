#include "./inc/Psat.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include "sm_20_atomic_functions.h"
#include "sm_13_double_functions.h"
#define atomic_double
#define pi 3.141592653589793
void MyMemcpy(double *dst,double *src,int n,int HTD){
	if(HTD==1)
		cudaMemcpy(dst,src,n*sizeof(double),cudaMemcpyHostToDevice);
	else if(HTD==0)
		cudaMemcpy(dst,src,n*sizeof(double),cudaMemcpyDeviceToHost);
	else 
		printf("wrong in MyMemcpy\n");
}
__global__ void MyMemcpyD2D_kernel(double *dst,double *src,int n,double c){
	int tid=threadIdx.x;
	while(tid<n){
		dst[tid]=src[tid]*c;
		tid+=1024;
	}
}
void MyMemcpyD2D(double *dst,double *src,int n,double c){
	 MyMemcpyD2D_kernel<<<1,1024>>>(dst,src,n,c);
}
#ifdef atomic_double
__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, 
                        __double_as_longlong(val + 
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif
void Psat::formDAEX_gpu(){
	int k=0;
	for ( int i = 0; i < dae.n; i += 1 ) {
		dae.X[k++]=dae.x[i];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.X[k++]=dae.a[i];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.X[k++]=dae.V[i];
	}
}
void Psat::DAEXtox_gpu(){
	int k=0;
	for ( int i = 0; i < dae.n; i += 1 ) {
		dae.x[i]=dae.X[k++];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.a[i]=dae.X[k++];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.V[i]=dae.X[k++];
	}
}
void Psat::DAEXtox_2_gpu(double *x_in){
	int k=0;
	for ( int i = 0; i < dae.n; i += 1 ) {
		dae.x[i]=x_in[k++];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.a[i]=x_in[k++];
	}
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.V[i]=x_in[k++];
	}
}

void Psat::dyn_f_store_gpu(int iFlag){
	for ( int iStep = 0; iStep < simu.nSteps; iStep += 1 ) {
		int iStart=(iStep)*(dae.n+2*bus.n);
		int iEnd=(iStep+1)*(dae.n+2*bus.n);
		record.t[record.nhis]=simu.t_cur+(iStep+1)*simu.tStep;
		int k=0;
		for ( int j = iStart; j < iEnd; j += 1 ) {
			record.x[record.nhis*(dae.n+2*bus.n)+k++]=dae.X[j];
		}
		iStart=(iStep)*dae.n;
		iEnd=(iStep+1)*dae.n;
		k=0;
		for ( int j = iStart; j < iEnd; j += 1 ) {
			record.f[record.nhis*(dae.n)+k++]=dae.f[j];
		}
		record.nhis++;
	}
}
void Psat::dyn_f_integration_gpu(int iFlag){
	if (simu.multiSteps==1){
		simu.nSteps=settings.dyn_MulStep_nSteps;
		simu.multiSteps=2;
	}
	if(settings.dyn_isPredict==1)
		dyn_f_prediction_gpu(iFlag);
	if(iFlag==1){
		dae.X=solver_jfng_gpu(dae.X);
		debug((char*)"daeX",dae.n+2*bus.n,dae.X);
	}
}
void Psat::dyn_f_prediction_gpu(int iFlag){//to do multsteps,...
	int ord=settings.dyn_predict_model;
	double *x=new double [dae.n+2*bus.n]; 
	for(int i = 0;i<dae.n+2*bus.n;++i){
		x[i]=dae.X[i];
	}
	int isPre=0;
	if(ord==1){
		if(record.t[record.nhis-1]>fault.tFaultEnd+settings.dyn_tStep*3){
			isPre=1;
			for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
				x[i]=2*record.x[i+(record.nhis-1)*(dae.n+2*bus.n)]-record.x[i+(record.nhis-2)*(dae.n+2*bus.n)];
			}
		}
	}
	else if(ord==2){
		if(record.t[record.nhis-1]>fault.tFaultEnd+settings.dyn_tStep*8){
			isPre=1;
			for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
				x[i]=3*record.x[i+(record.nhis-1)*(dae.n+2*bus.n)]-3*record.x[i+(record.nhis-2)*(dae.n+2*bus.n)]+record.x[i+(record.nhis-3)*(dae.n+2*bus.n)];
			}
		}
	}
	else if(ord==3){
		if(record.t[record.nhis-1]>fault.tFaultEnd+settings.dyn_tStep*12){
			isPre=1;
			for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
				x[i]=4*record.x[i+(record.nhis-1)*(dae.n+2*bus.n)]-6*record.x[i+(record.nhis-2)*(dae.n+2*bus.n)]+4*record.x[i+(record.nhis-3)*(dae.n+2*bus.n)]-record.x[i+(record.nhis-4)*(dae.n+2*bus.n)];
			}
		}
	}
	if(isPre==0){
		simu.nSteps=1;
		simu.multiSteps=0;
	}
	for(int i = 0;i<dae.n+2*bus.n;++i){
		dae.X[i]=x[i];
	}
	delete []x;
}
double * Psat::solver_jfng_gpu(double *x){
	int n=dae.n+2*bus.n;
	double *sol=new double [n];
	double *sol_dev;
	cudaMalloc((void**)&sol_dev, n * sizeof(double));
	double gamma=0.9;
	// double ierr=0;
	int maxit=solver.jfng.newton.maxit;
	
	double rat;
	double stop_tol=solver.jfng.newton.tol;
	solver.jfng.gmres.tol=solver.jfng.gmres.tol_fixed;
	double *f0=new double [n];
	double *f0_dev;
	cudaMalloc((void**)&f0_dev, n * sizeof(double));
	double *fold=new double [n];
	double *fold_dev;
	cudaMalloc((void**)&f0_dev, n * sizeof(double));
	int itc=0;
	f0=dyn_f_dae_gpu(x,simu.t_next);
	cudaMemcpy(f0_dev,f0,n*sizeof(double),cudaMemcpyHostToDevice);
	double *f0_temp;
	cudaMalloc((void**)&f0_temp, n * sizeof(double));
	cudaMemcpy(f0_temp,f0_dev,n*sizeof(double),cudaMemcpyDeviceToHost);
	double fnrm=norm(n,f0)/sqrt((double)n); 
	double fnrmo=1;
	simu.neval=0;
	while(fnrm > stop_tol&& itc < maxit){
		rat=fnrm/fnrmo;
		fnrmo=fnrm;
		itc=itc+1;
		if(solver.jfng.newton.method==1&&(solver.jfng.precond.inner==1||solver.jfng.precond.outer==1))
		{
			updatePreconditioner_gpu(1);
		}
		if(solver.jfng.newton.method==1){
			pre_gmres_gpu(f0,x,step);
			pre_gmres(f0,x,step);
		}
		for ( int i = 0; i < n; i += 1 ) {
			fold[i]=f0[i];
		}
		for ( int i = 0; i < n; i += 1 ) {
			x[i]+=step[i];
		}
		f0=dyn_f_dae_gpu(x,simu.t_next);
		fnrm=norm(n,f0)/sqrt((double)n); 
		rat=fnrm/fnrmo;
		if(solver.jfng.precond.outer==1 && solver.jfng.precond.outerStop == 2 && solver.update.num < solver.update.maxnum){
			int i=solver.update.num;
			for(int j=0;j<n;++j){
				solver.deltx[j+i*n]=step[j];
				solver.delty[j+i*n]=f0[j]-fold[j];
			}
			solver.update.num++;
		}
		if(solver.jfng.gmres.tol>0){
			double etaold=solver.jfng.gmres.tol;
			double etanew=gamma*rat*rat;
			if(gamma*etaold*etaold>0.1)
				etanew=max(etanew,gamma*etaold*etaold);
			solver.jfng.gmres.tol=min(etanew,solver.jfng.gmres.tol);
			solver.jfng.gmres.tol=max(solver.jfng.gmres.tol,0.5*stop_tol/fnrm);
		}
	}
	for ( int i = 0; i < n; i += 1 ) {
		sol[i]=x[i];
	}
	delete []fold;
	delete []f0;
	return sol;
}
double * Psat::dyn_f_dae_gpu(double *x_in,double t0)//todo multiSteps
{
	int pos=-1;
	int n=dae.n+2*bus.n;
	double *x_rec=new double [n];
	double *f_rec=new double [dae.n];
	double *f_out=new double [n];
	DAEXtox_2_gpu(x_in);
	if(settings.dyn_lf!=1){
		for ( int i = 0; i < record.nhis; i += 1 ) {
			if(abs(record.t[i]-t0)<1e-8)
				pos=i;
		}
		if(pos==-1){
			printf("no such time moment in his record");
			return NULL;
		}
		for ( int i = 0; i < n; i += 1 ) {
			x_rec[i]=record.x[i+pos*n];
			if(i<dae.n){
				f_rec[i]=record.f[i+pos*dae.n];
			}
		}
	}
	// for ( int i = 0; i < bus.n; i += 1 ) {
	//   printf("daeV\t%.16lf\n",dae.V[i]);
	//   printf("daeV\t%.16lf\n",dae.X[i+dae.n+bus.n]);
	// }
	// getchar();
	fm_lf_1();
	fm_mn_1();
	fm_syn(1);
	fm_sw_1();
	for ( int i = 0; i < bus.n; i += 1 ) {
		dae.g[i]=dae.gp[i];
		dae.g[i+bus.n]=dae.gq[i];
	}
	for ( int i = 0; i < dae.n; i += 1 ) {
		dae.f[i]=0;
	}
	fm_syn(3);
	for(int i=0;i<n;++i)
		f_out[i]=0;
	if(settings.dyn_lf!=1){
		for ( int i = 0; i < dae.n; i += 1 ) {
			f_out[i]=x_in[i]-x_rec[i]-simu.tStep*0.5*(dae.f[i]+f_rec[i]);
		}
	}
	for ( int i = dae.n; i < n; i += 1 ) {
		f_out[i]=dae.g[i-dae.n];
	}
	delete []x_rec;
	delete []f_rec;
	return f_out;
}
__global__ void updatePreconditioner_kernel1(int n,int j,double dest,double *p,double *q,double *deltx,double *delty){
	int tid=threadIdx.x;
	
	while(tid<n){
		p[tid]=dest*deltx[tid+j*n];
		q[tid]=delty[tid+j*n];
		tid+=1024;
	}
	__syncthreads();
}
__global__ void updatePreconditioner_kernel2(int n,int j,double dest,double *p,double *q,double *deltx,double *delty){
	__syncthreads();
#ifdef atomic_double
	__shared__  double temp[1];
	double a;
#else
	__shared__  float temp[1];
	float a;
#endif
	temp[0]=0;
	__syncthreads();
	int tid=threadIdx.x;

	while(tid<n){
		a=deltx[tid+j*n]*q[tid];
		atomicAdd( &temp[0],a);
		tid+=1024;
	}
	
	__syncthreads();
	//printf("%d\n",count);
	//printf("%d\n",count);
	tid=threadIdx.x;
	//if(tid==0)
	//	printf("0 %f\n",temp);
	//else if(tid==1)
	//	printf("1 %f\n",temp);
	while(tid<n){
		q[tid]+=temp[0]*delty[tid+j*n];
		tid+=1024;
	}
	__syncthreads();
}
__global__ void updatePreconditioner_kernel3(int n,int j,double dest,double *p,double *q,double *deltx,double *delty){

#ifdef atomic_double
	__shared__  double temp[1];
	double a;
#else
	__shared__  float temp[1];
	float a;
#endif
	temp[0]=0;
	__syncthreads();
	int tid=threadIdx.x;
	//printf("%d %d\n",j,tid);
	while(tid<n){
		a=p[tid]*q[tid];
		atomicAdd( &temp[0],a);
		tid+=1024;
	}
	__syncthreads();
	////__syncthreads();
	////printf("temp\t%f\n",temp);
	tid=threadIdx.x;
	
	//printf("%d\n",tid);
	while(tid<n){
		deltx[tid+j*n]=p[tid]/temp[0];
		delty[tid+j*n]=p[tid]-q[tid];
		//printf("%lf\n",deltx[tid+j*n]);
		tid+=1024;
	}
	__syncthreads();
}
void Psat::updatePreconditioner_gpu(int iFlag){
	int n=dae.n+2*bus.n;
	double *p=new double[n]; 
	double *q=new double [n];
	double *p_dev;
	double *q_dev;
	cudaMalloc((void**)&p_dev,n* sizeof(double));
	cudaMalloc((void**)&q_dev,n* sizeof(double));
	cudaMemset(p_dev,0,n*sizeof(double));
	cudaMemset(q_dev,0,n*sizeof(double));
	double temp=0;
	cudaMemcpy(solver.deltx_dev,solver.deltx,n*solver.update.maxnum * sizeof(double),cudaMemcpyHostToDevice);
	cudaMemcpy(solver.delty_dev,solver.delty,n*solver.update.maxnum * sizeof(double),cudaMemcpyHostToDevice);
	while(solver.update.last_num+1<solver.update.num){
		updatePreconditioner_kernel1<<<1,n>>>(n,solver.update.last_num,solver.update.dest,p_dev,q_dev,solver.deltx_dev,solver.delty_dev);
		for ( int j = 0; j < solver.update.last_num; j += 1 ) {
			updatePreconditioner_kernel2<<<1,n>>>(n,j,solver.update.dest,p_dev,q_dev,solver.deltx_dev,solver.delty_dev);
		}
		updatePreconditioner_kernel3<<<1,n>>>(n,solver.update.last_num,solver.update.dest,p_dev,q_dev,solver.deltx_dev,solver.delty_dev);
		solver.update.last_num++;
	}//end of while
	cudaMemcpy(solver.deltx,solver.deltx_dev,n*solver.update.maxnum * sizeof(double),cudaMemcpyDeviceToHost);
	cudaMemcpy(solver.delty,solver.delty_dev,n*solver.update.maxnum * sizeof(double),cudaMemcpyDeviceToHost);
	delete []p; 
	delete []q;
}
__global__ void pre_gmres_kernel1(double *g,double rho){
	int tid=threadIdx.x;
	if(tid==0)
		g[tid]=rho;
}
void Psat::pre_gmres_gpu(double *f0,double *xc,double *x){
	double errtol=solver.jfng.gmres.tol;
	int kmax=solver.jfng.gmres.maxit;
	int reorth=solver.jfng.gmres.reorth;
	int n=(dae.n+2*bus.n)*simu.nSteps;
	double rho;
	double *xc_dev;
	cudaMalloc((void**)&xc_dev,n* sizeof(double));
	double *x_dev;
	cudaMalloc((void**)&x_dev,n* sizeof(double));
	double *f0_dev;
	cudaMalloc((void**)&x_dev,n* sizeof(double));
	double *r_dev;
	cudaMalloc((void**)&r_dev,n* sizeof(double));
	double *z_dev;
	cudaMalloc((void**)&z_dev,n*kmax* sizeof(double));
	double *v_dev;
	cudaMalloc((void**)&v_dev,n*kmax* sizeof(double));
	double *h_dev;
	cudaMalloc((void**)&h_dev,kmax*kmax* sizeof(double));
	double *g_dev;
	cudaMalloc((void**)&g_dev,(kmax+1)* sizeof(double));
	double *s_dev;
	cudaMalloc((void**)&s_dev,(kmax+1)* sizeof(double));
	double *c_dev;
	cudaMalloc((void**)&c_dev,(kmax+1)* sizeof(double));
	cudaMemset(g_dev,0,(kmax+1)* sizeof(double));
	cudaMemset(s_dev,0,(kmax+1)* sizeof(double));
	cudaMemset(c_dev,0,(kmax+1)* sizeof(double));
	cudaMemset(x_dev,0,n* sizeof(double));
	double *r=new double [n];
	//double *test=new double [kmax+1];
	double *testn=new double [n];
	MyMemcpy(r_dev,f0,n,1);
	MyMemcpyD2D(r_dev,r_dev,n,-1.0);
	MyMemcpy(r,r_dev,n,0);
	
	rho=norm(n,r);
	//cudaMemset(g_dev,1,(kmax+1)*sizeof(double));
	pre_gmres_kernel1<<<1,1>>>(g_dev,rho);
	//MyMemcpy(test,g_dev,kmax+1,0);
	//debug("g",kmax+1,test);
	errtol=errtol*rho;
	errstep=rho;
	if(rho<errtol)
		return;
	MyMemcpyD2D(v_dev,r_dev,n,1/rho);
	//MyMemcpy(testn,v_dev+n-1,1,0);
	//debug("test n",1,testn);
	int k=0;
	//while(0){
	while((rho > errtol)&&(k<kmax)){
		k++;
		int kk=k-1;
		MyMemcpyD2D(z_dev+kk*n,v_dev+kk*n,n,1.0);
		if(solver.jfng.precond.inner==1||solver.jfng.precond.outer==1){
			precondition_gpu(z_dev+kk*n);
		}
		//MyMemcpy(testn,z_dev+kk*n,n,0);
		//debug("z_dev",n,testn);
	}
//	while((rho > errtol)&&(k<kmax)){
//		for ( int i = 0; i < n; i += 1 ) {
//			z_temp[i]=solver.jfng.gmres.z[i+kk*n];
//		}
//		dirder(xc,z_temp,f0,v_temp);
//		if(solver.jfng.precond.inner==1 && solver.jfng.precond.innerStop == 2 && solver.update.num < solver.update.maxnum){
//			int i=solver.update.num;
//			for(int j=0;j<n;++j){
//				solver.deltx[j+i*n]=solver.jfng.gmres.z[j+kk*n];
//				solver.delty[j+i*n]=v_temp[j];
//			}
//			solver.update.num++;
//		}
//		double normav=norm(n,v_temp);
//
//		for ( int j = 0; j < k; j += 1 ) {
//			double temp=0;
//			for (int i=0;i<n;++i)
//				temp+=solver.jfng.gmres.v[i+j*n]*v_temp[i];
//			solver.jfng.gmres.h[kk+j*kmax]=temp;
//			for (int i=0;i<n;++i)
//				v_temp[i]=v_temp[i]-temp*solver.jfng.gmres.v[i+j*n];
//		}
//		solver.jfng.gmres.h[kk+(kk+1)*kmax]=norm(n,v_temp);
//		double normav2=solver.jfng.gmres.h[kk+(kk+1)*kmax];
//		if((reorth==1&&abs(normav+0.001*normav2-normav)<1e-8)||reorth==3)
//		{
//			for ( int j = 0; j < k; j += 1 ) {
//				double hr=0;
//				for (int i=0;i<n;++i)
//					hr+=solver.jfng.gmres.v[i+j*n]*v_temp[i];
//				solver.jfng.gmres.h[kk+j*kmax]+=hr;
//				for (int i=0;i<n;++i)
//					v_temp[i]=v_temp[i]-hr*solver.jfng.gmres.v[i+j*n];
//			}
//			solver.jfng.gmres.h[kk+(kk+1)*kmax]=norm(n,v_temp);
//		}
//		if(solver.jfng.gmres.h[kk+(kk+1)*kmax]!=0){
//			for ( int i = 0; i < n; i += 1 ) {
//				v_temp[i]=v_temp[i]/solver.jfng.gmres.h[kk+(kk+1)*kmax];
//			}
//		}
//		if(k>1)
//			givapp(solver.jfng.gmres.c,solver.jfng.gmres.s,solver.jfng.gmres.h,k-1);
//		double nu=sqrt(solver.jfng.gmres.h[kk+(kk+1)*kmax]*solver.jfng.gmres.h[kk+(kk+1)*kmax]+solver.jfng.gmres.h[kk+(kk)*kmax]*solver.jfng.gmres.h[kk+(kk)*kmax]);
//		if(nu!=0)
//		{
//			solver.jfng.gmres.c[kk]=solver.jfng.gmres.h[kk+kk*kmax]/nu;
//			solver.jfng.gmres.s[kk]=-solver.jfng.gmres.h[kk+(kk+1)*kmax]/nu;
//			solver.jfng.gmres.h[kk+kk*kmax]=solver.jfng.gmres.c[kk]*solver.jfng.gmres.h[kk+kk*kmax]-solver.jfng.gmres.s[kk]*solver.jfng.gmres.h[kk+(kk+1)*kmax];
//			solver.jfng.gmres.h[kk+(kk+1)*kmax]=0;
//			double w1=solver.jfng.gmres.c[kk]*g[kk]-solver.jfng.gmres.s[kk]*g[kk+1];
//			double w2=solver.jfng.gmres.s[kk]*g[kk]+solver.jfng.gmres.c[kk]*g[kk+1];
//			g[kk]=w1;
//			g[kk+1]=w2;
//		}
//		rho=abs(g[kk+1]);
//		for (int i=0;i<n;++i)
//			solver.jfng.gmres.v[i+(kk+1)*n]=v_temp[i];
//	}
//	double *h_temp=new double [k*k];
//	double *z_temp2=new double [n*k];
//	int *ipiv=new int[k];
//	for ( int i = 0; i < k; i += 1 ) {
//		for ( int j = 0; j < k; j += 1 ) {
//			h_temp[i+j*k]=solver.jfng.gmres.h[j+i*kmax];
//		}
//	}
//	for ( int i = 0; i < k; i += 1 ) {
//		for ( int j = 0; j < n; j += 1 ) {
//			z_temp2[j+i*n]=solver.jfng.gmres.z[j+i*n];
//		}
//	}
//	double alpha=1;
//	double beta=0;
//	
//#ifndef gpu
//	MyDgesv(0,k,1,h_temp,k,ipiv,g,k);
//#else
//	MyDeviceDgesv(0,k,1,h_temp,k,ipiv,g,k);
//#endif
//#ifndef gpu
//	MyDgemv(0,n,k,alpha,z_temp2,n,g,1,beta,x,1);
//#else
//	MyDeviceDgemv(0,n,k,alpha,z_temp2,n,g,1,beta,x,1);
//#endif
//	// debug("x",n,x);
//	inner_it_count=k; 
//	simu.neval+=k;
//	delete []r;
//	delete []g;
//	delete []v_temp;
//	delete []z_temp;
//	delete []h_temp;
//	delete []ipiv;
//	delete []z_temp2;
} 
void Psat::givapp_gpu(double *c,double *s,double *vin,int k){
	int kmax=solver.jfng.gmres.maxit;
	for ( int i = 0; i < k; i += 1 ) {
		double w1=c[i]*vin[k+i*kmax]-s[i]*vin[k+(i+1)*kmax];
		double w2=s[i]*vin[k+i*kmax]+c[i]*vin[k+(i+1)*kmax];
		vin[k+i*kmax]=w1;
		vin[k+(i+1)*kmax]=w2;
		// printf("c\t%lf\ts\t%lf\n",c[i],s[i]);
		// printf("w1\t%lf\tw2\t%lf\n",w1,w2);
	}
}
__global__ void precondition_kernel(int n,int j,double *x,double *deltx,double *delty){

#ifdef atomic_double
	__shared__  double temp[1];
	double a;
#else
	__shared__  float temp[1];
	float a;
#endif
	int tid=threadIdx.x;
	temp[0]=0;
	__syncthreads();
	while(tid<n){
		a=deltx[tid+j*n]*x[tid];
		atomicAdd( &temp[0],a);
		tid+=1024;
	}
	__syncthreads();
	tid=threadIdx.x;
	while(tid<n){
		x[tid]+=temp[0]*delty[tid+j*n];
		tid+=1024;
	}
	__syncthreads();
}
void Psat::precondition_gpu(double *x){//todo multiSteps
	int kmax=solver.jfng.gmres.maxit;
	int n=dae.n+2*bus.n;
	MyMemcpyD2D(x,x,n,(double)solver.update.dest);
	for ( int i = 0; i <solver.update.last_num; i += 1 ) {
		precondition_kernel<<<1,n>>>(n,i,x,solver.deltx_dev,solver.delty_dev);
	}
}
void Psat::dirder_gpu(double *x,double *w,double *f0,double *z){
	int n=dae.n+2*bus.n;
	double temp=norm(n,w);
	double epsnew=solver.jfng.findiff;
	if(temp<1e-10)
	{
		for ( int i = 0; i < n; i += 1 ) {
			z[i]=0;
		}
	}
	epsnew=epsnew/temp;
	temp=norm(n,x);
	if(temp>0)
		epsnew=epsnew*temp;
	double *del=new double [n];
	for ( int i = 0; i < n; i += 1 ) {
		del[i]=x[i]+epsnew*w[i];
	}
	double *f1=new double [n];
	f1=dyn_f_dae(del,simu.t_next);
	for ( int i = 0; i < n; i += 1 ) {
		z[i]=(f1[i]-f0[i])/epsnew;
	}
	delete []del;
	delete []f1;
}
void Psat::dyn_f_increaseTimeSteps_gpu(int iFlag){
	switch(simu.converged){
	case 1:
		if(simu.newton_iteration>=15)
			simu.tStep=max(simu.tStep*0.9,settings.dyn_tStep_min);
		if(simu.newton_iteration<=10)
			simu.tStep=min(simu.tStep*1.3,settings.dyn_tStep_max);
		if(settings.fixt)
			simu.tStep=settings.dyn_tStep;
		break;
	case 0:
		simu.tStep=settings.dyn_tStep*0.5;
		if(simu.tStep<settings.dyn_tStep_min)
			simu.tStep=settings.dyn_tStep_min;
		break;
	}
	simu.t_cur=simu.t_next;
	simu.t_next=simu.t_cur+simu.nSteps*simu.tStep;
	double tempo_min=simu.t_next;
	for ( int i = 0; i < 4*fault.n; i += 1 ) {
		if(simu.t_switch[i]-simu.t_cur>1e-6&&simu.t_next-simu.t_switch[i]>1e-6)
			if(tempo_min>simu.t_switch[i])
				tempo_min=simu.t_switch[i];
	}
	simu.t_next=tempo_min;
	simu.tStep=simu.t_next-simu.t_cur;
	// for ( int i = 0; i < 4*fault.n; i += 1 ) {
	//   if(simu.t_next<simu.t_switch[4*fault.n-1]&&simu.t_switch[i]<simu.t_next){
	//     printf("I am here!!!\n");
	//     printf("%lf\n%lf\n",simu.t_next,simu.t_switch[i]);
	//     simu.t_next=simu.t_switch[i];
	//     simu.tStep=simu.t_next-simu.t_cur;
	//     simu.t_switch[i]=simu.t_switch[4*fault.n-1]*2;
	//     break;
	//   }
	// }
	if(simu.t_next>fault.tFaultStart+simu.tStep){
		solver.jfng.precond.inner=1;
		solver.jfng.precond.outer=1;
		if(simu.neval>0&&simu.neval<8){
			solver.jfng.precond.innerStop=1;
			solver.jfng.precond.outerStop=1;
			if(settings.dyn_isMulStep==1)
				if(simu.multiSteps==0)
					simu.multiSteps=1;
		}
		else if (simu.neval>7){
			solver.jfng.precond.innerStop=2;
			solver.jfng.precond.outerStop=2;
		}
	}
}
void Psat::dyn_f_dealFaults_gpu(int iFalg){
	formDAEX_gpu();
	for ( int i = 0; i < fault.n; i += 1 ) {
		int h=fault.bus[i];
		if(abs(simu.t_next-fault.con[i][4])<1e-8){
			for ( int j = 0; j < bus.n; j += 1 ) {
				V_bak[j]=dae.V[j];
			}
			printf("applying fault at t = %lf s\n",simu.t_next);
			for (int j=0;j<bus.n;++j){

				dae.V[j]=0.6;
			}
			dae.V[h]=0.01;
			for(int j=0;j<syn.n;++j){
				int k=syn.bus[j];
				dae.V[k]=dae.X[dae.n+k+bus.n];
			}
			shunt.g[h]=fault.dat[2+i*5]+fault.dat[0+i*5];
			shunt.b[h]=fault.dat[3+i*5]+fault.dat[1+i*5];
			formDAEX_gpu();
			fm_y();
			dyn_f_iniSolver(2);
			settings.dyn_lf=1;
			dae.X=solver_jfng_gpu(dae.X);
			settings.dyn_lf=2;
		}// fault intervention
		else if(abs(simu.t_next-fault.con[i][5])<1e-8){
			printf("Clearing fault at t = %lf s\n",simu.t_next);
			//getchar();
			shunt.g[h]=fault.dat[2+i*5];
			shunt.b[h]=fault.dat[3+i*5];
			fm_y();
			printf("LY done\n");
			dyn_f_iniSolver(2);
			solver.update.dest=1;
			settings.dyn_lf=1;
			for ( int i = 0; i < bus.n; i += 1 ) {
				dae.V[i]=V_bak[i];
			}
			if(syn.n>0){
				double mean_delta=0;
				for ( int i = 0; i < syn.n; i += 1 ) {
					int k=syn.delta_idx[i];
					mean_delta+=dae.x[k];
				}
				mean_delta=mean_delta/syn.n;
				for ( int i = 0; i < bus.n-boundarynode.n; i += 1 ) {
					int k=boundarynode.indexG[i];
					dae.a[k]=mean_delta-fault.delta+fault.ang[k];
				}
			}
			else{
				for ( int i = 0; i < bus.n; i += 1 ) {
					dae.a[i]=fault.ang[i];
				}
			}
			formDAEX_gpu();
			dae.X=solver_jfng_gpu(dae.X);
			settings.dyn_lf=2;
		}//end of else if

	}//end of for
	formDAEX_gpu();
}
