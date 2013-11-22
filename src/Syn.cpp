#include "../inc/Syn.h"
Syn::Syn()
{
	n=0;
	is2_n=0;
	is3_n=0;
	is4_n=0;
	is51_n=0;
	is52_n=0;
	is53_n=0;
	is6_n=0;
	is8_n=0;
	jay=complex<double>(0,1);
}
Syn::~Syn(){

}
void Syn::synDelete()
{
	if(n==0)
		return;
	for (int i=0;i<n;++i)
		delete []con[i];
	delete []con;
	delete []bus;
	delete []Pg;
	delete []Qg;
	delete []J11;
	delete []J12;
	delete []J21;
	delete []J22;
	delete []delta_idx;
	delete []omega_idx;
	delete []e1q_idx;
	delete []e1d_idx;
	delete []e2q_idx;
	delete []e2d_idx;
	delete []psiq_idx;
	delete []psid_idx;
	delete []delta;
	delete []omega;
	delete []e1q;
	delete []e1d;
	delete []e2q;
	delete []e2d;
	delete []psiq;
	delete []psid;
	delete []Gp;
	delete []Gq;
	delete []vf;
	delete []pm;
	delete []ord;
	delete []is2;
	delete []is3;
	delete []is4;
	delete []is51;
	delete []is52;
	delete []is53;
	delete []c1;
	delete []c2;
	delete []c3;

	delete []is6;
	delete []is8;
	delete []bs2;
	delete []bs3;
	delete []bs4;
	delete []bs51;
	delete []bs52;
	delete []bs53;
	delete []bs6;
	delete []bs8;

}
void Syn::init()
{
	if(n==0)
		return;
	con=new double*[n];
	for (int i=0;i<n;i++)
	{
		con[i]=new double[24];
	}
	for (int i=0;i<n;++i)
		for (int j=0;j<24;++j)
			con[i][j]=0;
	Pg=new double [n];
	Id=new double [n];
	Iq=new double [n];
	Qg=new double [n];
	J11=new double [n];
	J12=new double [n];
	J21=new double [n];
	J22=new double [n];
	delta_idx=new int [n];
	omega_idx=new int [n];
	e1q_idx=new int [n];
	e1d_idx=new int [n];
	e2q_idx=new int [n];
	e2d_idx=new int [n];
	psiq_idx=new int [n];
	psid_idx=new int [n];
	delta=new double [n];
	omega=new double [n];
	e1q=new double [n];
	e1d=new double [n];
	e2q=new double [n];
	e2d=new double [n];
	psiq=new double [n];
	psid=new double [n];
	vf=new double [n];
	pm=new double [n];
	ord=new int [n];
	is2=new int [n];
	is3=new int [n];
	is4=new int [n];
	is51=new int [n];
	is52=new int [n];
	is53=new int [n];
	is6=new int [n];
	is8=new int [n];
	c1=new double [n];
	c2=new double [n];
	c3=new double [n];

	bs2=new int [n];
	bs3=new int [n];
	bs4=new int [n];
	bs51=new int [n];
	bs52=new int [n];
	bs53=new int [n];
	bs6=new int [n];
	bs8=new int [n];
	Gp=new double [8*n];
	Gq=new double [8*n];
	for (int i=0;i<n;++i)
	{
		Pg[i]=0;
		Qg[i]=0;
		J11[i]=0;
		J12[i]=0;
		J21[i]=0;
		J22[i]=0;
	}
	for (int i=0;i<8*n;++i)
	{
		Gp[i]=0;
		Gq[i]=0;
	}
}
int Syn::check(Bus _bus)
{
	if(n==0)
		return 1;
	bus=new int [n];
	for (int i=0;i<n;++i){
		bus[i]=_bus.internel[(int)con[i][0]-1];
		//printf("%d\n",bus[i]);
	}
	return 1;
}
int Syn::fm_dynidx(DAE dae)
{
	int syn_ord;
	for (int i=0;i<n;++i){
		delta_idx[i]=dae.n+1-1;//-1 for the array start with 0 not 1,so do this
		omega_idx[i]=dae.n+2-1;
		syn_ord=(int)con[i][4]*10;//because 5.1 can't use switch,so X10 to make it int.
		switch(syn_ord){
		case 20:
			dae.n+=2;
			break;
		case 30:
			e1q_idx[i]=dae.n+3-1;
			dae.n+=3;
			break;
		case 40:
			e1q_idx[i]=dae.n+3-1;
			e1d_idx[i]=dae.n+4-1;
			dae.n+=4;
			break;
		case 51:
			e1q_idx[i]=dae.n+3-1;
			e1d_idx[i]=dae.n+4-1;
			e2d_idx[i]=dae.n+5-1;
			dae.n+=5;
			break;
		case 52:
			e1q_idx[i]=dae.n+3-1;
			e2q_idx[i]=dae.n+4-1;
			e2d_idx[i]=dae.n+5-1;
			dae.n+=5;
			break;
		case 53:
			e1q_idx[i]=dae.n+3-1;
			psid_idx[i]=dae.n+4-1;
			psiq_idx[i]=dae.n+5-1;
			dae.n+=5;
			break;
		case 60:
			e1q_idx[i]=dae.n+3-1;
			e1d_idx[i]=dae.n+4-1;
			e2q_idx[i]=dae.n+5-1;
			e2d_idx[i]=dae.n+6-1;
			dae.n+=6;
			break;
		case 80:
			e1q_idx[i]=dae.n+3-1;
			e1d_idx[i]=dae.n+4-1;
			e2q_idx[i]=dae.n+5-1;
			e2d_idx[i]=dae.n+6-1;
			psid_idx[i]=dae.n+7-1;
			psiq_idx[i]=dae.n+8-1;
			dae.n+=8;
			break;
		}
		vf[i]=2;
		pm[i]=1;
		//printf("%lf\t%lf\n",e1q[i],e1d[i]);
	}// end of for
	return dae.n;
}
void Syn::fm_synit(Settings settings,DAE dae,Bus _bus){
	double *ra=new double [n];
	double *xl=new double [n];
	double *xd=new double [n];
	double *xq=new double [n];
	double *xd1=new double [n];
	double *xq1=new double [n];
	double *xd2=new double [n];
	double *xq2=new double [n];
	double *Td10=new double [n];
	double *Tq10=new double [n];
	double *Td20=new double [n];
	double *Tq20=new double [n];
	double *K0=new double [n];
	double *Kp=new double [n];
	double *Taa=new double [n];
	double *Vg=new double [n];
	double *ag=new double [n];
	complex<double> *V=new complex<double> [n];
	complex<double> *S=new complex<double> [n];
	complex<double> *I=new complex<double> [n];
	complex<double> *Vdq=new complex<double>[n];
	complex<double> *Idq=new complex<double>[n];
	double *vd=new double [n];
	double *vq=new double [n];
	double *_delta=new double [n];
	if (n==0)
		return;
	int *idx=new int [n];
	int idx_n=0;
	for (int i=0;i<n;++i){
		ord[i]=(int)con[i][4]*10;
		switch(ord[i]){
		case 20:
			is2[is2_n]=i;
			bs2[is2_n++]=bus[i];
			break;
		case 30:
			is3[is3_n]=i;
			bs3[is3_n++]=bus[i];
			break;
		case 40:
			is4[is4_n]=i;
			bs4[is4_n++]=bus[i];
			break;
		case 51:
			is51[is51_n]=i;
			bs51[is51_n++]=bus[i];
			break;
		case 52:
			bs52[is52_n]=bus[i];
			is52[is52_n++]=i;
			break;
		case 53:
			bs53[is53_n]=bus[i];
			is53[is53_n++]=i;
			break;
		case 60:
			bs6[is6_n]=bus[i];
			is6[is6_n++]=i;
			break;
		case 80:
			bs8[is8_n]=bus[i];
			is8[is8_n++]=i;
			break;
		}//end of switch
	}// end of for
	//for (int i=0;i<is4_n;++i)
	//	printf("%d\n",bs4[i]);
	for (int i=0;i<n;++i){
		if(i<is2_n){
			con[is2[i]][19]=0;
			con[is2[i]][20]=0;
		}
		if(i<is53_n){
			con[is53[i]][19]=0;
			con[is53[i]][20]=0;
		}
		if(i<is8_n){
			con[is8[i]][19]=0;
			con[is8[i]][20]=0;
		}
	}
	for(int i=0;i<n;++i){
		if(con[i][17]<=0){
			if (settings.conv)
				con[i][17]=con[i][17]*con[i][1]/settings.mva;
			printf("Inertia cannot be <=0.M=10 will be used\n");
		}
		if(con[i][8]<=0){
			con[i][8]=0.302;
			if (settings.conv)
				con[i][8]=con[i][8]/con[i][1]*settings.mva;
			printf("X'd cannot be <=0.X'd=0.302 will be used\n");
		}
		if(ord[i]!=20){
			if(con[i][7]<=0){
				con[i][7]=1.90;
				if(settings.conv)
					con[i][7]=con[i][7]/con[i][1]*settings.mva;
				printf("xd cannot be <=0.xd=1.90 will be used\n");
			}
			if(con[i][10]<=0){
				con[i][10]=8.00;
				printf("T'd0 cannot be <=0.T d0=8.00 will be used\n");
			}
		}
		if(ord[i]==52||ord[i]==60||ord[i]==80){
			//todo
		}
		if(ord[i]==40||ord[i]==51||ord[i]==60||ord[i]==80){
			if(con[i][15]<=0){
				con[i][15]=0.80;
				printf("T q0 cannot be <=0.T q0=0.80 will be used\n");
			}
		}
		if(ord[i]!=20&&ord[i]!=52){
			if(con[i][12]<=0){
				con[i][12]=1.70;
				if(settings.conv)
					con[i][12]=con[i][12]/con[i][1]*settings.mva;
				printf("xq cannot be <=0.xq=1.70 will be used\n");
			}
		}
	}
	int *synbus=new int [n];
	synbus=sort(bus,n);
	int n_old=-1;//this code section is a little complex.....
	int n_new=0;
	for (int i=0;i<n;++i){
		n_new=synbus[i];
		if(n_new!=n_old){
			idx_n=0;
			for(int j=0;j<n;++j)
				if(bus[j]==n_new)
					idx[idx_n++]=j;
			if(idx_n==1){
				if(con[idx[0]][21]!=1){
					printf("Warning:active power ratio of generator #%d must be 1\n",idx[0]);
					con[idx[0]][21]=1;
				}
				if(con[idx[0]][22]!=1){
					printf("Warning:Reactive power ratio of generator #%d must be 1\n",idx[0]);
					con[idx[0]][22]=1;
				}
			}
			else if(idx_n>1){
				double ratiop=0;
				double ratioq=0;
				for (int j=0;j<idx_n;++j){
					ratiop+=con[idx[j]][21];
					ratioq+=con[idx[j]][22];
				}
				if(abs(ratiop-1)>1e-5){
					printf("warning:the sum of active power ratios of generators  must be 1");
					for (int j=0;j<idx_n;++j){
						con[idx[j]][21]=1/idx_n;
					}
				}
				else
					con[idx[0]][21]=con[idx[0]][21]-ratiop+1;
				if(abs(ratioq-1)>1e-5){
					printf("warning:the sum of reactive power ratios of generators  must be 1");
					for (int j=0;j<idx_n;++j){
						con[idx[j]][22]=1/idx_n;
					}
				}
				else
					con[idx[0]][22]=con[idx[0]][22]-ratiop+1;

			}//end of else if
			n_old=n_new;
		}
	}
	for (int i=0;i<n;++i){
		ra[i]=con[i][6];
		xl[i]=con[i][5];
		xd[i]=con[i][7];
		xq[i]=con[i][12];
		xd1[i]=con[i][8];
		xq1[i]=con[i][13];
		xd2[i]=con[i][9];
		xq2[i]=con[i][14];
		Td10[i]=con[i][10];
		Tq10[i]=con[i][15];
		Td20[i]=con[i][11];
		Tq20[i]=con[i][16];
		K0[i]=con[i][19];
		Kp[i]=con[i][20];
		Taa[i]=con[i][23];
	}
	for(int i=0;i<is2_n;++i)
		xq[is2[i]]=xd1[is2[i]];
	//for(int i=0;i<n;++i)
	//	printf("%d\n",omega[i]);
	for (int i=0;i<n;++i){
		dae.x[omega_idx[i]]=1;
		Pg[i]=-_bus.Pg[bus[i]]*con[i][21];
		Qg[i]=-_bus.Qg[bus[i]]*con[i][22];
		Vg[i]=dae.V[bus[i]];
		ag[i]=dae.a[bus[i]];
		V[i]=Vg[i]*exp(jay*ag[i]);
		S[i]=-Pg[i]+jay*Qg[i];
		I[i]=S[i]/conj(V[i]);
		_delta[i]=std::arg(V[i]+(ra[i]+jay*(xq[i]+xl[i])*I[i]));
		dae.x[delta_idx[i]]=_delta[i];
		Vdq[i]=V[i]*exp(-jay*(_delta[i]-pi/2));
		Idq[i]=I[i]*exp(-jay*(_delta[i]-pi/2));
		vd[i]=Vdq[i].real();
		vq[i]=Vdq[i].imag();
		Id[i]=Idq[i].real();
		Iq[i]=Idq[i].imag();
		pm[i]=(vq[i]+ra[i]*Iq[i])*Iq[i]+(vd[i]+ra[i]*Id[i])*Id[i];
		//printf("%lf\n",pm[i]);
		//printf("%lf\t%lf\t%lf\t%lf\n",vd[i],vq[i],Id[i],Iq[i]);
	}
	if(is2_n>0)
		fm_synit_mod2(is2_n,is2,vd,vq);
	if(is3_n>0)
		fm_synit_mod3(is3_n,is3,vd,vq,dae);
	if(is4_n>0)
		fm_synit_mod4(is4_n,is4,vd,vq,dae);
	if(is51_n>0)
		fm_synit_mod51(is51_n,is51,vd,vq);
	if(is52_n>0)
		fm_synit_mod52(is52_n,is2,vd,vq);
	if(is53_n>0)
		fm_synit_mod53(is53_n,is53,vd,vq);
	if(is6_n>0)
		fm_synit_mod6(is6_n,is8,vd,vq);
	if(is8_n>0)
		fm_synit_mod8(is8_n,is8,vd,vq);
	for (int i=0;i<n;++i){
		vf[i]+=Kp[i]*(pm[i]+Pg[i]);
		//printf("%lf\n",vf[i]);
	}

	//printf("%d\t%d\n",pv.n-pv_n,sw.n-sw_n);
	//for (int i=0;i<dae.n;++i)
	//	printf("%lf\n",dae.x[i]);
	//for(int i=0;i<n;++i)
	//	printf("%d\n",synbus[i]);
	delete []vd;
	delete []vq;
	delete []Vdq;
	delete []Idq;
	delete []ag;
	delete []Vg;
	delete []V;
	delete []S;
	delete []I;
	delete []_delta;
	delete []idx;
	delete []synbus;
	delete []ra   ;
	delete []xl   ;
	delete []xd   ;
	delete []xq   ;
	delete []xd1  ;
	delete []xq1  ;
	delete []xd2  ;
	delete []xq2  ;
	delete []Td10 ;
	delete []Tq10 ;
	delete []Td20 ;
	delete []Tq20 ;
	delete []K0   ;
	delete []Kp   ;
	delete []Taa  ;
}
void Syn::fm_synit_mod2(int is2_n,int *is2,double *vd,double *vq){
	for(int i=0;i<is2_n;++i){
		double ra=con[is2[i]][6];
		double xl=con[is2[i]][5];
//		double xd=con[is2[i]][7];
//		double xq=con[is2[i]][12];
		double xd1=con[is2[i]][8];
		double K;
		K=1/(ra*ra+xd1*xd1+xd1*xl+xl*xd1+xl*xl);
		c1[is2[i]]=ra*K;
		c2[is2[i]]=(xd1+xl)*K;
		c3[is2[i]]=(xd1+xl)*K;
		vf[is2[i]]=vq[is2[i]]+ra*Iq[is2[i]]+(xd1+xl)*Id[is2[i]];
	}
}
void Syn::fm_synit_mod3(int is3_n,int *is3,double *vd,double *vq,DAE dae){
	for(int i=0;i<is3_n;++i){
		int k=is3[i];
		double ra=con[k][6];
		double xl=con[k][5];
		double xd=con[k][7];
		double xq=con[k][12];
		double xd1=con[k][8];
		double K;
		K=1/(ra*ra+xq*xd1+xq*xl+xl*xd1+xl*xl);
		c1[k]=ra*K;
		c2[k]=(xd1+xl)*K;
		c3[k]=(xq+xl)*K;
		dae.x[e1q_idx[k]]=vq[k]+ra*Iq[k]+(xd1+xl)*Id[k];
		vf[k]=dae.x[e1q_idx[k]]+(xd-xd1)*Id[k];
	}
}
void Syn::fm_synit_mod4(int is4_n,int *is4,double *vd,double *vq,DAE dae){
	for(int i=0;i<is4_n;++i){
		int k=is4[i];
		double ra=con[k][6];
		double xl=con[k][5];
		double xd=con[k][7];
//		double xq=con[k][12];
		double xd1=con[k][8];
		double xq1=con[k][13];
		double K;
		K=1/(ra*ra+xq1*xd1+xq1*xl+xl*xd1+xl*xl);
		c1[k]=ra*K;
		c2[k]=(xd1+xl)*K;
		c3[k]=(xq1+xl)*K;
		dae.x[e1q_idx[k]]=vq[k]+ra*Iq[k]+(xd1+xl)*Id[k];
		dae.x[e1d_idx[k]]=vd[k]+ra*Id[k]-(xq1+xl)*Iq[k];
		vf[k]=dae.x[e1q_idx[k]]+(xd-xd1)*Id[k];
		//printf("%lf\n",dae.x[e1d[k]]);
	}
}
void Syn::fm_synit_mod51(int is51_n,int *is51,double *vd,double *vq){
	getchar();
}
void Syn::fm_synit_mod52(int is52_n,int *is52,double *vd,double *vq){
	getchar();
}
void Syn::fm_synit_mod53(int is53_n,int *is53,double *vd,double *vq){
	getchar();
}
void Syn::fm_synit_mod6(int is6_n,int *is6,double *vd,double *vq){
	getchar();
}
void Syn::fm_synit_mod8(int is8_n,int *is8,double *vd,double *vq){
	getchar();
}
