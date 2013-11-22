#include "../inc/DAE.h"
DAE::DAE()
{
	kg=0;
	n=1;
	npf=0;
}
DAE::~DAE()
{
	
	
}
void DAE::daeDelete()
{
	if (n==0)
		return ;
	delete []a;
	delete []V;
	delete []J11;
	delete []J21;
	delete []J12;
	delete []J22;
	delete []Jlf;
	delete []Jlfv;
	delete []Jlfd;
	delete []gp;
	delete []glfp;
	delete []gq;
	delete []glfq;
	delete []tn;
}
int DAE::check(Bus bus)
{
	if(bus.n==0)
		return 0;
	Fx=new double [n*n];
	f=new double [n];
	x=new double [n];
	Fy=new double [n*2*bus.n];
	Gx=new double [2*bus.n*n];
	a=new double [bus.n];
	V=new double [bus.n];
	J11=new double [bus.n*bus.n];
	J21=new double [bus.n*bus.n];
	J12=new double [bus.n*bus.n];
	J22=new double [bus.n*bus.n];
	Jlf=new double [2*bus.n*2*bus.n];
	Jlfv=new double [2*bus.n*2*bus.n];
	Jlfd=new double [2*bus.n*2*bus.n];
	gp=new double [bus.n];
	glfp=new double [bus.n];
	gq=new double [bus.n];
	glfq=new double [bus.n];
	g=new double [2*bus.n];
	for (int i=0;i<n*n;++i)
	{
		Fx[i]=0;
	}
	for (int i=0;i<n*2*bus.n;++i)
	{
		Gx[i]=0;
		Fy[i]=0;
	}
	for (int i=0;i<bus.n*bus.n;++i)
	{
		J11[i]=0;
		J21[i]=0;
		J12[i]=0;
		J22[i]=0;
	}
	for (int i=0;i<2*bus.n*2*bus.n;++i)
	{
		Jlf[i]=0;
		Jlfv[i]=0;
		Jlfd[i]=0;
	}
	for (int i=0;i<bus.n;++i)
	{
		a[i]=0;
		V[i]=0;
		gp[i]=0;
		glfp[i]=0;
		gq[i]=0;
		glfq[i]=0;
	}
	for (int i=0;i<bus.n;++i)
	{
		if(bus.con[i][2]<0.5)
			printf("Warning: some initial guess voltage amplitudes are too low.\n");
		if(bus.con[i][2]>1.5)
			printf("Warning: some initial guess voltage amplitudes are too high.\n");
		V[i]=bus.con[i][2];
		a[i]=bus.con[i][3];
		//printf("%lf\n",V[i]);
	}
	return 1;
}
int DAE::check_2(SW sw)
{
	for (int i=0;i<sw.n;++i)
	{
		V[sw.bus[i]]=sw.con[i][3];
		a[sw.bus[i]]=sw.con[i][4];
	}
	return 1;
	//for(int i=0;i<n;++i)
	//{
	//	printf("%lf\t%lf\n",V[i],a[i]);
	//}
}
