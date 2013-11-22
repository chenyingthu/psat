#include "../inc/PV.h"
PV::PV()
{
	n=0;
}
PV::~PV()
{

}
void PV::pvDelete()
{
	if(n==0)
		return ;
	delete []bus;
	for (int i=0;i<n;++i){
		delete []con[i];
		con[i]=NULL;
	}
	delete []con;
	con=NULL;
	for (int i=0;i<n;++i){
		delete []store[i];
		store[i]=NULL;
	}
	delete []store;
	store=NULL;
}
void PV::init()
{
	con=new double*[n];
	for (int i=0;i<n;i++)
	{
		con[i]=new double[24];
	}
	for (int i=0;i<n;++i)
		for (int j=0;j<24;++j)
			con[i][j]=0;
	store=new double*[n];
	for (int i=0;i<n;i++)
	{
		store[i]=new double[24];
	}
	for (int i=0;i<n;++i)
		for (int j=0;j<24;++j)
			store[i][j]=0;
}
int PV::check(Bus _bus,Clpsat clpsat,DAE dae){
	if (n==0)
		return 1;
	if(clpsat.init){
		for (int i=0;i<n;++i)
			for (int j=0;j<24;++j)
				store[i][j]=con[i][j];
	}
	bus = new int [n];
	for (int i=0;i<n;++i){
		bus[i]=_bus.internel[(int)con[i][0]-1];
		//printf("%d\n",bus[i]);
	}
	for (int i=0;i<n;++i)
	{
		int state=0;
		for (int j=0;j<n;++j)
		{
			if (bus[i]==bus[j])
				state++;
			if(state>1)
			{
				printf("Eror:More than One PV conentted one bus\n");
				return 0;
			}
		}
	}
	for (int i=0;i<n;++i){
		dae.V[bus[i]]=con[i][4];
		//printf("%lf\n",dae.V[bus[i]]);
	}
	return 1;
}
void PV::fm_synit(int _n,int *_bus){
	int temp=0;
	for(int i=0;i<_n;++i){
		for(int j=0;j<n;++j){
			if(_bus[i]==bus[j]){
				bus[j]=0;
				temp++;
			}
		}
	}
	n=n-temp;
}
