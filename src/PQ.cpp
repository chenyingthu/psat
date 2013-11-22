#include "../inc/PQ.h"
PQ::PQ()
{
	n=0;
}
PQ::~PQ()
{

}
void PQ::pqDelete()
{
	if(n==0)
		return;
	for (int i=0;i<n;++i){
		delete []con[i];
		con[i]=NULL;
	}
	delete []con;
	con=NULL;
	//for (int i=0;i<n;++i){
	//	delete [24]store[i];
	//	store[i]=NULL;
	//}
	//delete [n]store;
	//store=NULL;
	delete []bus;
	delete []P0;
	delete []Q0;
}
void PQ::init()
{
	if(n==0)
		return ;
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
	P0=new double [n];
	Q0=new double [n];
}
int PQ::check(Bus _bus ,Clpsat clpsat)
{
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
				printf("Eror:More than One PQ conentted one bus\n");
				return 0;
			}
		}
	}
	for (int i=0;i<n;++i)
	{
		P0[i]=con[i][3];
		Q0[i]=con[i][4];
		//printf("%lf\t%lf\n",P0[i],Q0[i]);
	}

	return 1;
}
