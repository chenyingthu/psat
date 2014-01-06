#include "../inc/SW.h"
SW::SW()
{
	n=0;
}
SW::~SW()
{
	
}
void SW::swDelete()
{
	if(n==0)
		return ;
	for (int i=0;i<n;++i)
		delete []con[i];
	delete []con;

	delete []bus;
}
void SW::init()
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
}
int SW::check(Bus _bus){
	if(n==0){
		printf("Error: No slack bus found.\n");
		return 0;
	}
	bus=new int [n];
	for (int i=0;i<n;++i){
		bus[i]=_bus.internel[(int)con[i][0]-1];
		//printf("%d\n",bus[i]);
	}
	//_bus.internel[0]=10;
	for (int i=0;i<n;++i)
	{
		int state=0;
		for (int j=0;j<n;++j)
		{
			if (bus[i]==bus[j])
				state++;
			if(state>1)
			{
				printf("Eror:More than One slack conentted one bus\n");
				return 0;
			}
		}
	}
	return 1;
	
}
void SW::fm_synit(int _n,int *_bus){
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
