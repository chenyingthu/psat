#include "../inc/Fault.h"
Fault::Fault()
{
	n=0;
	delta=0;
}
Fault::~Fault(){
	
}
void Fault::faultDelete()
{
	if(n==0)
		return;
	for (int i=0;i<n;++i)
		delete []con[i];
	delete []con;
	delete []bus;
}
void Fault::init()
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
int Fault::check(Bus _bus)
{
	if(n==0)
		return 1;
	bus=new int [n];
	for (int i=0;i<n;++i){
		bus[i]=_bus.internel[(int)con[i][0]-1];
		//printf("%d\n",bus[i]);
	}
}