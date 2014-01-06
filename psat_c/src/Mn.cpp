#include "../inc/Mn.h"
Mn::Mn(){
	n=0;
	init_n=0;
}
void Mn::Mndelete(int _n){
	for (int i=0;i<(n+_n);++i)
		delete []con[i];
	delete []con;
	for (int i=0;i<(n+_n);++i)
		delete []store[i];
	delete []store;

	if(init_n!=0){
		delete []init;
	}
		delete []bus;
    
}
void Mn::Mninit(int _n){
	if((n+_n)==0)
		return;
	con=new double*[n+_n];
	bus=new int [n+_n];
	for (int i=0;i<(n+_n);i++)
	{
		con[i]=new double[24];
	}
	for (int i=0;i<(n+_n);++i)
		for (int j=0;j<24;++j)
			con[i][j]=0;
	store=new double*[n+_n];
	for (int i=0;i<(n+_n);i++)
	{
		store[i]=new double[24];
	}
	for (int i=0;i<(n+_n);++i)
		for (int j=0;j<24;++j)
			store[i][j]=0;
}
void Mn::check(Bus _bus){
    if(n==0)
        return;
    for (int i=0;i<n;++i){
        bus[i]=_bus.internel[(int)con[i][0]];
        if((int)con[i][7]==0)
            init[init_n++]=i;
    }

}
