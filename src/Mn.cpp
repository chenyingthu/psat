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
	delete []bus;
	if(init_n!=0)
		delete []init;
}
void Mn::Mninit(int _n){
	if((n+_n)==0)
		return;
	con=new double*[n+_n];
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