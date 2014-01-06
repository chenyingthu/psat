#include "../inc/Shunt.h"
#include<stdlib.h>
#include<stdio.h>
#include <string.h>
using namespace std;
Shunt::Shunt()
{
	n=0;
}
void Shunt::init(int _n){
	con=new double*[_n];
	for (int i=0;i<_n;i++)
	{
		con[i]=new double[24];
	}
	bus=new int [_n];
	g=new double [_n];
	b=new double [_n];
}
void Shunt::Shuntdelete(int _n)
{
	for (int i=0;i<_n;++i)
		delete []con[i];
	delete []con;
	delete []bus;
	delete []g;
	delete []b;
}
void Shunt::check(Bus _bus,Settings settings)
{
	memset(g,0,_bus.n*sizeof(double));
	memset(b,0,_bus.n*sizeof(double));
	memset(bus,-1,_bus.n*sizeof(int));
	if(n>0){
		for (int i=0;i<n;++i){
			bus[i]=_bus.internel[(int)con[i][0]-1];
			con[i][4]=settings.mva*con[i][4]/con[i][1];
			con[i][5]=settings.mva*con[i][5]/con[i][1];
		}
		for (int i=0;i<n;++i){
			g[bus[i]]=con[i][4];
			b[bus[i]]=con[i][5];
		}
	}
}
