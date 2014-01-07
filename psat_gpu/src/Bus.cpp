#include "../inc/Bus.h"
#include <stdio.h>
Bus::Bus()
{
	n=0;
	island_n=0;
}
Bus::~Bus()
{
	//for (int i=0;i<n;++i)
	//	delete []con[i];
	//delete []con;
	//delete []internel;
	//delete []Pl;
	//delete []Ql;
	//delete []Pg;
	//delete []Qg;
}
void Bus::init()
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
int Bus::check()
{
	Pl=new double [n];
	Ql=new double [n];
	Pg=new double [n];
	Qg=new double [n];
	for (int i=0;i<n;++i)
	{
		Pl[i]=0;
		Ql[i]=0;
		Pg[i]=0;
		Qg[i]=0;
	}
	if(con==NULL)
	{
		printf("The data file does not seem to be in a valid format: no bus found.\n");
		return 0;
	}
	internel=new int [n];
	for (int i=0;i<n;++i)
		internel[i]=i;
	return 1;
}
void Bus::busDelete(){
	if(n==0)
		return ;
	for (int i=0;i<n;++i)
		delete []con[i];
	delete []con;
	delete []internel;
	delete []Pl;
	delete []Ql;
	delete []Pg;
	delete []Qg;
}
