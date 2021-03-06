#include "../inc/Line.h"
Line::Line()
{
	n=0;
}
Line::~Line()
{
	
}
void Line::lineDelete()
{
	if(n==0)
		return;
	for (int i=0;i<n;++i)
		delete []con[i];
	delete []con;
	delete []to;
	delete []from;
	delete []Y;
}
void Line::init(int nodes)
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
	Y=new complex<double>[nodes*nodes];
	for (int i=0;i<nodes *nodes;++i)
		Y[i]=complex<double>(0,0);
}
int Line::check(Settings settings){
	if(n>0)
	{
		from=new int [n];
		to=new int [n];
		for (int i=0;i<n;++i)
		{
			from[i]=(int)con[i][0]-1;
			to[i]=(int)con[i][1]-1;
//			printf("%d\n",from[i]);
		}
		for (int i=0;i<n;++i)
		{
			if(con[i][10]<0.0000005&&con[i][10]>-0.00000005)
				con[i][10]=1;
		}
		int *idx=new int [n];
		int num=0;
		double *temp=new double[n];
		for (int i=0;i<n;++i){
			temp[i]=con[i][5];
			//printf("%d\t",temp[i]);
		}

		num=find(temp,n,0.0,idx);

		if(num!=0)
		{
			for (int i=0;i<num;++i)
			{
				double XB=0;
				XB=con[idx[i]][3]*con[idx[i]][3]/con[idx[i]][2];
				con[idx[i]][7]=con[idx[i]][5]*con[idx[i]][7]/XB;
				con[idx[i]][8]=con[idx[i]][5]*con[idx[i]][8]/XB*settings.rad;
				con[idx[i]][9]=con[idx[i]][5]*con[idx[i]][9]*XB/settings.rad;
			}
			for (int i=0;i<num;++i)
			{
				if(con[idx[i]][9]>10){
					printf("Warning: Some line susceptances are too high...\n");
					con[idx[i]][9]=con[idx[i]][9]/1e6;
				}
			}
		}
		delete []temp;
		

	}

	return 1;
}
void Line::fm_y_line(int nodes,Shunt shunt)
{

	double pi=3.1416;
	complex<double> jay=complex<double>(0,1);//jay=i
	
	double *chrg=new double [n];
	complex<double> *y=new complex<double>[n]; 
	complex<double> *ts=new complex<double>[n]; 
	double *ts2=new double[n]; 

	for ( int i = 0; i < nodes*nodes; i += 1 ) {
	  Y[i]=complex<double>(0,0);
	}
	for (int i=0;i<n;++i){
		chrg[i]=con[i][9];
		//printf("%lf\n",chrg[i]);
		y[i]=1.0/(con[i][7]+jay*con[i][8]);
		ts[i]=con[i][10]*exp(jay*con[i][11]*pi/180.0);
		ts2[i]=(ts[i]*conj(ts[i])).real();
	//	cout<<y[i]<<endl;
//	        cout<<y[i]+jay*chrg[i]/2.0<<endl;
		Y[to[i]+from[i]*nodes]+=-y[i]/std::conj(ts[i]);
		Y[from[i]+to[i]*nodes]+=-y[i]/ts[i];
		Y[from[i]+from[i]*nodes]+=(y[i]+jay*chrg[i]/2.0)/ts2[i];
		Y[to[i]+to[i]*nodes]+=y[i]+jay*chrg[i]/2.0;
		
	}
	for(int i=0;i<nodes;++i){
		Y[i+i*nodes]+=shunt.g[i]+jay*shunt.b[i];
		//printf("%lf\t%lf\t",shunt.g[i],shunt.b[i]);
//		cout<<Y[i+i*nodes]<<endl;
	}
	delete []chrg;
	delete []y;
	delete []ts;
	delete []ts2;

      
	//getchar();
//FILE *fp=fopen("Y.txt","w");
//for (int i=0;i<nodes;++i){
//	for (int j=0;j<nodes;++j)
//	{
//		fprintf(fp,"%lf+%lfi\t",Y[j+i*nodes].real(),Y[j+i*nodes].imag());
//	}
//	fprintf(fp,"\n");
//}
//fclose(fp);
}
