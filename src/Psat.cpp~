#include "../inc/Psat.h"
#define pi 3.141592653589793
Psat::Psat(){
    jay=Complex(0,1);
    nBus=0;
    nFault=0;
    nLine=0;
    nPQ=0;
    nPV=0;
    nSW=0;
    nSyn=0;
    nMn=0;
    nShunt=0;
}
Psat::~Psat(){
    bus.busDelete();
    line.lineDelete();
    sw.swDelete();
    pv.pvDelete();
    mn.Mndelete(pq.n);
    pq.pqDelete();
    fault.faultDelete();
    syn.synDelete();
    shunt.Shuntdelete(bus.n);
    varout.deleteVarout();
    delete []fault.dat;
    delete []fault.V;
    delete []fault.ang;
    delete []settings.tempi;
}
void Psat::specifySystem(){
	double type=0;
	int i=0;
	FILE *fp=fopen("/home/yuzhitong/workspace/bishe/psat/caseDfn/d_039ieee.dat","r");
	while(fscanf(fp,"%lf",&type)!=EOF){
        if(i%25==0){
            switch((int)type){
                case 1:bus.n++;break;
                case 2:line.n++;break;
                case 3:sw.n++;break;
                case 4:pv.n++;break;
                case 5:pq.n++;break;
                case 6:fault.n++;break;
                case 7:syn.n++;break;
                case 8:mn.n++;break;
                case 9:shunt.n++;break;
                default:break;
            }
		}
		i++;
	}
	fclose(fp);
}
void Psat::init(){
    bus.init();
    line.init(bus.n);
    sw.init();
    pv.init();
    pq.init();
    fault.init();
    syn.init();
    mn.Mninit(pq.n);
    shunt.init(bus.n);
}
void Psat::formConMatrix(){
    double type=0;
    int i=0;
    FILE *fp=fopen("/home/yuzhitong/workspace/bishe/psat/caseDfn/d_039ieee.dat","r");
	while(fscanf(fp,"%lf",&type)!=EOF){
        if(i%25==0){
            switch((int)type){
                case 1:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&bus.con[nBus][j]);
                        i++;
                    }
                    nBus++;
                    break;
                case 2:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&line.con[nLine][j]);
                        i++;
                    }
                    nLine++;
                    break;
                case 3:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&sw.con[nSW][j]);
                        i++;
                    }
                    nSW++;
                    break;
                case 4:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&pv.con[nPV][j]);
                        i++;
                    }
                    nPV++;
                    break;
                case 5:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&pq.con[nPQ][j]);
                        i++;
                    }
                    nPQ++;
                    break;
                case 6:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&fault.con[nFault][j]);
                        i++;
                    }
                    nFault++;
                    break;
                case 7:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&syn.con[nSyn][j]);
                        i++;
                    }
                    nSyn++;
                    break;
                case 8:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&mn.con[nMn][j]);
                        i++;
                    }
                    nMn++;
                    break;
                case 9:
                    for(int j=0;j<24;++j){
                        fscanf(fp,"%lf",&shunt.con[nShunt][j]);
                        i++;
                    }
                    nShunt++;
                    break;
                default:break;
            }
		}
		i++;
	}
	if((bus.n==nBus)&&(line.n==nLine)&&(sw.n==nSW)&&(pv.n==nPV)&&(pq.n==nPQ)&&(fault.n==nFault)&&(syn.n==nSyn)&&(mn.n==nMn)&&(shunt.n==nShunt))
        printf("read successfully!\n");
	fclose(fp);
}
void Psat::initialLF(){
    settings.show=0;
    settings.init=0;
    if(!settings.init){
        int check=fm_ncomp();
        if(!check){
            //printf("i am here\n");
            return;
        }
        if(settings.conv){
            fm_base_line();
            fm_base_pq();
            fm_base_pv();
            fm_base_sw();
	    fm_base_mn();
            fm_base_syn();
        }
        fm_y();
        settings.refbus=sw.bus[0];
        fm_dynlf();
    }
}
int Psat::fm_ncomp(){
    if(!bus.check())
        return 0;
    if(!line.check(settings))
        return 0;
    shunt.check(bus,settings);
    if(!dae.check(bus))
        return 0;
    if(!sw.check(bus))
        return 0;
    if(!dae.check_2(sw))
        return 0;
    if(!pv.check(bus,clpsat,dae))
        return 0;
    if(!pq.check(bus,clpsat))
        return 0;
    mn.check(bus);
    if(!syn.check(bus))
        return 0;
    if(!fault.check(bus))
        return 0;
    return 1;
}
void Psat::fm_base_line(){
    int idx_n=0;
    int i=0;
    double *VL1=new double [line.n];
    double *V1=new double [line.n];
    double *V2=new double [line.n];
    int *idx=new int [line.n];
    for ( i = 0; i < line.n; i += 1 ) {
      VL1[i]=line.con[i][3];
      V1[i]=bus.con[line.from[i]][1];
      V2[i]=bus.con[line.to[i]][1];
      if((int)line.con[i][6]!=0)
      {
        idx[idx_n++]=i;
      }
    }

    for ( i = 0; i < idx_n; i += 1 ) {
      double kt=line.con[idx[i]][6];
      double KT=V1[idx[i]]/V2[idx[i]];
      double corr=abs(kt-KT)/KT;
      if(corr>0.1)
	printf("Tap ratio of transformer #%d differs more than 10%% from the bases difined\n",idx[i]);
    }

    for ( i = 0; i < line.n; i += 1 ) {
      if(abs(VL1[i]-V1[i])>0.1)
	printf("Voltage of Line # %d differs more than 10%% froem the defined",i);
      line.con[i][7]=line.con[i][7]/line.con[i][2]*settings.mva;
      line.con[i][8]=line.con[i][8]/line.con[i][2]*settings.mva;
      line.con[i][9]=line.con[i][9]*line.con[i][2]/settings.mva;
    }
    delete []idx;
    delete []VL1;
    delete []V1;
    delete []V2;
}
void Psat::fm_base_sw(){


  for (int i = 0; i < sw.n; i += 1 ) {
    sw.con[i][5]=sw.con[i][5]*sw.con[i][1]/settings.mva;
    sw.con[i][6]=sw.con[i][6]*sw.con[i][1]/settings.mva;
    sw.con[i][9]=sw.con[i][9]*sw.con[i][1]/settings.mva;
  }
}
void Psat::fm_base_pq(){


  for (int i = 0; i < pq.n; i += 1 ) {
    pq.con[i][3]=pq.con[i][3]*pq.con[i][1]/settings.mva;
    pq.con[i][4]=pq.con[i][4]*pq.con[i][1]/settings.mva;
    pq.P0[i]=pq.con[i][3];
    pq.Q0[i]=pq.con[i][4];
    //printf("%lf\t %lf\n",pq.P0[i],pq.Q0[i]);
  }
}
void Psat::fm_base_pv(){
  for (int i = 0; i < pv.n; i += 1 ) {
    pv.con[i][3]=pv.con[i][3]*pv.con[i][1]/settings.mva;
    pv.con[i][5]=pv.con[i][5]*pv.con[i][1]/settings.mva;
    pv.con[i][6]=pv.con[i][6]*pv.con[i][1]/settings.mva;
  }
}
void Psat::fm_base_mn(){
  for (int i = 0; i < mn.n; i += 1 ) {
    mn.con[i][3]=mn.con[i][3]*mn.con[i][1]/settings.mva;
    mn.con[i][4]=mn.con[i][4]*mn.con[i][1]/settings.mva;
  }
}
void Psat::fm_base_syn(){
  int j[8]={5,6,7,8,9,12,13,14};
  for (int k = 0; k < syn.n; k += 1 ) {
    for (int i = 0; i < 8; i += 1 ) {
       syn.con[k][j[i]]=(syn.con[k][j[i]]/syn.con[i][1])*settings.mva;
//       printf("%lf\t",syn.con[k][j[i]]);
    }
//    printf("\n");
    syn.con[k][17]=syn.con[k][17]*syn.con[k][1]/settings.mva;
  }

}
void Psat::fm_y(){
  line.fm_y_line(bus.n,shunt);
}
void Psat::fm_dynlf(){}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_spf
 *  Description:  solve power flow
 * =====================================================================================
 */
void Psat::fm_spf(){
  printf("Newton-Raphson Method for Power Flow Computation\n");
  int nodyn=1;
  dae.f[0]=0;
  dae.x[0]=0;
  dae.Fx[0]=1;
  resetBoundNode();
  int FDPF=settings.pfsolver;
  switch(FDPF){
    case 1:
      printf("PF solver Newton-Raphson method\n");
      break;
    case 2:
      printf("PF solver :XB fast decoupled power flow\n");
      break;
    case 3:
      printf("PF solver:BX fast decoupled power flow\n");
      break;
    default:
      break;
  }
  switch(settings.distrsw){
    case 1:
      printf("Distributed slack bus model\n");
      break;
    case 0:
      printf("Single slack bus model\n");
      break;
    default:
      break;
  }
  int iter_max=settings.lfmit;
  if(iter_max<2)
    iter_max=2;
  int iteration=0;
  double tol=settings.lftol;
  double err_max=tol+1;
   
  double *deltf_nosw=new double [dae.n+2*bus.n-2*boundarynode.n];
   double *tmp_nosw=new double [(dae.n+2*bus.n-2*boundarynode.n)*(dae.n+2*bus.n-2*boundarynode.n)];

   double *tmp=new double [(dae.n+2*bus.n)*(dae.n+2*bus.n)];
   double *deltf=new double [dae.n+2*bus.n];
  while ( err_max>tol&&iteration<=iter_max ) {
   fm_lf_1();
   fm_pq_1(); 
   fm_pv_1();
   fm_sw_1();
   for ( int i = 0; i < bus.n; i += 1 ) {
     dae.g[i]=dae.gp[i];
     dae.g[i+bus.n]=dae.gq[i];
   }
   fm_lf_2();
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlf);

   fm_pq_2();
   fm_pv_2();
   fm_sw_2();
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlfv);
   if(nodyn==1)
     dae.Fx[0]=1;


   for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
     if(i<dae.n)
       deltf[i]=dae.f[i];
     else 
       deltf[i]=dae.g[i-dae.n];
   }
   cat4matrix(dae.Fx,dae.n,dae.n,dae.Gx,dae.n,2*bus.n,dae.Fy,2*bus.n,dae.n,dae.Jlfv,2*bus.n,2*bus.n,tmp);
   for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
     int k=boundarynode.indexAll[i];
     deltf_nosw[i]=deltf[k];
   }

   for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
     int row=boundarynode.indexAll[i];
     for ( int j = 0; j < dae.n+2*bus.n-2*boundarynode.n; j += 1 ) {
       int col=boundarynode.indexAll[j];
       tmp_nosw[j+i*(dae.n+2*bus.n-2*boundarynode.n)]=tmp[col+row*(dae.n+2*bus.n)];
     }
   }

   ipiv=new int[dae.n+2*bus.n-2*boundarynode.n];
   LAPACKE_dgesv(LAPACK_COL_MAJOR,dae.n+2*bus.n-2*boundarynode.n,1,tmp_nosw,dae.n+2*bus.n-2*boundarynode.n,ipiv,deltf_nosw,dae.n+2*bus.n-2*boundarynode.n);
   
   for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
     deltf_nosw[i]=-deltf_nosw[i];
   }
   
   for ( int i = 0; i < dae.n; i += 1 ) {
     dae.x[i]+=deltf_nosw[i];
   }
   for ( int i = 0; i < bus.n-boundarynode.n; i += 1 ) {
    int h=boundarynode.indexG[i]; 
    dae.a[h]+=deltf_nosw[dae.n+i];
    dae.V[h]+=deltf_nosw[dae.n+i+bus.n-boundarynode.n];
   }
   
   



   iteration++;
   err_max=0;
   for(int i=0;i<dae.n+2*bus.n-2*boundarynode.n;++i)
   {
     if(abs(deltf_nosw[i])>err_max)
       err_max=abs(deltf_nosw[i]);
   }
  }
   memset(dae.gp,0,bus.n*sizeof(double));
   memset(dae.gq,0,bus.n*sizeof(double));
   fm_pq_1();
   memcpy(bus.Pl,dae.gp,bus.n*sizeof(double));
   memcpy(bus.Ql,dae.gq,bus.n*sizeof(double));
   if(line.n>0)
     fm_lf_1();
   memset(dae.gp,0,bus.n*sizeof(double));
   memset(dae.gq,0,bus.n*sizeof(double));
   fm_pq_1();
   for(int i=0;i<bus.n;++i){
     bus.Pg[i]=dae.gp[i]+dae.glfp[i];
     bus.Qg[i]=dae.gq[i]+dae.glfq[i];
//     printf("%lf\t%lf\n",bus.Pg[i],bus.Qg[i]);
   }
   if(nodyn==1)
     dae.n=0;
   int dynordold=dae.n;
   dae.npf=dae.n;
   fm_dynidx();
   if(dae.n!=0){
     dae.f=new double [dae.n-dynordold];
     dae.x=new double [dae.n-dynordold];
     dae.Fx=new double [dae.n*dae.n];
     dae.Fy=new double [dae.n*2*bus.n];
     dae.Gx=new double [2*bus.n*dae.n];
     memset(dae.x,1,(dae.n-dynordold)*sizeof(double));
     memset(dae.Fx,0,dae.n*dae.n*sizeof(double));
     memset(dae.Fy,0,dae.n*2*bus.n*sizeof(double));
     memset(dae.Gx,0,2*bus.n*dae.n*sizeof(double));
   }
   for ( int i = 0; i < dae.n; i += 1 ) {
     dae.f[i]=1;
     dae.x[i]=1;

   }
   settings.init=1;
   fm_synit();
   fm_excin();
   fm_syn(0);

   memset(dae.gp,0,bus.n*sizeof(double));
   memset(dae.gq,0,bus.n*sizeof(double));
   memset(dae.J11,0,bus.n*bus.n*sizeof(double));
   memset(dae.J12,0,bus.n*bus.n*sizeof(double));
   memset(dae.J21,0,bus.n*bus.n*sizeof(double));
   memset(dae.J22,0,bus.n*bus.n*sizeof(double));

   fm_lf_1();
   fm_pq_1();
   fm_syn(1);
   fm_pv_1();
   fm_sw_1();
   for ( int i = 0; i < bus.n; i += 1 ) {
     dae.g[i]=dae.gp[i];
     dae.g[i+bus.n]=dae.gq[i];
   }
   fm_lf_2();
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlf);

   fm_pq_2();

   fm_syn(2);
   fm_pv_2();
   fm_sw_2();

   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlfv);

   fm_syn(3);
   if(dae.n>0){
     dae.Fx=new double [dae.n*dae.n];
     dae.Fy=new double [dae.n*2*bus.n];
     dae.Gx=new double [2*bus.n*dae.n];
     for (int i=0;i<dae.n*dae.n;++i)
       dae.Fx[i]=0;
     for(int i=0;i<dae.n*2*bus.n;++i){
       dae.Fy[i]=0;
       dae.Gx[i]=0;
     }
   }
   fm_syn(4);
//   FILE	*fp;										/* output-file pointer */
//
//   fp	= fopen( "tmp", "w" );
//   if ( fp == NULL ) {
//     exit (EXIT_FAILURE);
//   }
//
//   for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
//     for ( int j = 0; j <dae.n+2*bus.n; j += 1 ) {
//       fprintf(fp,"%lf\t",tmp[j+i*(dae.n+2*bus.n)]);
//     }
//     fprintf(fp,"\n");
//   }
//   if( fclose(fp) == EOF ) {			/* close output file   */
//     exit (EXIT_FAILURE);
//   }
   delete []deltf_nosw;
   delete []tmp_nosw;
   delete []deltf;
   delete []tmp;


}
void Psat::resetBoundNode(){
 int *indexAll=new int [dae.n+2*bus.n]; 
 int *indexG=new int [bus.n];
 
 boundarynode.indexAll=new int [dae.n+2*bus.n-2*boundarynode.n];
 boundarynode.indexG=new int[bus.n-boundarynode.n];
 for (int i = 0; i < dae.n+2*bus.n; i += 1 ) {
   indexAll[i]=i;
   if(i<bus.n)
     indexG[i]=i;
 }
 for (int i = 0; i < boundarynode.n; i += 1 ) {
   int nbus=boundarynode.bnode[i];
   for (int j = 0; j < dae.n+2*bus.n; j += 1 ) {
    if(indexAll[j]==nbus+dae.n) 
      indexAll[j]=-1;
    if(indexAll[j]==nbus+dae.n+bus.n)
      indexAll[j]=-1;
    }
   for (int j = 0; j < bus.n; j += 1 ) {
     if(indexG[j]==nbus)
       indexG[j]=-1;
   }
   dae.V[nbus]=boundarynode.Voltage[i];
   dae.a[nbus]=boundarynode.Angle[i];
 }
 int temp=0;
 boundarynode.indexAll=new int [dae.n+2*bus.n-2*boundarynode.n];
 boundarynode.indexG=new int [bus.n-boundarynode.n];
 for (int i = 0; i < dae.n+2*bus.n; i += 1 ) {
   if(indexAll[i]!=-1)
     boundarynode.indexAll[temp++]=indexAll[i];
 }
// printf("%d\t%d\n",temp,dae.n+2*bus.n-2*boundarynode.n);
 temp=0;
 for (int i = 0; i < bus.n; i += 1 ) {
   if(indexG[i]!=-1)
     boundarynode.indexG[temp++]=i;
 }
//
// for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
//   printf("%d\n",boundarynode.indexAll[i]);
// }
// printf("%d\t%d\n",temp,bus.n-boundarynode.n);
 delete []indexAll;
 delete []indexG;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_lf_1
 *  Description:  solve load flow when flag =1
 * =====================================================================================
 */
void Psat::fm_lf_1(){
  alpha=Complex(1,0);
  beta=Complex(0,0);
  Complex *Vc=new Complex[bus.n];
  Complex *S=new Complex[bus.n];
  for ( int i = 0; i < bus.n; i += 1 ) {
    Vc[i]=dae.V[i]*exp(jay*dae.a[i]);
//    cout<<Vc[i]<<endl;
  }
  cblas_zgemv(CblasColMajor,CblasTrans,bus.n,bus.n,&alpha,line.Y,bus.n,Vc,1,&beta,S,1);
  for ( int i = 0; i < bus.n; i += 1 ) {
    S[i]=Vc[i]*conj(S[i]);
    dae.gp[i]=S[i].real();
    dae.gq[i]=S[i].imag();
    dae.glfp[i]=dae.gp[i];
    dae.glfq[i]=dae.gq[i];
//    printf ( "%lf\t%lf\n",dae.glfp[i],dae.glfq[i] );
  }
//  cout<<"test"<<endl;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_pq_1
 *  Description:  for dae.gp dae.gq with pq
 * =====================================================================================
 */
void Psat::fm_pq_1(){
  
  for ( int i = 0; i < pq.n; i += 1 ) {
    int k=pq.bus[i];
    dae.gp[k]=pq.con[i][3]+dae.gp[k];
    dae.gq[k]=pq.con[i][4]+dae.gq[k];
    if(dae.V[k]<pq.con[i][6]&&(int)pq.con[i][7]==1)
    {
      dae.gp[k]=pq.con[i][3]*dae.V[k]*dae.V[k]/pq.con[i][6]/pq.con[i][6]+dae.gp[k]-pq.con[i][3];
      dae.gq[k]=pq.con[i][4]*dae.V[k]*dae.V[k]/pq.con[i][6]/pq.con[i][6]+dae.gq[k]-pq.con[i][4];
    }
    else if(dae.V[k]>pq.con[i][5]&&(int)pq.con[i][7]==1)
    {
      dae.gp[k]=pq.con[i][3]*dae.V[k]*dae.V[k]/pq.con[i][5]/pq.con[i][5]+dae.gp[k]-pq.con[i][3];
      dae.gq[k]=pq.con[i][4]*dae.V[k]*dae.V[k]/pq.con[i][5]/pq.con[i][5]+dae.gq[k]-pq.con[i][4];
    }
//    printf ( "%lf\t%lf\n",dae.gp[k],dae.gq[k] );
  }
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_pv_1
 *  Description:  solve dae.gp dae.gq with pv
 * =====================================================================================
 */
void Psat::fm_pv_1(){
  for ( int i = 0; i < pv.n; i += 1 ) {
    double K=1+dae.kg*pv.con[i][9];
    int k=pv.bus[i];
    dae.gp[k]=dae.gp[k]-K*pv.con[i][3];
    dae.gq[k]=0;
  }
//  for ( int i = 0; i < bus.n; i += 1 ) {
//    printf ( "%lf\t%lf\n",dae.gp[i],dae.gq[i] );
//  }
  /*-----------------------------------------------------------------------------
   *  todo:pv2pq
   *-----------------------------------------------------------------------------*/
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_sw_1
 *  Description:  solve pf for sw with flag 1
 * =====================================================================================
 */
void Psat::fm_sw_1(){
  
  for ( int i = 0; i < sw.n; i += 1 ) {
    int k=sw.bus[i];
    dae.gp[k]=0;
    dae.gq[k]=0;
  }
}
void Psat::fm_lf_2(){
  Complex *U=new Complex [bus.n];
  Complex *V=new Complex [bus.n];
  Complex *I=new Complex [bus.n];
  Complex *Vc=new Complex [bus.n*bus.n];
  Complex *Vn=new Complex [bus.n*bus.n];
  Complex *Ic=new Complex [bus.n*bus.n];
  Complex *temp=new Complex [bus.n*bus.n];
  Complex *dS=new Complex[bus.n*bus.n];
  for ( int i = 0; i < bus.n; i += 1 ) {
    U[i]=exp(jay*dae.a[i]);
    V[i]=dae.V[i]*U[i];
  }
  cblas_zgemv(CblasColMajor,CblasTrans,bus.n,bus.n,&alpha,line.Y,bus.n,V,1,&beta,I,1);
  for ( int i = 0; i < bus.n; i += 1 ) {
    Vc[i+i*bus.n]=V[i];
    Vn[i+i*bus.n]=U[i];
    Ic[i+i*bus.n]=I[i];
  }
    cblas_zgemm(CblasColMajor,CblasNoTrans, CblasNoTrans, bus.n, bus.n, bus.n, &alpha,  line.Y, bus.n, Vn,bus.n, &beta, temp, bus.n);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,bus.n,bus.n,bus.n,&alpha,Vc,bus.n,conj_(bus.n*bus.n,temp),bus.n,&beta,dS,bus.n);
    cblas_zgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,bus.n,bus.n,bus.n,&alpha,conj_(bus.n*bus.n,Ic),bus.n,Vn,bus.n,&beta,temp,bus.n);
    for ( int i = 0; i < bus.n*bus.n; i += 1 ) {
     dS[i]=dS[i]+temp[i]; 
     dae.J12[i]=dS[i].real();
     dae.J22[i]=dS[i].imag();
    }
    /*-----------------------------------------------------------------------------
     *  with cblasColMajor ,output has been trans,be Carefull!!!
     *-----------------------------------------------------------------------------*/
    cblas_zgemm(CblasColMajor,CblasNoTrans, CblasNoTrans, bus.n, bus.n, bus.n, &alpha,  line.Y, bus.n, Vc,bus.n, &beta, temp, bus.n);

    for ( int i = 0; i < bus.n*bus.n; i += 1 ) {
      temp[i]=conj(Ic[i]-temp[i]);
    }
    cblas_zgemm(CblasColMajor,CblasNoTrans, CblasNoTrans, bus.n, bus.n, bus.n, &jay,  Vc, bus.n, temp,bus.n, &beta, dS, bus.n);
//    for ( int i = 0; i < bus.n; i += 1 ) {
//      for ( int j = 0; j < bus.n; j += 1 ) {
//	cout<<dS[i+j*bus.n];
//	printf("\t");
//      }
//      printf("\n");
//    }
    for ( int i = 0; i < bus.n*bus.n; i += 1 ) {
      dae.J11[i]=dS[i].real();
      dae.J21[i]=dS[i].imag();
    }

    for ( int i = 0; i < bus.n; i += 1 ) {
      dae.J11[i+i*bus.n]+=1e-6;
      dae.J22[i+i*bus.n]+=1e-6;
    }
    delete []dS;
    delete []temp;
    delete []U;
    delete []V;
    delete []I;
    delete []Vc;
    delete []Vn;
    delete []Ic;
}

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_pq_2
 *  Description:  algebraic Jacobian matrices
 * =====================================================================================
 */
void Psat::fm_pq_2(){
  for ( int i = 0; i < pq.n; i += 1 ) {
    int h=pq.bus[i];
    if(dae.V[h]<pq.con[i][6]&&(int)pq.con[i][7]==1){
      dae.J12[h+h*bus.n]+=2*pq.con[i][3]*dae.V[h]/pq.con[i][6];
      dae.J22[h+h*bus.n]+=2*pq.con[i][4]*dae.V[h]/pq.con[i][6];
    }
    else if(dae.V[h]>pq.con[i][5]&&(int)pq.con[i][7]==1){
      dae.J12[h+h*bus.n]+=2*pq.con[i][3]*dae.V[h]/pq.con[i][5];
      dae.J22[h+h*bus.n]+=2*pq.con[i][4]*dae.V[h]/pq.con[i][5];
    }
  }
}
void Psat::fm_pv_2(){
  for ( int i = 0; i < pv.n; i += 1 ) {
    int h=pv.bus[i];
    for ( int j = 0; j < bus.n; j += 1 ) {
      dae.J21[h+j*bus.n]=0; 
      dae.J22[h+j*bus.n]=0;
      dae.J12[j+h*bus.n]=0;
      dae.J22[j+h*bus.n]=0;
    }
  }
  for ( int i = 0; i < pv.n; i += 1 ) {
    int h=pv.bus[i];
    dae.J22[h+h*bus.n]+=1;
  }
}
void Psat::fm_sw_2(){
  for ( int i = 0; i < sw.n; i += 1 ) {
    int h=sw.bus[i];
    for ( int j = 0; j < bus.n; j += 1 ) {
      dae.J11[h+j*bus.n]=0;
      dae.J12[h+j*bus.n]=0;
      dae.J21[h+j*bus.n]=0;
      dae.J22[h+j*bus.n]=0;
      dae.J11[j+h*bus.n]=0;
      dae.J12[j+h*bus.n]=0;
      dae.J21[j+h*bus.n]=0;
      dae.J22[j+h*bus.n]=0;
    }
  }
  for ( int i = 0; i < sw.n; i += 1 ) {
    int h=sw.bus[i];
    dae.J11[h+h*bus.n]+=1;
    dae.J22[h+h*bus.n]+=1;
  }
}
void Psat::fm_dynidx(){
  dae.n=syn.fm_dynidx(dae);
//  printf("%d\n",dae.n);
}
void Psat::fm_synit(){

  syn.fm_synit(settings,dae,bus);
  pv.fm_synit(syn.n,syn.bus);
  sw.fm_synit(syn.n,syn.bus);
}
void Psat::fm_excin(){}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_syn
 *  Description:  flag = 1 algebraic equations
 *                flag = 2 algebraic jacobians
 *                flag = 3 differential equations;
 *                flag = 4 state matrix
 * =====================================================================================
 */
void Psat::fm_syn(int flag){
  switch ( flag ) {
    case 1:
      fm_syn1(1);
      break;
    case 2:      
      fm_syn2(1);
      break;
    case 3:	
      fm_syn3(1);
      break;
    case 4:
      fm_syn4(1);
      break;
    default:	
      break;
  }				/* -----  end switch  ----- */
}
void Psat::fm_syn1(int flag){
  double ag,ss,cc;
  for ( int i = 0; i < syn.n; i += 1 ) {
    int k1=syn.delta_idx[i];
    int k2=syn.omega_idx[i];
    syn.delta[i]=dae.x[k1];
    syn.omega[i]=dae.x[k2];
    if(i<syn.is3_n){
      int k=syn.is3[i];
      syn.e1q[k]=dae.x[syn.e1q_idx[k]];
    }
    if(i<syn.is4_n){
      int k=syn.is4[i];
      syn.e1d[k]=dae.x[syn.e1d_idx[k]];
      syn.e1q[k]=dae.x[syn.e1q_idx[k]];
    }
    if(i<syn.is51_n){
    }
    if(i<syn.is52_n){
    }
    if(i<syn.is53_n){
    }
    if(i<syn.is6_n){
    }
    if(i<syn.is8_n){
    }

  }
  for ( int i = 0; i < syn.n; i += 1 ) {
    int k=syn.bus[i];
    ag=dae.a[k];
    ss=sin(syn.delta[i]-ag);
    cc=cos(syn.delta[i]-ag);
    syn.Id[i]=-syn.c1[i]*dae.V[k]*ss-syn.c3[i]*dae.V[k]*cc;
    syn.Iq[i]=syn.c2[i]*dae.V[k]*ss-syn.c1[i]*dae.V[k]*cc;
//    printf("%lf\t%lf\n",syn.Id[i],syn.Iq[i]);
  }
  for ( int i = 0; i < syn.n; i += 1 ) {
    if(i<syn.is2_n){
      int k=syn.is2[i];
      syn.Id[k]+=syn.c3[k]*syn.vf[k];
      syn.Iq[k]+=syn.c1[k]*syn.vf[k];
    }
    if(i<syn.is3_n){
      int k=syn.is3[i];
      syn.Id[k]+=syn.c3[k]*syn.e1q[k];
      syn.Iq[k]+=syn.c1[k]*syn.e1q[k];
    }
    if(i<syn.is4_n){
      int k=syn.is4[i];
      syn.Id[k]+=syn.c1[k]*syn.e1d[k]+syn.c3[k]*syn.e1q[k];
      syn.Iq[k]+=syn.c2[k]*syn.e1d[k]+syn.c1[k]*syn.e1q[k];
    }
    if(i<syn.is51_n){
    }
    if(i<syn.is52_n){
    }
    if(i<syn.is53_n){
    }
    if(i<syn.is6_n){
    }
    if(i<syn.is8_n){
    }
  }

  for ( int i = 0; i < syn.n; i += 1 ) {
    int k=syn.bus[i];
    ag=dae.a[k];
    ss=sin(syn.delta[i]-ag);
    cc=cos(syn.delta[i]-ag);
    syn.Pg[i]=-dae.V[k]*(syn.Id[i]*ss+syn.Iq[i]*cc);
    syn.Qg[i]=-dae.V[k]*(syn.Id[i]*cc-syn.Iq[i]*ss);
    dae.gp[k]+=syn.Pg[i];
    dae.gq[k]+=syn.Qg[i];
//    printf("%lf\t%lf\n",syn.Pg[i],syn.Qg[i]);
  }
  
//  for ( int i = 0; i < bus.n; i += 1 ) {
//    printf("%lf\t%lf\n",dae.gp[i],dae.gq[i]);
//  }
}
void Psat::fm_syn2(int flag){
  double ag,ss,cc;
  double M1,M2,M3,M4;
  for ( int i = 0; i < syn.n; i += 1 ) {
    int k=syn.bus[i];
    ag=dae.a[k];
    ss=sin(syn.delta[i]-ag);
    cc=cos(syn.delta[i]-ag);

    M1=dae.V[k]*(syn.c1[i]*cc-syn.c3[i]*ss);
    M2=-dae.V[k]*(syn.c2[i]*cc+syn.c1[i]*ss);
    M3=-(syn.c1[i]*ss+syn.c3[i]*cc);
    M4=syn.c2[i]*ss-syn.c1[i]*cc;

    syn.J11[i]=dae.V[k]*(-M1*ss+syn.Id[i]*cc-M2*cc-syn.Iq[i]*ss);
    syn.J12[i]=-syn.Id[i]*ss-syn.Iq[i]*cc-dae.V[k]*(M3*ss+M4*cc);
    syn.J21[i]=dae.V[k]*(-M1*cc-syn.Id[i]*ss+M2*ss-syn.Iq[i]*cc);
    syn.J22[i]=-syn.Id[i]*cc+syn.Iq[i]*ss-dae.V[k]*(M3*cc-M4*ss);

    dae.J11[k+k*bus.n]+=syn.J11[i];
    dae.J12[k+k*bus.n]+=syn.J12[i];
    dae.J21[k+k*bus.n]+=syn.J21[i];
    dae.J22[k+k*bus.n]+=syn.J22[i];

//    printf("%lf\t%lf\t%lf\t%lf\n",syn.J11[i],syn.J12[i],syn.J21[i],syn.J22[i]);

  }
}
void Psat::fm_syn3(int flag){
  double iM,D;
  double ra,xd,xq,xd1,xq1,Td10,Tq10,K0,Kp;
  double a34,a35,a45,a44,b43,b44;
  double *Vf=new double [syn.n];

  for ( int i = 0; i < syn.n; i += 1 ) {
    ra=syn.con[i][6];
    D=syn.con[i][18];
    iM=1/syn.con[i][17];
    K0=syn.con[i][19];
    Kp=syn.con[i][20];
    dae.f[syn.delta_idx[i]]=settings.rad*(syn.omega[i]-1);
    dae.f[syn.omega_idx[i]]=(syn.pm[i]+syn.Pg[i]-ra*(syn.Id[i]*syn.Id[i]+syn.Iq[i]*syn.Iq[i])-D*(syn.omega[i]-1))*iM;
    Vf[i]=syn.vf[i]+K0*(syn.omega[i]-1)-Kp*(syn.pm[i]+syn.Pg[i]);
  }
//  for ( int i = 0; i < dae.n; i += 1 ) {
//    printf("%lf\n",dae.f[i]);
//  }
  for ( int i = 0; i < syn.n; i += 1 ) {
    if(i<syn.is3_n){
      int k=syn.is3[i];
      Td10=syn.con[k][10];
      xd=syn.con[k][7];
      xd1=syn.con[k][8];
      a34=1/Td10;
      a35=a34*(xd-xd1);
      dae.f[syn.e1q_idx[k]]=-a34*syn.e1q[k]-a35*syn.Id[k]+a34*Vf[k];
    }
    if(i<syn.is4_n){
      int k=syn.is4[i];
      Td10=syn.con[k][10];
      xd=syn.con[k][7];
      xd1=syn.con[k][8];
      xq=syn.con[k][12];
      xq1=syn.con[k][13];
      Tq10=syn.con[k][15];
      a44=1/Td10;
      a45=a44*(xd-xd1);
      b43=1/Tq10;
      b44=b43*(xq-xq1);
      dae.f[syn.e1q_idx[k]]=-a44*syn.e1q[k]-a45*syn.Id[k]+a44*Vf[k];
      dae.f[syn.e1d_idx[k]]=-b43*syn.e1d[k]+b44*syn.Iq[k];
    }
    if(i<syn.is51_n){
    }
    if(i<syn.is52_n){
    }
    if(i<syn.is53_n){
    }
    if(i<syn.is6_n){
    }
    if(i<syn.is8_n){
    }
  }
//  
//  for ( int i = 0; i < dae.n; i += 1 ) {
//    printf("%lf\n",dae.f[i]);
//  }
  delete []Vf;
}
void Psat::fm_syn4(int flag){
  double ag,ss,cc,iM,D;
  double *M1=new double [syn.n];
  double *M2=new double [syn.n];
  double *M3=new double [syn.n];
  double *M4=new double [syn.n];

  double *N1=new double [syn.n];
  double *N2=new double [syn.n];
  for ( int i = 0; i < syn.n; i += 1 ) {
    int k=syn.bus[i];
    ag=dae.a[k];
    ss=sin(syn.delta[i]-ag);
    cc=cos(syn.delta[i]-ag);
    iM=1/syn.con[i][17];
    D=syn.con[i][18];
    M1[i]=dae.V[k]*(syn.c1[i]*cc-syn.c3[i]*ss);
    M2[i]=-dae.V[k]*(syn.c2[i]*cc+syn.c1[i]*ss);
    M3[i]=-(syn.c1[i]*ss+syn.c3[i]*cc); 
    M4[i]=syn.c2[i]*ss-syn.c1[i]*cc;
    double Wn=settings.rad;
    double ra=syn.con[i][6];
//    common jacobians

    int delta=syn.delta_idx[i];
    int omega=syn.omega_idx[i];
    dae.Fx[omega+delta*dae.n]+=Wn;
    dae.Fx[omega+omega*dae.n]+=-iM*D;
    dae.Fx[delta+omega*dae.n]+=(-syn.J11[i]+2*ra*(syn.Id[i]*M2[i]+syn.Iq[i]*M1[i]))*iM;

    dae.Fy[k+omega*2*bus.n]+=(syn.J11[i]-2*ra*(syn.Id[i]*M2[i]+syn.Iq[i]*M1[i]))*iM;
    dae.Fy[k+bus.n+omega*2*bus.n]+=(syn.J12[i]-2*ra*(syn.Id[i]*M3[i]+syn.Iq[i]*M4[i]))*iM;

    syn.Gp[0+i*8]=-syn.J11[i];
    syn.Gq[0+i*8]=-syn.J21[i];
    dae.Gx[delta+k*dae.n]+=-syn.J11[i];
    dae.Gx[delta+(k+bus.n)*dae.n]+=-syn.J21[i];

    syn.Gp[2+8*i]=dae.V[k]*(-syn.c3[i]*ss-syn.c1[i]*cc);
    syn.Gp[3+8*i]=dae.V[k]*(-syn.c1[i]*ss+syn.c2[i]*cc);
    syn.Gq[2+8*i]=dae.V[k]*(-syn.c3[i]*cc+syn.c1[i]*ss);
    syn.Gq[3+8*i]=dae.V[k]*(-syn.c1[i]*cc-syn.c2[i]*ss);
    N1[i]=(syn.Gp[2+8*i]-2*ra*(syn.Id[i]*syn.c3[i]+syn.Iq[i]*syn.c1[i]))*iM;
    N2[i]=(syn.Gp[3+8*i]-2*ra*(syn.Id[i]*syn.c1[i]+syn.Iq[i]*syn.c2[i]))*iM;
  }

  for ( int i = 0; i < syn.n; i += 1 ) {
    if(i<syn.is3_n){
      int k=syn.is3[i];
      syn.e1q[k]=dae.x[syn.e1q_idx[k]];
      double Td10=syn.con[k][10];
      double xd=syn.con[k][7];
      double xd1=syn.con[k][8];
      double a34=1/Td10;
      double a35=a34*(xd-xd1);
      double Kp=syn.con[k][20];
      double K0=syn.con[k][19];
      int o3=syn.omega_idx[k];
      int e1q3=syn.e1q_idx[k];
      dae.Fx[e1q3+o3*dae.n]+=N1[k];
      dae.Fx[syn.delta_idx[k]+e1q3*dae.n]+=a35*M1[k]+a34*Kp*syn.J11[k];
      dae.Fx[o3+e1q3*dae.n]+=a34*K0;
      dae.Fx[e1q3+e1q3*dae.n]+=-a34-a35*syn.c3[k]+a34*Kp*syn.Gp[2+k*8];

      dae.Gx[e1q3+syn.bs3[i]*dae.n]+=syn.Gp[2+k*8];
      dae.Gx[e1q3+(syn.bs3[i]+bus.n)*dae.n]+=syn.Gq[2+k*8];

      dae.Fy[syn.bs3[i]+e1q3*2*bus.n]+=-a35*M1[k]-a34*Kp*syn.J11[k];
      dae.Fy[syn.bs3[i]+bus.n+e1q3*2*bus.n]+=-a35*M3[k]-a34*Kp*syn.J12[k];
    }

    if(i<syn.is4_n){
      int k=syn.is4[i];
      syn.e1q[k]=dae.x[syn.e1q_idx[k]];
      syn.e1d[k]=dae.x[syn.e1d_idx[k]];
      double Td10=syn.con[k][10];
      double xd=syn.con[k][7];
      double xd1=syn.con[k][8];
      double Tq10=syn.con[k][15];
      double xq=syn.con[k][12];
      double xq1=syn.con[k][13];
      double a44=1/Td10;
      double a45=a44*(xd-xd1);
      double b43=1/Tq10;
      double b44=b43*(xq-xq1);
      double Kp=syn.con[k][20];
      double K0=syn.con[k][19];
      int o4=syn.omega_idx[k];
      int e1q4=syn.e1q_idx[k];
      int e1d4=syn.e1d_idx[k];

      dae.Fx[e1q4+o4*dae.n]+=N1[k];
      dae.Fx[e1d4+o4*dae.n]+=N2[k];

      dae.Fx[syn.delta_idx[k]+e1q4*dae.n]+=a45*M1[k]+a44*Kp*syn.J11[k];
      dae.Fx[o4+e1q4*dae.n]+=a44*K0;
      dae.Fx[e1q4+e1q4*dae.n]+=-a44-a45*syn.c3[k]+a44*Kp*syn.Gp[2+k*8];

      dae.Fx[e1d4+e1q4*dae.n]+=-a45*syn.c1[k]+a44*Kp*syn.Gp[3+k*8];
      dae.Fx[syn.delta_idx[k]+e1d4*dae.n]+=b44*M2[k];
      dae.Fx[e1q4+e1d4*dae.n]+=b44*syn.c1[k];
      dae.Fx[e1d4+e1d4*dae.n]+=-b43-b44*syn.c2[k];


      dae.Gx[e1q4+syn.bs4[i]*dae.n]+=syn.Gp[2+k*8];
      dae.Gx[e1d4+syn.bs4[i]*dae.n]+=syn.Gp[3+k*8];
      dae.Gx[e1q4+(syn.bs4[i]+bus.n)*dae.n]+=syn.Gq[2+k*8];
      dae.Gx[e1d4+(syn.bs4[i]+bus.n)*dae.n]+=syn.Gq[3+k*8];

      dae.Fy[syn.bs4[i]+e1q4*2*bus.n]+=-a45*M1[k]-a44*Kp*syn.J11[k];
      dae.Fy[syn.bs4[i]+bus.n+e1q4*2*bus.n]+=-a45*M3[k]-a44*Kp*syn.J12[k];
      dae.Fy[syn.bs4[i]+e1d4*2*bus.n]+=b44*M2[k];
      dae.Fy[syn.bs4[i]+bus.n+e1d4*2*bus.n]+=b44*M4[k];
    }
    if(i<syn.is51_n){
    }
    if(i<syn.is52_n){
    }
    if(i<syn.is53_n){
    }
    if(i<syn.is6_n){
    }
    if(i<syn.is8_n){
    }
  }
  delete []N1;
  delete []N2;
  delete []M1;
  delete []M2;
  delete []M3;
  delete []M4;
}
void Psat::fm_int_intial(){

  if(dae.n==0&&clpsat.init==0)
    printf("Time domain simulation aborted.\n");
  switch(settings.method){
    case 1:
      printf("Implicit Euler integration method\n");
      break;
    case 2:
      printf("Trapezoidal integration method\n");
      break;
    default:
      break;
  }
  if(clpsat.pq2z)
    settings.pq2z=1;
  if(settings.pq2z&&settings.init<2&&pq.n>0){
    for ( int i = 0; i < pq.n; i += 1 ) {
      pq.con[i][3]=pq.con[i][3]/(dae.V[bus.internel[(int)pq.con[i][0]-1]]*dae.V[bus.internel[(int)pq.con[i][0]-1]]);
      pq.con[i][4]=pq.con[i][4]/(dae.V[bus.internel[(int)pq.con[i][0]-1]]*dae.V[bus.internel[(int)pq.con[i][0]-1]]);
      for ( int j = 0; j < 8; j += 1 ) {
	if(j<5){
	  mn.con[mn.n+i][j]=pq.con[i][j];
	}
	else if(j<7){
	  mn.con[mn.n+i][j]=2;
	}
	else if(j<8){
	  mn.con[mn.n+i][j]=1;
	}
      }
    }
    mn.n+=pq.n;
    pq.n=0;
    for ( int i = 0; i < mn.n; i += 1 ) {
      mn.bus[i]=bus.internel[(int)mn.con[i][0]-1];
//      printf("%d\n",mn.bus[i]);
    }
  }
  for(int i=0;i<bus.n*bus.n;++i){
    dae.J11[i]=0;
    dae.J12[i]=0;
    dae.J21[i]=0;
    dae.J22[i]=0;
  }
  dae.t=settings.t0;
  fm_lf_1();
  fm_mn_1();
  fm_syn(1);
  for ( int i = 0; i < boundarynode.n; i += 1 ) {
    int k=boundarynode.bnode[i];
    boundarynode.P[i]=dae.gp[k];
    boundarynode.Q[i]=dae.gq[k];
  }
  fm_sw_1();
  for ( int i = 0; i < bus.n; i += 1 ) {
     dae.g[i]=dae.gp[i];
     dae.g[i+bus.n]=dae.gq[i];
//     printf("%lf\t%lf\n",dae.g[i],dae.g[i+bus.n]);
  }
  fm_lf_2();
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlf);

   fm_mn_2();
   fm_syn(2);
   fm_sw_2();

   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlfv);
   fm_syn(3);
   if(dae.n>0){
     dae.Fx=new double [dae.n*dae.n];
     dae.Fy=new double [dae.n*2*bus.n];
     dae.Gx=new double [2*bus.n*dae.n];
     for ( int i = 0; i < dae.n*dae.n; i += 1 ) {
       dae.Fx[i]=0;
     }
     for ( int i = 0; i < dae.n*2*bus.n; i += 1 ) {
       dae.Fy[i]=0;
       dae.Gx[i]=0;
     }
   }
   fm_syn(4);

     
   dae.tn=new double [dae.n];
   for ( int i = 0; i < dae.n; i += 1 ) {
     dae.tn[i]=dae.f[i];
   }
   for ( int i = 0; i < bus.n; i += 1 ) {
     dae.g[i]=dae.gp[i];
     dae.g[i+bus.n]=dae.gq[i];
  }

   fm_out_0(0,0);
   fm_out_2(settings.t0,1);
   if(fault.n>0)
     fm_fault_0(0);
   if(fault.n>0){
     double *tempiguasto=new double [2*fault.n];
     settings.tempi=new double [4*fault.n];
     for ( int i = 0; i < fault.n; i += 1 ) {
       tempiguasto[i]=fault.con[i][4];
       tempiguasto[i+fault.n]=fault.con[i][5];
     }
     tempiguasto=sortD(tempiguasto,2*fault.n);
     for ( int i = 0; i < 2*fault.n; i += 1 ) {
       settings.tempi[2*i]=tempiguasto[i]-1e-6;
       settings.tempi[2*i+1]=tempiguasto[i];
     }
//     for ( int i = 0; i < 4*fault.n; i += 1 ) {
//       printf("%lf\n",settings.tempi[i]);
//     }
     delete []tempiguasto;
   }
//   FILE	*fp;										/* output-file pointer */
//
//   fp	= fopen( "Jlfv", "w" );
//   if ( fp == NULL ) {
//     exit (EXIT_FAILURE);
//   }
//
//   for ( int i = 0; i < 2*bus.n; i += 1 ) {
//     for ( int j = 0; j <2*bus.n; j += 1 ) {
//       fprintf(fp,"%lf\t",dae.Jlfv[j+i*(2*bus.n)]);
//     }
//     fprintf(fp,"\n");
//   }
//   if( fclose(fp) == EOF ) {			/* close output file   */
//     exit (EXIT_FAILURE);
//   }
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_mn
 *  Description:  flag = 0 initialization
 *                flag = 1 algebraic equations
 *                flag =2 algebraic jacobians
 * =====================================================================================
 */
void Psat::fm_mn_1(){
  if(settings.init){
    for ( int i = 0; i < mn.n; i += 1 ) {
       int k=mn.bus[i];
       dae.gp[k]+=mn.con[i][3]*pow(dae.V[k],(int)mn.con[i][5]);
       dae.gq[k]+=mn.con[i][4]*pow(dae.V[k],(int)mn.con[i][6]);
    }
  }
  else if(mn.init_n>0){
    for ( int i = 0; i < mn.init_n; i += 1 ) {
      int k=mn.init[i];
      int h=mn.bus[k];
      dae.gp[h]+=mn.con[k][3]*pow(dae.V[h],(int)mn.con[k][5]);
      dae.gq[h]+=mn.con[k][4]*pow(dae.V[h],(int)mn.con[k][6]);
    }
  }
}
void Psat::fm_mn_2(){
  if(settings.init){
    for ( int i = 0; i < mn.n; i += 1 ) {
      int k=mn.bus[i];
      dae.J12[k+k*bus.n]+=mn.con[i][3]*mn.con[i][5]*pow(dae.V[k],(int)mn.con[i][5]-1);
      dae.J22[k+k*bus.n]+=mn.con[i][4]*mn.con[i][6]*pow(dae.V[k],(int)mn.con[i][6]-1);
    }
  }
  else if(mn.init_n>0){
    for ( int i = 0; i < mn.init_n; i += 1 ) {
      int k=mn.init[i];
      int h=mn.bus[k];
      dae.J12[h+h*bus.n]+=mn.con[k][3]*mn.con[k][5]*pow(dae.V[h],(int)mn.con[i][5]-1);
      dae.J22[h+h*bus.n]+=mn.con[k][4]*mn.con[k][6]*pow(dae.V[k],(int)mn.con[i][6]-1);
    }
  }
}
void Psat::fm_mn_0(){}
void Psat::fm_fault(int flag,double t){}
void Psat::fm_fault_0(double t){
  fault.dat=new double [fault.n*5];
  fault.V=new double [bus.n];
  fault.ang=new double [bus.n];
  for ( int i = 0; i < bus.n; i += 1 ) {
    fault.ang[i]=dae.a[i];
  }
  if(syn.n>0){
    for ( int i = 0; i < syn.n; i += 1 ) {
      int k=syn.delta_idx[i];
      fault.delta+=dae.x[k];
    }
    fault.delta/=syn.n;
  }
  else
    fault.delta=0;
  Complex *x=new Complex [fault.n];
  Complex *y=new Complex [fault.n];

  for ( int i = 0; i < fault.n; i += 1 ) {
    int k=fault.bus[i];
    x[i]=fault.con[i][6]+jay*fault.con[i][7];
    if(abs(x[i])==0)
      x[i]=jay*1e-6;
    y[i]=1.0/x[i];
    fault.dat[0+i*5]=y[i].real();
    fault.dat[1+i*5]=y[i].imag();
    fault.dat[2+i*5]=shunt.g[k];
    fault.dat[3+i*5]=shunt.b[k];
  }
  delete []x;
  delete []y;
}
void Psat::fm_fault_1(double t){
  for ( int i = 0; i < fault.n; i += 1 ) {
    int h=fault.bus[i];
    if(abs(t-fault.con[i][4])<1e-8){
      printf("applying fault at t = %lf s\n",t);
      shunt.g[h]=fault.dat[2+i*5]+fault.dat[0+i*5];
      shunt.b[h]=fault.dat[3+i*5]+fault.dat[1+i*5];
      fm_y();
      int conv=fm_nrlf(40,1e-4);
      if(conv)
	printf("i am here\n");
    }// fault intervention
    else if(abs(t-fault.con[i][5])<1e-8){
      printf("Clearing fault at t = %lf s\n",t);
      shunt.g[h]=fault.dat[2+i*5];
      shunt.b[h]=fault.dat[3+i*5];
      fm_y();
      printf("LY done\n");
      if(syn.n>0){
	double mean_delta=0;
	for ( int i = 0; i < syn.n; i += 1 ) {
	  int k=syn.delta_idx[i];
	  mean_delta+=dae.x[k];
	}
	mean_delta=mean_delta/syn.n;
	for ( int i = 0; i < bus.n-boundarynode.n; i += 1 ) {
	  int k=boundarynode.indexG[i];
	  dae.a[k]=mean_delta-fault.delta+fault.ang[k];
	}
      }
      else{
	for ( int i = 0; i < bus.n; i += 1 ) {
	  dae.a[i]=fault.ang[i];
	}
      }
      int conv=fm_nrlf(40,1e-4);
      if(conv){
	printf("i am here too\n");
      }
    }//end of else if

  }//end of for
}
double Psat:: fm_tstep(int flag,int convergency,int iteration,double t){
  switch(flag){
    case 1:
      fm_tstep_1(convergency,iteration,t);
      break;
    case 2:
      fm_tstep_2(convergency,iteration,t);
      break;
    default:
      break;
  }
  return settings.deltat;
}
void Psat::fm_tstep_1(int convergency,int iteration,double t){
  double *tempJlfv=new double [2*bus.n*2*bus.n];
  double *tempGx=new double [2*bus.n*dae.n];
  double *tempFy=new double [dae.n*2*bus.n];
  double *As=new double [dae.n*dae.n];
  double *wr=new double [dae.n];
  double *wi=new double [dae.n];
  double *vr=new double [dae.n*dae.n];
  double *vl=new double [dae.n*dae.n];
  for ( int i = 0; i < 2*bus.n; i += 1 ) {
    for ( int j = 0; j < 2*bus.n; j += 1 ) {
      tempJlfv[j+i*2*bus.n]=dae.Jlfv[i+j*2*bus.n];
    }
  }
  for ( int i = 0; i < 2*bus.n; i += 1 ) {
    for ( int j = 0; j < dae.n; j += 1 ) {
      tempGx[i+j*2*bus.n]=dae.Gx[j+i*dae.n];
    }
  }


  for ( int i = 0; i < dae.n; i += 1 ) {
    for ( int j = 0; j < 2*bus.n; j += 1 ) {
      tempFy[i+j*dae.n]=dae.Fy[j+i*2*bus.n];
    }
  }
  ipiv=new int[2*bus.n];
  double alpha1=1;
  double beta1=0;
  LAPACKE_dgesv(LAPACK_COL_MAJOR,2*bus.n,dae.n,tempJlfv,2*bus.n,ipiv,tempGx,2*bus.n);
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,dae.n,dae.n,2*bus.n,alpha1,tempFy,dae.n,tempGx,2*bus.n,beta1,As,dae.n);

  for ( int i = 0; i < dae.n; i += 1 ) {
    for ( int j = 0; j < dae.n; j += 1 ) {
      As[j+i*dae.n]=dae.Fx[i+j*dae.n]-As[j+i*dae.n];
    }
  }
  LAPACKE_dgeev(LAPACK_COL_MAJOR,'N','N',dae.n,As,dae.n,wr,wi,vl,dae.n,vr,dae.n);
  double freq=0;
  int freq_max=0;
  for(int i=0;i<dae.n;++i){
    if(freq<abs(wi[i])){
      freq=abs(wi[i]);
      freq_max=i;
    }
  }
  freq=sqrt(wr[freq_max]*wr[freq_max]+wi[freq_max]*wi[freq_max]);
  if(freq>settings.freq)
    freq=settings.freq;
  double deltaT=abs(settings.tf-settings.t0);
  double Tstep=1/freq; 
  if(5*Tstep<deltaT/100)
    settings.deltamax=5*Tstep;
  else
    settings.deltamax=deltaT/100;
  settings.chunk=100;
  settings.deltat=Tstep;
  settings.deltamin=Tstep/64;
  if(settings.fixt==1){
    if(settings.tstep<0){
      printf("Error:fixed time step is negative or zero");
      settings.fixt=0;
    }
    else if(settings.tstep<settings.deltamin)
      printf("Warning:fixed time stemp is less than estimated minimum time step");
    else
      settings.deltat=settings.tstep;
  }
  delete []wr;
  delete []wi;
  delete []vr;
  delete []vl;
  delete []tempJlfv;
  delete []tempGx;
  delete []tempFy;
  delete []As;
}
void Psat::fm_tstep_2(int convergency,int iteration,double t){}
void Psat::fm_out_0(double t,int k){
  settings.chunk=200;
  int chunk=settings.chunk;
  varout.t=new double [chunk];
  varout.numOfStep=0;
  if(dae.n>0){
    varout.x=new double [chunk*dae.n];
    varout.f=new double [chunk*dae.n];

    for ( int i = 0; i < chunk*dae.n; i += 1 ) {
      varout.x[i]=0;
      varout.f[i]=0;
    }
  }
  varout.V=new double [chunk*bus.n];
  varout.ang=new double [chunk*bus.n];
  for ( int i = 0; i < chunk*bus.n; i += 1 ) {
    varout.V[i]=0;
    varout.ang[i]=0;
  }
  if(syn.n>0){
    varout.Pm=new double [chunk*syn.n];
    varout.Vf=new double [chunk*syn.n];
    for ( int i = 0; i < chunk*syn.n; i += 1 ) {
      varout.Pm[i]=0;
      varout.Vf[i]=0;
    }
  }
}
void Psat::fm_out_1(double t,int k){}
void Psat::fm_out_2(double t,int k){
  int kk=k-1;//for array start with zero
  if((int)varout.t[kk]==0){
    varout.numOfStep++;
    if(varout.numOfStep>40)
      varout.numOfStep=1;
  }
  varout.t[kk]=t;
  for ( int i = 0; i < dae.n; i += 1 ) {
    varout.x[i+kk*dae.n]=dae.x[i];
    varout.f[i+kk*dae.n]=dae.f[i];
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    varout.V[i+kk*bus.n]=dae.V[i];
    varout.ang[i+kk*bus.n]=dae.a[i];
    printf("out %lf\t%lf\n",dae.V[i],dae.a[i]);
  }
  for ( int i = 0; i < syn.n; i += 1 ) {
    varout.Pm[i+kk*bus.n]=syn.pm[i];
    varout.Vf[i+kk*bus.n]=syn.vf[i];
  }
}
void Psat::fm_out_3(double t,int k){}
void Psat::fm_int_dyn(double t0,double tf,double h){
  int k=1;
  for ( int i = 0; i < settings.chunk; i += 1 ) {
    if(abs(varout.t[i]-settings.t0)<1e-10){
      k=i+1;
      break;
    }
  }
  double t=t0;
  dae.Ac=new double [(dae.n+2*bus.n)*(dae.n+2*bus.n)];
  dae.tn=new double [dae.n];
  t_end=t;
  nrecord=k;
  nexttstep=h;
  while(t_end<tf&&(t_end+nexttstep)>t){
    fm_int_step(t_end,nexttstep,settings.tempi,nrecord);
    printf("k:\t%d\n",nrecord);
  }
  printf("i am out\n");
  delete []dae.Ac;
  delete []dae.tn;
}
void Psat::fm_int_step(double t,double h,double *tempi,int k){
  double h_old=h;
  if((t+h)>settings.tf)
    h=settings.tf-t;
  double tempo=t+h;
  double tempo_min=tempo;
  for ( int i = 0; i < 4*fault.n; i += 1 ) {
    if(tempi[i]-t>1e-6&&tempo-tempi[i]>1e-6)
      if(tempo_min>tempi[i])
	tempo_min=tempi[i];
  }
  tempo=tempo_min;
  h=tempo-t;
  dae.t=tempo;
  int kk=k-1; 
  double *xa=new double [dae.n];
  double *anga=new double [bus.n];
  double *Va=new double [bus.n];
  double *fn=new double [dae.n];
  for ( int i = 0; i < dae.n; i += 1 ) {
    dae.x[i]=varout.x[i+kk*dae.n];
    dae.f[i]=varout.f[i+kk*dae.n];
    xa[i]=dae.x[i];
    fn[i]=dae.f[i];
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    dae.V[i]=varout.V[i+kk*bus.n];
    dae.a[i]=varout.ang[i+kk*bus.n];
    anga[i]=dae.a[i];
    Va[i]=dae.V[i];
  }
  fault.tFaultStart=fault.con[0][4];
  fault.tFaultEnd=fault.con[0][5];
  if(k<settings.chunk&&varout.t[kk+1]==tempo){
    for ( int i = 0; i < dae.n; i += 1 ) {
      dae.x[i]=varout.x[i+k*dae.n];
    }
    for ( int i = 0; i < bus.n; i += 1 ) {
      dae.V[i]=varout.V[i+k*bus.n];
      dae.a[i]=varout.ang[i+k*bus.n];
    }
  }
  else{
    if(settings.t0>fault.tFaultEnd+3*h){
      for ( int i = 0; i < dae.n; i += 1 ) {
	dae.x[i]=3*varout.x[i+kk*dae.n]-3*varout.x[i+(kk-1)*dae.n]+varout.x[i+(kk-2)*dae.n];
      }
      for ( int i = 0; i < bus.n; i += 1 ) {
	dae.V[i]=3*varout.V[i+kk*bus.n]-3*varout.V[i+(kk-1)*bus.n]+varout.V[i+(kk-2)*bus.n];
	dae.a[i]=3*varout.ang[i+kk*bus.n]-3*varout.ang[i+(kk-1)*bus.n]+varout.ang[i+(kk-2)*bus.n];
      }
    }
  }
  if(abs(tempo-fault.tFaultEnd)<1e-8){
    for ( int i = 0; i < settings.chunk; i += 1 ) {
      if(abs(varout.t[i]-fault.tFaultStart)<1e-10){
	kk=i;
	for ( int i = 0; i < bus.n; i += 1 ) {
	  dae.V[i]=varout.V[i+kk*bus.n];
	  dae.a[i]=varout.ang[i+kk*bus.n];
	}
	break;
      }
    }
  }
  resetBoundNode();

  for ( int i = 0; i < 4*fault.n; i += 1 ) {
    if(abs(tempi[i]-tempo)<1e-5){
      if(fault.n>0){
	printf("fault\t%lf\n",tempo);
	fm_fault_1(tempo);
	break;
      }
    }
  }
  double err_max=1;
  int iterazione=1;
  int AllN=dae.n+2*bus.n-2*boundarynode.n;
  int GN=bus.n-boundarynode.n;
  double *inc=new double [AllN];
  double *tempAc=new double [AllN*AllN];
  double *deltf=new double [dae.n+2*bus.n];
  int *ipiv=new int[AllN];
  while(err_max>settings.dyntol&&iterazione<=settings.dynmit){
    fm_lf_1();
    fm_mn_1();
    fm_syn(1);

    for ( int i = 0; i < boundarynode.n; i += 1 ) {
      int k=boundarynode.bnode[i];
      boundarynode.P[i]=dae.gp[k];
      boundarynode.Q[i]=dae.gq[k];
    }
    fm_sw_1();
    for ( int i = 0; i < bus.n; i += 1 ) {
      dae.g[i]=dae.gp[i];
      dae.g[i+bus.n]=dae.gq[i];
    }
    fm_syn(3);
    if(iterazione%3==1){
      fm_lf_2();
      cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlf);

      fm_mn_2();
      fm_syn(2);
      fm_sw_2();

      cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlfv);
      for ( int i = 0; i < dae.n*dae.n; i += 1 ) {
	dae.Fx[i]=0;
      }
      for ( int i = 0; i < dae.n*2*bus.n; i += 1 ) {
	dae.Fy[i]=0;
	dae.Gx[i]=0;
      }
      fm_syn(4);
      for ( int i = 0; i < sw.n; i += 1 ) {
	int k=sw.bus[i];
	for ( int j = 0; j < dae.n; j += 1 ) {
	  dae.Fy[k+j*2*bus.n]=0;
	  dae.Fy[k+bus.n+j*2*bus.n]=0;
	  dae.Gx[j+k*dae.n]=0;
	  dae.Gx[j+(k+bus.n)*dae.n]=0;
	}
      }
      for ( int i = 0; i >= pv.n; i -= 1 ) {
	int k=pv.bus[i];
	for ( int j = 0; j < dae.n; j += 1 ) {
	  dae.Fy[k+bus.n+j*2*bus.n]=0;
	  dae.Gx[j+(k+bus.n)*dae.n]=0;
	}
      }
      switch(settings.method){
	case 1:
	  cat4matrix_Ac(dae.Fx,dae.n,dae.n,dae.Gx,dae.n,2*bus.n,dae.Fy,2*bus.n,dae.n,dae.Jlfv,2*bus.n,2*bus.n,dae.Ac,h);
	  break;
	case 2:
	  cat4matrix_Ac(dae.Fx,dae.n,dae.n,dae.Gx,dae.n,2*bus.n,dae.Fy,2*bus.n,dae.n,dae.Jlfv,2*bus.n,2*bus.n,dae.Ac,0.5*h);
	  break;
	default:
	  break;
      }
    }//den of if iterazione
    switch(settings.method){
      case 1:
	for ( int i = 0; i < dae.n; i += 1 ) {
	  dae.tn[i]=dae.x[i]-xa[i]-h*dae.f[i];
	}
	break;
      case 2:
	for ( int i = 0; i < dae.n; i += 1 ) {
	  dae.tn[i]=dae.x[i]-xa[i]-h*0.5*dae.f[i];
	}
	break;
      default:
	break;
    }
    for ( int i = 0; i < dae.n+2*bus.n; i += 1 ) {
      if(i<dae.n)
	deltf[i]=dae.tn[i];
      else
	deltf[i]=dae.g[i-dae.n];
    }
    for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
      int k=boundarynode.indexAll[i];
      inc[i]=deltf[k];
    }

    for ( int i = 0; i < dae.n+2*bus.n-2*boundarynode.n; i += 1 ) {
      int row=boundarynode.indexAll[i];
      for ( int j = 0; j < dae.n+2*bus.n-2*boundarynode.n; j += 1 ) {
	int col=boundarynode.indexAll[j];
	tempAc[j+i*(dae.n+2*bus.n-2*boundarynode.n)]=dae.Ac[col+row*(dae.n+2*bus.n)];
      }
    }
   LAPACKE_dgesv(LAPACK_COL_MAJOR,AllN,1,tempAc,AllN,ipiv,inc,AllN);
   err_max=0;
   for ( int i = 0; i <AllN; i += 1 ) {
     inc[i]=-inc[i];
     if(abs(inc[i])>err_max)
       err_max=abs(inc[i]);
     if(i<dae.n){
       dae.x[i]+=inc[i];
     }
     else if(i<dae.n+GN){
       int k=boundarynode.indexG[i-dae.n];
       dae.a[k]+=inc[i];
     }
     else{
       int k=boundarynode.indexG[i-dae.n-GN];
       dae.V[k]+=inc[i];
     }
   }
   iterazione++;
  }//end of while
  if(iterazione>settings.dynmit){
    h=fm_tstep(2,0,iterazione,t);
    for ( int i = 0; i < dae.n; i += 1 ) {
      dae.x[i]=xa[i];
      dae.f[i]=fn[i];
    }
    for ( int i = 0; i < bus.n; i += 1 ) {
      dae.a[i]=anga[i];
      dae.V[i]=Va[i];
    }
    printf("Error::new ton iteration for dynamic not converged\n");
  }
  else{
    h=h_old;
    t=tempo;
    k=k+1;
    if(k>settings.chunk)
      fm_out_1(t,k);
    fm_out_2(t,k);
  }
  printf("%lf\n",tempo);
  t_end=t;
  nrecord=k;
  nexttstep=h;
  delete []deltf;
  delete []ipiv;
  delete []fn;
  delete []inc;
  delete []xa;
  delete []anga;
  delete []Va;
}
/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  fm_nrlf
 *  Description:  solve power flow with locked ste variables
 *                itermax=max number of iteraions
 *                tol= convergence tolerance
 *                conv =1 if convergence reached ,0 otherwise
 *
 * =====================================================================================
 */
int Psat::fm_nrlf(int iter_max,double tol){
  int conv=1;
  int iteration=0;
  double *Vguess=new double [bus.n];
  double *y2=new double [2*bus.n]; 
  int *index=new int [2*bus.n-2*boundarynode.n];
  double *tempJlfv=new double [(2*bus.n-2*boundarynode.n)*(2*bus.n-2*boundarynode.n)];
  double *inc=new double [2*bus.n-2*boundarynode.n];
  for ( int i = 0; i < bus.n*bus.n; i += 1 ) {
    dae.J11[i]=0;
    dae.J12[i]=0;
    dae.J21[i]=0;
    dae.J22[i]=0;

  }
  for ( int i = 0; i < bus.n-boundarynode.n; i += 1 ) {
    index[i]=boundarynode.indexG[i];
    index[i+bus.n-boundarynode.n]=boundarynode.indexG[i]+bus.n;
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    inc[2*i]=1;
    inc[2*i+1]=1;
    Vguess[i]=1;
  }
  for ( int i = 0; i < sw.n; i += 1 ) {
    int k=sw.bus[i];
    Vguess[k]=dae.V[k];
  }
  for ( int i = 0; i < pv.n; i += 1 ) {
    int k=pv.bus[i];
    Vguess[k]=dae.V[k];
  }
  
  for ( int i = 0; i < bus.n-boundarynode.n; i += 1 ) {
    int k=boundarynode.indexG[i];
    dae.V[k]=Vguess[k];
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    y2[i]=dae.a[i];
    y2[i+bus.n]=dae.V[i];
  }
  double err_max=1;
  int *ipiv=new int[2*(bus.n-boundarynode.n)];
  while(err_max>tol&&iteration<iter_max){
    for(int i=0;i<bus.n;++i){
      dae.gp[i]=0;
      dae.gq[i]=0;
    }
    fm_lf_1();
    fm_mn_1();
    fm_syn(1);
    for ( int i = 0; i < bus.n; i += 1 ) {
      dae.g[i]=dae.gp[i];
      dae.g[i+bus.n]=dae.gq[i];
    }
    fm_lf_2();
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlf);
   fm_mn_2();
   fm_syn(2);
   cat4matrix(dae.J11,bus.n,bus.n,dae.J12,bus.n,bus.n,dae.J21,bus.n,bus.n,dae.J22,bus.n,bus.n,dae.Jlfv);
   getchar();
   for ( int i = 0; i < bus.island_n; i += 1 ) {
     int k=bus.island[i];
     for ( int j = 0; j < 2*bus.n; j += 1 ) {
       dae.Jlfv[j+k*2*bus.n]=0;
       dae.Jlfv[k+j*2*bus.n]=0;
       dae.Jlfv[j+(k+bus.n)*2*bus.n]=0;
       dae.Jlfv[k+bus.n+j*2*bus.n]=0;
       dae.Jlfv[k+k*2*bus.n]=0;
       dae.Jlfv[k+bus.n+(k+bus.n)*2*bus.n]=0;
       dae.g[k]=0;
       dae.g[k+bus.n]=0;
       dae.V[k]=0;
       dae.a[k]=0;
     }
   }
   for ( int i = 0; i < 2*(bus.n-boundarynode.n); i += 1 ) {
     int row=index[i];
     for ( int j = 0; j < 2*(bus.n-boundarynode.n); j += 1 ) {
       int col=index[j];
       tempJlfv[i+j*2*(bus.n-boundarynode.n)]=dae.Jlfv[col+row*2*(bus.n)];
     }
   }
   for (int i=0;i < 2*(bus.n-boundarynode.n); i +=1 ){
     int k=index[i];
     inc[i]=dae.g[k];
   }
   LAPACKE_dgesv(LAPACK_COL_MAJOR,2*bus.n-2*boundarynode.n,1,tempJlfv,2*bus.n-2*boundarynode.n,ipiv,inc,2*bus.n-2*boundarynode.n);
   for (int i=0;i < 2*(bus.n-boundarynode.n); i +=1 ){
     int k=index[i];
     y2[k]+=-inc[i];
   }
   for ( int i = 0; i < bus.n; i += 1 ) {
     dae.a[i]=y2[i];
     dae.V[i]=y2[i+bus.n];
   }
   err_max=0;
   for(int i=0;i<2*(bus.n-boundarynode.n);++i){
     if(abs(inc[i])>err_max)
       err_max=abs(inc[i]);
   }
   iteration++;
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    if(dae.V[i]<=1e-6)
      dae.a[i]=0;
  }
  for ( int i = 0; i < bus.n; i += 1 ) {
    if(abs(dae.V[i])<1e-4&&(abs(dae.a[i]))>2*pi){
      while(abs(dae.a[i])>2*pi){
	if( dae.a[i]<0)
	  dae.a[i]+=2*pi;
	else
	  dae.a[i]=dae.a[i]-2*pi;
      }
    }
  }
  if(iteration>=iter_max){
    printf("Solution of algebraic equations failed.\n");
    conv=0;
  }
  else{
    printf("Solution of algebraic equations completed in %d iterations",iteration);
  }
  fm_syn(3);
  delete []index;
  delete []Vguess;
  delete []inc;
  delete []ipiv;
  return conv;
}
