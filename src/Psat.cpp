#include "../inc/Psat.h"
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
    //mn.Mndelete(pq.n);
   // varout.deleteVarout();
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
   
   

//   FILE	*fp;										/* output-file pointer */
//
//   fp	= fopen( "Jlfv", "w" );
//   if ( fp == NULL ) {
//     exit (EXIT_FAILURE);
//   }
//
//   for ( int i = 0; i < 2*bus.n; i += 1 ) {
//     for ( int j = 0; j < 2*bus.n; j += 1 ) {
//       fprintf(fp,"%lf\t",dae.Jlfv[j+i*2*bus.n]);
//     }
//     fprintf(fp,"\n");
//   }
//   if( fclose(fp) == EOF ) {			/* close output file   */
//     exit (EXIT_FAILURE);
//   }

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
     memset(dae.f,1,(dae.n-dynordold)*sizeof(double));
     memset(dae.x,1,(dae.n-dynordold)*sizeof(double));
     memset(dae.Fx,0,dae.n*dae.n*sizeof(double));
     memset(dae.Fy,0,dae.n*2*bus.n*sizeof(double));
     memset(dae.Gx,0,2*bus.n*dae.n*sizeof(double));
   }
   settings.init=1;

   getchar();
   delete []deltf_nosw;
   delete []tmp_nosw;
   delete []deltf;
   delete []tmp;


}
void Psat::resetBoundNode(){
 int *indexAll=new int [dae.n+2*bus.n]; 
 int *indexG=new int [bus.n];
 
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
     dae.J22[i]=dS[i].imag()+1e-6;
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
      dae.J11[i]=dS[i].real()+1e-6;
      dae.J21[i]=dS[i].imag();
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
  printf("%d\n",dae.n);
}
void Psat::fm_synit(){}
void Psat::fm_excin(){}
void Psat::fm_syn(int flag){}
void Psat::fm_syn1(int flag){}
void Psat::fm_syn2(int flag){}
void Psat::fm_syn3(int flag){}
void Psat::fm_syn4(int flag){}
void Psat::fm_int_intial(){}
void Psat::fm_mn_1(){}
void Psat::fm_mn_2(){}
void Psat::fm_mn_0(){}
void Psat::fm_fault(int flag,double t){}
void Psat::fm_fault_0(double t){}
void Psat::fm_fault_1(double t){}
double Psat:: fm_tstep(int flag,int convergency,int iteration,double t){return 0.0;}
void Psat::fm_tstep_1(int convergency,int iteration,double t){}
void Psat::fm_tstep_2(int convergency,int iteration,double t){}
void Psat::fm_out_0(double t,int k){}
void Psat::fm_out_1(double t,int k){}
void Psat::fm_out_2(double t,int k){}
void Psat::fm_out_3(double t,int k){}
void Psat::fm_int_dyn(double t0,double tf,double h){}
void Psat::fm_int_step(double t,double h,double *tempi,int k){}
