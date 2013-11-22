#include "../inc/Psat.h"
Psat::Psat(){
    jay=Complex(1,0);
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
    if(!dae.check_2(sw))
        return 0;
    if(!sw.check(bus))
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
    idx_n=1;
    double *VL1=new double [line.n];
    double *kt=new double [line.n];
    double *V1=new double [bus.n];
    double *V2=new double [bus.n];
    delete []VL1;
    delete []kt;
    delete []V1;
    delete []V2;
}
void Psat::fm_base_sw(){}
void Psat::fm_base_pq(){}
void Psat::fm_base_pv(){}
void Psat::fm_base_syn(){}
void Psat::fm_y(){}
void Psat::fm_dynlf(){}
void Psat::fm_spf(){}
void Psat::resetBoundNode(){}
void Psat::fm_lf_1(){}
void Psat::fm_pq_1(){}
void Psat::fm_pv_1(){}
void Psat::fm_sw_1(){}
void Psat::fm_lf_2(){}
void Psat::fm_pq_2(){}
void Psat::fm_pv_2(){}
void Psat::fm_sw_2(){}
void Psat::fm_dynidx(){}
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
