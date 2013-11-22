#include "../inc/Psat.h"
Psat::Psat(){
    jay=Complex(1,0);
}
Psat::~Psat(){}
void Psat::specifySystem(){}
void Psat::init(){}
void Psat::formConMatrix(){}
void Psat::initialLF(){}
int Psat::fm_ncomp(){return 0;}
void Psat::fm_base_line(){}
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
