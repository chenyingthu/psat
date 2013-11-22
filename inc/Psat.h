#ifndef _PSAT_H
#define _PSAT_H
//extern "C"
#include <stdio.h>
#include <set>
#include "Bus.h"
#include "DAE.h"
#include "Fault.h"
#include "Line.h"
#include "PQ.h"
#include "PV.h"
#include "settings.h"
#include "SW.h"
#include "Syn.h"
#include "Clpsat.h"
#include "BoundaryNode.h"
#include "Mn.h"
#include "Varout.h"
#include "Shunt.h"
#include <cmath>
#include <complex>
//#include <cstdlib>
#include <iostream>
#include <lapacke.h>
typedef  complex<double> Complex;
using namespace std;
class Psat{
public:
	Bus bus;
	DAE dae;
	Fault fault;
	Line line;
	PQ pq;
	PV pv;
	Settings settings;
	SW sw;
	Syn syn;
	Mn mn;
	Shunt shunt;
	Varout varout;
	Clpsat clpsat;
	BoundaryNode boundarynode;
	int *ipiv;
	int nBus,nFault,nLine,nPQ,nPV,nSW,nSyn,nMn,nShunt;
	complex<double> jay;
	complex<double>alpha;
	complex<double>beta;
	Psat();
	~Psat();
	void specifySystem();
	void init();
	void formConMatrix();
	void initialLF();
	int fm_ncomp();
	void fm_base();
	void fm_base_line();
	void fm_base_sw();
	void fm_base_pq();
	void fm_base_pv();
	void fm_base_syn();
	void fm_base_mn();
	void fm_y();
	void fm_dynlf();
	void fm_spf();
	void resetBoundNode();
	void fm_lf_1();
	void fm_pq_1();
	void fm_pv_1();
	void fm_sw_1();
	void fm_lf_2();
	void fm_pq_2();
	void fm_pv_2();
	void fm_sw_2();
	void fm_dynidx();
	void fm_synit();
	void fm_excin();
	void fm_syn(int flag);
	void fm_syn1(int flag);
	void fm_syn2(int flag);
	void fm_syn3(int flag);
	void fm_syn4(int flag);
	void fm_int_intial();
	void fm_mn_1();
	void fm_mn_2();
	void fm_mn_0();
	void fm_fault(int flag,double t);
	void fm_fault_0(double t);
	void fm_fault_1(double t);
	double fm_tstep(int flag,int convergency,int iteration,double t);
	void fm_tstep_1(int convergency,int iteration,double t);
	void fm_tstep_2(int convergency,int iteration,double t);
	void fm_out_0(double t,int k);
	void fm_out_1(double t,int k);
	void fm_out_2(double t,int k);
	void fm_out_3(double t,int k);
	void fm_int_dyn(double t0,double tf,double h);
	void fm_int_step(double t,double h,double *tempi,int k);
};
#endif
