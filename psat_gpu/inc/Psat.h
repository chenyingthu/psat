#ifndef _PSAT_H
#define _PSAT_H
#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))

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
#include "Solver.h"
#include "Simu.h"
#include "Record.h"
#include <math.h>
#include <complex>
//#include <cstdlib>
#include <iostream>
// #include <lapacke.h>
// #include <cblas.h>
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
	Solver solver;
	Simu simu;
	Record record;
	int *ipiv;
	int nBus,nFault,nLine,nPQ,nPV,nSW,nSyn,nMn,nShunt;
	double t_end,nexttstep;
	int nrecord;
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
	int fm_nrlf(int iter_max,double tol);

	void dyn_f_iniNet(int);
	void dyn_f_iniSyn(int);
	void dyn_f_iniDae(int);
	void dyn_f_iniSimu(int);
	void dyn_f_iniSolver(int);
	void dyn_f_store(int);
	void formDAEX();
	void DAEXtox();
	void DAEXtox_2(double *);
	void dyn_f_integration(int iFlag);
	void dyn_f_increaseTimeSteps(int);
	void dyn_f_dealFaults(int iFlag);
	void dyn_f_prediction(int iFlag);
	double * solver_jfng(double *);
	double *dyn_f_dae(double *,double );
	void updatePreconditioner(int);
	void pre_gmres(double *,double *,double *);
	double *step;
	double errstep;
	int inner_it_count;
	void precondition(double *,int );
	void dirder(double *,double * ,double *,double *);
	void givapp(double *,double *,double *,int);

	void debug(char* str,int n,double *);
	double *V_bak;
};
#endif
