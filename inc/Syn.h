#ifndef _SYN_H
#define _SYN_H
#include <stdio.h>
#include "Bus.h"
#include "DAE.h"
#include "settings.h"
#include "PV.h"
#include "SW.h"
#include "../lib/myLib.h"
#include <complex>
class Syn
{
public:
	Syn();
	~Syn();

	/* data */
	complex<double> jay;
	double **con;
	int n;
	int *bus;
	double *Id;
	double *Iq;
	double *Pg;
	double *Qg;
	double *J11;
	double *J12;
	double *J21;
	double *J22;
	int *delta_idx;
	int *omega_idx;
	int *e1q_idx;
	int *e2q_idx;
	int *e1d_idx;
	int *e2d_idx;
	int *psiq_idx;
	int *psid_idx;
	double *delta;
	double *omega;
	double *e1q;
	double *e2q;
	double *e1d;
	double *e2d;
	double *psiq;
	double *psid;
	double *pm;
	double *vf;
	double *Gp;
	double *Gq;
	double *c1;
	double *c2;
	double *c3;
	int *ord;
	int *is2;
	int is2_n;
	int *is3;
	int is3_n;
	int *is4;
	int is4_n;
	int *is51;
	int is51_n;
	int *is52;
	int is52_n;
	int *is53;
	int is53_n;
	int *is6;
	int is6_n;
	int *is8;
	int is8_n;

	int *bs2;
	int *bs3;
	int *bs4;
	int *bs51;
	int *bs52;
	int *bs53;
	int *bs6;
	int *bs8;
	void init();
	void synDelete();
	int check(Bus bus);
	int fm_dynidx(DAE dae);
	void fm_synit(Settings settings,DAE dae,Bus bus);
	void fm_synit_mod2(int is2_n,int *is2,double *vd,double *vq);
	void fm_synit_mod3(int is3_n,int *is3,double *vd,double *vq,DAE dae);
	void fm_synit_mod4(int is4_n,int *is4,double *vd,double *vq,DAE dae);
	void fm_synit_mod51(int is51_n,int *is51,double *vd,double *vq);
	void fm_synit_mod52(int is52_n,int *is52,double *vd,double *vq);
	void fm_synit_mod53(int is53_n,int *is53,double *vd,double *vq);
	void fm_synit_mod6(int is6_n,int *is6,double *vd,double *vq);
	void fm_synit_mod8(int is8_n,int *is8,double *vd,double *vq);
};
#endif
