#ifndef _DAE_H
#define _DAE_H
#include "Bus.h"
#include "SW.h"
#include <stdio.h>
class DAE{
public:

	double *a;
	double *V;
	int kg;
	int n;
	int npf;
	double *g;
	double *gp;
	double *gq;
	double *f;
	double *x;
	double *glfp;
	double *glfq;
	double *J11;
	double *J12;
	double *J21;
	double *J22;
	double *Jlf;
	double *Jlfv;
	double *Jlfd;
	double *Fx;
	double *Fy;
	double *Gx;
	double *Ac;
	double *tn;
	double t;

	DAE();
	~DAE();
	void init();
	int check(Bus bus);
	void daeDelete();
	int check_2(SW sw);
};
#endif