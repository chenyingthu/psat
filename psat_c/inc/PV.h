#ifndef _PV_H
#define _PV_H
#include "Bus.h"
#include "Clpsat.h"
#include "DAE.h"
class PV
{
public:
	
	PV();
	~PV();
	/* data */
	double **con;
	int n;
	int *bus;
	double *pq;
	double **store;

	void init();
	void pvDelete();
	int check(Bus bus,Clpsat clpsat,DAE dae);
	void fm_synit(int _n,int *_bus);
};
#endif