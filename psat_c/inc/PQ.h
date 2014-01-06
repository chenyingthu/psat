#ifndef _PQ_H
#define _PQ_H
#include "Bus.h"
#include "Clpsat.h"
#include <stdio.h>
class PQ
{
public:
	PQ();
	~PQ();

	/* data */
	double **con;
	int n;
	int *bus;
	double *P0;
	double *Q0;
	double **store;

	void init();
	void pqDelete();
	int check(Bus bus,Clpsat clpsats);
};
#endif