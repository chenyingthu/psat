#ifndef _FAULT_H
#define _FAULT_H
#include "Bus.h"
#include <stdio.h>
class Fault
{
public:
	Fault();
	~Fault();

	/* data */

	double **con;
	int n;
	int *bus;
	double *dat;
	double *V;
	double *ang;
	double delta;
	double tFaultStart;
	double tFaultEnd;
	
	void init();
	void faultDelete();
	int check(Bus bus);
};
#endif