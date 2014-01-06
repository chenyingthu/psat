#ifndef _BOUNDARYNODE_H
#define _BOUNDARYNODE_H
#include <stdio.h>
class BoundaryNode{
public:
	int n;
	int *bnode;
	double *Voltage;
	double *Angle;
	double *P;
	double *Q;
	int *indexG;
	int *indexAll;

	BoundaryNode();
	~BoundaryNode();
	void init(int bus_n,int dae_n);
	void boundaryNodeDelete();
};
#endif