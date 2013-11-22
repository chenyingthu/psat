#ifndef _SHUNT_H
#define _SHUNT_H
#include "Bus.h"
#include "settings.h"
#include <set>
class Shunt{
public:
	double **con;
	int *bus;
	double *g;
	double *b;
	int n;

	Shunt();
	void init(int n);
	void Shuntdelete(int n);
	void check(Bus bus,Settings settings);
};
#endif