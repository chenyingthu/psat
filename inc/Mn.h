#ifndef _MN_H
#define _MN_H
#include "Bus.h"
class Mn{
	public:
	double **con;
	int n;
	int *bus;
	int init_n;
	int *init;
	double **store;
	Mn();
	void Mninit(int n);
	void Mndelete(int n);
    void check(Bus _bus);
};
#endif
