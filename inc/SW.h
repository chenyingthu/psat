#ifndef _SW_H
#define _SW_H
#include "Bus.h"
#include <stdio.h>
class SW{
public:
	double **con;
	int n;
	int *bus;

	SW();
	~SW();

	void init();
	void swDelete();
	int check(Bus bus);
	void fm_synit(int _n,int *_bus);
};
#endif