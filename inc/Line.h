#ifndef _LINE_H
#define _LINE_H
#include "../lib/myLib.h"
#include "settings.h"
#include "Shunt.h"
#include <complex>
#include <iostream>
using namespace std;
class Line{
public:
	double **con;
	int n;
	complex<double> *Y;
	double *Bp;
	double *Bpp;
	int *from;
	int *to;

	Line();
	~Line();

	void init(int n);
	int check(Settings settings);
	void lineDelete();
	void fm_y_line(int n,Shunt shunt);
};
#endif