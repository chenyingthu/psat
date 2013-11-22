#ifndef _VAROUT_H
#define _VAROUT_H
class Varout
{
public:
	Varout();
	~Varout();
	void deleteVarout();
	/* data */
	int numOfStep;
	double *t;
	double *f;
	double *x;
	double *V;
	double *ang;
	double *Pm;
	double *Vf;
};


#endif // VAROUT_H