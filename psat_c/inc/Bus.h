#ifndef _BUS_H
#define _BUS_H
class Bus{
public:
	double **con;
	int n;
	int island_n;
	int *internel;
	double *Pg;
	double *Qg;
	double *Pl;
	double *Ql;
	double *island;
	

	Bus();
	~Bus();

	void init();
	int check();
	void busDelete();
};
#endif
