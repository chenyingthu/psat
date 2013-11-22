#ifndef _MN_H
#define _MN_H
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

};
#endif 