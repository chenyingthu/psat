#ifndef _SETTINGS_H
#define _SETTINGS_H
class Settings {
public:
	int dac;
	int pq2z;
	int pv2pq;
	int showlf;
	int init;
	int conv;
	int red;
	int plot;
	int plottype;
	int method;
	int show;
	int vs;
	int ok;
	int pfsolver;
	double deltat;
	double deltamax;
	double deltamin;
	int chunk;
	int mv;
	int iter;
	int isStatic;
	int zoom;
	double freq;
	int beep;
	double dyntol;
	int dynmit;
	double lftol;
	double lfmit;
	double lftime;
	double t0;
	double tf;
	double mva;
	double rad;
	int distrsw;
	int refbus;
	int color;
	int fixt;
	double tstep;
	int locksnap;
	int octave;
	int local;
	int absvalues;
	int shuntvalues;
	int violations;
	double *tempi;
	int dyn_isPredict;
	int dyn_predict_model;
	int dyn_isMulStep;
	int dyn_MulStep_nSteps;
	int dyn_lf;
	double dyn_tStep;
	double dyn_tStep_min;
	double dyn_tStep_max;
	char exportType[20];
	Settings();
	~Settings();
	
};
#endif
