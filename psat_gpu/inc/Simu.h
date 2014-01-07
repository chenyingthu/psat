//=====================================================================
// 
//=====================================================================
#ifndef SIMU_HEADER
#define SIMU_HEADER
class Simu{
  public:
  int multiSteps;
  int nSteps;
  double t_start;
  double t_end;
  double t_cur;
  double t_next;
  double tStep;
  int converged;
  int newton_iteration;
  int neval;
  double *t_switch;
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------
#endif
//=====================================================================
