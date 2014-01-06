//=====================================================================
// 
//=====================================================================
#ifndef SOLVER_HEADER
#define SOLVER_HEADER
//---------------------------------------------------------------------
//---------------------------------------------------------------------
class Newton{
  public:
  double tol;
  int maxit;
  int debug;
  int method;
};
class Gmres{
  public:
  double tol_fixed;
  double tol;
  int maxit;
  int reorth;
  double *h;
  double *z;
  double *v;
  double *c;
  double *s;
};
class Precond{
  public:
  int inner;
  int outer;
  int innerStop;
  int outerStop;
};
class Jfng{
  public:
  Newton newton;
  Gmres gmres;
  Precond precond;

  double findiff;

};
class Update{
  public:
  int num;
  int last_num;
  int maxnum;
  int dest;
};
class Solver{
  public:
  Jfng jfng;
  Update update;

  double *P;
  double *deltx;
  double *delty;
  void dyn_f_iniSolver(int);
};
#endif
//=====================================================================
