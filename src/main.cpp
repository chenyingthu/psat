#include <stdio.h>
//#include "../inc/psat.h"
#include "../lib/myLib.h"
#include "../inc/Psat.h"
#include <lapacke.h>
#include <iostream>
#include <complex>
#include <stdio.h>
//typedef culaFloatComplex cuComplex;
//class cuComplex{
//public:
//	double r,i;
//	cuComplex(){}
//	cuComplex(double a,double b):r(a),i(b){}
//	double magnitude2(void){
//		return r*r+i*i;
//	}
//	cuComplex operator*(const cuComplex&a){
//		return cuComplex(r*a.r-i*a.i, i*a.r+r*a.i);
//	}
//	cuComplex operator+(const cuComplex &a) {
//		return cuComplex(r+a.r,i+a.i);
//	}
//	friend ostream &operator<<(ostream &output,cuComplex &d){
//	{
//		output<<"("<<d.r<<","<<d.i<<") ";
//		output<<endl;
//	}
//	return output;
//	}
//};
	//cuComplex &operator*(const cuComplex&a){
	//	return cuComplex(x*a.x-y*a.y, y*a.x+x*a.y);
	//}
	//cuComplex &operator+(const cuComplex &a) {
	//	return cuComplex(x+a.x,y+a.y);
	//}
	//ostream &operator<<(ostream &output,cuComplex &d){
	//	output<<"("<<d.x<<","<<d.y<<") ";
	//	output<<endl;
	//}
int dyn_main(int iFlag,int sys,int isPredict,int model,int isMultStep,int nSteps,double tStep);
int main()
{
//	culaInitialize();
	// Psat psat;
	// psat.specifySystem();
	// psat.init();
	// psat.formConMatrix();
	// psat.initialLF();
	// psat.fm_spf();
	// psat.fm_int_intial();
	// psat.resetBoundNode();
	// double h;
	// h=psat.fm_tstep(1,1,0,psat.settings.t0);
	// psat.fm_int_dyn(psat.settings.t0,psat.settings.tf,h);

	// return 0;

 dyn_main(1,39,0,2,0,2,0.005);
}
int dyn_main(int iFlag,int sys,int isPredict,int model,int isMultStep,int nSteps,double tStep){
  Psat psat;
  psat.dyn_f_iniNet(1);
  psat.dyn_f_iniSyn(1);
  psat.dyn_f_iniDae(1);
  psat.dyn_f_iniSolver(1);
  psat.dyn_f_iniSimu(1);
  return 0;
}
//int main()
//{
//     char JOBU;
//     char JOBVT;
//     int i;
//     //数据类型integer是fortran里的。这里在C++下可以使用的原因是f2c.h文件中已经作了定义
//     integer M = SIZE;
//     integer N = SIZE;
//     integer LDA = M;
//     integer LDU = M;
//     integer LDVT = N;
//     integer LWORK;
//     integer INFO;
//
//     integer mn = min( M, N );
//
//     integer MN = max( M, N );
//
//     double a[SIZE*SIZE] = { 16.0, 5.0, 9.0 , 4.0, 2.0, 11.0, 7.0 , 14.0, 3.0, 10.0, 6.0, 15.0, 13.0, 8.0, 12.0, 1.0};
//     double s[SIZE];
//     double wk[201];
//     double uu[SIZE*SIZE];
//     double vt[SIZE*SIZE];
//
//       JOBU = 'A';
//
//       JOBVT = 'A';
//
//    LWORK = 201;
//
///* Subroutine int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n,
//        doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
//        ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
//        integer *info)
//*/
//    dgesvd_( &JOBU, &JOBVT, &M, &N, a, &LDA, s, uu, &LDU, vt, &LDVT, wk, &LWORK, &INFO);
//
//    printf("INFO=%d /n", INFO );
//    for ( i= 0; i< SIZE; i++ ) {
//        printf("s[ %d ] = %f/n", i, s[ i ] );
//    }
//   // system("pause");
//    return 0;
//}
