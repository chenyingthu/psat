#include "../inc/Varout.h"
Varout::Varout(){
}
Varout::~Varout(){
}
void Varout::deleteVarout(){
	delete []t;
	delete []f;
	delete []x;
	delete []V;
	delete []ang;
	delete []Pm;
	delete []Vf;
}