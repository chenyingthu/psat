#include "../inc/BoundaryNode.h"
BoundaryNode::BoundaryNode()
{
	n=0;
}
BoundaryNode::~BoundaryNode()
{
}
void BoundaryNode::boundaryNodeDelete()
{
	if(indexAll!=NULL)
		;//delete []indexAll;
	if(indexG!=NULL)
		delete []indexG;
	if (n==0)
		return ;
}
void BoundaryNode::init(int bus_n,int dae_n)
{
	
	if(n==0)
		return ;
	bnode=new int [n];
	Voltage=new double [n];
	Angle=new double [n];
	P=new double [n];
	Q=new double [n];
	for (int i=0;i<n;++i){
		Voltage[i]=0;
		Angle[i]=0;
		P[i]=0;
		Q[i]=0;
	}

}