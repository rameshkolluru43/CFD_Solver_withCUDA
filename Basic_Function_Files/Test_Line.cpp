#include "../Basic_Function_Files/Geometry_Header.h"


int main()
{

	Point P1,P2;
	Line L1,L2;
	int Nop=31;
	bool bothsides =true;
	
	P1(0.0,0.0,0.0);
	P2(0.0,0.0,1.0);
	
	L1.generate(P1,P2,Nop);
	L1.Print();
	
	L2.generate(P1,P2,Nop,bothsides);
	L2.Print();

	return 0;

}