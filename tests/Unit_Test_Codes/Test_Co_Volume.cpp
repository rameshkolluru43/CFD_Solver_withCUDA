#include "definitions.h"


int main()
{
	V_D P1(3,0.0),P2(3,0.0),P3(3,0.0),P4(3,0.0);
	double A1=0.0,A2=0.0,A3=0.0, Sum_of_normals = 0.0, Sum_of_areas = 0.0;
	
//Cell covering all quadrants
        
    P1[0] = 0.0; P1[1] = -1.0;P1[2] = 0.0;
	P2[0] = 1.0; P2[1] = 0.0;P2[2] = 0.0;
	P3[0] = 0.0; P3[1] = 1.0;P3[2] = 0.0;
	P4[0] = -1.0; P4[1] = 0.0;P4[2] = 0.0;
	
	Construct_Cell(P1,P2,P3,P4);

        
return 0;
}
