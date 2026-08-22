#include"../Basic_Function_Files/headers.hpp"


vector<double>& Face::Flux_From_Face_CUSP(Cell * CC, Cell * NC)
{
	Flux[0]=0.0;Flux[1]=0.0;Flux[2]=0.0;Flux[3]=0.0;Flux[4]=0.0;
	return Flux;
}