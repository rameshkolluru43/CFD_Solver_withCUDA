#ifndef MPI_UTILS_H
#define MPI_UTILS_H

#include <vector>

void CFD_MPI_Initialize(int *argc, char ***argv);
void CFD_MPI_Finalize();
bool CFD_MPI_Is_Enabled();
bool CFD_MPI_Is_Root();
int CFD_MPI_Rank();
int CFD_MPI_Size();

bool CFD_MPI_Owns_Cell(int cellIndex);
void CFD_MPI_Build_Local_Leaf_Cell_List(std::vector<int> &leafCells);

double CFD_MPI_Global_Min(double value);
double CFD_MPI_Global_Max(double value);
void CFD_MPI_Global_Sum4(double values[4]);

void CFD_MPI_Synchronize_Solution_State();
void CFD_MPI_Barrier();

#endif // MPI_UTILS_H
