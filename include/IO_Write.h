#ifndef IO_WRITE_H
#define IO_WRITE_H
#include "definitions.h"
#include "Globals.h"

/*-------------------Initialize,Output and Printing Functions--------------------*/
void Write_Cell_Info(const int &);
void Write_Solution(const string &, const int &);
/// void Write_Solution(const string &);
void Read_Write_Grid(const string &, const string &);
void Append_Solution(const string &, const string &);
void Write_VTK_File(const string &, const string &);
void Write_VTK_File(const string &);
void Write_Error_File(const string &);
void Write_Limiter_File(const string &);
void Write_CF_File(const string &);
void Write_A_MatrixToFile(vector<V_D> &, const string &);
void Write_b_VectorToFile(V_D &, const string &);
void createOutputDirectories();
void findGridFiles(string &, int &, int &, vector<string> &);
void searchGridFiles();
/*------------------------------------------------------------*/
#endif // #ifndef IO_WRITE_H