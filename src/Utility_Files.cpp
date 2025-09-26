#include "definitions.h"
#include "Globals.h"
#include "Utilities.h"

void Maximum(double &a, double &b, double &c, double &Max_V)
{
  if ((a >= b) && (a >= c))
    Max_V = a;
  else if ((b >= a) && (b >= c))
    Max_V = b;
  else
    Max_V = c;
}

void Minimum(double &a, double &b, double &c, double &Min_E)
{
  if ((a <= b) && (a <= c))
    Min_E = a;
  else if ((b <= a) && (b <= c))
    Min_E = b;
  else
    Min_E = c;
}

void Maximum(double &a, double &b, double &Max)
{
  if (a >= b)
    Max = a;
  else
    Max = b;
}

void Minimum(double &a, double &b, double &Min_E)
{
  if ((a <= b))
    Min_E = a;
  else
    Min_E = b;
}

void Print(V_I &vect)
{
  for (unsigned int index = 0; index < vect.size(); index++)
    cout << vect[index] << "\t";
  cout << endl;
}

void Print(V_D &Vect)
{
  int index = 0, size = Vect.size();
  //      cout<<size<<endl;
  //	cout<<"Iam being called"<<endl;
  switch (size)
  {
  case 25:
    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 5; j++)
      {
        index = i * 5 + j;
        cout << Vect[index] << "\t";
      }
      cout << endl;
    }
    break;
  case 15:
    for (int i = 0; i < 5; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        index = i * 3 + j;
        cout << Vect[index] << "\t";
      }
      cout << endl;
    }
    break;
  case 4:
    for (index = 0; index < size; index++)
    {
      cout << Vect[index] << "\t";
    }
    cout << endl;
    break;
  case 8:
    for (index = 0; index < size; index++)
    {
      cout << Vect[index] << "\t";
    }
    cout << endl;
    break;

  default:
    for (index = 0; index < size; index++)
      cout << Vect[index] << "\t";
    cout << endl;
    break;
  }
  // 			cout<<"--------------------------------------------------------\n";
}

void Print(vector<V_D> &Vect)
{
  //	cout<<"I am being called"<<endl;
  for (unsigned int i = 0; i < Vect.size(); i++)
  {
    //		cout<<i<<"\t";
    for (unsigned int index = 0; index < Vect[i].size(); index++)
      cout << Vect[i][index] << "\t";
    cout << endl;
  }
  //	cout<<"--------------------------------------------------------\n";
}

void Print(vector<bool> &Vect)
{
  for (unsigned int index = 0; index < Vect.size(); index++)
    cout << Vect[index] << "\t";
  cout << endl;
  cout << "--------------------------------------------------------\n";
}

// Resets a given vector to zero
void Vector_Reset(V_D &Vect)
{
  for (unsigned int i = 0; i < Vect.size(); i++)
    Vect[i] = 0.0;
  //	mu=0.0;
}

void Vector_Reset(vector<V_D> &Vect)
{
  for (unsigned int i = 0; i < Vect.size(); i++)
  {
    for (unsigned int j = 0; j < Vect[i].size(); j++)
    {
      Vect[i][j] = 0.0;
    }
  }

  //	mu=0.0;
}

void Print(Cell &cell)
{
  cout << "Cell ID: " << cell.cellID << "\t";
  cout << "Cell Type: " << cell.cellType << "\t";
  cout << "Cell Dimension: " << cell.Dimension << "\t";
  cout << "No of Boundary Faces: " << cell.NoBoundaryFaces << "\t";
  cout << "No of Faces: " << cell.numFaces << "\t";
  cout << "No of Nodes: " << cell.numNodes << "\t";
  cout << "No of Neighbours: " << cell.Neighbours.size() << "\n";
  cout << "Has Boundary Face: " << cell.hasBoundaryface << "\t";
  cout << "Left Face: " << cell.Left_Face << "\t";
  cout << "Right Face: " << cell.Right_Face << "\t";
  cout << "Top Face: " << cell.Top_Face << "\t";
  cout << "Bottom Face: " << cell.Bottom_Face << "\t";
  cout << "Interior Face: " << cell.Interior_Face << "\n";
  cout << "Node Indices: \t";
  for (int node : cell.nodeIndices)
  {
    cout << node << "\t";
  }
  cout << endl;
  cout << "   Cell Verticies: \t";
  if (cell.Cell_Vertices.size() % 3 != 0)
  {
    std::cerr << "Warning: Number of vertex entries is not a multiple of 3." << std::endl;
  }

  for (size_t i = 0; i < cell.Cell_Vertices.size(); i += 3)
  {
    cout << "("
         << cell.Cell_Vertices[i] << ", "
         << cell.Cell_Vertices[i + 1] << ", "
         << cell.Cell_Vertices[i + 2] << ") \t";
  }
  cout << endl;
  cout << "Cell Neighbours: ";
  for (int node : cell.Neighbours)
  {
    cout << node << "\t";
  }
  cout << endl;
  cout << "Face Normals: ";
  for (int i = 0; i < cell.Face_Normals.size(); i += 2)
  {
    cout << "( " << cell.Face_Normals[i] << " " << cell.Face_Normals[i + 1] << " )\t";
  }
  cout << endl;
  cout << "Face Areas: ";
  for (int i = 0; i < cell.Face_Areas.size(); i++)
  {
    cout << cell.Face_Areas[i] << "\t";
  }
  cout << endl;
  cout << ".......................................................";
  cout << endl;
}