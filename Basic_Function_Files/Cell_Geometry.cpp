#include "../Basic_Files/definitions.h"
/*
Terminology used in this code
Cell_Q ------- represents the vector for storing rho, rhou, rhov,rhow and rhoEt
Cell_Flux------ Total Flux for a given cell
*/

// Constructor for cell, no points => Cell fromed by this constructor will be used for ghost cells
Cell::Cell()
{
	Cell_Q.resize(5,0.0);
	Cell_Flux.resize(5,0.0);
	Viscous_Flux.resize(5,0.0);
	Q_new.resize(5,0.0);
	Del2_Q_1.resize(5,0.0);	Del2_Q_2.resize(5,0.0);	Del2_Q_3.resize(5,0.0);
	Del3_Q_1.resize(5,0.0);	Del3_Q_2.resize(5,0.0);	Del3_Q_3.resize(5,0.0);
	Del4_Q_1.resize(5,0.0);	Del4_Q_2.resize(5,0.0);	Del4_Q_3.resize(5,0.0);
	Del5_Q_1.resize(5,0.0);	Del5_Q_2.resize(5,0.0);	Del5_Q_3.resize(5,0.0);
	Del6_Q_1.resize(5,0.0);	Del6_Q_2.resize(5,0.0);	Del6_Q_3.resize(5,0.0);
	No_of_Faces=6;
	Cell_Rho=0.0;Cell_Pressure=0.0;Cell_Et=0.0;Cell_Temparature=0.0;
	is_Ghost_Cell = true;
}
//Destructor for cell
Cell::~Cell(){ }
/*
overloaded () operator for creating a cell. It takes 6 points and constructs a triangular prism with 5 faces. used along
the axis of the pipe so the bottom face will not be present. Front face corresponds to inlet, back face corresponds to
exit .  Upon calling this function with required number of points it creates five Face objects with the points
calculates their area and unit normals on the faces gets area magnitude, calculates cell center, cell volume.
*/
void Cell::operator()(const Point& o,const Point& a,const Point& c,const Point& o1,const Point& a1,const Point& c1)
{
//  	cout<<"Cell with 5 faces called\t";
	cell_nop = 6;
	No_of_Faces=5;
	Front_Face(o,c,a,0);	Back_Face(o1,a1,c1,1); 	Left_Face(o,o1,c1,c,2); 	Right_Face(o,a,a1,o1,3); 	Top_Face(a,c,c1,a1,4);
	Volume();
	FF_Area = Front_Face.Get_Face_Area();BaF_Area=Back_Face.Get_Face_Area();
	TF_Area = Top_Face.Get_Face_Area();LF_Area=Left_Face.Get_Face_Area();RF_Area=Right_Face.Get_Face_Area();
	Average_total_area = (FF_Area+BaF_Area+TF_Area+LF_Area+RF_Area+BoF_Area)/No_of_Faces;
	Cell_Avg_Length=volume/Average_total_area;
	inv_vol =1.0/volume;
	Cell_Center(o,a,c,o,o1,a1,c1,o1);
 	Check_Cell();
	is_Ghost_Cell = false;
	Diagonal_Vector(o,a1);
//  	cout<<"Cell Created with 5 faces\n";
}

/*
overloaded () operator for creating a cell. It takes 8 points and constructs a hexahedra prism with 6 faces.
Front face corresponds to inlet, back face corresponds to exit .
 Upon calling this function with 8 points it creates six Face objects  with their magnitude,  calculates cell center, cell volume
*/
void Cell::operator()(const Point& o,const Point& a,const Point& b,const Point& c,const Point& o1,const Point& a1,const Point& b1,const Point& c1)
{
//    	cout<<"Cell with 6 faces called\n";
	cell_nop=8;
	No_of_Faces=6;
// 	o.Print();a.Print();b.Print();c.Print(); o1.Print();a1.Print();b1.Print();c1.Print();
	Front_Face(o,c,b,a,0);	Back_Face(o1,a1,b1,c1,1);
	Top_Face(a,b,b1,a1,4);	Bottom_Face(o,o1,c1,c,5);
	Left_Face(c,c1,b1,b,2);	Right_Face(o,a,a1,o1,3);
	Volume();
	FF_Area = Front_Face.Get_Face_Area();		BaF_Area=Back_Face.Get_Face_Area();
	TF_Area = Top_Face.Get_Face_Area();			BoF_Area=Bottom_Face.Get_Face_Area();
	LF_Area=Left_Face.Get_Face_Area();			RF_Area=Right_Face.Get_Face_Area();
	Average_total_area = (FF_Area+BaF_Area+TF_Area+LF_Area+RF_Area+BoF_Area)/No_of_Faces;
	Cell_Avg_Length=volume/Average_total_area;
	inv_vol =1.0/volume;
	Cell_Center(o,a,b,c,o1,a1,b1,c1);
	Check_Cell();
	is_Ghost_Cell = false;
	Diagonal_Vector(o,b);
//  	cout<<"Cell Created with 6 faces\n";
}

//********************Indexes of Neighbouring Cells****************************************//
void Cell::operator()(const int& i0,const int & i1,const int & i2,const int & i3,const int & i4,const int & i5,const int & i6)
{
// 	cout<<"In cell neighbours assigning\n";
	self_index=i0;// corresponds to its own index
	index1=i1; // left cell
	index2=i2; // right cell
	index3=i3; // top cell
	index4=i4; // bottom cell
	index5=i5;// back cell
	index6=i6;// front cell
//  	cout<<self_index<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<index5<<"\t"<<index6<<endl;
}
//************Calculates Volume of the Cell based of Gauss Theorem***********************
void Cell::Volume()
{
        if(cell_nop==6)
        {
		volume = ((1.0/3.0)*(Front_Face.Get_RmiddotA()+Back_Face.Get_RmiddotA()
							+Top_Face.Get_RmiddotA()+Left_Face.Get_RmiddotA()
							+Right_Face.Get_RmiddotA()));
	}
        else
        {
		volume =((1.0/3.0)*(Front_Face.Get_RmiddotA()+Back_Face.Get_RmiddotA()
						+Top_Face.Get_RmiddotA()+Bottom_Face.Get_RmiddotA()
						+Left_Face.Get_RmiddotA()+Right_Face.Get_RmiddotA()));
        }
//       cout<<"cell volume\t"<<volume<<endl;
}

//***************************************Calculates Mid point of the Cell*********************************
void Cell::Cell_Center(const Point& p1,const Point & p2,const Point & p3,const Point & p4,const Point & p5,const Point & p6,const Point & p7,const Point & p8)
{
        if(cell_nop==6)
        {
		Cell_MidPoint = (p1+p2+p3+p5+p6+p7)/6.0;
        }
        else
        {
		Cell_MidPoint = (p1+p2+p3+p4+p5+p6+p7+p8)/8.0;
        }
      //  cout<<"Calculating Cell center"<<endl;
}

//returns the mid point of the cell
const Point& Cell::Get_Cell_Center() const
{
	return Cell_MidPoint;
}

const Vector & Cell::Get_Diagonal_Vector() const
{
	return Diagonal_Vector;
}
//this function returns area magnitude of the face required
const double& Cell::Cell_Face_Area(const int& i) const
{
	switch(i)
	{
		case 0:
			return FF_Area;
		case 1:
			return BaF_Area;
		case 2:
			return LF_Area;
		case 3:
			return RF_Area;
		case 4:
			return TF_Area;
		case 5:
			return BoF_Area;
		default :
			cout<<"Enter 0 -5 for fetching Area magnitude, default returing Front_Face Magnitude\n";
			return FF_Area;
	}
}

// returns number of faces for a given cell
const int& Cell::no_of_Faces() const
{
	return No_of_Faces;
}


//returns the correspoding face object
const Face& Cell::Get_Face(const int& i) const 
{
	switch(i)
	{
		case 0:
			return Front_Face;
		case 1:
			return Back_Face;
		case 2:
			return Left_Face;
		case 3:
			return Right_Face;
		case 4:
			return Top_Face;
		case 5:
			return Bottom_Face;
		default :
			cout<<"Enter 0 - 5 for fetching Faces returning Front Face as default Please check the face number\n";
			return Front_Face;
	}
}

// returns 1/Cell_Volume for a given cell
const double& Cell::Get_Volume() const
{
	return inv_vol;
}
//tells us whether a cell is a ghost cell or a physical cell
const bool& Cell::Is_Ghost_Cell() const
{
	return is_Ghost_Cell;
}
// prints indexs of neighbouring cells on std out -- terminal
void Cell::Print_Neighbours() const
{
	cout<<self_index<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<index5<<"\t"<<index6<<endl;
}

const double& Cell::Get_Celldx() const
{
	return Cell_Avg_Length;
}
//returns the neighbouring cell indicies based on facenumber
const int& Cell::Get_Neighbours(const int& faceno) const
{
	switch(faceno)
	{
		case 0:
			return index6;	// front face neighbour
		case 1:
			return index5;	//back face neighbour
		case 2:
			return index1;	// left face neighbour
		case 3:
			return index2;	// right face neighbour
		case 4:
			return index3;	// top face neighbour
		case 5:
			return index4;	// bottom face neighbour
		default :
			return self_index; // self
	}
}

//This function when called  sets points to neighbouring cells 
void Cell::Set_NeighbourCells(const vector<Cell*>* const Cell_List)
{
// 	cout<<self_index<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<index5<<"\t"<<index6<<"\n";
	Left_Cell=Cell_List->at(index1);
	Right_Cell=Cell_List->at(index2);
	Top_Cell=Cell_List->at(index3);
	Bottom_Cell=Cell_List->at(index4);
	Back_Cell=Cell_List->at(index5);
	Front_Cell=Cell_List->at(index6);
}

Cell* Cell::Get_Cell(const int & i) 
{
	switch(i)
	{
		case 0:
			return Front_Cell;	// front face neighbour
		case 1:
			return Back_Cell;	//back face neighbour
		case 2:
			return Left_Cell;	// left face neighbour
		case 3:
			return Right_Cell;	// right face neighbour
		case 4:
			return Top_Cell;	// top face neighbour
		case 5:
			return Bottom_Cell;	// bottom face neighbour
		default :
			return this; // self
	}
}

//writes the geometrical features to a file
void Cell::Write_Cell_Info(const string& ipfile) const
{
	ofstream fileout(ipfile.c_str(),ios::out|ios::app);
	Vector Temp;
	if(fileout.is_open())
	{
                fileout<<"Cell index\t\t"<<self_index<<"\tNumber of Faces\t"<<No_of_Faces<<endl;
		fileout<<"Mid point of the cell\t"<<Cell_MidPoint.Get_x()<<"\t"<<Cell_MidPoint.Get_y()<<"\t"<<Cell_MidPoint.Get_z()<<endl;
		fileout<<"Volume of the cell\t"<<volume<<endl;
		fileout<<"Cell Area Vectors\n";
		Temp = Front_Face.Get_FaceArea_Vector();
		fileout<<"Front Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp = Back_Face.Get_FaceArea_Vector();
		fileout<<"Back Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp =Top_Face.Get_FaceArea_Vector();
		fileout<<"Top Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		if(No_of_Faces==6)
		{
		    Temp = Bottom_Face.Get_FaceArea_Vector();
		    fileout<<"Bottom Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		}
		Temp = Left_Face.Get_FaceArea_Vector();
		fileout<<"Left Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp = Right_Face.Get_FaceArea_Vector();
		fileout<<"Right Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		fileout<<"Cell Normals\n";
		Temp = Front_Face.Get_FaceNormal();
		fileout<<"Front Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp = Back_Face.Get_FaceNormal();
		fileout<<"Back Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp =Top_Face.Get_FaceNormal();
		fileout<<"Top Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		if(No_of_Faces==6)
		{
		    Temp = Bottom_Face.Get_FaceNormal();
		    fileout<<"Bottom Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		}
		Temp = Left_Face.Get_FaceNormal();
		fileout<<"Left Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		Temp = Right_Face.Get_FaceNormal();
		fileout<<"Right Face\t"<< Temp(1)<<"\t"<<Temp(2)<<"\t"<<Temp(3)<<endl;
		fileout<<self_index<<"\t"<<index1<<"\t"<<index2<<"\t"<<index3<<"\t"<<index4<<"\t"<<index5<<"\t"<<index6<<endl;
		fileout<<"------------------------------------------------------------------------------------\n";
	}
	else
	{
		cout<<"Could not open file to write\n";
	}
}

void Cell::Check_Cell() const
{
	Vector Area_Vector;
  	Area_Vector = Front_Face.Get_FaceArea_Vector() + Back_Face.Get_FaceArea_Vector()
			 		+Left_Face.Get_FaceArea_Vector()  + Right_Face.Get_FaceArea_Vector()
					+ Top_Face.Get_FaceArea_Vector();
	if(No_of_Faces==6)
	    Area_Vector += Bottom_Face.Get_FaceArea_Vector();
// 	cout<<"Total Area\t";Area_Vector.Print();
	if(volume<0.0)
		cout<<"Index\t"<<self_index<<"\t"<<volume<<"\tCell with Negative Volume\n";
}
