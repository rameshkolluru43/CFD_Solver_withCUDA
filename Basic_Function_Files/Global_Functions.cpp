#include "../Basic_Function_Files/headers.hpp"

vector<double>  Q_Mod(5,0.0);

vector<double> & operator*(vector<double>& Q,const double & alpha)
{
	for(unsigned int i=0;i<Q.size();i++)
		Q_Mod[i] = Q[i]*alpha;
	return Q_Mod;
}

vector<double> & operator+(vector<double>&Q1,vector<double>&Q2)
{
	if(Q1.size()==Q2.size())
	{
		for(unsigned int i=0;i<Q1.size();i++)
			Q_Mod[i] = Q1[i]+Q2[i];
		return Q_Mod;
	}
	else
	{
		cout<<"Sizes of the passed STL vectors are different cant add these two STL vectors, returning Zero by default\n";
		return Q_Mod;
	}
}

vector<double> & operator-(vector<double>&Q1,vector<double>&Q2)
{
	if(Q1.size()==Q2.size())
	{
		for(unsigned int i=0;i<Q1.size();i++)
			Q_Mod[i] = Q1[i]-Q2[i];
		return Q_Mod;
	}
	else
	{
		cout<<"Sizes of the passed STL vectors are different cant subtract these two STL vectors, returning Zero by default\n";
		return Q_Mod;
	}
}