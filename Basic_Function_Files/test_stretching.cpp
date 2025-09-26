#include<iostream>
#include<vector>
#include<cmath>

using namespace std;

int main()
{
vector<double> x,x1;
int n=21;
double eta,zi,beta =-2.0;
eta=(1.0/(n-1));
zi=(1.0/(n-1));
x.resize(n,0.0);
x1.resize(n,0.0);
for(int i=0;i<n;i++)
{
	x[i]=9.5*zi*i;
	cout<<i<<"\t"<<x[i]<<"\t";
	x[i]=9.5*((exp(beta*zi*i)-1)/(exp(beta)-1));
	cout<<x[i]<<endl;
}
return 0;
}