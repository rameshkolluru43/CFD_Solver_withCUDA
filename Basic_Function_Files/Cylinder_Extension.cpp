void Cylinder::Identify_Cells_3h()
{
	int p=0;
	for(int k=0;k<nlength-1;k++)
	{
		for(int j=0;j<nr-2;j++)
		{
			for(int i=0;i<ntheta-1;i++)
			{
				cout<<i+j*(ntheta-1)+k*(ntheta-1)*(nr-1)<<"\t";
				a = p*8;
				if(i==0)
				{
					if(j==0)
					{
						Cell_Point_Data_3h[a+0]=(i-1)+(ntheta-1)+j*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+1]=(i-1)+(ntheta-1)+(j+3)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+4]=(i-1)+(ntheta-1)+j*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+5]=(i-1)+(ntheta-1)+(j+3)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
					else
					{
						Cell_Point_Data_3h[a+0]=(i-1)+(ntheta-1)+(j-1)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+1]=(i-1)+(ntheta-1)+(j+2)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+4]=(i-1)+(ntheta-1)+(j-1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+5]=(i-1)+(ntheta-1)+(j+2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
				}
				else
				{
					if(j==0)
					{
						Cell_Point_Data_3h[a+0]=(i-1)+j*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+1]=(i-1)+(j+3)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+4]=(i-1)+j*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+5]=(i-1)+(j+3)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
					else
					{
						Cell_Point_Data_3h[a+0]=(i-1)+(j-1)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+1]=(i-1)+(j+2)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+4]=(i-1)+(j-1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+5]=(i-1)+(j+2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
				}
				if(i==(ntheta-3)||(i==(ntheta-2)))
				{
					if(j==0)
					{
						Cell_Point_Data_3h[a+3]=(i+2)-(ntheta-1)+j*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+2]=(i+2)-(ntheta-1)+(j+3)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+7]=(i+2)-(ntheta-1)+j*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+6]=(i+2)-(ntheta-1)+(j+3)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
					else
					{
						Cell_Point_Data_3h[a+3]=(i+2)-(ntheta-1)+(j-1)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+2]=(i+2)-(ntheta-1)+(j+2)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+7]=(i+2)-(ntheta-1)+(j-1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+6]=(i+2)-(ntheta-1)+(j+2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
				}
				else
				{
					if(j==0)
					{
						Cell_Point_Data_3h[a+3]=(i+2)+j*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+2]=(i+2)+(j+3)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+7]=(i+2)+j*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+6]=(i+2)+(j+3)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
					else
					{
						Cell_Point_Data_3h[a+3]=(i+2)+(j-1)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+2]=(i+2)+(j+2)*(ntheta-1)+k*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+7]=(i+2)+(j-1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
						Cell_Point_Data_3h[a+6]=(i+2)+(j+2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
					}
				}
// cout<<Cell_Point_Data_3h[a+0]<<"\t"<<Cell_Point_Data_3h[a+1]<<"\t"<<Cell_Point_Data_3h[a+2]<<"\t"<<Cell_Point_Data_3h[a+3]<<"\t"<<Cell_Point_Data_3h[a+4]<<"\t"<<Cell_Point_Data_3h[a+5]<<"\t"<<Cell_Point_Data_3h[a+6]<<"\t"<<Cell_Point_Data_3h[a+7]<<endl;
			p++;
			}
// cout<<"--------------------------------------\n";
		}
		int j=0;
		for(int i=0;i<(ntheta-1);i++)
		{
			j=nr-2;
			cout<<i+j*(ntheta-1)+k*(ntheta-1)*(nr-1)<<"\t";
			if(i==0)
			{
				Cell_Point_Data_3h[a+0]=(i-1)+(ntheta-1)+(j-2)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+1]=(i-1)+(ntheta-1)+(j+1)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+4]=(i-1)+(ntheta-1)+(j-2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+5]=(i-1)+(ntheta-1)+(j+1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
			}
			else
			{
				Cell_Point_Data_3h[a+0]=(i-1)+(j-2)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+1]=(i-1)+(j+1)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+4]=(i-1)+(j-2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+5]=(i-1)+(j+1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
			}
			if((i==(ntheta-2))||(i==(ntheta-3)))
			{
				Cell_Point_Data_3h[a+3]=(i+2)-(ntheta-1)+(j+1)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+2]=(i+2)-(ntheta-1)+(j-2)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+6]=(i+2)-(ntheta-1)+(j+1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+7]=(i+2)-(ntheta-1)+(j-2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
			}
			else
			{
				Cell_Point_Data_3h[a+3]=(i+2)+(j-2)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+2]=(i+2)+(j+1)*(ntheta-1)+k*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+6]=(i+2)+(j+1)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
				Cell_Point_Data_3h[a+7]=(i+2)+(j-2)*(ntheta-1)+(k+1)*(ntheta-1)*nr;
			}
// cout<<Cell_Point_Data_3h[a+0]<<"\t"<<Cell_Point_Data_3h[a+1]<<"\t"<<Cell_Point_Data_3h[a+2]<<"\t"<<Cell_Point_Data_3h[a+3]<<"\t"<<Cell_Point_Data_3h[a+4]<<"\t"<<Cell_Point_Data_3h[a+5]<<"\t"<<Cell_Point_Data_3h[a+6]<<"\t"<<Cell_Point_Data_3h[a+7]<<endl;
			p++
			}
// 			cout<<"--------------------------------------\n";
}
cout<<"Identifying 3h cells done"<<endl;
}