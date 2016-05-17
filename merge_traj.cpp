#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cstdlib>
#include<string>
#include<time.h>
#define PI 3.1415926	//math Pi
#define MOL_BEAD 300
#define BOX_NUM 200
using namespace std;
class coords
{	public:
	double c[3];			//{x, y, z} = {c[0], c[1], c[2]}
	coords()				//constructor without parameters
	{
		for(int d=0;d<3;d++)
			c[d]=0;
	}
	coords(double radius)	//constructor with parameters, creates uniformly distributed points on sphere
	{
		for(int d=0;d<3;d++)
			c[d]=2*(double)rand()/RAND_MAX-1;
		double r=sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]);
		for(int d=0;d<3;d++)	//http://mathworld.wolfram.com/SpherePointPicking.html, equation (16)
			c[d]=radius*c[d]/r;
	}
	void transform(double phi,double theta,double *shift)
	{
		double temp[3];
		double ang=0;
		double RMx[3][3]={{1,       0,		  0},	//rotation matrix for x axis
						  {0,cos(ang),-sin(ang)},
						  {0,sin(ang), cos(ang)}};
		ang=theta;
		double RMy[3][3]={{ cos(ang),0,sin(ang)},	//rotation matrix for rotate theta about y axis
						  {		   0,1,		  0},
						  {-sin(ang),0,cos(ang)}};
		ang=phi;
		double RMz[3][3]={{ cos(ang),-sin(ang),0},	//rotation matrix for rotate phi about z axis
						  { sin(ang), cos(ang),0},
						  {		   0,		 0,1}};

		temp[0]=c[0];temp[1]=c[1];temp[2]=c[2];
		for(int i=0;i<3;i++)
		{
			c[i]=0;
			for(int j=0;j<3;j++)
				c[i]+=RMy[i][j]*temp[j];			//rotate theta about y axis
		}

		temp[0]=c[0];temp[1]=c[1];temp[2]=c[2];
		for(int i=0;i<3;i++)
		{
			c[i]=0;
			for(int j=0;j<3;j++)
				c[i]+=RMz[i][j]*temp[j];			//rotate phi about z axis
		}

		for(int d=0;d<3;d++)						//shift along a vector
			c[d]+=shift[d];

	}
};

int main(int argc, char *argv[])
{
	if(argc!=3)
	{	
		cerr<<"Format: "<<argv[0]
			<<" <trajectory common title> <division n on theta(0~Pi)>"<<endl;
		exit(-1);
	}
	ifstream np_coords(argv[1],ios::in);
	int div=atoi(argv[2]);
	double 	NpR 	=200.0;		//nano particle radius
	double ang_conf[BOX_NUM][2];//(phi,theta)
	double xyz_conf[BOX_NUM][3];//(x,y,z) also represents the shift vector
	int box_num=0;
	double d_theta=(double)PI/div;
	double theta,phi;
	for(double theta=0;theta<=PI+0.1;theta+=d_theta)// +0.1 fix the rounding error 
	{							//generate uniformed grids on sphere
		double d_phi;
		if(theta==0||theta==PI)
		{
			d_phi=2*PI;
		}	
		else
		{	
			d_phi=2*PI/round(2*PI*sin(theta)/d_theta);
		}
			//cerr <<d_theta<<'\t'<<d_phi<<endl;
		for(double phi=0;phi<=PI+0.1;phi+=d_phi)
		{
			ang_conf[box_num][0]=phi;
			ang_conf[box_num][1]=theta;
			xyz_conf[box_num][0]=NpR*sin(theta)*cos(phi);
			xyz_conf[box_num][1]=NpR*sin(theta)*sin(phi);
			xyz_conf[box_num][2]=NpR*cos(theta);
			//cerr<<phi<<"\t"<<theta<<endl;
			box_num++;
		}
	}
	cerr<<"Need "<<box_num<<" simulation boxes"<<endl;
	coords xyz;
	for(int n=0;n<box_num;n++)
	{	
		cerr<<"processing box "<<n+1<<endl; 
		char traj_orig_fname[50];
		char traj_tran_fname[50];
		sprintf(traj_orig_fname,"%s_orig_%d",argv[1],n);
		sprintf(traj_tran_fname,"%s_tran_%d",argv[1],n);//change the file name for real case

		//I'm using "./mol_coords.xyz" as original trajectory
		sprintf(traj_orig_fname,"%s",argv[1]);
		sprintf(traj_tran_fname,"TranCoords/%s_tran_%d.xyz",argv[1],n+1);
		cerr<<" Opened "<<traj_orig_fname<<". Write to "<<traj_tran_fname<<endl;
		ifstream in(traj_orig_fname,ios::in);
		ofstream out(traj_tran_fname,ios::out);
		while(in.good())
		{
			in>>xyz.c[0]>>xyz.c[1]>>xyz.c[2];
			xyz.c[2]-=NpR;
			xyz.transform(ang_conf[n][0],ang_conf[n][1],xyz_conf[n]);
			out<<xyz.c[0]<<'\t'<<xyz.c[1]<<'\t'<<xyz.c[2]<<'\t'<<endl;
		}
		in.close();
		out.close();
	}	

}