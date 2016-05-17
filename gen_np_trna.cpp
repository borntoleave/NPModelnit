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
using namespace std;
//extern double rm[3][3];	//the best way to pass the matrix. but failed 
//void rotmat(char axis,double ang)	//define rotation matrix
//{
//switch (axis)
//	{
//		case 'x':
//		{
//			double rm[3][3]={{1,       0,		 0},	//rotation matrix for x axis
//							 {0,cos(ang),-sin(ang)},
//							 {0,sin(ang), cos(ang)}};
//			break;
//		}
//		case 'y':
//		{
//			double rm[3][3]={{ cos(ang),0,sin(ang)},	//rotation matrix for y axis
//							 {		  0,1,		 0},
//							 {-sin(ang),0,cos(ang)}};
//			break;
//		}
//		case 'z':
//		{
//			double rm[3][3]={{cos(ang),-sin(ang),0},	//rotation matrix for z axis
//							 {sin(ang), cos(ang),0},
//							 {		 0,		   0,1}};
//			break;
//		}
//	}
//}

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
	if (argc!=3)
	{
		cerr<<"Format: "<<argv[0]<<" <layer> <molecules per layer>"<<endl;
		exit(-1);
	}
	int layer_num	=atoi(argv[1]);
	int m_per_layer	=atoi(argv[2]);
	double 	AuR		=1.46;		//gold atom radius
	double 	NpR 	=200.0;		//nano particle radius
	double 	MoL 	=60.0;		//pdb's longest length
	double 	BoxL 	=4.0*MoL; 	//simulation box size
	double 	NpSurfA	=4.0*PI*NpR*NpR; 	//nano particle surface area
	int 	AuN		=round(NpSurfA/AuR/AuR/50);	//gold atom number
	double 	ThetaRange	=asin(BoxL/sqrt(2)/NpR);//not used

	ofstream surf_coords("./surf_coords.xyz",ios::out);				//save the coordinates of Au atoms on the surface
	ofstream surf_inbox_coords("./surf_inbox_coords.xyz",ios::out);	//save the coordinates of Au atoms within simulation box
	
	srand (time(NULL));			//change random seed
	//cerr<<time(NULL)<<endl;	//show random seed

	//vvvvvvvvvvvv  generate nano particle surface  vvvvvvvvvvvvvvvvvvv 
	int surf_inbox_num=0;
	for (int a=0;a<AuN;a++)		//generate the NP coordinates in simulation box
	{	
		coords NanoP(NpR);		//create uniformly distributed points on sphere
		surf_coords<<setprecision(4)<<NanoP.c[0]<<'\t'<<NanoP.c[01]<<'\t'<<NanoP.c[2]<<endl;
		if(fabs(NanoP.c[0])<=BoxL/2&&fabs(NanoP.c[1])<=BoxL/2&&NanoP.c[2]>0)	//pick out coordinates within box
		{
			surf_inbox_coords<<setprecision(4)<<NanoP.c[0]<<'\t'<<NanoP.c[1]<<'\t'<<NanoP.c[2]<<endl;
			surf_inbox_num++;
		}
	}
	surf_coords.close();
	surf_inbox_coords.close();
	//^^^^^^^^^^^^^  generate nano particle surface  ^^^^^^^^^^^^^^^^^^^^

	//vvvvvvvvvvvvv  generate pdb  vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	coords pdb[MOL_BEAD];
	int bead_num=0;
	double mean[3]={0,0,0};
	ifstream mol_pdb("./fmet_coords.xyz");
	while(mol_pdb.good())
	{
		mol_pdb>>pdb[bead_num].c[0]>>pdb[bead_num].c[1]>>pdb[bead_num].c[2];
		for(int d=0;d<3;d++)
			mean[d]+=pdb[bead_num].c[d];
		bead_num++;
	}
	--bead_num;
	mol_pdb.close();
	for(int d=0;d<3;d++)
		mean[d]=(double)mean[d]/bead_num;
	for(int b=0;b<bead_num;b++)
	{	
		for(int d=0;d<3;d++)
			pdb[bead_num].c[d]-=mean[d];
	}

	ofstream mol_coords("./mol_coords.xyz",ios::out);
	for (int l=1;l<=layer_num;l++)
	{
		for (int m=0;m<m_per_layer;m++)
		{	
			coords *Molecule=new coords[bead_num];
			for(int b=0;b<bead_num;b++)
				Molecule[b]=pdb[b];
			double rand_theta=(double)PI*rand()/RAND_MAX;
			double rand_phi=(double)2*PI*rand()/RAND_MAX;
			double shift[3]=
			{
				MoL*cos(2*PI*m/m_per_layer),
				MoL*sin(2*PI*m/m_per_layer),
				NpR+1.5*l*MoL
			};
			for(int b=0;b<bead_num;b++)
			{
				Molecule[b].transform(rand_phi,rand_theta,shift);
				mol_coords<<setprecision(4)<<Molecule[b].c[0]<<'\t'<<Molecule[b].c[1]<<'\t'<<Molecule[b].c[2]<<endl;
			}
		}
	}
	mol_coords.close();
	
	cerr<<"total beads on surface:\t"<<AuN<<endl
		<<"=> beads of inbox surface:\t"<<surf_inbox_num<<endl
		<<"total molecule chain:\t"<<layer_num*m_per_layer<<endl
		<<"beads of one molecule:\t"<<bead_num<<endl
		<<"=> beads of all molecules:\t"<<bead_num*layer_num*m_per_layer<<endl;
}
