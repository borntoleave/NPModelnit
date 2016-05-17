#include<cstdlib>
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<string>
#define RADIUS 100
#define PI 3.14159
#define DIST_RATIO 1.5
using namespace std;
class coords
{public:
	double x;
	double y;
	double z;
};


int main(int argc, char *argv[])
{
if (argc!=3)
	{
	cerr<<"Format: "<<argv[0]<<" <surface bead #> <molecule #>"<<endl;
	exit(-1);}
int surf_num=atoi(argv[1]);
int mole_in_num=atoi(argv[2]);
int mole_out_num=mole_in_num*DIST_RATIO*DIST_RATIO;
double u,v,the,phi;//http://mathworld.wolfram.com/SpherePointPicking.html
coords *surf_coord=new coords [surf_num];
coords *mole_in_coord=new coords [mole_in_num];
coords *mole_out_coord=new coords [mole_out_num];
cerr<<"NP: "<<surf_num<<"\tinner: "<<mole_in_num<<"\touter: "<<mole_out_num<<endl;
for(int sn=0;sn<surf_num;sn++)
	{
	u=(double)rand()/RAND_MAX;
	v=(double)rand()/RAND_MAX;
	phi=2*PI*u;
	the=acos(2*v-1);
	surf_coord[sn].x=RADIUS*sin(the)*cos(phi);
	surf_coord[sn].y=RADIUS*sin(the)*sin(phi);
	surf_coord[sn].z=RADIUS*cos(the);
	cout<<surf_coord[sn].x<<'\t'<<surf_coord[sn].y<<'\t'<<surf_coord[sn].z<<endl;
	}
for(int mn=0;mn<mole_in_num;mn++)
	{
	u=(double)rand()/RAND_MAX;
	v=(double)rand()/RAND_MAX;
	//cout<<u<<endl;
	phi=2*PI*u;
	the=acos(2*v-1);
	mole_in_coord[mn].x=DIST_RATIO*RADIUS*sin(the)*cos(phi);
	mole_in_coord[mn].y=DIST_RATIO*RADIUS*sin(the)*sin(phi);
	mole_in_coord[mn].z=DIST_RATIO*RADIUS*cos(the);
	cout<<mole_in_coord[mn].x<<'\t'<<mole_in_coord[mn].y<<'\t'<<mole_in_coord[mn].z<<endl;
	}
for(int mn=0;mn<mole_out_num;mn++)
	{
	u=(double)rand()/RAND_MAX;
	v=(double)rand()/RAND_MAX;
	//cout<<u<<endl;
	phi=2*PI*u;
	the=acos(2*v-1);
	mole_out_coord[mn].x=DIST_RATIO*DIST_RATIO*RADIUS*sin(the)*cos(phi);
	mole_out_coord[mn].y=DIST_RATIO*DIST_RATIO*RADIUS*sin(the)*sin(phi);
	mole_out_coord[mn].z=DIST_RATIO*DIST_RATIO*RADIUS*cos(the);
	cout<<mole_out_coord[mn].x<<'\t'<<mole_out_coord[mn].y<<'\t'<<mole_out_coord[mn].z<<endl;
	}

}
