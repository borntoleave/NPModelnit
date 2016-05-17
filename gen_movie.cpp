//rna_length is not the lines of rna beads!!!!!!!
#include<iostream>
#include<fstream>
#include<cstdlib>
#include<string>
#include<cstring>
#include<iomanip>
#define RNA_LENGTH 2000		
#define PROTEIN_LENGTH 2000	//for C_ALPHA
#define TIS_LENGTH RNA_LENGTH*3		//for P,B,S representation
#define CGM_LENGTH PROTEIN_LENGTH+TIS_LENGTH
using namespace std;

struct coords
{
double x;
double y;
double z;
}coords[CGM_LENGTH];

class pdbline
{
public:
char line[2048],buffer[10];
char atom[5],         atom_type[5],atom_name[4],chain_id[2],                             tail[20];
int          atom_num,                                      resi_num;
double                                                                x, y, z, a_occ,b_f;
//PDB<| 4       7    1     4      1     3      1    1          4    4 8  8  8     6   6     19 |>
//idx:      4       12           17           21          22        30 38 46 54    60  66  
pdbline()
	{
	atom_num=resi_num=x=y=z=a_occ=b_f=0;
	}
int read(std::ifstream &raw_pdb)
	{
	 raw_pdb.getline(line,2048);
	 	strncpy(atom,		line,	4);atom[4]		='\0';														//ATOM
		if(strcmp(atom,"ATOM"))return 1;
		if(!strcmp(atom,"END ")||!strcmp(atom,"TER "))return 2;
  	 if(!strcmp(atom,"ATOM"))
  	 {
		strncpy(buffer,		line+4,	7);	buffer[7]	='\0';atom_num	=atoi(buffer);								//atom_num	
   		strncpy(atom_type,	line+12,4);	atom_type[4]='\0';														//atom_type
   		strncpy(atom_name,	line+17,3);	atom_name[3]='\0';														//atom_name
		strncpy(chain_id,	line+21,1);	chain_id[1]	='\0';														//chain_id
		strncpy(buffer,		line+22,4);	buffer[4]	='\0';resi_num	=atoi(buffer);								//resi_num
		strncpy(buffer,		line+30,8);	buffer[8]	='\0';x		=atof(buffer);									//x
		strncpy(buffer,		line+38,8);	buffer[8]	='\0';y		=atof(buffer);									//y
		strncpy(buffer,		line+46,8);	buffer[8]	='\0';z		=atof(buffer);									//z
		strncpy(buffer,		line+54,6);	buffer[6]	='\0';a_occ	=atof(buffer);									//atom_occupancy
		strncpy(buffer,		line+60,6);	buffer[6]	='\0';b_f	=atof(buffer);									//beta_factor
		strncpy(tail,			line+66,19);tail[19]	='\0';														//tail info
	}
	return 0;
	}
void write()
	{
	cout	<<atom<<setw(7)<<atom_num<<setw(5)<<atom_type<<setw(4)<<atom_name<<setw(2)<<chain_id<<setw(4)<<resi_num<<setiosflags(ios::fixed)<<setprecision(3)
			<<setw(12)<<x<<setw(8)<<y<<setw(8)<<z//<<setw(6)<<a_occ<<setw(6)<<b_f<<tail
			<<endl;
	}
void fwrite(std::ofstream &out_pdb)
	{
	out_pdb	<<atom<<setw(7)<<atom_num<<setw(5)<<atom_type<<setw(4)<<atom_name<<setw(2)<<chain_id<<setw(4)<<resi_num<<setiosflags(ios::fixed)<<setprecision(3)
			<<setw(12)<<x<<setw(8)<<y<<setw(8)<<z//<<setw(6)<<a_occ<<setw(6)<<b_f<<tail
			<<endl;
	}
};

int rna_length,protein_length,cgm_length;
pdbline cgm[CGM_LENGTH];


void print_psf(char domain);

int main(int argn, char* argv[])
{
if(argn!=4&&argn!=5){cerr<<"Format: "<<argv[0]<<" <CGM pdb> <trajectory coordinates> <output pdb> <region[p,r,a]>"<<endl;exit(-1);}
int frame=0;double x,y,z,x0,y0,z0;char domain;
if (argn==4)domain='a';else domain=argv[4][0];

string buffer;
ifstream tis(argv[1],ios::in);
ifstream traj(argv[2],ios::in);
ofstream movie(argv[3],ios::out);
while (tis.good())
{int token=cgm[cgm_length].read(tis);
if(token==2)break;
if(token==1)continue;
if(!strcmp(cgm[cgm_length].atom_type," CA "))protein_length++;
cgm_length++;
}
rna_length=(cgm_length-protein_length)/3+1;
cout <<protein_length<<" amino acids and "<<rna_length<<" nucleotides. "<<cgm_length<<" beads in total."<<endl;
while (traj.good())
{
x0=y0=z0=0;
for(int n=0;n<cgm_length;n++)
	{
	traj>>coords[n].x>>coords[n].y>>coords[n].z;
	if(domain=='r'&&n<protein_length)continue;
	else if(domain=='p'&&n>=protein_length)continue;
	else {x0+=coords[n].x;y0+=coords[n].y;z0+=coords[n].z;}
	}
frame++;
if(!(frame%1))
	{printf("\b\b\b\b\b\b%d",frame);
		if(domain=='r'){x0/=(cgm_length-protein_length);y0/=(cgm_length-protein_length);z0/=(cgm_length-protein_length);}
		else if(domain=='p'){x0/=protein_length;y0/=protein_length;z0/=protein_length;}
		else {x0/=cgm_length;y0/=cgm_length;z0/=cgm_length;}
	for(int n=0;n<cgm_length;n++)
		{
		if(domain=='r'&&n<protein_length)continue;
		else if(domain=='p'&&n>=protein_length)break;
		else {cgm[n].x=coords[n].x-x0;cgm[n].y=coords[n].y-y0;cgm[n].z=coords[n].z-z0;
		cgm[n].fwrite(movie);}
		}
	movie <<"TER\nEND"<<endl;
	}
		
}
cout<<endl;
cerr <<"Movie file is "<<argv[3]<<endl;
print_psf(domain);//print corresponding psf file
tis.clear();tis.close();
traj.clear();traj.close();
movie.clear();movie.close();
}

void print_psf(char domain)	//print out the PSF info
{
char outname[60];int bond=0,bonds[CGM_LENGTH][2];
int cgm_l,start,end,shift;
sprintf(outname,"./movie.psf");
ofstream out(outname,ios::out);
out	<<"PSF CMAP"<<endl<<endl
		<<setw(8)<<0<<" !NTITLE"<<endl<<endl;
if(domain=='r'){start=protein_length;end=cgm_length;}
else if(domain=='p'){start=0;end=protein_length;}
else {start=0;end=cgm_length;}
//--------------atoms----------------------------------
out	<<setw(8)<<end-start<<"!NATOM"<<endl;		
for(int cgm_l=start;cgm_l<end;cgm_l++)
out <<setw(8)<<cgm_l+1-start<<setw(5)<<cgm[cgm_l].chain_id<<' '<<setiosflags(ios::left)<<setw(4)<<cgm[cgm_l].resi_num
		<<resetiosflags(ios::left)<<setiosflags(ios::right)<<setw(4)<<cgm[cgm_l].atom_name<<setw(5)<<cgm[cgm_l].atom_type<<setw(5)<<cgm[cgm_l].atom_type
		<<setiosflags(ios::fixed)<<setprecision(6)<<setw(12)<<0.000000<<setprecision(4)<<setw(14)<<12.0110<<setw(12)<<0<<endl<<resetiosflags(ios::right);
out <<endl;
//---------------bonds-------------------------------
if(domain!='r')//p or a
	for(int l=0;l<protein_length-1;l++)
		{bonds[bond][0]=l+1;bonds[bond][1]=l+2;bond++;}
if(domain!='p')//r or a
	{if(domain=='r')shift=0;else shift=protein_length;
	for(int l=0;l<rna_length;l++)
		{
		if(l)
		{bonds[bond][0]=l*3-2+shift;bonds[bond][1]=l*3+shift;bond++;
		bonds[bond][0]=l*3+shift;bonds[bond][1]=l*3+1+shift;bond++;}
		bonds[bond][0]=l*3+1+shift;bonds[bond][1]=l*3+2+shift;bond++;
		}
	}
out	<<setw(8)<<bond<<" !NBOND: bonds"<<endl;
for(int b=0;b<bond;b++)
{if(!(b%4)&&b)out<<endl;
out <<setw(8)<<bonds[b][0]<<setw(8)<<bonds[b][1];
}
cerr <<"PSF file is ./movie.psf"<<endl;
out.clear();out.close();
}

