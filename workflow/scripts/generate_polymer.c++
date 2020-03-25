#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<fstream>
#include<sstream>
#include<vector>
#include <unistd.h>
#include<ctime>


int main(int argc, char* argv[]){
std::cerr << "TYPE: 1.filename 2.N" << std::endl;

int N = atoi(argv[1]);
int mode = atoi(argv[2]);

double pos[N][3];

std::cout <<"LAMMPS data file from restart file: timestep = 0, procs = 1" << std::endl;
std::cout << std::endl;

std::cout << 1*N << " atoms"<< std::endl;
std::cout << 1*N-1 << " bonds"<< std::endl;
std::cout << 1*N-2 << " angles"<< std::endl;
std::cout << std::endl;

std::cout << 4 << " atom types" << std::endl;
std::cout << 1 << " bond types" << std::endl;
std::cout << 1 << " angle types" << std::endl;
std::cout << std::endl;

std::cout << "-50 50"  << " xlo xhi "<< std::endl;
std::cout << "-50 50"  << " ylo yhi "<< std::endl;
std::cout << "-50 50"  << " zlo zhi "<< std::endl;
std::cout << std::endl;

std::cout << "Masses "<< std::endl;
std::cout << std::endl;
std::cout << " 1 1 " << std::endl;
std::cout << " 2 2 " << std::endl;
std::cout << " 3 1 " << std::endl;
std::cout << " 4 2 " << std::endl;
std::cout << std::endl;

std::cout << "Atoms "<< std::endl;
std::cout << std::endl;

int n;

// Build array of atom type
int type_array[N];
int type;

for (int i=0; i<N; i++) {
	// Mode 0 = Constant atom type (1)
	// Mode 1 = Random atom type (1 - 4)
	// Other mode = Clustered atom type - distribute 4 types equally
	if (mode == 0) {
		type = 1;
	} else if (mode == 1) {
		type = (rand() % 4) + 1;
	} else {
		if (i < N * 0.25) {
			type = 1;
		} else if (i < N * 0.5) {
			type = 2;
		} else if (i < N * 0.75) {
			type = 3;
		} else {
			type = 4;
		}
	}
	type_array[i] = type;
}

for(n=0;n<N;n++){
double r=1.1;

	if(n==0){
	pos[n][0]=rand()*1.0/RAND_MAX;
	pos[n][1]=rand()*1.0/RAND_MAX;
	pos[n][2]=rand()*1.0/RAND_MAX;
    }

    if(n>0){
    double phi=rand()*1.0/RAND_MAX*2.0*M_PI;
    double theta=rand()*1.0/RAND_MAX*M_PI;
    pos[n][0]=pos[n-1][0]+r*sin(theta)*cos(phi);
    pos[n][1]=pos[n-1][1]+r*sin(theta)*sin(phi);
    pos[n][2]=pos[n-1][2]+r*cos(theta);
    }
    //NEED TO std::cout :
    //INDEX MOLECULE TYPE X Y Z IX IY IZ
    std::cout << n+1 << " 1 " << type_array[n] << " " << pos[n][0] << " " <<  pos[n][1] << " "  << pos[n][2] << " 0 0 0 " << std::endl;
}

//CREATE BONDS BETWEEN BEADS
//THIS IS A LINEAR POLYMER SO N-1 BONDS
int nbonds=1;
std::cout << "Bonds"<< std::endl;
std::cout << std::endl;
for(int n=1;n<N;n++){
std::cout << nbonds << " " << 1  << " " << n << " " << n%N+1 << std::endl;
nbonds++;
}

std::cout << std::endl;

//CREATE ANGLES BETWEEN BEADS
//THIS IS A LINEAR POLYMER SO N-2 ANGLES
int nangles=1;
std::cout << "Angles"<< std::endl;
std::cout << std::endl;
for(int n=0;n<N-2;n++){
if(n<N-2)std::cout << nangles << " " << 1  << " " << n+1<< " " << n+2 << " " << n+3 << std::endl;
if(n==N-2)std::cout << nangles << " " << 1  << " " << n+1 << " " << n+2 << " " << 1 << std::endl;
if(n==N-1)std::cout << nangles << " " << 1  << " " << n+1 << " " << 1 << " " << 2 << std::endl;
nangles++;
}

return 0;
}
