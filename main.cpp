#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <complex>
#define PI 3.141592653589793
#include <cmath>
using namespace std;
#include "myvector.h"
#include "mylist.h"
#include "matrixAlgebra.h"
#include "atom.h"
#include <libconfig.h++>
using namespace libconfig;


#define bufsize 1000
FILE *out;
FILE *in;
FILE *print_h_out;

double qe = 1.6021766209e-19; //Coulomb
double epsilon0 = 8.854187817e-12; //farads per metre

const double a = 0.5431; //1
const int n = 1; //20;
const double volumeSize = n*a;
const int N = n*n*n*8; //n*n*n*8;
#define BOUNDARYCOEFF 1.0   //li

V<double> k(0,0,0);
complex<double> li(0,1.0);

#define PRINT_H
#define CSR
#define HERMITIAN_TEST
#define SO
#define BOYKIN2002
//#define BOYKIN2010
#define UI
#define CSRBINARY
#define CSRTEXT

#define BANDSTRUCTURE
#define DEFORMX 0.0
#define DEFORMY 0.0
#define DEFORMZ 0.0

//{ //init TB parameters
double SPS		= 0.0;
double PDS		= 0.0;
double S_PS		= 0.0;
double SDS		= 0.0;
double PDP		= 0.0;
double S_DS		= 0.0;
double SSS		= 0.0;
double S_SS		= 0.0;
double SS_S		= 0.0;
double S_S_S	= 0.0;
double PPS		= 0.0;
double PPP		= 0.0;
double DDS		= 0.0;
double DDP		= 0.0;
double DDD		= 0.0;

double a_lattice = 0.0;
double bond_length = 0.0;
double Es = 0.0;
double Ep = 0.0;
double Es_ = 0.0;
double Ed = 0.0;
double delta = 0.0;

double strain_exponent_V_S_S_Sigma = 0.0;
double strain_exponent_V_S_P_Sigma = 0.0;
double strain_exponent_V_Sstar_P_Sigma = 0.0;
double strain_exponent_V_S_Sstar_Sigma = 0.0;
double strain_exponent_V_Sstar_S_Sigma = 0.0;
double strain_exponent_V_Sstar_Sstar_Sigma = 0.0;
double strain_exponent_V_Sstar_D_Sigma = 0.0;
double strain_exponent_V_S_D_Sigma = 0.0;
double strain_exponent_V_P_P_Sigma = 0.0;
double strain_exponent_V_P_P_Pi = 0.0;
double strain_exponent_V_P_D_Sigma = 0.0;
double strain_exponent_V_P_D_Pi = 0.0;
double strain_exponent_V_D_D_Sigma = 0.0;
double strain_exponent_V_D_D_Pi = 0.0;
double strain_exponent_V_D_D_Delta = 0.0;

double strain_constant_K_S_S = 0.0;
double strain_constant_K_S_P = 0.0;
double strain_constant_K_Sstar_P = 0.0;
double strain_constant_K_Sstar_S = 0.0;
double strain_constant_K_S_Sstar = 0.0;
double strain_constant_K_Sstar_D = 0.0;
double strain_constant_K_Sstar_Sstar = 0.0;
double strain_constant_K_S_D = 0.0;
double strain_constant_K_P_D = 0.0;
double strain_constant_K_P_P = 0.0;
double strain_constant_K_D_D = 0.0;

double Z_eff = 0.0;
double r_P_D = 0.0;
double r2_P_P = 0.0;
double r2_D_D = 0.0;
double r3_P_D = 0.0;
double r4_D_D = 0.0;

double Energy_shift = 0.0;
//}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

int alex_fgets (char* s, int n, FILE* in) {
	int i,f;
	char c;
	if (n<2) return 0;
	i=0;
	f=0;
	while (1) {
		if (i>=n-1) break;
		c=fgetc(in);
		if (c==EOF) break;
		if (c=='\r') continue;
		f=1;
		if (c=='\n') break;
		s[i]=c;
		i++;
	}
	s[i]=0;
	return f;
}
/*
complex<double> slaterkoster(int alpha, List<Atom>* i, int beta, List<Atom>* j) {
	V<double> v1(i->d->x, i->d->y, i->d->z);
	V<double> v2;
	v2.x = j->d->x;
	v2.y = j->d->y;
	v2.z = j->d->z;

	bool ifSIcurrentNeighbour = i->d->sostav;  //ifSI(v2);
	bool ifSIcurrentPsi = j->d->sostav;
	if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == false)) {
		//GE-GE
		SPS = 2.73135;
		PDS = -2.00115;
		S_PS = 2.68638;
		SDS = -2.64779;
		PDP = 2.10953;
		S_DS = -1.12312;
		SSS = -1.39456;
		S_SS = -2.0183;
		SS_S = -2.0183;
		S_S_S = -3.5668;
		PPS = 4.28921;
		PPP = -1.73707;
		DDS = -1.32941;
		DDP = 2.56261;
		DDD = -1.9512;
	}
	if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == true)) {
		//SI-SI
		SPS = 3.02562;
		PDS = -1.35554;
		S_PS = 3.15565;
		SDS = -2.28485;
		PDP = 2.38479;
		S_DS = -0.80993;
		SSS = -1.95933;
		S_SS = -1.5223;
		SS_S = -1.5223;
		S_S_S = -4.24135;
		PPS = 4.10364;
		PPP = -1.51801;
		DDS = -1.68136;
		DDP = 2.5888;
		DDD = -1.814;
	}
	if (ifSIcurrentNeighbour != ifSIcurrentPsi) {
		//GE-SI
		SPS   = 2.878485;	// (2.73135+3.02562)/2.0;
		PDS   = -1.678345;	// (-2.00115-1.35554)/2.0;
		S_PS  = 2.921015;	// (2.68638+3.15565)/2.0;
		SDS   = -2.46632;	// (-2.64779-2.28485)/2.0;
		PDP   = 2.24716;	// (2.10953+2.38479)/2.0;
		S_DS  = -0.966525;	// (-1.12312-0.80993)/2.0;
		SSS   = -1.676945;	// (-1.39456-1.95933)/2.0;
		S_SS  = -1.7703;	// (-2.0183-1.5223)/2.0;
		SS_S  = -1.7703;	// (-2.0183-1.5223)/2.0;
		S_S_S = -3.904075;	// (-3.5668-4.24135)/2.0;
		PPS   = 4.196425;	// (4.28921+4.10364)/2.0;
		PPP   = -1.62754;	// (-1.73707-1.51801)/2.0;
		DDS   = -1.505385;	// (-1.32941-1.68136)/2.0;
		DDP   = 2.575705;	// (2.56261+2.5888)/2.0;
		DDD   = -1.8826;	// (-1.9512-1.814)/2.0;
	}

	complex<double> boundaryCoeff(1.0, 0.0);
	V<double> v12;
	v12 = v2 - v1;
	if (fabs(v12.x) > volumeSize/2.0) {
		v12.x -= volumeSize*v12.x/fabs(v12.x);
		//boundaryCoeff = li;
	}
	if (fabs(v12.y) > volumeSize/2.0)
		v12.y -= volumeSize*v12.y/fabs(v12.y);
	if (fabs(v12.z) > volumeSize/2.0)
		v12.z -= volumeSize*v12.z/fabs(v12.z);
	v12 = (1.0/sqrt(v12*v12))*v12;
	double l = v12.x;
	double m = v12.y;
	double n = v12.z;

	switch(alpha) {
		case 0: switch(beta) {
					case 0: return   boundaryCoeff * SSS;
					case 1: return - boundaryCoeff * l * SPS;
					case 2: return - boundaryCoeff * m * SPS;
					case 3: return - boundaryCoeff * n * SPS;
					case 4: return   boundaryCoeff * SS_S;
					case 5: return   boundaryCoeff * sqrt(3.0)*l*m*SDS;			//xy
					case 6: return   boundaryCoeff * sqrt(3.0)*m*n*SDS;			//yz
					case 7: return   boundaryCoeff * sqrt(3.0)*n*l*SDS;			//zx
					case 8: return   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
					case 9: return   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2
				}
		case 1: switch(beta) {
					case 0: return   boundaryCoeff * l*SPS;
					case 1: return   boundaryCoeff * (l*l*(PPS-PPP) + PPP);
					case 2: return   boundaryCoeff * l*m*(PPS-PPP);
					case 3: return   boundaryCoeff * l*n*(PPS-PPP);
					case 4: return   boundaryCoeff * l*S_PS;
					case 5: return - boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
					case 6: return - boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
					case 7: return - boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
					case 8: return - boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
					case 9: return - boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				}
		case 2: switch(beta) {
					case 0: return   boundaryCoeff * m*SPS;
					case 1: return   boundaryCoeff * m*l*(PPS-PPP);
					case 2: return   boundaryCoeff * (m*m*(PPS-PPP) + PPP);
					case 3: return   boundaryCoeff * m*n*(PPS-PPP);
					case 4: return   boundaryCoeff * m*S_PS;
					case 5: return - boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
					case 6: return - boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
					case 7: return - boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
					case 8: return - boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
					case 9: return - boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				}
		case 3: switch(beta) {
					case 0: return   boundaryCoeff * n*SPS;
					case 1: return   boundaryCoeff * n*l*(PPS-PPP);
					case 2: return   boundaryCoeff * n*m*(PPS-PPP);
					case 3: return   boundaryCoeff * (n*n*(PPS-PPP) + PPP);
					case 4: return   boundaryCoeff * n*S_PS;
					case 5: return - boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
					case 6: return - boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
					case 7: return - boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
					case 8: return - boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
					case 9: return - boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				}
		case 4: switch(beta) {
					case 0: return   boundaryCoeff * SS_S;
					case 1: return - boundaryCoeff * l * S_PS;
					case 2: return - boundaryCoeff * m * S_PS;
					case 3: return - boundaryCoeff * n * S_PS;
					case 4: return   boundaryCoeff * S_S_S;
					case 5: return   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
					case 6: return   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
					case 7: return   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
					case 8: return   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
					case 9: return   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
				}
		case 5: switch(beta) {
					case 0: return   boundaryCoeff * sqrt(3.0)*l*m*SDS;
					case 1: return   boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
					case 2: return   boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
					case 3: return   boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
					case 4: return   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
					case 5: return   boundaryCoeff * (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
					case 6: return   boundaryCoeff * (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
					case 7: return   boundaryCoeff * (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
					case 8: return   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
					case 9: return   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				}
		case 6: switch(beta) {
					case 0: return   boundaryCoeff * sqrt(3.0)*m*n*SDS;
					case 1: return   boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
					case 2: return   boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
					case 3: return   boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
					case 4: return   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
					case 5: return   boundaryCoeff * (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
					case 6: return   boundaryCoeff * (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
					case 7: return   boundaryCoeff * (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
					case 8: return   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
					case 9: return   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				}
		case 7: switch(beta) {
					case 0: return   boundaryCoeff * sqrt(3.0)*n*l*SDS;
					case 1: return   boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
					case 2: return   boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
					case 3: return   boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
					case 4: return   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
					case 5: return   boundaryCoeff * (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
					case 6: return   boundaryCoeff * (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
					case 7: return   boundaryCoeff * (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
					case 8: return   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
					case 9: return   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				}
		case 8: switch(beta) {
					case 0: return   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;
					case 1: return   boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
					case 2: return   boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
					case 3: return   boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
					case 4: return   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
					case 5: return   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
					case 6: return   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
					case 7: return   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
					case 8: return   boundaryCoeff * (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
					case 9: return   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				}
		case 9: switch(beta) {
					case 0: return   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;
					case 1: return   boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
					case 2: return   boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
					case 3: return   boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
					case 4: return   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
					case 5: return   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
					case 6: return   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
					case 7: return   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
					case 8: return   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
					case 9: return   boundaryCoeff * ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				}
	}

}
*/
complex<double> operator *(double f, complex<double> c) {
	complex<double> tmp(f*c.real(), f*c.imag());
	return tmp;
}


int main(int argc, char * argv[]) {

	Config cfg;
	try {
		cfg.readFile("init.cfg");
	}
	catch(const FileIOException &fioex) {
		std::cerr << "I/O error while reading file." << std::endl;
		return(EXIT_FAILURE);
	}
	catch(const ParseException &pex) {
		std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine() << " - " << pex.getError() << std::endl;
		return(EXIT_FAILURE);
	}

	if (argc==4) {
		k.x = stof(argv[1], NULL);
		k.y = stof(argv[2], NULL);
		k.z = stof(argv[3], NULL);
	}
	else if (argc!=1) { cout << "wrong argc, argv" << endl; exit(1); }


#ifdef BANDSTRUCTURE
	for(k.x=-3.3; k.x<=6.3; k.x+=0.05) {
//	for(k.x=4.76-0.1; k.x<=4.76+0.1; k.x+=0.2) { /// delta x-mass
		k.y = k.x;
		k.z = k.x;	
		if(k.x > 0) {
			k.y = 0.0;
			k.z = 0.0;
		}


/*       ///// delta y-mass
	k.x=4.76;
	for(k.y=-1.0; k.y<=1.0; k.y+=0.2) {
*/
#endif


#ifdef PRINT_H
	try {
		const char* print_h_filename = cfg.lookup("print_h_filename");
		print_h_out = fopen(print_h_filename, "w");  if (print_h_out==NULL) fprintf(stderr,"Error!\ncan't open the file %s for writing.\n",print_h_filename),exit(0);
	}
	catch(const SettingNotFoundException &nfex) {
		cerr << "No print_h_filename" << endl;
	}
#endif

//////////////////////////  INIT            /////////////////////////////////////////////////
		//{
		double* E[10];
		E[0] = &Es;
		E[1] = &Ep;
		E[2] = &Ep;
		E[3] = &Ep;
		E[4] = &Es_;
		E[5] = &Ed;
		E[6] = &Ed;
		E[7] = &Ed;
		E[8] = &Ed;
		E[9] = &Ed;

		double* strain_constant_K[10][10];
		for (int i=0; i<10; i++)
			for (int j=0; j<10; j++) {
				strain_constant_K[i][j] = NULL;
			}
		strain_constant_K[0][0] = &strain_constant_K_S_S;
		strain_constant_K[0][1] = &strain_constant_K_S_P;
		strain_constant_K[0][2] = &strain_constant_K_S_P;
		strain_constant_K[0][3] = &strain_constant_K_S_P;
		strain_constant_K[0][4] = &strain_constant_K_Sstar_S;
		strain_constant_K[0][5] = &strain_constant_K_S_D;
		strain_constant_K[0][6] = &strain_constant_K_S_D;
		strain_constant_K[0][7] = &strain_constant_K_S_D;
		strain_constant_K[0][8] = &strain_constant_K_S_D;
		strain_constant_K[0][9] = &strain_constant_K_S_D;

		strain_constant_K[1][0] = &strain_constant_K_S_P;
		strain_constant_K[1][1] = &strain_constant_K_P_P;
		strain_constant_K[1][2] = &strain_constant_K_P_P;
		strain_constant_K[1][3] = &strain_constant_K_P_P;
		strain_constant_K[1][4] = &strain_constant_K_Sstar_P;
		strain_constant_K[1][5] = &strain_constant_K_P_D;
		strain_constant_K[1][6] = &strain_constant_K_P_D;
		strain_constant_K[1][7] = &strain_constant_K_P_D;
		strain_constant_K[1][8] = &strain_constant_K_P_D;
		strain_constant_K[1][9] = &strain_constant_K_P_D;

		strain_constant_K[2][0] = &strain_constant_K_S_P;
		strain_constant_K[2][1] = &strain_constant_K_P_P;
		strain_constant_K[2][2] = &strain_constant_K_P_P;
		strain_constant_K[2][3] = &strain_constant_K_P_P;
		strain_constant_K[2][4] = &strain_constant_K_Sstar_P;
		strain_constant_K[2][5] = &strain_constant_K_P_D;
		strain_constant_K[2][6] = &strain_constant_K_P_D;
		strain_constant_K[2][7] = &strain_constant_K_P_D;
		strain_constant_K[2][8] = &strain_constant_K_P_D;
		strain_constant_K[2][9] = &strain_constant_K_P_D;

		strain_constant_K[3][0] = &strain_constant_K_S_P;
		strain_constant_K[3][1] = &strain_constant_K_P_P;
		strain_constant_K[3][2] = &strain_constant_K_P_P;
		strain_constant_K[3][3] = &strain_constant_K_P_P;
		strain_constant_K[3][4] = &strain_constant_K_Sstar_P;
		strain_constant_K[3][5] = &strain_constant_K_P_D;
		strain_constant_K[3][6] = &strain_constant_K_P_D;
		strain_constant_K[3][7] = &strain_constant_K_P_D;
		strain_constant_K[3][8] = &strain_constant_K_P_D;
		strain_constant_K[3][9] = &strain_constant_K_P_D;

		strain_constant_K[4][0] = &strain_constant_K_S_Sstar;
		strain_constant_K[4][1] = &strain_constant_K_Sstar_P;
		strain_constant_K[4][2] = &strain_constant_K_Sstar_P;
		strain_constant_K[4][3] = &strain_constant_K_Sstar_P;
		strain_constant_K[4][4] = &strain_constant_K_Sstar_Sstar;
		strain_constant_K[4][5] = &strain_constant_K_Sstar_D;
		strain_constant_K[4][6] = &strain_constant_K_Sstar_D;
		strain_constant_K[4][7] = &strain_constant_K_Sstar_D;
		strain_constant_K[4][8] = &strain_constant_K_Sstar_D;
		strain_constant_K[4][9] = &strain_constant_K_Sstar_D;



		strain_constant_K[5][0] = &strain_constant_K_S_D;
		strain_constant_K[5][1] = &strain_constant_K_P_D;
		strain_constant_K[5][2] = &strain_constant_K_P_D;
		strain_constant_K[5][3] = &strain_constant_K_P_D;
		strain_constant_K[5][4] = &strain_constant_K_Sstar_D;
		strain_constant_K[5][5] = &strain_constant_K_D_D;
		strain_constant_K[5][6] = &strain_constant_K_D_D;
		strain_constant_K[5][7] = &strain_constant_K_D_D;
		strain_constant_K[5][8] = &strain_constant_K_D_D;
		strain_constant_K[5][9] = &strain_constant_K_D_D;

		strain_constant_K[6][0] = &strain_constant_K_S_D;
		strain_constant_K[6][1] = &strain_constant_K_P_D;
		strain_constant_K[6][2] = &strain_constant_K_P_D;
		strain_constant_K[6][3] = &strain_constant_K_P_D;
		strain_constant_K[6][4] = &strain_constant_K_Sstar_D;
		strain_constant_K[6][5] = &strain_constant_K_D_D;
		strain_constant_K[6][6] = &strain_constant_K_D_D;
		strain_constant_K[6][7] = &strain_constant_K_D_D;
		strain_constant_K[6][8] = &strain_constant_K_D_D;
		strain_constant_K[6][9] = &strain_constant_K_D_D;

		strain_constant_K[7][0] = &strain_constant_K_S_D;
		strain_constant_K[7][1] = &strain_constant_K_P_D;
		strain_constant_K[7][2] = &strain_constant_K_P_D;
		strain_constant_K[7][3] = &strain_constant_K_P_D;
		strain_constant_K[7][4] = &strain_constant_K_Sstar_D;
		strain_constant_K[7][5] = &strain_constant_K_D_D;
		strain_constant_K[7][6] = &strain_constant_K_D_D;
		strain_constant_K[7][7] = &strain_constant_K_D_D;
		strain_constant_K[7][8] = &strain_constant_K_D_D;
		strain_constant_K[7][9] = &strain_constant_K_D_D;

		strain_constant_K[8][0] = &strain_constant_K_S_D;
		strain_constant_K[8][1] = &strain_constant_K_P_D;
		strain_constant_K[8][2] = &strain_constant_K_P_D;
		strain_constant_K[8][3] = &strain_constant_K_P_D;
		strain_constant_K[8][4] = &strain_constant_K_Sstar_D;
		strain_constant_K[8][5] = &strain_constant_K_D_D;
		strain_constant_K[8][6] = &strain_constant_K_D_D;
		strain_constant_K[8][7] = &strain_constant_K_D_D;
		strain_constant_K[8][8] = &strain_constant_K_D_D;
		strain_constant_K[8][9] = &strain_constant_K_D_D;

		strain_constant_K[9][0] = &strain_constant_K_S_D;
		strain_constant_K[9][1] = &strain_constant_K_P_D;
		strain_constant_K[9][2] = &strain_constant_K_P_D;
		strain_constant_K[9][3] = &strain_constant_K_P_D;
		strain_constant_K[9][4] = &strain_constant_K_Sstar_D;
		strain_constant_K[9][5] = &strain_constant_K_D_D;
		strain_constant_K[9][6] = &strain_constant_K_D_D;
		strain_constant_K[9][7] = &strain_constant_K_D_D;
		strain_constant_K[9][8] = &strain_constant_K_D_D;
		strain_constant_K[9][9] = &strain_constant_K_D_D;

		//}
//////////////////////////  ATOMS INIT      /////////////////////////////////////////////////
		//{
		List<Atom> *current = NULL;

		WaveFunction psi;

		#ifdef UI
			cout << "initializing " << N << " atoms" << endl;
		#endif
		
		double currInit = 0;
		for(double x=0; x<=volumeSize-a/2.0; x+=a) {
			for(double y=0; y<=volumeSize-a/2.0; y+=a) {
				for(double z=0; z<=volumeSize-a/2.0; z+=a) {
					Atom* newAtom = new Atom();
					newAtom->x = x;
					newAtom->y = y;
					newAtom->z = z;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a/4.0;
					newAtom->y = y + a/4.0;
					newAtom->z = z + a/4.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a/2.0;
					newAtom->y = y + a/2.0;
					newAtom->z = z;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a*3.0/4.0;
					newAtom->y = y + a*3.0/4.0;
					newAtom->z = z + a/4.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a/2.0;
					newAtom->y = y;
					newAtom->z = z + a/2.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a*3.0/4.0;
					newAtom->y = y + a/4.0;
					newAtom->z = z + a*3.0/4.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x;
					newAtom->y = y + a/2.0;
					newAtom->z = z + a/2.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					newAtom = new Atom();
					newAtom->x = x + a/4.0;
					newAtom->y = y + a*3.0/4.0;
					newAtom->z = z + a*3.0/4.0;
					Add_tail(newAtom, &(psi.atoms), &(psi.atoms_tail));
					newAtom = NULL;

					currInit += 8;
				}
			}
		//cout << currInit << endl;
		}
		//}
//////////////////////////  DEFORMATION     /////////////////////////////////////////////////
		//{
		current = psi.atoms;
		while (current) {
			current->d->dx = DEFORMX*current->d->x;
			current->d->dy = DEFORMY*current->d->y;
			current->d->dz = DEFORMZ*current->d->z;
			current = current->next;
		}


/*		char* strain_filename = "./../QD_strain.xyz";
		in = fopen(strain_filename, "r");
		if (in==NULL) fprintf(stderr, "error: can't open the file %s for reading.\n", strain_filename), exit(1);

		current = psi.atoms;
		while (current) {
			double strained_x, strained_y, strained_z;
			double d, dx, dy, dz;
			int sostav;
			double Exx, Eyy, Ezz, Exy, Exz, Eyz;
			char buf[bufsize];

			alex_fgets(buf, bufsize, in);
			int scanerr = sscanf(buf, "%i %lf %lf %lf %lf %lf %lf %lf %lf %lf", &sostav, &strained_x, &strained_y, &strained_z, &d, &dx, &dy, &dz, &Exx, &Exy, &Exz, &Eyy, &Eyz, &Ezz);
			if (scanerr == char_traits<char>::eof())
			if ((fabs(current->d->x - (strained_x-dx)) > 1e-6) || (fabs(current->d->y - (strained_y-dy)) > 1e-6) || (fabs(current->d->z - (strained_z-dz)) > 1e-6)) {
				fprintf(stderr, "error: the position of sigma(x,y,z) does not equal the position of atom(x,y,z)\n");
				exit(1);
			}

			current->d->dx = dx;
			current->d->dy = dy;
			current->d->dz = dz;


			current->d->sostav = sostav;
			
			current = current->next;
		}

		fclose(in);

*/
		//}
//////////////////////////  NEIGHBOURS      /////////////////////////////////////////////////
		//{
		#ifdef UI
			cout << "linking atom neighbours" << endl;
		#endif
		current = psi.atoms;
		int AtomCounter = 1;
		int neighboursProgress = N/10;
		double dist  = (11.0+4.0*sqrt(6.0))*a*a/8.0/8.0;
		double dist2 = (2.0*sqrt(2.0)+sqrt(11.0))*(2.0*sqrt(2.0)+sqrt(11.0))*a*a/8.0/8.0;
		while (current) {
			current->d->n = AtomCounter;
			AtomCounter++;
			if (AtomCounter > neighboursProgress) {
				#ifdef UI
				cout << "." << endl;
				#endif
				neighboursProgress += N/10;
			}
			int i = 0;
			List<Atom> *current2 = psi.atoms;
			while (current2) {
				double xx = current->d->x - current2->d->x;
				if (fabs(xx) > volumeSize/2.0) xx -= volumeSize*xx/fabs(xx);
				double yy = current->d->y - current2->d->y;
				if (fabs(yy) > volumeSize/2.0) yy -= volumeSize*yy/fabs(yy);
				double zz = current->d->z - current2->d->z;
				if (fabs(zz) > volumeSize/2.0) zz -= volumeSize*zz/fabs(zz);
				double rSquared = xx*xx + yy*yy + zz*zz;
				if ((current != current2) && (rSquared < dist)) {
					Add(current2->d, &(current->d->neighbours));
					i++;
				}
//				if ((dist < rSquared) && (rSquared < dist2)) {
//					Add(current2->d, &(current->d->second_neighbours));
//				}
				current2 = current2->next;
			}
			if (i!=4)
				fprintf(stderr, "error: atom %f has %f neighbours.\n", current->d->n, i), exit(1);
			current = current->next;
		}
		#ifdef UI
		cout << endl;
		#endif
		//}
//////////////////////////  CSR init        /////////////////////////////////////////////////
		//{
		#ifdef UI
		cout << "initializing CSR" << endl;
		#endif
#ifdef SO
	#ifdef UCSR
		//     onsite: diag  nondiag           first neighbour    SO
		//long CSRsize = 20*N + 45*N             + 4*(100)*N   +   6*N;     //UCSR
		long CSRsize = 20*N                      + (400-19)*N   +   6*N;
	#else
		#ifdef BOYKIN2002
			long CSRsize = 2*100*N                   + 4*(100)*2*N   +   6*2*N; //CSR
		#else
			long CSRsize = 20*N                      + 4*(100)*2*N   +   6*2*N; //CSR
		#endif
	#endif
#else
		long CSRsize = 20*N + 45*N             + 4*(100)*N;		//UCSR
#endif

		#ifdef UI
		cout << "init " << CSRsize << " nonzero elements" << endl;
		#endif
		long currentCSRval = 0;

		complex<double>		* CSRval = (complex<double>*)malloc(CSRsize*sizeof(complex<double>));
		long				* CSRrows = (long*)malloc((20*N+1)*sizeof(long));
		long				* CSRcols = (long*)malloc(CSRsize*sizeof(long));
		if (!((CSRval) && (CSRrows) && (CSRcols)))
			cerr << "error: not enough heap memory" << endl;



		/////////////////////////////////////////////////////////////////////////////////////
		complex<double> *Htmp = (complex<double>*)malloc((10 * 20*N)*sizeof(complex<double>));
		//}
//////////////////////////  spin plus       /////////////////////////////////////////////////
		//{
		#ifdef UI
		cout << "spin plus" << endl;
		#endif
		int csrProgress = N/10;
		List<Atom> *currentPsi = psi.atoms;
		while (currentPsi) {
			if (currentPsi->d->n == csrProgress) {
				#ifdef UI
				cout << "." << endl;
				#endif
				csrProgress += N/10;
			}
			for(long i=0; i<10 * 20*N; i++)
				Htmp[i]=0.0;
			V<double> v1(currentPsi->d->x + currentPsi->d->dx, currentPsi->d->y + currentPsi->d->dy, currentPsi->d->z + currentPsi->d->dz);




			bool ifSIcurrentPsi = currentPsi->d->sostav; //ifSI(v1);
			if (ifSIcurrentPsi == true) {
				a_lattice = 0.54310;
				bond_length = sqrt(3.0)/4.0 * a_lattice;
				Es = -2.15168;
				Ep = 4.22925;
				Es_ = 19.1165;
				Ed = 13.7895;
				delta = 0.03978;
				Energy_shift = 27;
			}
			else {
				a_lattice = 0.56579060;
				bond_length = sqrt(3.0)/4.0 * a_lattice;
				Es = -1.95617;
				Ep = 5.3097;
				Es_ = 19.296;
				Ed = 13.5806;
				delta = 0.20264;
				Energy_shift = 27.77;
			}
			

			Htmp[0*20*N + (currentPsi->d->n-1)*10    ] = Es;
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1] = Ep;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2] = Ep;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3] = Ep;
			Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4] = Es_;
			Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5] = Ed;
			Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6] = Ed;
			Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7] = Ed;
			Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8] = Ed;
			Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9] = Ed;

#ifdef SO
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2       ] = - delta/3.0*li;	// px+ py+
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1       ] =   delta/3.0*li;	// py+ px+
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] =   delta/3.0;	// px+ pz-
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] = - delta/3.0*li;	// py+ pz-
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] = - delta/3.0;	// pz+ px-
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] =   delta/3.0*li; // pz+ py-
#endif


			List<Atom> *currentNeighbour = currentPsi->d->neighbours;
			while (currentNeighbour) {
				V<double> v2;
				v2.x = currentNeighbour->d->x + currentNeighbour->d->dx;
				v2.y = currentNeighbour->d->y + currentNeighbour->d->dy;
				v2.z = currentNeighbour->d->z + currentNeighbour->d->dz;

				bool ifSIcurrentNeighbour = currentNeighbour->d->sostav;  //ifSI(v2);
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == false)) {
					//GE-GE
					SPS = 2.73135;
					PDS = -2.00115;
					S_PS = 2.68638;
					SDS = -2.64779;
					PDP = 2.10953;
					S_DS = -1.12312;
					SSS = -1.39456;
					S_SS = -2.0183;
					SS_S = -2.0183;
					S_S_S = -3.5668;
					PPS = 4.28921;
					PPP = -1.73707;
					DDS = -1.32941;
					DDP = 2.56261;
					DDD = -1.9512;
					strain_exponent_V_S_S_Sigma = 1.99551;
					strain_exponent_V_S_P_Sigma = 1.29303;
					strain_exponent_V_Sstar_P_Sigma = 5.0000;
					strain_exponent_V_S_Sstar_Sigma = 0.0000;
					strain_exponent_V_Sstar_S_Sigma = 0.0000;
					strain_exponent_V_Sstar_Sstar_Sigma = 2.38823;
					strain_exponent_V_Sstar_D_Sigma = 0.75134;
					strain_exponent_V_S_D_Sigma = 2.79244;
					strain_exponent_V_P_P_Sigma = 1.13641;
					strain_exponent_V_P_P_Pi = 1.74803;
					strain_exponent_V_P_D_Sigma = 2.68784;
					strain_exponent_V_P_D_Pi = 4.36921;
					strain_exponent_V_D_D_Sigma = 5.00000;
					strain_exponent_V_D_D_Pi = 0.69769;
					strain_exponent_V_D_D_Delta = 3.06253;

					strain_constant_K_S_S = 0.0000;
					strain_constant_K_S_P = 1.2507;
					strain_constant_K_Sstar_P = 2.6841;
					strain_constant_K_Sstar_S = 1.1766;
					strain_constant_K_S_Sstar = 1.1766;
					strain_constant_K_Sstar_D = 1.2276;
					strain_constant_K_Sstar_Sstar = 2.6841;
					strain_constant_K_S_D = 0.1524;
					strain_constant_K_P_D = 0.1143;
					strain_constant_K_P_P = 0.3626;
					strain_constant_K_D_D = 1.9527;
					Z_eff = 0.0;
					r_P_D = 0.0;
					r2_P_P = 0.0;
					r2_D_D = 0.0;
					r3_P_D = 0.0;
					r4_D_D = 0.0;
				}
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == true)) {
					//SI-SI
					SPS = 3.02562;
					PDS = -1.35554;
					S_PS = 3.15565;
					SDS = -2.28485;
					PDP = 2.38479;
					S_DS = -0.80993;
					SSS = -1.95933;
					S_SS = -1.5223;
					SS_S = -1.5223;
					S_S_S = -4.24135;
					PPS = 4.10364;
					PPP = -1.51801;
					DDS = -1.68136;
					DDP = 2.5888;
					DDD = -1.814;
					strain_exponent_V_S_S_Sigma	= 0.56247;
					strain_exponent_V_S_P_Sigma = 2.36548;
					strain_exponent_V_Sstar_P_Sigma = 0.34492;
					strain_exponent_V_S_Sstar_Sigma = 0.13203;
					strain_exponent_V_Sstar_S_Sigma = 0.13203;
					strain_exponent_V_Sstar_Sstar_Sigma = 0.19237;
					strain_exponent_V_Sstar_D_Sigma = 1.08601;
					strain_exponent_V_S_D_Sigma = 2.56720;
					strain_exponent_V_P_P_Sigma = 0.2;
					strain_exponent_V_P_P_Pi = 1.6777;
					strain_exponent_V_P_D_Sigma = 0.2;
					strain_exponent_V_P_D_Pi = 4.4325;
					strain_exponent_V_D_D_Sigma = 0.1;
					strain_exponent_V_D_D_Pi = 6;
					strain_exponent_V_D_D_Delta = 5.9997;

					strain_constant_K_S_S = 1.0875;
					strain_constant_K_S_P = 0.4001;
					strain_constant_K_Sstar_P = 1.8597;
					strain_constant_K_Sstar_S = 1.1374;
					strain_constant_K_S_Sstar = 1.1374;
					strain_constant_K_Sstar_D = 0.2978;
					strain_constant_K_Sstar_Sstar = 0.5986;
					strain_constant_K_S_D = 0;
					strain_constant_K_P_D = 0.0564;
					strain_constant_K_P_P = 0;
					strain_constant_K_D_D = 2.1071;
					Z_eff = 3.0;
					r_P_D = 0.1249;
					r2_P_P = 0.1227;
					r2_D_D = 0.5147;
					r3_P_D = 0.0;
					r4_D_D = 0.7661;
				}
				if (ifSIcurrentNeighbour != ifSIcurrentPsi) {
					//GE-SI
					SPS   = 2.878485;	// (2.73135+3.02562)/2.0;
					PDS   = -1.678345;	// (-2.00115-1.35554)/2.0;
					S_PS  = 2.921015;	// (2.68638+3.15565)/2.0;
					SDS   = -2.46632;	// (-2.64779-2.28485)/2.0;
					PDP   = 2.24716;	// (2.10953+2.38479)/2.0;
					S_DS  = -0.966525;	// (-1.12312-0.80993)/2.0;
					SSS   = -1.676945;	// (-1.39456-1.95933)/2.0;
					S_SS  = -1.7703;	// (-2.0183-1.5223)/2.0;
					SS_S  = -1.7703;	// (-2.0183-1.5223)/2.0;
					S_S_S = -3.904075;	// (-3.5668-4.24135)/2.0;
					PPS   = 4.196425;	// (4.28921+4.10364)/2.0;
					PPP   = -1.62754;	// (-1.73707-1.51801)/2.0;
					DDS   = -1.505385;	// (-1.32941-1.68136)/2.0;
					DDP   = 2.575705;	// (2.56261+2.5888)/2.0;
					DDD   = -1.8826;	// (-1.9512-1.814)/2.0;

					strain_exponent_V_S_S_Sigma = 1.27899;			//(0.56247+1.99551)/2.0;
					strain_exponent_V_S_P_Sigma = 1.829255;			//(1.29303+2.36548)/2.0;
					strain_exponent_V_Sstar_P_Sigma = 2.67246;		//(5.0000+0.34492)/2.0;
					strain_exponent_V_S_Sstar_Sigma = 0.13203;		//(0.0000+0.13203)/2.0;
					strain_exponent_V_Sstar_S_Sigma = 0.13203;		//(0.0000+0.13203)/2.0;
					strain_exponent_V_Sstar_Sstar_Sigma = 1.2903;	//(2.38823+0.19237)/2.0;
					strain_exponent_V_Sstar_D_Sigma = 0.918675;		//(0.75134+1.08601)/2.0;
					strain_exponent_V_S_D_Sigma = 2.67982;			//(2.79244+2.56720)/2.0;
					strain_exponent_V_P_P_Sigma = 0.668205;			//(1.13641+0.2)/2.0;
					strain_exponent_V_P_P_Pi = 1.712865;			//(1.74803+1.6777)/2.0;
					strain_exponent_V_P_D_Sigma = 1.44392;			//(2.68784+0.2)/2.0;
					strain_exponent_V_P_D_Pi = 4.400855;			//(4.36921+4.4325)/2.0;
					strain_exponent_V_D_D_Sigma = 2.55;				//(5.00000+0.1)/2.0;
					strain_exponent_V_D_D_Pi = 3.348845;			//(0.69769+6.0)/2.0;
					strain_exponent_V_D_D_Delta = 4.531115;			//(3.06253+5.9997)/2.0;

					strain_constant_K_S_S = 1.0875;				//(0.0+1.0875)/2.0;
					strain_constant_K_S_P = 0.8254;				//(1.2507+0.4001)/2.0;
					strain_constant_K_Sstar_P = 2.2719;			//(2.6841+1.8597)/2.0;
					strain_constant_K_Sstar_S = 1.157;			//(1.1766+1.1374)/2.0;
					strain_constant_K_S_Sstar = 1.157;			//(1.1766+1.1374)/2.0;
					strain_constant_K_Sstar_D = 0.7627;			//(1.2276+0.2978)/2.0;
					strain_constant_K_Sstar_Sstar = 1.64135;	//(2.6841+0.5986)/2.0;
					strain_constant_K_S_D = 0.1524;				//(0.1524+0.0)/2.0;
					strain_constant_K_P_D = 0.08535;			//(0.1143+0.0564)/2.0;
					strain_constant_K_P_P = 0.3626;				//(0.3626+0.0)/2.0;
					strain_constant_K_D_D = 2.0299;				//(1.9527+2.1071)/2.0;
					Z_eff = 0.0;
					r_P_D = 0.0;
					r2_P_P = 0.0;
					r2_D_D = 0.0;
					r3_P_D = 0.0;
					r4_D_D = 0.0;
				}

				complex<double> boundaryCoeff(1.0, 0.0);
				V<double> v12;
				v12 = v2 - v1;
				double volumeSizeX = volumeSize + DEFORMX*volumeSize;
				double volumeSizeY = volumeSize + DEFORMY*volumeSize;
				double volumeSizeZ = volumeSize + DEFORMZ*volumeSize;

				#ifdef BANDSTRUCTURE
					double kL = 0;
					if (fabs(v12.x) > volumeSizeX/2.0) {
						v12.x -= volumeSizeX*v12.x/fabs(v12.x);
						kL += k.x*volumeSizeX*(v12.x/fabs(v12.x));
					}
					if (fabs(v12.y) > volumeSizeY/2.0) {
						v12.y -= volumeSizeY*v12.y/fabs(v12.y);
						kL += k.y*volumeSizeY*(v12.y/fabs(v12.y));
					}
					if (fabs(v12.z) > volumeSizeZ/2.0) {
						v12.z -= volumeSizeZ*v12.z/fabs(v12.z);
						kL += k.z*volumeSizeZ*(v12.z/fabs(v12.z));
					}
					complex<double> ikL(0,kL);
					boundaryCoeff = exp(ikL);
				#else
					if (fabs(v12.x) > volumeSizeX/2.0) {
						v12.x -= volumeSizeX*v12.x/fabs(v12.x);
						boundaryCoeff = BOUNDARYCOEFF;
					}
					if (fabs(v12.y) > volumeSizeY/2.0)
						v12.y -= volumeSizeY*v12.y/fabs(v12.y);
					if (fabs(v12.z) > volumeSizeZ/2.0)
						v12.z -= volumeSizeZ*v12.z/fabs(v12.z);
				#endif
				V<double> v12norm = (1.0/sqrt(v12*v12))*v12;
				double l = v12norm.x;
				double m = v12norm.y;
				double n = v12norm.z;
#ifdef BOYKIN2002
				V<double> v1unstrained;
				v1unstrained.x = currentPsi->d->x;
				v1unstrained.y = currentPsi->d->y;
				v1unstrained.z = currentPsi->d->z;
				V<double> v2unstrained;
				v2unstrained.x = currentNeighbour->d->x;
				v2unstrained.y = currentNeighbour->d->y;
				v2unstrained.z = currentNeighbour->d->z;
				V<double> v12unstrained;
				v12unstrained = v2unstrained - v1unstrained;

				// calculate d0
				if (fabs(v12unstrained.x) > volumeSize/2.0)
					v12unstrained.x -= volumeSize*v12unstrained.x/fabs(v12unstrained.x);
				if (fabs(v12unstrained.y) > volumeSize/2.0)
					v12unstrained.y -= volumeSize*v12unstrained.y/fabs(v12unstrained.y);
				if (fabs(v12unstrained.z) > volumeSize/2.0)
					v12unstrained.z -= volumeSize*v12unstrained.z/fabs(v12unstrained.z);
				V<double> v12unstrainedNorm = (1.0/sqrt(v12unstrained*v12unstrained))*v12unstrained;

				double d  = sqrt(v12*v12);
				double d0 = sqrt(v12unstrained*v12unstrained);


				SPS		*= pow(d0/d, strain_exponent_V_S_P_Sigma);
				PDS		*= pow(d0/d, strain_exponent_V_P_D_Sigma);
				S_PS	*= pow(d0/d, strain_exponent_V_Sstar_P_Sigma);
				SDS		*= pow(d0/d, strain_exponent_V_S_D_Sigma);
				PDP		*= pow(d0/d, strain_exponent_V_P_D_Pi);
				S_DS	*= pow(d0/d, strain_exponent_V_Sstar_D_Sigma);
				SSS		*= pow(d0/d, strain_exponent_V_S_S_Sigma);
				S_SS	*= pow(d0/d, strain_exponent_V_Sstar_S_Sigma);
				SS_S	*= pow(d0/d, strain_exponent_V_S_Sstar_Sigma);
				S_S_S	*= pow(d0/d, strain_exponent_V_Sstar_Sstar_Sigma);
				PPS		*= pow(d0/d, strain_exponent_V_P_P_Sigma);
				PPP		*= pow(d0/d, strain_exponent_V_P_P_Pi);
				DDS		*= pow(d0/d, strain_exponent_V_D_D_Sigma);
				DDP		*= pow(d0/d, strain_exponent_V_D_D_Pi);
				DDD		*= pow(d0/d, strain_exponent_V_D_D_Delta);
#endif

				//{
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * SSS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 1] = - boundaryCoeff * l * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 2] = - boundaryCoeff * m * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 3] = - boundaryCoeff * n * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * SS_S;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;			//xy
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;			//yz
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;			//zx
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				Htmp[1*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * l*SPS;
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * l*m*(PPS-PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * l*n*(PPS-PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * l*S_PS;
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 5] = - boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 6] = - boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 7] = - boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 8] = - boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 9] = - boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				Htmp[2*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * m*SPS;
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * m*l*(PPS-PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * m*n*(PPS-PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * m*S_PS;
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 5] = - boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 6] = - boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 7] = - boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 8] = - boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 9] = - boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				Htmp[3*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * n*SPS;
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * n*l*(PPS-PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * n*m*(PPS-PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (n*n*(PPS-PPP) + PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * n*S_PS;
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 5] = - boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 6] = - boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 7] = - boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 8] = - boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 9] = - boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				Htmp[4*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * SS_S;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 1] = - boundaryCoeff * l * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 2] = - boundaryCoeff * m * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 3] = - boundaryCoeff * n * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * S_S_S;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;



				Htmp[5*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				Htmp[6*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				Htmp[7*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				Htmp[8*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				Htmp[9*20*N + (currentNeighbour->d->n-1)*10    ] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 1] =   boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 2] =   boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 3] =   boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 4] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 5] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 6] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 7] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 8] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 9] =   boundaryCoeff * ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}




#ifdef BOYKIN2002
				v12unstrained = v2unstrained - v1unstrained;
				#ifdef BANDSTRUCTURE
					kL = 0;
					if (fabs(v12unstrained.x) > volumeSize/2.0) {
						v12unstrained.x -= volumeSize*v12unstrained.x/fabs(v12unstrained.x);
						kL += k.x*volumeSize*(v12unstrained.x/fabs(v12unstrained.x));
					}
					if (fabs(v12unstrained.y) > volumeSize/2.0) {
						v12unstrained.y -= volumeSize*v12unstrained.y/fabs(v12unstrained.y);
						kL += k.y*volumeSize*(v12unstrained.y/fabs(v12unstrained.y));
					}
					if (fabs(v12unstrained.z) > volumeSize/2.0) {
						v12unstrained.z -= volumeSize*v12unstrained.z/fabs(v12unstrained.z);
						kL += k.z*volumeSize*(v12unstrained.z/fabs(v12unstrained.z));
					}
					ikL = li*kL;
					boundaryCoeff = exp(ikL);
				#else
					boundaryCoeff = 1.0;
					if (fabs(v12unstrained.x) > volumeSize/2.0)
						boundaryCoeff = BOUNDARYCOEFF;
				#endif
				// v12unstrainedNorm already calculated (but v12unstrained is incorrect here)
				l = v12unstrainedNorm.x;
				m = v12unstrainedNorm.y;
				n = v12unstrainedNorm.z;

				complex<double> slaterkoster[10][10];
				for (int i=0; i<10; i++)
					for (int j=0; j<10; j++)
						slaterkoster[i][j] = 0.0;


				//{
				slaterkoster[0][0] =   boundaryCoeff * SSS;
				slaterkoster[0][1] = - boundaryCoeff * l * SPS;
				slaterkoster[0][2] = - boundaryCoeff * m * SPS;
				slaterkoster[0][3] = - boundaryCoeff * n * SPS;
				slaterkoster[0][4] =   boundaryCoeff * SS_S;
				slaterkoster[0][5] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;			//xy
				slaterkoster[0][6] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;			//yz
				slaterkoster[0][7] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;			//zx
				slaterkoster[0][8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				slaterkoster[0][9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				slaterkoster[1][0] =   boundaryCoeff * l*SPS;
				slaterkoster[1][1] =   boundaryCoeff * (l*l*(PPS-PPP) + PPP);
				slaterkoster[1][2] =   boundaryCoeff * l*m*(PPS-PPP);
				slaterkoster[1][3] =   boundaryCoeff * l*n*(PPS-PPP);
				slaterkoster[1][4] =   boundaryCoeff * l*S_PS;
				slaterkoster[1][5] = - boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				slaterkoster[1][6] = - boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				slaterkoster[1][7] = - boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				slaterkoster[1][8] = - boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				slaterkoster[1][9] = - boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				slaterkoster[2][0] =   boundaryCoeff * m*SPS;
				slaterkoster[2][1] =   boundaryCoeff * m*l*(PPS-PPP);
				slaterkoster[2][2] =   boundaryCoeff * (m*m*(PPS-PPP) + PPP);
				slaterkoster[2][3] =   boundaryCoeff * m*n*(PPS-PPP);
				slaterkoster[2][4] =   boundaryCoeff * m*S_PS;
				slaterkoster[2][5] = - boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				slaterkoster[2][6] = - boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				slaterkoster[2][7] = - boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				slaterkoster[2][8] = - boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				slaterkoster[2][9] = - boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				slaterkoster[3][0] =   boundaryCoeff * n*SPS;
				slaterkoster[3][1] =   boundaryCoeff * n*l*(PPS-PPP);
				slaterkoster[3][2] =   boundaryCoeff * n*m*(PPS-PPP);
				slaterkoster[3][3] =   boundaryCoeff * (n*n*(PPS-PPP) + PPP);
				slaterkoster[3][4] =   boundaryCoeff * n*S_PS;
				slaterkoster[3][5] = - boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				slaterkoster[3][6] = - boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				slaterkoster[3][7] = - boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				slaterkoster[3][8] = - boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				slaterkoster[3][9] = - boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				slaterkoster[4][0] =   boundaryCoeff * SS_S;
				slaterkoster[4][1] = - boundaryCoeff * l * S_PS;
				slaterkoster[4][2] = - boundaryCoeff * m * S_PS;
				slaterkoster[4][3] = - boundaryCoeff * n * S_PS;
				slaterkoster[4][4] =   boundaryCoeff * S_S_S;
				slaterkoster[4][5] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				slaterkoster[4][6] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				slaterkoster[4][7] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				slaterkoster[4][8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				slaterkoster[4][9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;



				slaterkoster[5][0] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;
				slaterkoster[5][1] =   boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				slaterkoster[5][2] =   boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				slaterkoster[5][3] =   boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				slaterkoster[5][4] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				slaterkoster[5][5] =   boundaryCoeff * (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				slaterkoster[5][6] =   boundaryCoeff * (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				slaterkoster[5][7] =   boundaryCoeff * (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				slaterkoster[5][8] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				slaterkoster[5][9] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				slaterkoster[6][0] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;
				slaterkoster[6][1] =   boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				slaterkoster[6][2] =   boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				slaterkoster[6][3] =   boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				slaterkoster[6][4] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				slaterkoster[6][5] =   boundaryCoeff * (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				slaterkoster[6][6] =   boundaryCoeff * (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				slaterkoster[6][7] =   boundaryCoeff * (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				slaterkoster[6][8] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				slaterkoster[6][9] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				slaterkoster[7][0] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;
				slaterkoster[7][1] =   boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				slaterkoster[7][2] =   boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				slaterkoster[7][3] =   boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				slaterkoster[7][4] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				slaterkoster[7][5] =   boundaryCoeff * (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				slaterkoster[7][6] =   boundaryCoeff * (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				slaterkoster[7][7] =   boundaryCoeff * (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				slaterkoster[7][8] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				slaterkoster[7][9] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				slaterkoster[8][0] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				slaterkoster[8][1] =   boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				slaterkoster[8][2] =   boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				slaterkoster[8][3] =   boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				slaterkoster[8][4] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				slaterkoster[8][5] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				slaterkoster[8][6] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				slaterkoster[8][7] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				slaterkoster[8][8] =   boundaryCoeff * (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				slaterkoster[8][9] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				slaterkoster[9][0] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;
				slaterkoster[9][1] =   boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				slaterkoster[9][2] =   boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				slaterkoster[9][3] =   boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				slaterkoster[9][4] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
				slaterkoster[9][5] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				slaterkoster[9][6] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				slaterkoster[9][7] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				slaterkoster[9][8] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				slaterkoster[9][9] =   boundaryCoeff * ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}

				for (int gamma=0; gamma<10; gamma++)
					for (int alpha=0; alpha<10; alpha++) {
						for (int beta=0; beta<10; beta++) {
							double Kjbia = *strain_constant_K[beta][alpha];
							double Kigjb = *strain_constant_K[gamma][beta];
							double Ejb = *E[beta] - Energy_shift;
							double Eia = *E[alpha] - Energy_shift;
							double Eig = *E[gamma] - Energy_shift;
							
							Htmp[gamma*20*N + (currentPsi->d->n-1)*10 + alpha] += (Htmp[gamma*20*N + (currentNeighbour->d->n-1)*10 + beta]*Htmp[beta*20*N + (currentNeighbour->d->n-1)*10 + alpha] - slaterkoster[gamma][beta]*slaterkoster[beta][alpha]) * (-Kjbia/(Ejb+Eia)/2.0 -Kigjb/(Eig+Ejb)/2.0 -Kigjb*Kjbia*(1.0/(Ejb+Eia)+1.0/(Eig+Ejb))/4.0);
						}
					}
#endif



#ifdef BOYKIN2010
				//complex<double> modified_slaterkoster[10][10];
				//for (int i=0; i<10; i++)
				//	for (int j=0; j<10; j++)
				//		modified_slaterkoster[i][j] = 0.0;

				double coeff = - qe*Zeff/4.0/PI/epsilon0/d; //eV

				double r2pp = r2_P_P*d0*d0/d/d;
				double r_pd = r_P_D*d0/d;
				double r3pd = r3_P_D*d0*d0*d0/d/d/d;
				double r2dd = r2_D_D*d0*d0/d/d;
				double r4dd = r4_D_D*d0*d0*d0*d0/d/d/d/d;

				PPS = coeff            * (2.0/5.0*r2pp);
				PPP = coeff            * (-r2pp/5.0);
				PDS = coeff/sqrt(15.0) * (2.0*r_pd + 9.0/7.0*r3pd);
				PDP = coeff/sqrt(5.0)  * (    r_pd - 3.0/7.0*r3pd);
				DDS = coeff            * 2.0/7.0*(r2dd  +  r4dd);
				DDP = coeff            * (r2dd/7.0   -   4.0/21.0*r4dd);
				DDD = coeff            * (-2.0/7.0*r2dd   +   r4dd/21.0);

				SSS		= 0.0;
				SS_S 	= 0.0;
				SPS		= 0.0;
				S_PS	= 0.0;
				SDS		= 0.0;
				S_DS	= 0.0;
				S_S_S	= 0.0;

				l = v12norm.x;
				m = v12norm.y;
				n = v12norm.z;


/*
				//{
				modified_slaterkoster[0][0] =   SSS;
				modified_slaterkoster[0][1] = - l * SPS;
				modified_slaterkoster[0][2] = - m * SPS;
				modified_slaterkoster[0][3] = - n * SPS;
				modified_slaterkoster[0][4] =   SS_S;
				modified_slaterkoster[0][5] =   sqrt(3.0)*l*m*SDS;			//xy
				modified_slaterkoster[0][6] =   sqrt(3.0)*m*n*SDS;			//yz
				modified_slaterkoster[0][7] =   sqrt(3.0)*n*l*SDS;			//zx
				modified_slaterkoster[0][8] =   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				modified_slaterkoster[0][9] =   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				modified_slaterkoster[1][0] =   l*SPS;
				modified_slaterkoster[1][1] =   (l*l*(PPS-PPP) + PPP);
				modified_slaterkoster[1][2] =   l*m*(PPS-PPP);
				modified_slaterkoster[1][3] =   l*n*(PPS-PPP);
				modified_slaterkoster[1][4] =   l*S_PS;
				modified_slaterkoster[1][5] = - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[1][6] = - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[1][7] = - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[1][8] = - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				modified_slaterkoster[1][9] = - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				modified_slaterkoster[2][0] =   m*SPS;
				modified_slaterkoster[2][1] =   m*l*(PPS-PPP);
				modified_slaterkoster[2][2] =   (m*m*(PPS-PPP) + PPP);
				modified_slaterkoster[2][3] =   m*n*(PPS-PPP);
				modified_slaterkoster[2][4] =   m*S_PS;
				modified_slaterkoster[2][5] = - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[2][6] = - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[2][7] = - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[2][8] = - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				modified_slaterkoster[2][9] = - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				modified_slaterkoster[3][0] =   n*SPS;
				modified_slaterkoster[3][1] =   n*l*(PPS-PPP);
				modified_slaterkoster[3][2] =   n*m*(PPS-PPP);
				modified_slaterkoster[3][3] =   (n*n*(PPS-PPP) + PPP);
				modified_slaterkoster[3][4] =   n*S_PS;
				modified_slaterkoster[3][5] = - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[3][6] = - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[3][7] = - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[3][8] = - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				modified_slaterkoster[3][9] = - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				modified_slaterkoster[4][0] =   SS_S;
				modified_slaterkoster[4][1] = - l * S_PS;
				modified_slaterkoster[4][2] = - m * S_PS;
				modified_slaterkoster[4][3] = - n * S_PS;
				modified_slaterkoster[4][4] =   S_S_S;
				modified_slaterkoster[4][5] =   sqrt(3.0)*l*m*S_DS;
				modified_slaterkoster[4][6] =   sqrt(3.0)*m*n*S_DS;
				modified_slaterkoster[4][7] =   sqrt(3.0)*n*l*S_DS;
				modified_slaterkoster[4][8] =   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				modified_slaterkoster[4][9] =   (n*n-0.5*(l*l+m*m))*S_DS;



				modified_slaterkoster[5][0] =   sqrt(3.0)*l*m*SDS;
				modified_slaterkoster[5][1] =   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[5][2] =   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[5][3] =   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[5][4] =   sqrt(3.0)*l*m*S_DS;
				modified_slaterkoster[5][5] =   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				modified_slaterkoster[5][6] =   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				modified_slaterkoster[5][7] =   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				modified_slaterkoster[5][8] =   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				modified_slaterkoster[5][9] =   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				modified_slaterkoster[6][0] =   sqrt(3.0)*m*n*SDS;
				modified_slaterkoster[6][1] =   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[6][2] =   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[6][3] =   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[6][4] =   sqrt(3.0)*m*n*S_DS;
				modified_slaterkoster[6][5] =   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				modified_slaterkoster[6][6] =   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				modified_slaterkoster[6][7] =   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				modified_slaterkoster[6][8] =   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[6][9] =   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				modified_slaterkoster[7][0] =   sqrt(3.0)*n*l*SDS;
				modified_slaterkoster[7][1] =   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[7][2] =   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[7][3] =   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[7][4] =   sqrt(3.0)*n*l*S_DS;
				modified_slaterkoster[7][5] =   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				modified_slaterkoster[7][6] =   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				modified_slaterkoster[7][7] =   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				modified_slaterkoster[7][8] =   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[7][9] =   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				modified_slaterkoster[8][0] =   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				modified_slaterkoster[8][1] =   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				modified_slaterkoster[8][2] =   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				modified_slaterkoster[8][3] =   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				modified_slaterkoster[8][4] =   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				modified_slaterkoster[8][5] =   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				modified_slaterkoster[8][6] =   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[8][7] =   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[8][8] =   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				modified_slaterkoster[8][9] =   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				modified_slaterkoster[9][0] =   (n*n-0.5*(l*l+m*m))*SDS;
				modified_slaterkoster[9][1] =   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				modified_slaterkoster[9][2] =   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				modified_slaterkoster[9][3] =   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				modified_slaterkoster[9][4] =   (n*n-0.5*(l*l+m*m))*S_DS;
				modified_slaterkoster[9][5] =   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				modified_slaterkoster[9][6] =   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				modified_slaterkoster[9][7] =   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				modified_slaterkoster[9][8] =   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				modified_slaterkoster[9][9] =   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}


				for (int i=0; i<10; i++)
					for (int j=0; j<10; j++)
						Htmp[i*20*N + (currentPsi->d->n-1)*10 + j] += modified_slaterkoster[i][j];
*/

				//{
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 0] +=   SSS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 1] += - l * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 2] += - m * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 3] += - n * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 4] +=   SS_S;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 5] +=   sqrt(3.0)*l*m*SDS;			//xy
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 6] +=   sqrt(3.0)*m*n*SDS;			//yz
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 7] +=   sqrt(3.0)*n*l*SDS;			//zx
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 8] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 9] +=   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 0] +=   l*SPS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1] +=   (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2] +=   l*m*(PPS-PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3] +=   l*n*(PPS-PPP);
				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 4] +=   l*S_PS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 9] += - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 0] +=   m*SPS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1] +=   m*l*(PPS-PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2] +=   (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3] +=   m*n*(PPS-PPP);
				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 4] +=   m*S_PS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 9] += - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 0] +=   n*SPS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1] +=   n*l*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2] +=   n*m*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3] +=   (n*n*(PPS-PPP) + PPP);
				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 4] +=   n*S_PS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 9] += - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 0] +=   SS_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 1] += - l * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 2] += - m * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 3] += - n * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4] +=   S_S_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 5] +=   sqrt(3.0)*l*m*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 6] +=   sqrt(3.0)*m*n*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 7] +=   sqrt(3.0)*n*l*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 8] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 9] +=   (n*n-0.5*(l*l+m*m))*S_DS;



				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 5] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 6] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 7] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8] +=   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 0] +=   (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 1] +=   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 2] +=   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 3] +=   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 4] +=   (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 5] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 6] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 7] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 8] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9] +=   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}



				// unstrained
				coeff = qe*Zeff/4.0/PI/epsilon0/d0; //eV

				PPS = coeff            * (2.0/5.0*r2_P_P);
				PPP = coeff            * (-r2_P_P/5.0);
				PDS = coeff/sqrt(15.0) * (2.0*r_P_D + 9.0/7.0*r3_P_D);
				PDP = coeff/sqrt(5.0)  * (    r_P_D - 3.0/7.0*r3_P_D);
				DDS = coeff            * 2.0/7.0*(r2_D_D  +  r4_D_D);
				DDP = coeff            * (r2_D_D/7.0   -   4.0/21.0*r4_D_D);
				DDD = coeff            * (-2.0/7.0*r2_D_D   +   r4_D_D/21.0);

				l = v12unstrainedNorm.x;
				m = v12unstrainedNorm.y;
				n = v12unstrainedNorm.z;

				//{
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 0] +=   SSS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 1] += - l * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 2] += - m * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 3] += - n * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 4] +=   SS_S;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 5] +=   sqrt(3.0)*l*m*SDS;			//xy
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 6] +=   sqrt(3.0)*m*n*SDS;			//yz
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 7] +=   sqrt(3.0)*n*l*SDS;			//zx
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 8] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 9] +=   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 0] +=   l*SPS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1] +=   (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2] +=   l*m*(PPS-PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3] +=   l*n*(PPS-PPP);
				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 4] +=   l*S_PS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 9] += - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 0] +=   m*SPS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1] +=   m*l*(PPS-PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2] +=   (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3] +=   m*n*(PPS-PPP);
				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 4] +=   m*S_PS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 9] += - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 0] +=   n*SPS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1] +=   n*l*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2] +=   n*m*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3] +=   (n*n*(PPS-PPP) + PPP);
				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 4] +=   n*S_PS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 5] += - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 6] += - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 7] += - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 8] += - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 9] += - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 0] +=   SS_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 1] += - l * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 2] += - m * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 3] += - n * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4] +=   S_S_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 5] +=   sqrt(3.0)*l*m*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 6] +=   sqrt(3.0)*m*n*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 7] +=   sqrt(3.0)*n*l*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 8] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 9] +=   (n*n-0.5*(l*l+m*m))*S_DS;



				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 5] +=   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 6] +=   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7] +=   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 8] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 0] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 1] +=   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 2] +=   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 3] +=   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 4] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 5] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 6] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 7] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8] +=   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 9] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 0] +=   (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 1] +=   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 2] +=   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 3] +=   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 4] +=   (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 5] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 6] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 7] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 8] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9] +=   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}

#endif


				currentNeighbour = currentNeighbour->next;
			}
			

			for(long i=0; i<10; i++) {
				bool ifTheFirst = false;
				for(long j=0; j<20*N; j++) {
#ifdef UCSR
					if ( !((fabs(Htmp[i*20*N + j].real()) < 1e-14) && (fabs(Htmp[i*20*N + j].imag()) < 1e-14)) && (j >= i+10*(currentPsi->d->n-1)) )
#endif
#ifdef CSR
					if ( !((fabs(Htmp[i*20*N + j].real()) < 1e-14) && (fabs(Htmp[i*20*N + j].imag()) < 1e-14))                                    )
#endif
					{
						CSRval[currentCSRval] = Htmp[i*20*N + j];
						if (!ifTheFirst) {
							CSRrows[i+10*(currentPsi->d->n-1)] = currentCSRval+1;
							ifTheFirst = true;
						}
						CSRcols[currentCSRval] = j+1;
						currentCSRval++;
					}
				}
			}


#ifdef PRINT_H
			// Hamiltonean printer
			for(int i=0; i<10; i++) {
				for(int j=0; j<20*N; j++) {
					fprintf(print_h_out, "(");
					if (fabs(Htmp[i*20*N+j].real()) < 1e-14)
						fprintf(print_h_out, "        ");
					else {
						if (Htmp[i*20*N+j].real() >  0.) fprintf(print_h_out, "+");
						if (Htmp[i*20*N+j].real() == 0.) fprintf(print_h_out, " ");
						if (Htmp[i*20*N+j].real() <  0.) fprintf(print_h_out, "-");
						fprintf(print_h_out, "%6.3f ", fabs(Htmp[i*20*N+j].real()));
					}
					if (fabs(Htmp[i*20*N+j].imag()) < 1e-14)
						fprintf(print_h_out, "       )\t");
					else {
						if (Htmp[i*20*N+j].imag() >  0.) fprintf(print_h_out, "+");
						if (Htmp[i*20*N+j].imag() == 0.) fprintf(print_h_out, " ");
						if (Htmp[i*20*N+j].imag() <  0.) fprintf(print_h_out, "-");
						fprintf(print_h_out, "%6.3f)\t", fabs(Htmp[i*20*N+j].imag()));
					}
				}
				fprintf(print_h_out, "\n");
			}
#endif



			currentPsi = currentPsi->next;
		}
		#ifdef UI
		cout << endl;
		#endif
		//}
//////////////////////////  spin minus      /////////////////////////////////////////////////
		//{
		currentPsi = psi.atoms;
		#ifdef UI
		cout << "spin minus" << endl;
		#endif
		csrProgress = N/10;
		while (currentPsi) {
			if (currentPsi->d->n == csrProgress) {
				#ifdef UI
				cout << "." << endl;
				#endif
				csrProgress += N/10;
			}
			for(long i=0; i<10 * 20*N; i++)
				Htmp[i]=0.0;
			V<double> v1(currentPsi->d->x + currentPsi->d->dx, currentPsi->d->y + currentPsi->d->dy, currentPsi->d->z + currentPsi->d->dz);




			bool ifSIcurrentPsi = currentPsi->d->sostav; //ifSI(v1);
			if (ifSIcurrentPsi == true) {
				a_lattice = 0.54310;
				bond_length = sqrt(3.0)/4.0 * a_lattice;
				Es = -2.15168;
				Ep = 4.22925;
				Es_ = 19.1165;
				Ed = 13.7895;
				delta = 0.03978;
				Energy_shift = 27;
			}
			else {
				a_lattice = 0.56579060;
				bond_length = sqrt(3.0)/4.0 * a_lattice;
				Es = -1.95617;
				Ep = 5.3097;
				Es_ = 19.296;
				Ed = 13.5806;
				delta = 0.20264;
				Energy_shift = 27.77;
			}
			

			Htmp[0*20*N + (currentPsi->d->n-1)*10     + 10*N] = Es;
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] = Ep;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] = Ep;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] = Ep;
			Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] = Es_;
			Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] = Ed;
			Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] = Ed;
			Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] = Ed;
			Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] = Ed;
			Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] = Ed;

#ifdef SO
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] =   delta/3.0*li;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] = - delta/3.0*li;
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3       ] = - delta/3.0;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3       ] = - delta/3.0*li;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1       ] =   delta/3.0;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2       ] =   delta/3.0*li;
#endif


			List<Atom> *currentNeighbour = currentPsi->d->neighbours;
			while (currentNeighbour) {
				V<double> v2;
				v2.x = currentNeighbour->d->x + currentNeighbour->d->dx;
				v2.y = currentNeighbour->d->y + currentNeighbour->d->dy;
				v2.z = currentNeighbour->d->z + currentNeighbour->d->dz;

				bool ifSIcurrentNeighbour = currentNeighbour->d->sostav; //ifSI(v2);
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == false)) {
					//GE-GE
					SPS = 2.73135;
					PDS = -2.00115;
					S_PS = 2.68638;
					SDS = -2.64779;
					PDP = 2.10953;
					S_DS = -1.12312;
					SSS = -1.39456;
					S_SS = -2.0183;
					SS_S = -2.0183;
					S_S_S = -3.5668;
					PPS = 4.28921;
					PPP = -1.73707;
					DDS = -1.32941;
					DDP = 2.56261;
					DDD = -1.9512;
					strain_exponent_V_S_S_Sigma = 1.99551;
					strain_exponent_V_S_P_Sigma = 1.29303;
					strain_exponent_V_Sstar_P_Sigma = 5.0000;
					strain_exponent_V_S_Sstar_Sigma = 0.0000;
					strain_exponent_V_Sstar_S_Sigma = 0.0000;
					strain_exponent_V_Sstar_Sstar_Sigma = 2.38823;
					strain_exponent_V_Sstar_D_Sigma = 0.75134;
					strain_exponent_V_S_D_Sigma = 2.79244;
					strain_exponent_V_P_P_Sigma = 1.13641;
					strain_exponent_V_P_P_Pi = 1.74803;
					strain_exponent_V_P_D_Sigma = 2.68784;
					strain_exponent_V_P_D_Pi = 4.36921;
					strain_exponent_V_D_D_Sigma = 5.00000;
					strain_exponent_V_D_D_Pi = 0.69769;
					strain_exponent_V_D_D_Delta = 3.06253;

					strain_constant_K_S_S = 0.0000;
					strain_constant_K_S_P = 1.2507;
					strain_constant_K_Sstar_P = 2.6841;
					strain_constant_K_Sstar_S = 1.1766;
					strain_constant_K_S_Sstar = 1.1766;
					strain_constant_K_Sstar_D = 1.2276;
					strain_constant_K_Sstar_Sstar = 2.6841;
					strain_constant_K_S_D = 0.1524;
					strain_constant_K_P_D = 0.1143;
					strain_constant_K_P_P = 0.3626;
					strain_constant_K_D_D = 1.9527;
					Z_eff = 0.0;
					r_P_D = 0.0;
					r2_P_P = 0.0;
					r2_D_D = 0.0;
					r3_P_D = 0.0;
					r4_D_D = 0.0;
				}
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == true)) {
					//SI-SI
					SPS = 3.02562;
					PDS = -1.35554;
					S_PS = 3.15565;
					SDS = -2.28485;
					PDP = 2.38479;
					S_DS = -0.80993;
					SSS = -1.95933;
					S_SS = -1.5223;
					SS_S = -1.5223;
					S_S_S = -4.24135;
					PPS = 4.10364;
					PPP = -1.51801;
					DDS = -1.68136;
					DDP = 2.5888;
					DDD = -1.814;
					strain_exponent_V_S_S_Sigma	= 0.56247;
					strain_exponent_V_S_P_Sigma = 2.36548;
					strain_exponent_V_Sstar_P_Sigma = 0.34492;
					strain_exponent_V_S_Sstar_Sigma = 0.13203;
					strain_exponent_V_Sstar_S_Sigma = 0.13203;
					strain_exponent_V_Sstar_Sstar_Sigma = 0.19237;
					strain_exponent_V_Sstar_D_Sigma = 1.08601;
					strain_exponent_V_S_D_Sigma = 2.56720;
					strain_exponent_V_P_P_Sigma = 0.2;
					strain_exponent_V_P_P_Pi = 1.6777;
					strain_exponent_V_P_D_Sigma = 0.2;
					strain_exponent_V_P_D_Pi = 4.4325;
					strain_exponent_V_D_D_Sigma = 0.1;
					strain_exponent_V_D_D_Pi = 6;
					strain_exponent_V_D_D_Delta = 5.9997;

					strain_constant_K_S_S = 1.0875;
					strain_constant_K_S_P = 0.4001;
					strain_constant_K_Sstar_P = 1.8597;
					strain_constant_K_Sstar_S = 1.1374;
					strain_constant_K_S_Sstar = 1.1374;
					strain_constant_K_Sstar_D = 0.2978;
					strain_constant_K_Sstar_Sstar = 0.5986;
					strain_constant_K_S_D = 0;
					strain_constant_K_P_D = 0.0564;
					strain_constant_K_P_P = 0;
					strain_constant_K_D_D = 2.1071;
					Z_eff = 3.0;
					r_P_D = 0.1249;
					r2_P_P = 0.1227;
					r2_D_D = 0.5147;
					r3_P_D = 0.0;
					r4_D_D = 0.7661;
				}
				if (ifSIcurrentNeighbour != ifSIcurrentPsi) {
					//GE-SI
					SPS   = 2.878485;	// (2.73135+3.02562)/2.0;
					PDS   = -1.678345;	// (-2.00115-1.35554)/2.0;
					S_PS  = 2.921015;	// (2.68638+3.15565)/2.0;
					SDS   = -2.46632;	// (-2.64779-2.28485)/2.0;
					PDP   = 2.24716;	// (2.10953+2.38479)/2.0;
					S_DS  = -0.966525;	// (-1.12312-0.80993)/2.0;
					SSS   = -1.676945;	// (-1.39456-1.95933)/2.0;
					S_SS  = -1.7703;	// (-2.0183-1.5223)/2.0;
					SS_S  = -1.7703;	// (-2.0183-1.5223)/2.0;
					S_S_S = -3.904075;	// (-3.5668-4.24135)/2.0;
					PPS   = 4.196425;	// (4.28921+4.10364)/2.0;
					PPP   = -1.62754;	// (-1.73707-1.51801)/2.0;
					DDS   = -1.505385;	// (-1.32941-1.68136)/2.0;
					DDP   = 2.575705;	// (2.56261+2.5888)/2.0;
					DDD   = -1.8826;	// (-1.9512-1.814)/2.0;

					strain_exponent_V_S_S_Sigma = 1.27899;			//(0.56247+1.99551)/2.0;
					strain_exponent_V_S_P_Sigma = 1.829255;			//(1.29303+2.36548)/2.0;
					strain_exponent_V_Sstar_P_Sigma = 2.67246;		//(5.0000+0.34492)/2.0;
					strain_exponent_V_S_Sstar_Sigma = 0.13203;		//(0.0000+0.13203)/2.0;
					strain_exponent_V_Sstar_S_Sigma = 0.13203;		//(0.0000+0.13203)/2.0;
					strain_exponent_V_Sstar_Sstar_Sigma = 1.2903;	//(2.38823+0.19237)/2.0;
					strain_exponent_V_Sstar_D_Sigma = 0.918675;		//(0.75134+1.08601)/2.0;
					strain_exponent_V_S_D_Sigma = 2.67982;			//(2.79244+2.56720)/2.0;
					strain_exponent_V_P_P_Sigma = 0.668205;			//(1.13641+0.2)/2.0;
					strain_exponent_V_P_P_Pi = 1.712865;			//(1.74803+1.6777)/2.0;
					strain_exponent_V_P_D_Sigma = 1.44392;			//(2.68784+0.2)/2.0;
					strain_exponent_V_P_D_Pi = 4.400855;			//(4.36921+4.4325)/2.0;
					strain_exponent_V_D_D_Sigma = 2.55;				//(5.00000+0.1)/2.0;
					strain_exponent_V_D_D_Pi = 3.348845;			//(0.69769+6.0)/2.0;
					strain_exponent_V_D_D_Delta = 4.531115;			//(3.06253+5.9997)/2.0;

					strain_constant_K_S_S = 1.0875;				//(0.0+1.0875)/2.0;
					strain_constant_K_S_P = 0.8254;				//(1.2507+0.4001)/2.0;
					strain_constant_K_Sstar_P = 2.2719;			//(2.6841+1.8597)/2.0;
					strain_constant_K_Sstar_S = 1.157;			//(1.1766+1.1374)/2.0;
					strain_constant_K_S_Sstar = 1.157;			//(1.1766+1.1374)/2.0;
					strain_constant_K_Sstar_D = 0.7627;			//(1.2276+0.2978)/2.0;
					strain_constant_K_Sstar_Sstar = 1.64135;	//(2.6841+0.5986)/2.0;
					strain_constant_K_S_D = 0.1524;				//(0.1524+0.0)/2.0;
					strain_constant_K_P_D = 0.08535;			//(0.1143+0.0564)/2.0;
					strain_constant_K_P_P = 0.3626;				//(0.3626+0.0)/2.0;
					strain_constant_K_D_D = 2.0299;				//(1.9527+2.1071)/2.0;
					Z_eff = 0.0;
					r_P_D = 0.0;
					r2_P_P = 0.0;
					r2_D_D = 0.0;
					r3_P_D = 0.0;
					r4_D_D = 0.0;
				}

				complex<double> boundaryCoeff(1.0, 0.0);
				V<double> v12;
				v12 = v2 - v1;
				double volumeSizeX = volumeSize + DEFORMX*volumeSize;
				double volumeSizeY = volumeSize + DEFORMY*volumeSize;
				double volumeSizeZ = volumeSize + DEFORMZ*volumeSize;

				#ifdef BANDSTRUCTURE
					double kL = 0;
					if (fabs(v12.x) > volumeSizeX/2.0) {
						v12.x -= volumeSizeX*v12.x/fabs(v12.x);
						kL += k.x*volumeSizeX*(v12.x/fabs(v12.x));
					}
					if (fabs(v12.y) > volumeSizeY/2.0) {
						v12.y -= volumeSizeY*v12.y/fabs(v12.y);
						kL += k.y*volumeSizeY*(v12.y/fabs(v12.y));
					}
					if (fabs(v12.z) > volumeSizeZ/2.0) {
						v12.z -= volumeSizeZ*v12.z/fabs(v12.z);
						kL += k.z*volumeSizeZ*(v12.z/fabs(v12.z));
					}
					complex<double> ikL(0,kL);
					boundaryCoeff = exp(ikL);
				#else
					if (fabs(v12.x) > volumeSizeX/2.0) {
						v12.x -= volumeSizeX*v12.x/fabs(v12.x);
						boundaryCoeff = BOUNDARYCOEFF;
					}
					if (fabs(v12.y) > volumeSizeY/2.0)
						v12.y -= volumeSizeY*v12.y/fabs(v12.y);
					if (fabs(v12.z) > volumeSizeZ/2.0)
						v12.z -= volumeSizeZ*v12.z/fabs(v12.z);
				#endif
				V<double> v12norm = (1.0/sqrt(v12*v12))*v12;
				double l = v12norm.x;
				double m = v12norm.y;
				double n = v12norm.z;
#ifdef BOYKIN2002
				V<double> v1unstrained;
				v1unstrained.x = currentPsi->d->x;
				v1unstrained.y = currentPsi->d->y;
				v1unstrained.z = currentPsi->d->z;
				V<double> v2unstrained;
				v2unstrained.x = currentNeighbour->d->x;
				v2unstrained.y = currentNeighbour->d->y;
				v2unstrained.z = currentNeighbour->d->z;
				V<double> v12unstrained;
				v12unstrained = v2unstrained - v1unstrained;

				// calculate d0
				if (fabs(v12unstrained.x) > volumeSize/2.0)
					v12unstrained.x -= volumeSize*v12unstrained.x/fabs(v12unstrained.x);
				if (fabs(v12unstrained.y) > volumeSize/2.0)
					v12unstrained.y -= volumeSize*v12unstrained.y/fabs(v12unstrained.y);
				if (fabs(v12unstrained.z) > volumeSize/2.0)
					v12unstrained.z -= volumeSize*v12unstrained.z/fabs(v12unstrained.z);
				V<double> v12unstrainedNorm = (1.0/sqrt(v12unstrained*v12unstrained))*v12unstrained;

				double d  = sqrt(v12*v12);
				double d0 = sqrt(v12unstrained*v12unstrained);


				SPS		*= pow(d0/d, strain_exponent_V_S_P_Sigma);
				PDS		*= pow(d0/d, strain_exponent_V_P_D_Sigma);
				S_PS	*= pow(d0/d, strain_exponent_V_Sstar_P_Sigma);
				SDS		*= pow(d0/d, strain_exponent_V_S_D_Sigma);
				PDP		*= pow(d0/d, strain_exponent_V_P_D_Pi);
				S_DS	*= pow(d0/d, strain_exponent_V_Sstar_D_Sigma);
				SSS		*= pow(d0/d, strain_exponent_V_S_S_Sigma);
				S_SS	*= pow(d0/d, strain_exponent_V_Sstar_S_Sigma);
				SS_S	*= pow(d0/d, strain_exponent_V_S_Sstar_Sigma);
				S_S_S	*= pow(d0/d, strain_exponent_V_Sstar_Sstar_Sigma);
				PPS		*= pow(d0/d, strain_exponent_V_P_P_Sigma);
				PPP		*= pow(d0/d, strain_exponent_V_P_P_Pi);
				DDS		*= pow(d0/d, strain_exponent_V_D_D_Sigma);
				DDP		*= pow(d0/d, strain_exponent_V_D_D_Pi);
				DDD		*= pow(d0/d, strain_exponent_V_D_D_Delta);
#endif

				//{
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * SSS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] = - boundaryCoeff * l * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] = - boundaryCoeff * m * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] = - boundaryCoeff * n * SPS;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * SS_S;
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;			//xy
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;			//yz
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;			//zx
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				Htmp[0*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				Htmp[1*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * l*SPS;
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * l*m*(PPS-PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * l*n*(PPS-PPP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * l*S_PS;
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] = - boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] = - boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] = - boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] = - boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] = - boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				Htmp[2*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * m*SPS;
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * m*l*(PPS-PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * m*n*(PPS-PPP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * m*S_PS;
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] = - boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] = - boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] = - boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] = - boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] = - boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				Htmp[3*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * n*SPS;
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * n*l*(PPS-PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * n*m*(PPS-PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (n*n*(PPS-PPP) + PPP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * n*S_PS;
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] = - boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] = - boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] = - boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] = - boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] = - boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				Htmp[4*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * SS_S;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] = - boundaryCoeff * l * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] = - boundaryCoeff * m * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] = - boundaryCoeff * n * S_PS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * S_S_S;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[4*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;



				Htmp[5*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				Htmp[6*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				Htmp[7*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				Htmp[8*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				Htmp[9*20*N + (currentNeighbour->d->n-1)*10     + 10*N] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 1 + 10*N] =   boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 2 + 10*N] =   boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 3 + 10*N] =   boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 4 + 10*N] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 5 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 6 + 10*N] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 7 + 10*N] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 8 + 10*N] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentNeighbour->d->n-1)*10 + 9 + 10*N] =   boundaryCoeff * ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}




#ifdef BOYKIN2002
				v12unstrained = v2unstrained - v1unstrained;
				#ifdef BANDSTRUCTURE
					kL = 0;
					if (fabs(v12unstrained.x) > volumeSize/2.0) {
						v12unstrained.x -= volumeSize*v12unstrained.x/fabs(v12unstrained.x);
						kL += k.x*volumeSize*(v12unstrained.x/fabs(v12unstrained.x));
					}
					if (fabs(v12unstrained.y) > volumeSize/2.0) {
						v12unstrained.y -= volumeSize*v12unstrained.y/fabs(v12unstrained.y);
						kL += k.y*volumeSize*(v12unstrained.y/fabs(v12unstrained.y));
					}
					if (fabs(v12unstrained.z) > volumeSize/2.0) {
						v12unstrained.z -= volumeSize*v12unstrained.z/fabs(v12unstrained.z);
						kL += k.z*volumeSize*(v12unstrained.z/fabs(v12unstrained.z));
					}
					ikL = li*kL;
					boundaryCoeff = exp(ikL);
				#else
					boundaryCoeff = 1.0;
					if (fabs(v12unstrained.x) > volumeSize/2.0)
						boundaryCoeff = BOUNDARYCOEFF;
				#endif
				// v12unstrainedNorm already calculated (but v12unstrained is incorrect here)
				l = v12unstrainedNorm.x;
				m = v12unstrainedNorm.y;
				n = v12unstrainedNorm.z;
				
				complex<double> slaterkoster[10][10];
				for (int i=0; i<10; i++)
					for (int j=0; j<10; j++)
						slaterkoster[i][j] = 0.0;


				//{
				slaterkoster[0][0] =   boundaryCoeff * SSS;
				slaterkoster[0][1] = - boundaryCoeff * l * SPS;
				slaterkoster[0][2] = - boundaryCoeff * m * SPS;
				slaterkoster[0][3] = - boundaryCoeff * n * SPS;
				slaterkoster[0][4] =   boundaryCoeff * SS_S;
				slaterkoster[0][5] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;			//xy
				slaterkoster[0][6] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;			//yz
				slaterkoster[0][7] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;			//zx
				slaterkoster[0][8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				slaterkoster[0][9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				slaterkoster[1][0] =   boundaryCoeff * l*SPS;
				slaterkoster[1][1] =   boundaryCoeff * (l*l*(PPS-PPP) + PPP);
				slaterkoster[1][2] =   boundaryCoeff * l*m*(PPS-PPP);
				slaterkoster[1][3] =   boundaryCoeff * l*n*(PPS-PPP);
				slaterkoster[1][4] =   boundaryCoeff * l*S_PS;
				slaterkoster[1][5] = - boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				slaterkoster[1][6] = - boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				slaterkoster[1][7] = - boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				slaterkoster[1][8] = - boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				slaterkoster[1][9] = - boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				slaterkoster[2][0] =   boundaryCoeff * m*SPS;
				slaterkoster[2][1] =   boundaryCoeff * m*l*(PPS-PPP);
				slaterkoster[2][2] =   boundaryCoeff * (m*m*(PPS-PPP) + PPP);
				slaterkoster[2][3] =   boundaryCoeff * m*n*(PPS-PPP);
				slaterkoster[2][4] =   boundaryCoeff * m*S_PS;
				slaterkoster[2][5] = - boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				slaterkoster[2][6] = - boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				slaterkoster[2][7] = - boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				slaterkoster[2][8] = - boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				slaterkoster[2][9] = - boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				slaterkoster[3][0] =   boundaryCoeff * n*SPS;
				slaterkoster[3][1] =   boundaryCoeff * n*l*(PPS-PPP);
				slaterkoster[3][2] =   boundaryCoeff * n*m*(PPS-PPP);
				slaterkoster[3][3] =   boundaryCoeff * (n*n*(PPS-PPP) + PPP);
				slaterkoster[3][4] =   boundaryCoeff * n*S_PS;
				slaterkoster[3][5] = - boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				slaterkoster[3][6] = - boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				slaterkoster[3][7] = - boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				slaterkoster[3][8] = - boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				slaterkoster[3][9] = - boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				slaterkoster[4][0] =   boundaryCoeff * SS_S;
				slaterkoster[4][1] = - boundaryCoeff * l * S_PS;
				slaterkoster[4][2] = - boundaryCoeff * m * S_PS;
				slaterkoster[4][3] = - boundaryCoeff * n * S_PS;
				slaterkoster[4][4] =   boundaryCoeff * S_S_S;
				slaterkoster[4][5] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				slaterkoster[4][6] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				slaterkoster[4][7] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				slaterkoster[4][8] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				slaterkoster[4][9] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;



				slaterkoster[5][0] =   boundaryCoeff * sqrt(3.0)*l*m*SDS;
				slaterkoster[5][1] =   boundaryCoeff * (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				slaterkoster[5][2] =   boundaryCoeff * (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				slaterkoster[5][3] =   boundaryCoeff * (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				slaterkoster[5][4] =   boundaryCoeff * sqrt(3.0)*l*m*S_DS;
				slaterkoster[5][5] =   boundaryCoeff * (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				slaterkoster[5][6] =   boundaryCoeff * (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				slaterkoster[5][7] =   boundaryCoeff * (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				slaterkoster[5][8] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				slaterkoster[5][9] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				slaterkoster[6][0] =   boundaryCoeff * sqrt(3.0)*m*n*SDS;
				slaterkoster[6][1] =   boundaryCoeff * (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				slaterkoster[6][2] =   boundaryCoeff * (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				slaterkoster[6][3] =   boundaryCoeff * (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				slaterkoster[6][4] =   boundaryCoeff * sqrt(3.0)*m*n*S_DS;
				slaterkoster[6][5] =   boundaryCoeff * (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				slaterkoster[6][6] =   boundaryCoeff * (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				slaterkoster[6][7] =   boundaryCoeff * (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				slaterkoster[6][8] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				slaterkoster[6][9] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				slaterkoster[7][0] =   boundaryCoeff * sqrt(3.0)*n*l*SDS;
				slaterkoster[7][1] =   boundaryCoeff * (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				slaterkoster[7][2] =   boundaryCoeff * (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				slaterkoster[7][3] =   boundaryCoeff * (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				slaterkoster[7][4] =   boundaryCoeff * sqrt(3.0)*n*l*S_DS;
				slaterkoster[7][5] =   boundaryCoeff * (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				slaterkoster[7][6] =   boundaryCoeff * (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				slaterkoster[7][7] =   boundaryCoeff * (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				slaterkoster[7][8] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				slaterkoster[7][9] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				slaterkoster[8][0] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				slaterkoster[8][1] =   boundaryCoeff * (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				slaterkoster[8][2] =   boundaryCoeff * (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				slaterkoster[8][3] =   boundaryCoeff * (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				slaterkoster[8][4] =   boundaryCoeff * sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				slaterkoster[8][5] =   boundaryCoeff * (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				slaterkoster[8][6] =   boundaryCoeff * (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				slaterkoster[8][7] =   boundaryCoeff * (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				slaterkoster[8][8] =   boundaryCoeff * (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				slaterkoster[8][9] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				slaterkoster[9][0] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*SDS;
				slaterkoster[9][1] =   boundaryCoeff * (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				slaterkoster[9][2] =   boundaryCoeff * (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				slaterkoster[9][3] =   boundaryCoeff * (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				slaterkoster[9][4] =   boundaryCoeff * (n*n-0.5*(l*l+m*m))*S_DS;
				slaterkoster[9][5] =   boundaryCoeff * (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				slaterkoster[9][6] =   boundaryCoeff * (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				slaterkoster[9][7] =   boundaryCoeff * (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				slaterkoster[9][8] =   boundaryCoeff * (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				slaterkoster[9][9] =   boundaryCoeff * ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}

				for (int gamma=0; gamma<10; gamma++)
					for (int alpha=0; alpha<10; alpha++) {
						for (int beta=0; beta<10; beta++) {
							double Kjbia = *strain_constant_K[beta][alpha];
							double Kigjb = *strain_constant_K[gamma][beta];
							double Ejb = *E[beta] - Energy_shift;
							double Eia = *E[alpha] - Energy_shift;
							double Eig = *E[gamma] - Energy_shift;
							
							Htmp[gamma*20*N + (currentPsi->d->n-1)*10 + alpha + 10*N] += (Htmp[gamma*20*N + (currentNeighbour->d->n-1)*10 + beta + 10*N]*Htmp[beta*20*N + (currentNeighbour->d->n-1)*10 + alpha + 10*N] - slaterkoster[gamma][beta]*slaterkoster[beta][alpha]) * (-Kjbia/(Ejb+Eia)/2.0 -Kigjb/(Eig+Ejb)/2.0 -Kigjb*Kjbia*(1.0/(Ejb+Eia)+1.0/(Eig+Ejb))/4.0);
						}
					}
#endif


#ifdef BOYKIN2010
				complex<double> modified_slaterkoster[10][10];
				for (int i=0; i<10; i++)
					for (int j=0; j<10; j++)
						modified_slaterkoster[i][j] = 0.0;

				double coeff = - qe*Zeff/4.0/PI/epsilon0/d; //eV

				double r2pp = r2_P_P*d0*d0/d/d;
				double r_pd = r_P_D*d0/d;
				double r3pd = r3_P_D*d0*d0*d0/d/d/d;
				double r2dd = r2_D_D*d0*d0/d/d;
				double r4dd = r4_D_D*d0*d0*d0*d0/d/d/d/d;

				PPS = coeff            * (2.0/5.0*r2pp);
				PPP = coeff            * (-r2pp/5.0);
				PDS = coeff/sqrt(15.0) * (2.0*r_pd + 9.0/7.0*r3pd);
				PDP = coeff/sqrt(5.0)  * (    r_pd - 3.0/7.0*r3pd);
				DDS = coeff            * 2.0/7.0*(r2dd  +  r4dd);
				DDP = coeff            * (r2dd/7.0   -   4.0/21.0*r4dd);
				DDD = coeff            * (-2.0/7.0*r2dd   +   r4dd/21.0);

				SSS		= 0.0;
				SS_S 	= 0.0;
				SPS		= 0.0;
				S_PS	= 0.0;
				SDS		= 0.0;
				S_DS	= 0.0;
				S_S_S	= 0.0;
/*
				//{
				modified_slaterkoster[0][0] =   SSS;
				modified_slaterkoster[0][1] = - l * SPS;
				modified_slaterkoster[0][2] = - m * SPS;
				modified_slaterkoster[0][3] = - n * SPS;
				modified_slaterkoster[0][4] =   SS_S;
				modified_slaterkoster[0][5] =   sqrt(3.0)*l*m*SDS;			//xy
				modified_slaterkoster[0][6] =   sqrt(3.0)*m*n*SDS;			//yz
				modified_slaterkoster[0][7] =   sqrt(3.0)*n*l*SDS;			//zx
				modified_slaterkoster[0][8] =   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				modified_slaterkoster[0][9] =   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				modified_slaterkoster[1][0] =   l*SPS;
				modified_slaterkoster[1][1] =   (l*l*(PPS-PPP) + PPP);
				modified_slaterkoster[1][2] =   l*m*(PPS-PPP);
				modified_slaterkoster[1][3] =   l*n*(PPS-PPP);
				modified_slaterkoster[1][4] =   l*S_PS;
				modified_slaterkoster[1][5] = - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[1][6] = - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[1][7] = - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[1][8] = - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				modified_slaterkoster[1][9] = - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				modified_slaterkoster[2][0] =   m*SPS;
				modified_slaterkoster[2][1] =   m*l*(PPS-PPP);
				modified_slaterkoster[2][2] =   (m*m*(PPS-PPP) + PPP);
				modified_slaterkoster[2][3] =   m*n*(PPS-PPP);
				modified_slaterkoster[2][4] =   m*S_PS;
				modified_slaterkoster[2][5] = - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[2][6] = - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[2][7] = - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[2][8] = - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				modified_slaterkoster[2][9] = - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				modified_slaterkoster[3][0] =   n*SPS;
				modified_slaterkoster[3][1] =   n*l*(PPS-PPP);
				modified_slaterkoster[3][2] =   n*m*(PPS-PPP);
				modified_slaterkoster[3][3] =   (n*n*(PPS-PPP) + PPP);
				modified_slaterkoster[3][4] =   n*S_PS;
				modified_slaterkoster[3][5] = - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[3][6] = - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[3][7] = - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[3][8] = - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				modified_slaterkoster[3][9] = - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				modified_slaterkoster[4][0] =   SS_S;
				modified_slaterkoster[4][1] = - l * S_PS;
				modified_slaterkoster[4][2] = - m * S_PS;
				modified_slaterkoster[4][3] = - n * S_PS;
				modified_slaterkoster[4][4] =   S_S_S;
				modified_slaterkoster[4][5] =   sqrt(3.0)*l*m*S_DS;
				modified_slaterkoster[4][6] =   sqrt(3.0)*m*n*S_DS;
				modified_slaterkoster[4][7] =   sqrt(3.0)*n*l*S_DS;
				modified_slaterkoster[4][8] =   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				modified_slaterkoster[4][9] =   (n*n-0.5*(l*l+m*m))*S_DS;



				modified_slaterkoster[5][0] =   sqrt(3.0)*l*m*SDS;
				modified_slaterkoster[5][1] =   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[5][2] =   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[5][3] =   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[5][4] =   sqrt(3.0)*l*m*S_DS;
				modified_slaterkoster[5][5] =   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				modified_slaterkoster[5][6] =   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				modified_slaterkoster[5][7] =   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				modified_slaterkoster[5][8] =   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				modified_slaterkoster[5][9] =   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				modified_slaterkoster[6][0] =   sqrt(3.0)*m*n*SDS;
				modified_slaterkoster[6][1] =   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[6][2] =   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				modified_slaterkoster[6][3] =   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[6][4] =   sqrt(3.0)*m*n*S_DS;
				modified_slaterkoster[6][5] =   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				modified_slaterkoster[6][6] =   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				modified_slaterkoster[6][7] =   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				modified_slaterkoster[6][8] =   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[6][9] =   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				modified_slaterkoster[7][0] =   sqrt(3.0)*n*l*SDS;
				modified_slaterkoster[7][1] =   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				modified_slaterkoster[7][2] =   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				modified_slaterkoster[7][3] =   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				modified_slaterkoster[7][4] =   sqrt(3.0)*n*l*S_DS;
				modified_slaterkoster[7][5] =   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				modified_slaterkoster[7][6] =   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				modified_slaterkoster[7][7] =   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				modified_slaterkoster[7][8] =   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[7][9] =   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				modified_slaterkoster[8][0] =   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				modified_slaterkoster[8][1] =   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				modified_slaterkoster[8][2] =   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				modified_slaterkoster[8][3] =   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				modified_slaterkoster[8][4] =   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				modified_slaterkoster[8][5] =   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				modified_slaterkoster[8][6] =   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[8][7] =   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				modified_slaterkoster[8][8] =   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				modified_slaterkoster[8][9] =   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				modified_slaterkoster[9][0] =   (n*n-0.5*(l*l+m*m))*SDS;
				modified_slaterkoster[9][1] =   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				modified_slaterkoster[9][2] =   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				modified_slaterkoster[9][3] =   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				modified_slaterkoster[9][4] =   (n*n-0.5*(l*l+m*m))*S_DS;
				modified_slaterkoster[9][5] =   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				modified_slaterkoster[9][6] =   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				modified_slaterkoster[9][7] =   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				modified_slaterkoster[9][8] =   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				modified_slaterkoster[9][9] =   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}


				for (int i=0; i<10; i++)
					for (int j=0; j<10; j++)
						Htmp[i*20*N + (currentPsi->d->n-1)*10 + j + 10*N] += modified_slaterkoster[i][j];
*/

				//{
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   SSS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] += - l * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] += - m * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] += - n * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   SS_S;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   sqrt(3.0)*l*m*SDS;			//xy
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   sqrt(3.0)*m*n*SDS;			//yz
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   sqrt(3.0)*n*l*SDS;			//zx
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   l*SPS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   l*m*(PPS-PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   l*n*(PPS-PPP);
				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   l*S_PS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   m*SPS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   m*l*(PPS-PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   m*n*(PPS-PPP);
				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   m*S_PS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   n*SPS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   n*l*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   n*m*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (n*n*(PPS-PPP) + PPP);
				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   n*S_PS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   SS_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] += - l * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] += - m * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] += - n * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   S_S_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   sqrt(3.0)*l*m*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   sqrt(3.0)*m*n*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   sqrt(3.0)*n*l*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (n*n-0.5*(l*l+m*m))*S_DS;



				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}

				// unstrained
				coeff = qe*Zeff/4.0/PI/epsilon0/d0; //eV

				PPS = coeff            * (2.0/5.0*r2_P_P);
				PPP = coeff            * (-r2_P_P/5.0);
				PDS = coeff/sqrt(15.0) * (2.0*r_P_D + 9.0/7.0*r3_P_D);
				PDP = coeff/sqrt(5.0)  * (    r_P_D - 3.0/7.0*r3_P_D);
				DDS = coeff            * 2.0/7.0*(r2_D_D  +  r4_D_D);
				DDP = coeff            * (r2_D_D/7.0   -   4.0/21.0*r4_D_D);
				DDD = coeff            * (-2.0/7.0*r2_D_D   +   r4_D_D/21.0);

				l = v12unstrainedNorm.x;
				m = v12unstrainedNorm.y;
				n = v12unstrainedNorm.z;

				//{
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   SSS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] += - l * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] += - m * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] += - n * SPS;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   SS_S;
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   sqrt(3.0)*l*m*SDS;			//xy
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   sqrt(3.0)*m*n*SDS;			//yz
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   sqrt(3.0)*n*l*SDS;			//zx
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;	// x^2 - y^2
				//Htmp[0*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (n*n-0.5*(l*l+m*m))*SDS;	//3z^2 - r^2

				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   l*SPS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (l*l*(PPS-PPP) + PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   l*m*(PPS-PPP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   l*n*(PPS-PPP);
				//Htmp[1*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   l*S_PS;
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[1*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);

				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   m*SPS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   m*l*(PPS-PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (m*m*(PPS-PPP) + PPP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   m*n*(PPS-PPP);
				//Htmp[2*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   m*S_PS;
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[2*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);

				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   n*SPS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   n*l*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   n*m*(PPS-PPP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (n*n*(PPS-PPP) + PPP);
				//Htmp[3*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   n*S_PS;
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] += - (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] += - (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] += - (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] += - (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				Htmp[3*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] += - (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);

				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   SS_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] += - l * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] += - m * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] += - n * S_PS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   S_S_S;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   sqrt(3.0)*l*m*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   sqrt(3.0)*m*n*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   sqrt(3.0)*n*l*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				//Htmp[4*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (n*n-0.5*(l*l+m*m))*S_DS;



				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*l*m*SDS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*l*m*PDS + m*(1.0-2.0*l*l)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*l*m*PDS + l*(1.0-2.0*m*m)*PDP);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*l*m*PDS - 2.0*l*m*n*PDP);
				//Htmp[5*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*l*m*S_DS;
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*l*m*l*m*DDS + (l*l+m*m-4.0*l*l*m*m)*DDP + (n*n+l*l*m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*l*m*m*n*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*l*m*n*l*DDS + m*n*(1.0-4.0*l*l)*DDP + m*n*(l*l-1.0)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[5*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);

				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*m*n*SDS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*m*n*PDS - 2.0*l*m*n*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*m*n*PDS + n*(1.0-2.0*m*m)*PDP);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*m*n*PDS + m*(1.0-2.0*n*n)*PDP);
				//Htmp[6*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*m*n*S_DS;
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*m*n*l*m*DDS + l*n*(1.0-4.0*m*m)*DDP + l*n*(m*m-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*m*n*m*n*DDS + (m*m+n*n-4.0*m*m*n*n)*DDP + (l*l+m*m*n*n)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*m*n*n*l*DDS + m*l*(1.0-4.0*n*n)*DDP + m*l*(n*n-1.0)*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[6*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);

				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)*n*l*SDS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)*l*n*l*PDS + n*(1.0-2.0*l*l)*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)*m*n*l*PDS - 2.0*l*m*n*PDP);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)*n*n*l*PDS + l*(1.0-2.0*n*n)*PDP);
				//Htmp[7*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)*n*l*S_DS;
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (3.0*n*l*l*m*DDS + n*m*(1.0-4.0*l*l)*DDP + n*m*(l*l-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (3.0*n*l*m*n*DDS + l*m*(1.0-4.0*n*n)*DDP + l*m*(n*n-1.0)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (3.0*n*l*n*l*DDS + (n*n+l*l-4.0*n*n*l*l)*DDP + (m*m+n*n*l*l)*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[7*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);

				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*SDS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (sqrt(3.0)/2.0*l*(l*l-m*m)*PDS + l*(1.0-l*l+m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (sqrt(3.0)/2.0*m*(l*l-m*m)*PDS - m*(1.0+l*l-m*m)*PDP);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (sqrt(3.0)/2.0*n*(l*l-m*m)*PDS - n*(l*l-m*m)*PDP);
				//Htmp[8*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   sqrt(3.0)/2.0*(l*l-m*m)*S_DS;
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (1.5*l*m*(l*l-m*m)*DDS + 2.0*l*m*(m*m-l*l)*DDP + 0.5*l*m*(l*l-m*m)*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (1.5*m*n*(l*l-m*m)*DDS - m*n*(1.0+2.0*(l*l-m*m))*DDP + m*n*(1.0+0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (1.5*n*l*(l*l-m*m)*DDS + n*l*(1.0-2.0*(l*l-m*m))*DDP - n*l*(1.0-0.5*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (0.75*(l*l-m*m)*(l*l-m*m)*DDS + (l*l+m*m-(l*l-m*m)*(l*l-m*m))*DDP + (n*n+0.25*(l*l-m*m)*(l*l-m*m))*DDD);
				Htmp[8*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);

				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 0 + 10*N] +=   (n*n-0.5*(l*l+m*m))*SDS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 1 + 10*N] +=   (l*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*l*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 2 + 10*N] +=   (m*(n*n-0.5*(l*l+m*m))*PDS - sqrt(3.0)*m*n*n*PDP);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 3 + 10*N] +=   (n*(n*n-0.5*(l*l+m*m))*PDS + sqrt(3.0)*n*(l*l+m*m)*PDP);
				//Htmp[9*20*N + (currentPsi->d->n-1)*10 + 4 + 10*N] +=   (n*n-0.5*(l*l+m*m))*S_DS;
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 5 + 10*N] +=   (sqrt(3.0)*l*m*(n*n-0.5*(l*l+m*m))*DDS - 2.0*sqrt(3.0)*l*m*n*n*DDP + sqrt(3.0)/2.0*l*m*(1.0+n*n)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 6 + 10*N] +=   (sqrt(3.0)*m*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*m*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*m*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 7 + 10*N] +=   (sqrt(3.0)*l*n*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*l*n*(l*l+m*m-n*n)*DDP - sqrt(3.0)/2.0*l*n*(l*l+m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 8 + 10*N] +=   (sqrt(3.0)/2.0*(l*l-m*m)*(n*n-0.5*(l*l+m*m))*DDS + sqrt(3.0)*n*n*(m*m-l*l)*DDP + sqrt(3.0)/4.0*(1.0+n*n)*(l*l-m*m)*DDD);
				Htmp[9*20*N + (currentPsi->d->n-1)*10 + 9 + 10*N] +=   ((n*n-0.5*(l*l+m*m))*(n*n-0.5*(l*l+m*m))*DDS + 3.0*n*n*(l*l+m*m)*DDP + 0.75*(l*l+m*m)*(l*l+m*m)*DDD);
				//}

#endif

				currentNeighbour = currentNeighbour->next;
			}


			for(long i=0; i<10; i++) {
				bool ifTheFirst = false;
				for(long j=0; j<20*N; j++) {
#ifdef UCSR
					if ( !((fabs(Htmp[i*20*N + j].real()) < 1e-14) && (fabs(Htmp[i*20*N + j].imag()) < 1e-14)) && (j >= i+10*(currentPsi->d->n-1)+10*N) )
#endif
#ifdef CSR
					if ( !((fabs(Htmp[i*20*N + j].real()) < 1e-14) && (fabs(Htmp[i*20*N + j].imag()) < 1e-14))                                        )
#endif
					{
						CSRval[currentCSRval] = Htmp[i*20*N + j];
						if (!ifTheFirst) {
							CSRrows[i+10*(currentPsi->d->n-1)+10*N] = currentCSRval+1;
							ifTheFirst = true;
						}
						CSRcols[currentCSRval] = j+1;
						currentCSRval++;
					}
				}
			}

	
#ifdef PRINT_H
			// Hamiltonean printer
			for(int i=0; i<10; i++) {
				for(int j=0; j<20*N; j++) {
					fprintf(print_h_out, "(");
					if (fabs(Htmp[i*20*N+j].real()) < 1e-14)
						fprintf(print_h_out, "        ");
					else {
						if (Htmp[i*20*N+j].real() >  0.) fprintf(print_h_out, "+");
						if (Htmp[i*20*N+j].real() == 0.) fprintf(print_h_out, " ");
						if (Htmp[i*20*N+j].real() <  0.) fprintf(print_h_out, "-");
						fprintf(print_h_out, "%6.3f ", fabs(Htmp[i*20*N+j].real()));
					}
					if (fabs(Htmp[i*20*N+j].imag()) < 1e-14)
						fprintf(print_h_out, "       )\t");
					else {
						if (Htmp[i*20*N+j].imag() >  0.) fprintf(print_h_out, "+");
						if (Htmp[i*20*N+j].imag() == 0.) fprintf(print_h_out, " ");
						if (Htmp[i*20*N+j].imag() <  0.) fprintf(print_h_out, "-");
						fprintf(print_h_out, "%6.3f)\t", fabs(Htmp[i*20*N+j].imag()));
					}
				}
				fprintf(print_h_out, "\n");
			}

#endif



			currentPsi = currentPsi->next;
		}
		#ifdef PRINT_H
			fclose(print_h_out);
		#endif
		#ifdef UI
		cout << endl;
		#endif
		//}
//////////////////////////  HERMITIAN   /////////////////////////////////////////////////////
#ifdef HERMITIAN_TEST
#ifdef CSR	
		CSRrows[20*N] = currentCSRval+1;
		long exactCSRsize = currentCSRval;
		long j = 0; //one-based
		for(long currentCSRindx=0; currentCSRindx<exactCSRsize; currentCSRindx++) {
			if(CSRrows[(j-1)+1] == currentCSRindx+1)
				j++;
			long i = CSRcols[currentCSRindx]; //one-based
			
			bool ifNotHermitian = false;
			// i->j; j->i;
			long transposed_indx = CSRrows[i-1]-1;
			while ( (CSRcols[transposed_indx] != j) && (transposed_indx+1 < CSRrows[(i-1)+1]) ) {
				transposed_indx++;
			}
			if (CSRcols[transposed_indx] != j)
				ifNotHermitian = true;
			if (fabs(CSRval[currentCSRindx].real() - CSRval[transposed_indx].real()) > 1e-14)
				ifNotHermitian = true;
			if (fabs(CSRval[currentCSRindx].imag() + CSRval[transposed_indx].imag()) > 1e-14)
				ifNotHermitian = true;
			
			if (ifNotHermitian) {
				fprintf(stderr, "error: Hamiltonean is not hermitian (i=%d, j=%d, CSRindx=%d)\n", i, j, currentCSRindx);
			}
		}
#endif
#endif
//////////////////////////  OUTPUT CSR  /////////////////////////////////////////////////////
		//{
		waitpid(-1, NULL, 0);
		char filename[100];
		//sprintf(filename, "../FCSR20x%i_%i", N, currentCSRval);
		sprintf(filename, "../FCSR");
		char filename_out[100];
#ifdef CSRTEXT		
		sprintf(filename_out, "%s.dat", filename);
		out = fopen(filename_out, "w");
		if (out==NULL) fprintf(stderr,"error: can't open the file %s for writing.\n",filename_out),exit(1);
		#ifdef UI
		printf("writing CSR to a text file %s\n", filename_out);
		#endif


		if (currentCSRval != CSRsize)
			fprintf(stderr, "error: currentCSRval == %i != CSRsize\n", currentCSRval);
		fprintf(out, "%i ", 20*N);
		fprintf(out, "%i\n", currentCSRval);
		for(long i=0; i<currentCSRval; i++)
			fprintf(out, "%f ", CSRval[i].real());
		fprintf(out, "\n");
		for(long i=0; i<currentCSRval; i++)
			fprintf(out, "%f ", CSRval[i].imag());
		fprintf(out, "\n");
		CSRrows[20*N] = currentCSRval+1;
		for(long i=0; i<20*N; i++)
			fprintf(out, "%i ", CSRrows[i]);
		fprintf(out, "\n");
		for(long i=0; i<currentCSRval; i++)
			fprintf(out, "%i ", CSRcols[i]);
		fclose(out);
#endif
#ifdef CSRBINARY
		sprintf(filename_out, "%s.bin", filename);
		out = fopen(filename_out, "wb"); if (out==NULL) fprintf(stderr,"Error!\ncan't open the file %s for writing.\n",filename_out),exit(1);
		#ifdef UI
		printf("writing CSR to a binary file %s\n", filename_out);
		#endif


		long CSRN = 20*N;
		fwrite(&CSRN, sizeof(CSRN), 1, out);
		fwrite(&currentCSRval, sizeof(currentCSRval), 1, out);
		CSRrows[20*N] = currentCSRval+1;
		fwrite(CSRval, sizeof(CSRval[0]), currentCSRval, out);
		fwrite(CSRrows, sizeof(CSRrows[0]), 20*N+1, out);
		fwrite(CSRcols, sizeof(CSRcols[0]), currentCSRval, out);
		fclose(out);
#endif
		//}
//////////////////////////  FINISH      /////////////////////////////////////////////////////
		//{
		#ifdef UI
		cout << "finished, freing memory" << endl;
		#endif
		for (int i=0; i<10; i++)
			for (int j=0; j<10; j++) {
				strain_constant_K[i][j] = NULL;
			}
		free(Htmp);
		free(CSRval);
		free(CSRrows);
		free(CSRcols);
		//}

#ifdef BANDSTRUCTURE
		#ifdef UI
		cout << "running solver" << endl;
		#endif
		char char_k_x[10];
		char char_k_y[10];
		char char_k_z[10];
		sprintf(char_k_x, "%f", k.x);
		sprintf(char_k_y, "%f", k.y);
		sprintf(char_k_z, "%f", k.z);
		char const * solver_path = cfg.lookup("solver_path");
		char const * solver_name = cfg.lookup("solver_name");
		//int execerrno = execl(solver_path, solver_name, char_k_x, char_k_y, char_k_z);
		pid_t solverpid = fork();
		if (solverpid == -1) { fprintf(stderr, "fork() failed\n"); exit(1); }
		if (solverpid ==  0) {
			int execerrno = execl("/mnt/storage/home/aakoshkarev/cpp/eigenmodel/slepc_solver/main", "main", char_k_x, char_k_y, char_k_z, (char*)0 );
			printf("error: can't execl %s\n", solver_path);
		}
		else {
			//int status;
			//waitpid(-1, &status, 0);
			//printf("solverpid = %i; status = %i\n", solverpid, status);
		}
	}          //k forcircle
#endif

	return 0;
}







































