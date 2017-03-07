#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;
#include "myvector.h"
#include "mylist.h"
#include "matrixAlgebra.h"
#include <complex>
#define PI 3.141592653589793
#include <cmath>

#define bufsize 1000
FILE *out;
FILE *in;

const double a = 0.5431; //1
const int n = 20;
const double volumeSize = n*a;
const int N = 64000; //n*n*n*8;



#include "atom.h"
#include <string>


V<double> k(0,0,0);
complex<double> li(0,1.0);

//#define PRINT_H
#define UCSR
#define SO
#define UI

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


////////////////////////////////////////////////

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


complex<double> operator *(double f, complex<double> c) {
	complex<double> tmp(f*c.real(), f*c.imag());
	return tmp;
}


int main() {




	char* out_filename = "./CSR.dat";
	

//	for(k.x=-3.3; k.x<=6.3; k.x+=0.2) {
/*	for(k.x=4.76-0.1; k.x<=4.76+0.1; k.x+=0.2) {
		k.y = k.x;
		k.z = k.x;	
		if(k.x > 0) {
			k.y = 0.0;
			k.z = 0.0;
		}
*/

/*
	k.x=4.76;
	for(k.y=-1.0; k.y<=1.0; k.y+=0.2) {
*/

		out = fopen(out_filename, "w");
		if (out==NULL) fprintf(stderr,"error: can't open the file %s for writing.\n",out_filename),exit(1);
		
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

//////////////////////////  DEFORMATION     /////////////////////////////////////////////////
		char* strain_filename = "./../QD_strain.xyz";
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


////////////////////////// NEIGHBOURS    ////////////////////////////////////////
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
				if ((dist < rSquared) && (rSquared < dist2)) {
					Add(current2->d, &(current->d->second_neighbours));
				}
				current2 = current2->next;
			}
			if (i!=4)
				fprintf(stderr, "error: atom %f has %f neighbours.\n", current->d->n, i), exit(1);
			current = current->next;
		}
		#ifdef UI
		cout << endl;
		#endif

///////////////////////////////////////// CSR init   /////////////////////////////////////

		#ifdef UI
		cout << "initializing CSR" << endl;
		#endif
#ifdef SO
		//     onsite: diag  nondiag           first neighbour    SO
		long CSRsize = 20*N + 45*N             + 4*(100)*N   +   6*N;     //UCSR
#else
		long CSRsize = 20*N + 45*N             + 4*(100)*N;		//UCSR
#endif
		#ifdef UI
		cout << CSRsize << " nonzero elements" << endl;
		#endif
		long currentCSRval = 0;

		complex<double>		* CSRval = (complex<double>*)malloc(CSRsize*sizeof(complex<double>));
		long				* CSRrows = (long*)malloc((10*N+1)*sizeof(long));
		long				* CSRcols = (long*)malloc(CSRsize*sizeof(long));
		if (!((CSRval) && (CSRrows) && (CSRcols)))
			cerr << "error: not enough heap memory" << endl;



///////////////////////////////////////////////////////////////////////////////////////////
		complex<double> *Htmp = (complex<double>*)malloc((10 * 20*N)*sizeof(complex<double>));

/////////////////////////////// spin plus ///////////////////////////////////
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
			V<double> v1(currentPsi->d->x, currentPsi->d->y, currentPsi->d->z);




			double a_lattice = 0.56579060;
			double bond_length = sqrt(3.0)/4.0 * a_lattice;
			double Es = -1.95617;
			double Ep = 5.3097;
			double Es_ = 19.296;
			double Ed = 13.5806;
			double delta = 0.20264;


			double strain_exponent_V_S_S_Sigma = 1.99551;
			double strain_exponent_V_S_P_Sigma = 1.29303;
			double strain_exponent_V_Sstar_P_Sigma = 5.0000;
			double strain_exponent_V_S_Sstar_Sigma = 0.0000;
			double strain_exponent_V_Sstar_S_Sigma = 0.0000;
			double strain_exponent_V_Sstar_Sstar_Sigma = 2.38823;
			double strain_exponent_V_Sstar_D_Sigma = 0.75134;
			double strain_exponent_V_S_D_Sigma = 2.79244;
			double strain_exponent_V_P_P_Sigma = 1.13641;
			double strain_exponent_V_P_P_Pi = 1.74803;
			double strain_exponent_V_P_D_Sigma = 2.68784;
			double strain_exponent_V_P_D_Pi = 4.36921;
			double strain_exponent_V_D_D_Sigma = 5.00000;
			double strain_exponent_V_D_D_Pi = 0.69769;
			double strain_exponent_V_D_D_Delta = 3.06253;


			double strain_constant_C_S_S = 0.0000;
			double strain_constant_C_S_P = 2.03278;
			double strain_constant_C_Sstar_P = 6.28624;
			double strain_constant_C_Sstar_S = 1.86887;
			double strain_constant_C_S_Sstar = 1.86887;
			double strain_constant_C_Sstar_D = 1.98112;
			double strain_constant_C_Sstar_Sstar = 6.28624;
			double strain_constant_C_S_D = 0.16396;
			double strain_constant_C_P_D = 0.12084;
			double strain_constant_C_P_P = 0.42830;
			double strain_constant_C_D_D = 3.85908;

			double Energy_shift = 27.77;

			// strain_constant_K_S_S = 0;
			// strain_constant_K_S_P = 0;
			// strain_constant_K_Sstar_P = 0;
			// strain_constant_K_Sstar_S = 0;
			// strain_constant_K_S_Sstar = 0;
			// strain_constant_K_Sstar_D = 0;
			// strain_constant_K_Sstar_Sstar = 0;
			// strain_constant_K_S_D = 0;
			// strain_constant_K_P_D = 0;
			// strain_constant_K_P_P = 0;
			// strain_constant_K_D_D = 0;

			bool ifSIcurrentPsi = currentPsi->d->sostav; //ifSI(v1);
			if (ifSIcurrentPsi == true) {
				a_lattice = 0.54310;
				bond_length = sqrt(3.0)/4.0 * a_lattice;
				Es = -2.15168;
				Ep = 4.22925;
				Es_ = 19.1165;
				Ed = 13.7895;
				delta = 0.03978;

				strain_exponent_V_S_S_Sigma = 0.56247;
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

				double strain_constant_K_S_S = 1.0875;
				double strain_constant_K_S_P = 0.4001;
				double strain_constant_K_Sstar_P = 1.8597;
				double strain_constant_K_Sstar_S = 1.1374;
				double strain_constant_K_S_Sstar = 1.1374;
				double strain_constant_K_Sstar_D = 0.2978;
				double strain_constant_K_Sstar_Sstar = 0.5986;
				double strain_constant_K_S_D = 0;
				double strain_constant_K_P_D = 0.0564;
				double strain_constant_K_P_P = 0;
				double strain_constant_K_D_D = 2.1071;
				double Z_eff = 3.0;
				double r_P_D = 0.1249;
				double r2_P_P = 0.1227;
				double r2_D_D = 0.5147;
				double r3_P_D = 0.0;
				double r4_D_D = 0.7661;
				Energy_shift = 27;
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
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 2       ] = - delta/3.0*li;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 1       ] =   delta/3.0*li;
			Htmp[1*20*N + (currentPsi->d->n-1)*10 + 3 + 20*N] =   delta/3.0;
			Htmp[2*20*N + (currentPsi->d->n-1)*10 + 3 + 20*N] = - delta/3.0*li;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 1 + 20*N] = - delta/3.0;
			Htmp[3*20*N + (currentPsi->d->n-1)*10 + 2 + 20*N] =   delta/3.0*li;
#endif


			/// first neighbour
			List<Atom> *currentNeighbour = currentPsi->d->neighbours;
			while (currentNeighbour) {
				V<double> v2;
				v2.x = currentNeighbour->d->x;
				v2.y = currentNeighbour->d->y;
				v2.z = currentNeighbour->d->z;

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

				currentNeighbour = currentNeighbour->next;
			}
			

			for(long i=0; i<5; i++) {
				bool ifTheFirst = false;
				for(long j=0; j<10*N; j++) {
#ifdef UCSR
					if ( !((Htmp[i*10*N + j].real() == 0) && (Htmp[i*10*N + j].imag() == 0)) && (j >= i+5*(currentPsi->d->n-1)) ) {
#endif
#ifdef CSR
					if ( !((Htmp[i*10*N + j].real() == 0) && (Htmp[i*10*N + j].imag() == 0))                                    ) {
#endif
						CSRval[currentCSRval] = Htmp[i*10*N + j];
						if (!ifTheFirst) {
							CSRrows[i+5*(currentPsi->d->n-1)] = currentCSRval+1;
							ifTheFirst = true;
						}
						CSRcols[currentCSRval] = j+1;
						currentCSRval++;
					}
				}
			}

#ifdef PRINT_H
			// Hamiltonean printer
			for(int i=0; i<5; i++) {
				for(int j=0; j<10*N; j++) {
					fprintf(out, "(");
					if (Htmp[i*10*N+j].real()>=0) fprintf(out, "+"); else fprintf(out, "-");
					fprintf(out, "%3.3f ", abs(Htmp[i*10*N+j].real()));
					if (Htmp[i*10*N+j].imag()>=0) fprintf(out, "+"); else fprintf(out, "-");
					fprintf(out, "%3.3f)\t", abs(Htmp[i*10*N+j].imag()));
				}
				fprintf(out, "\n");
			}
#endif

			currentPsi = currentPsi->next;
		}
		#ifdef UI
		cout << endl;
		#endif

///////////////////////////  spin   minus  ///////////////////////////////////
/*		currentPsi = psi.atoms;
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
			for(long i=0; i<5 * 10*N; i++)
				Htmp[i]=0.0;
			V<double> v1;
			v1.x = currentPsi->d->x;
			v1.y = currentPsi->d->y;
			v1.z = currentPsi->d->z;


			double delta		=	0.29;
			double Es		=	-5.88;
			double Ep		=	1.61;
			double Es_		=	6.39;
			double SSS		=	-6.78/4.0;
			double SPS		=	sqrt(3.0)/4.0*5.4649;
			double PPS		=	(1.61+2.0*4.9)/4.0;
			double PPP		=	(1.61-    4.9)/4.0;
			double S_S_S		=	0;
			double S_PS		=	sqrt(3.0)/4.0*5.2191;
			double SS_S		=	0;

			double A		=	1e-40;
			double B		=	1e-40;
			double C		=	1e-40;
			double B_		=	1e-40;
			double C_		=	1e-40;
			double D		=	1e-40;
			double E		=	1e-40;
			double F		=	1e-40;
			double G		=	1e-40;

			bool ifSIcurrentPsi = currentPsi->d->sostav; //ifSI(v1);
			if (ifSIcurrentPsi == true) {
				delta		=	0.04503;
				Es		=	-4.81341;
				Ep		=	1.77563;
				Es_		=	5.61342;
				SSS	=	-8.33255/4.0;
				SPS	=	sqrt(3.0)/4.0*5.86140;
				PPS	=	(1.69916+2.0*5.29091)/4.0;
				PPP	=	(1.69916-    5.29091)/4.0;
				S_PS	=	sqrt(3.0)/4.0*4.88308;
				A		=	0.01591/4.0;
				B		=	0.08002/4.0;
				C		=	-1.31699/4.0;		//-
				B_		=	-0.00579/4.0;
				C_		=	-0.50103/4.0;		//-
				D		=	0.00762/4.0;
				E		=	-0.10662/4.0;
				F		=	0.55067/4.0;
				G		=	-2.27784/4.0;		//-
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
				v2.x = currentNeighbour->d->x;
				v2.y = currentNeighbour->d->y;
				v2.z = currentNeighbour->d->z;

				bool ifSIcurrentNeighbour = currentNeighbour->d->sostav; //ifSI(v2);
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == false)) {
					SSS		=	-6.78/4.0;
					SPS		=	sqrt(3.0)/4.0*5.4649;
					PPS		=	(1.61+2.0*4.9)/4.0;
					PPP		=	(1.61-    4.9)/4.0;
					S_PS		=	sqrt(3.0)/4.0*5.2191;
				}
				if ((ifSIcurrentNeighbour == ifSIcurrentPsi) && (ifSIcurrentPsi == true)) {
					SSS	=	-8.33255/4.0;
					SPS	=	sqrt(3.0)/4.0*5.86140;
					PPS	=	(1.69916+2.0*5.29091)/4.0;
					PPP	=	(1.69916-    5.29091)/4.0;
					S_PS	=	sqrt(3.0)/4.0*4.88308;
				}
				if (ifSIcurrentNeighbour != ifSIcurrentPsi) {
					SSS		=	(-6.78-8.33255)/8.0;
					SPS		=	sqrt(3.0)/8.0*(5.86140+5.4649);
					PPS		=	(1.69916+2.0*5.29091)/8.0 + (1.61+2.0*4.9)/8.0;
					PPP		=	(1.69916-    5.29091)/8.0 + (1.61-    4.9)/8.0;
					S_PS	=	sqrt(3.0)/8.0*(4.88308+5.2191);
				}


				complex<double> boundaryCoeff(1.0, 0.0);
				V<double> v12;
				v12 = v2 - v1;
				if (fabs(v12.x) > volumeSize/2.0) {
					v12.x -= volumeSize*v12.x/fabs(v12.x);
					boundaryCoeff = li;
				}
				if (fabs(v12.y) > volumeSize/2.0)
					v12.y -= volumeSize*v12.y/fabs(v12.y);
				if (fabs(v12.z) > volumeSize/2.0)
					v12.z -= volumeSize*v12.z/fabs(v12.z);
				v12 = (1.0/sqrt(v12*v12))*v12;
				double l = v12.x;
				double m = v12.y;
				double n = v12.z;



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



				currentNeighbour = currentNeighbour->next;
			}



			for(long i=0; i<5; i++) {
				bool ifTheFirst = false;
				for(long j=0; j<10*N; j++) {
#ifdef UCSR
					if ( !((Htmp[i*10*N + j].real() == 0) && (Htmp[i*10*N + j].imag() == 0)) && (j >= i+5*(currentPsi->d->n-1)+5*N) ) {
#endif
#ifdef CSR
					if ( !((Htmp[i*10*N + j].real() == 0) && (Htmp[i*10*N + j].imag() == 0))                                        ) {
#endif
						CSRval[currentCSRval] = Htmp[i*10*N + j];
						if (!ifTheFirst) {
							CSRrows[i+5*(currentPsi->d->n-1)+5*N] = currentCSRval+1;
							ifTheFirst = true;
						}
						CSRcols[currentCSRval] = j+1;
						currentCSRval++;
					}
				}
			}

#ifdef PRINT_H
			// Hamiltonian printer
			for(int i=0; i<5; i++) {
				for(int j=0; j<10*N; j++) {
					fprintf(out, "(");
					if (Htmp[i*10*N+j].real()>=0) fprintf(out, "+"); else fprintf(out, "-");
					fprintf(out, "%3.3f ", abs(Htmp[i*10*N+j].real()));
					if (Htmp[i*10*N+j].imag()>=0) fprintf(out, "+"); else fprintf(out, "-");
					fprintf(out, "%3.3f)\t", abs(Htmp[i*10*N+j].imag()));
				}
				fprintf(out, "\n");
			}
#endif

			currentPsi = currentPsi->next;
		}
		#ifdef UI
		cout << endl;
		#endif
*/
////////////////////////////////////////////////////////////////////////////////////////////////////
		#ifdef UI
		cout << "writing CSR into file" << endl;
		#endif
		if (currentCSRval != CSRsize)
			cerr << "error: currentCSRval == " << currentCSRval << " != CSRsize" << endl;
		fprintf(out, "%i ", 10*N);
		fprintf(out, "%i\n", currentCSRval);
		for(long i=0; i<CSRsize; i++)
			fprintf(out, "%f ", CSRval[i].real());
		fprintf(out, "\n");
		for(long i=0; i<CSRsize; i++)
			fprintf(out, "%f ", CSRval[i].imag());
		fprintf(out, "\n");
		for(long i=0; i<10*N; i++)
			fprintf(out, "%i ", CSRrows[i]);
		fprintf(out, "\n");
		for(long i=0; i<CSRsize; i++)
			fprintf(out, "%i ", CSRcols[i]);
		fclose(out);


		#ifdef UI
		cout << "finished, freing memory" << endl;
		#endif
		free(Htmp);
		free(CSRval);
		free(CSRrows);
		free(CSRcols);

//	}          //k forcircle

	return 0;
}

