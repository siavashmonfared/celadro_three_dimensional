#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <stdlib.h>
#include <set>
#include <iterator>
#include <numeric> 
#include <chrono>
#include <random>
#include <complex>
#include <algorithm>

using namespace std;

int seed, mcsteps;
double d_max_mc;

double* xcf;
double* ycf;
double* zcf;
double* accepted_moves;






double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);

}

void write_simCard(const string& config, double gamma, double omegacc, double omegacw, double alpha, int nsteps, int ninfo, double Lx, double Ly, double Lz, int nsubsteps, int bc, double margin, int relax_time, int nphases, double mu, double lambda, double kappa, double rad, double xi, double wallThich, double wallKappa, double SPOL, double DPOL, double JPOL,double KPOL, double zetaS, double zetaQ,double SNEM,double KNEM,double JNEM,double WNEM, int count){
        
    FILE * sortie;
    char nomfic[256];
    sprintf(nomfic, "simCard_%d.dat", count);
    
    sortie = fopen(nomfic, "w");
    fprintf(sortie, "# Sample runcard\n");
    fprintf(sortie, "config = %s\n", config.c_str() );
	fprintf(sortie, "nsteps = %d\n", nsteps);
	fprintf(sortie, "ninfo = %d\n", ninfo);
	fprintf(sortie, "LX = %g\n", Lx);
	fprintf(sortie, "LY = %g\n", Ly);
	fprintf(sortie, "LZ = %g\n", Lz);
	fprintf(sortie, "nsubsteps = %d\n", nsubsteps);
	fprintf(sortie, "bc = %d\n", bc);
	fprintf(sortie, "margin = %g\n", margin);
	fprintf(sortie, "relax-time = %d\n", relax_time);
	fprintf(sortie, "nphases = %d\n", nphases);
	fprintf(sortie, "gamma = %g\n", gamma);
	fprintf(sortie, "mu = %g\n", mu);
	fprintf(sortie, "lambda = %g\n", lambda);
	fprintf(sortie, "kappa_cc = %g\n", kappa);
	fprintf(sortie, "R = %g\n", rad);
	fprintf(sortie, "xi = %g\n", xi);
	fprintf(sortie, "omega_cc = %g\n", omegacc);
	fprintf(sortie, "wall-thickness = %g\n", wallThich);
	fprintf(sortie, "kappa_cs = %g\n", wallKappa);
	fprintf(sortie, "omega_cs = %g\n", omegacw);
	fprintf(sortie, "alpha = %g\n", alpha);
	fprintf(sortie, "S-pol = %g\n", SPOL);
	fprintf(sortie, "D-pol = %g\n", DPOL);
	fprintf(sortie, "J-pol = %g\n", JPOL);
	fprintf(sortie, "K-pol = %g\n", KPOL);
	fprintf(sortie, "zetaS = %g\n", zetaS);
	fprintf(sortie, "zetaQ = %g\n", zetaQ);
	fprintf(sortie, "S-nem = %g\n", SNEM);
	fprintf(sortie, "K-nem = %g\n", KNEM);
	fprintf(sortie, "J-nem = %g\n", JNEM);
	fprintf(sortie, "W-nem = %g\n", WNEM);
		
    fclose(sortie);
     
}

void disorder_initial(double l, int np, double xcf[], double ycf[], double zcf[], int nx, int ny, int nz){
	
     double xc, yc, zc;
	int num;

	num = 0;
     for (int z = 0 ; z < nz ; z++){
     for (int y = 0 ; y < ny ; y++){
	for (int x = 0 ; x < nx ; x++){
	
     xc = round(l * x + l/2);
     yc = round(l * y + l/2);
	zc = round(l * z + l/2);
	xcf[num] = xc;
	ycf[num] = yc;
	zcf[num] = zc;
	
	num++;
	}
	}
	}
}

/*
void disorder_mc(size_t d_max, size_t rnd, double xcf[], double ycf[], double zcf[], int npart){
	
	double nxc,nyc,nzc,dmin,p1,p2,p3;
	int exp1, exp2, exp3;
	bool xlub, xllb, ylub, yllb, zlub, zllb, dcond;
	bool xgub, xglb, ygub, yglb, zgub, zglb;

	for (int i = 0 ; i < npart ; i++){

	exp1 = rand() %2;
	exp2 = rand() %2;
	exp3 = rand() %2;
		
	p1 = fRand(0,d_max);
	p2 = fRand(0,d_max);
	p3 = fRand(0,d_max);
		
    	xcf[i] = xcf[i] + pow(-1,exp1) * p1;
    	ycf[i] = ycf[i] + pow(-1,exp2) * p2;
	zcf[i] = zcf[i] + pow(-1,exp3) * p3;
		
	}
}
*/

void init_square_lattice(double l, int np, double xcf[], double ycf[], double zcf[], int nx, int ny, int nz){
	
     double xc, yc, zc;
	int num;

	num = 0;
     for (int z = 0 ; z < nz ; z++){
     for (int y = 0 ; y < ny ; y++){
	for (int x = 0 ; x < nx ; x++){
	
     xc = round(l * x + l/2);
     yc = round(l * y + l/2);
	zc = round(l * z + l/2);
	xcf[num] = xc;
	ycf[num] = yc;
	zcf[num] = zc;
	
	num++;
	}
	}
	}
}


void init_triangular_lattice(int nx, int ny, double domainWidth, double domainHeight, double zcoor, double xcf[], double ycf[], double zcf[], double R0) {
  //  double dx = domainWidth / (nx - 1); // Horizontal spacing
  //  double dy = domainHeight / (ny - 1); // Vertical spacing
    double dx = domainWidth / (nx); // Horizontal spacing
    double dy = domainHeight / (ny); // Vertical spacing
    
    double dx2 = dx * 0.75; // Horizontal spacing between alternate rows
    int num = 0;
    for (int j = 0; j < ny; j++) {
        for (int i = 0; i < nx; i++) {
            double x = i * dx;
            if (j % 2 == 1) {
                x += dx2;
            }
            double y = j * dy;
            xcf[num] = x;
            ycf[num] = y + R0;
            zcf[num] = zcoor;
            num++;
        }
    }
    cout<<"num :"<<num<<endl;
}

void disorder_mc(size_t d_max, size_t rnd, double xcf[], double ycf[], double zcf[], int npart){
	
	double nxc,nyc,nzc,dmin,p1,p2,p3;
	int exp1, exp2, exp3;
	bool xlub, xllb, ylub, yllb, zlub, zllb, dcond;
	bool xgub, xglb, ygub, yglb, zgub, zglb;

	for (int i = 0 ; i < npart ; i++){

	exp1 = rand() %2;
	exp2 = rand() %2;
	exp3 = rand() %2;
		
	p1 = fRand(0,d_max);
	p2 = fRand(0,d_max);
	p3 = fRand(0,d_max);
		
    	xcf[i] = xcf[i] + pow(-1,exp1) * p1;
    	ycf[i] = ycf[i] + pow(-1,exp2) * p2;
	zcf[i] = zcf[i] + pow(-1,exp3) * p3;
		
	}
}



void write_lattice(const string& _name,int np, double xcf[], double ycf[], double zcf[] , double zcoor){

	double xc,yc,zc,r;
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");
    	
	for(int j = 0 ; j < np ; j++){

    	xc = xcf[j];
	yc = ycf[j];
	zc = zcf[j];

	fprintf(sortie,"%g %g %g\n", xc,yc,zcoor);

	}
 	fclose(sortie);
}
 	

double find_min(vector<double> &vect){

	double small = vect[0];
	for (size_t i = 0 ; i < vect.size() ; i++){
	if(vect[i] < small) small = vect[i];
	}
	return small;
}



double compute_all_dist(int index, double ix, double iy, double iz, double txc[], double tyc[], double tzc[], int npart){

	vector<double> dist;

	for (int j = 0 ; j < npart ; j++){
	if(index != j){
	dist.push_back( sqrt(pow( ( txc[j]-ix ),2) + pow( ( tyc[j]-iy ),2) + pow( ( tzc[j]-iz ),2) ) );
	}
	}

	double dmin = find_min(dist);
	dist.clear();
	return dmin;
}















void write_posfile_mix_perc(int np, double xcf[], double ycf[], double zcf[] , double zcoor, int count){

	double xc,yc,zc,r;
	FILE * sortie;
	//char nomfic[256];
	//sprintf(nomfic, "input_str_%d.dat", count);
	sortie = fopen("input_str.dat","w+");
    	//sortie = fopen(nomfic, "w+");
    	
	for(int j = 0 ; j < np ; j++){

    	xc = xcf[j];
	yc = ycf[j];
	zc = zcf[j];

	fprintf(sortie,"%g %g %g\n",xc,yc,zcoor);
	}
	
 	fclose(sortie);
}


void write_summary( int count, double gam1, double gam2, double zetas1, double zetas2, double omegacc1,double omegacc2, double omegacw1,double omegacw2, double alpha, double xi, double zetaQ){

 	FILE * sortie; 
	sortie = fopen("simulation_parameter_summary.dat","a");
	fprintf(sortie,"%i %g %g %g %g %g %g %g %g %g %g %g\n",count,gam1,gam2,zetas1,zetas2,omegacc1,omegacc2,omegacw1,omegacw2,alpha,xi,zetaQ);
	fclose(sortie);
}





void  Export(const string& _name, double tmp[], int tot){
	
   const char * c = _name.c_str();
   FILE * sortie;
   sortie = fopen(c, "w+");	
   for (int k = 0 ; k < (tot) ; k++) fprintf(sortie,"%g\n",tmp[k]);
}





int main (){
	
	int nsteps = 1000;
	int ninfo = 10;
	//double Lx = 64.;
	//double Ly = 64.;
	double Lz = 40.;

	double ncells = 9.;
	double R0 = 8;
	double rad = R0;
	double lbox = 2. * R0;
	int nsqrt = sqrt(ncells);
	double Lx = nsqrt * lbox;
	double Ly = nsqrt * lbox;

	int nsubsteps = 10;
	int bc = 2;
	double margin = 18.;
	int relax_time = 100;
	double wall_thickness = 7.0;
	double zcoor = wall_thickness + R0/2.;


	double gamma = 0.008;
	double omegacc = 0.0008;
	double omegacw = 0.002;
	double alpha = 0.05;
	double xi = 1.;
	double zetaS = 0.;
	double zetaQ = 0.;

	double mu = 45;
	double lambda = 3;
	double kappa = 0.5;

	double wall_kappa = 0.15;
	double SPOL = 1.;
	double SNEM = 0.;
	double KNEM = 0.;
	double JNEM = 0.;
	double WNEM = 0.;
	double JPOL = 0.001;
	double KPOL = 0.001;
	double DPOL = 0.001;

	d_max_mc = 1.;
	int npx = nsqrt;
	int npy = nsqrt;
	int npz = 1;
	int nphases = ncells;

	cout<<"number of cells are : "<<nphases<<" with rad: "<<rad<<endl;

	double* xcf;
	double* ycf;
	double* zcf;
	xcf = new double[npx*npy*npz];
	ycf = new double[npx*npy*npz];
	zcf = new double[npx*npy*npz];

	init_square_lattice(lbox,nphases,xcf,ycf,zcf,npx,npy,npz);
	// disorder_mc(d_max_mc,seed,xcf,ycf,zcf,(npx*npy*npz) );
	// write_lattice("triangular_lattice.dat",nphases,xcf,ycf,zcf,zcoor);
	// init_triangular_lattice(nsqrt,nsqrt,Lx,Ly,zcoor,xcf,ycf,zcf,R0);

	string config = "input const";

	#define nZetaS 1
	const double zetaS_mix_a[nZetaS] = { 0. };//wt 
	const double zetaS_mix_b[nZetaS] = { 0. };// ecad 
	#define nGamma 1
	const double gamma_mix_a[nGamma] = { 0.008 };
	const double gamma_mix_b[nGamma] = { 0.008 };
	#define nOmegCC 1
	const double omegacc_mix_a[nOmegCC] = { 0.0008 };// wt 
	const double omegacc_mix_b[nOmegCC] = { 0.0008 };// ecad 
	#define nOmegCW 1
	const double omegacw_mix_a[nOmegCW] = { 0.0020  }; // wt 
	const double omegacw_mix_b[nOmegCW] = { 0.0020  };// ecad 
	#define nZetaQ 1
	const double zetaQ_mix_a[nZetaQ] = { 0. };// wt 
	const double zetaQ_mix_b[nZetaQ] = { 0. };// ecad 

	#define nAlpha 1
	const double alpha_mix[nAlpha] = { alpha };
	#define nXi 1
	const double xi_mix[nXi] = {xi};
	#define nKappa 1
	const double kappa_mix[nKappa] = {kappa};
	#define nMu 1
	const double mu_mix[nMu] = {mu};
	#define nRad 1
	const double rad_mix[nRad] = {rad};

	int count = 1;

	for (int i = 0 ; i < nGamma ; i++){
	for (int j = 0 ; j < nOmegCC ; j++){
	for (int k = 0 ; k < nOmegCW ; k++){
	for (int m = 0 ; m < nAlpha   ; m++){
	for (int ii = 0 ; ii < nZetaS ; ii++){
	for (int jj = 0 ; jj < nZetaQ ; jj++){
	for (int kk = 0 ; kk < nXi ; kk++){
	for (int mm = 0 ; mm < nKappa ; mm++){
	for (int ll = 0 ; ll < nMu ; ll++){
	for (int hh = 0 ; hh < nRad ; hh++){
 
		double GAMMA_A = gamma_mix_a[i];
		double GAMMA_B = gamma_mix_b[i];

		double ZETAS_A = zetaS_mix_a[ii];
		double ZETAS_B = zetaS_mix_b[ii];

		double OMEGACC_A = omegacc_mix_a[j];
		double OMEGACC_B = omegacc_mix_b[j];

		double OMEGACW_A = omegacw_mix_a[k];
		double OMEGACW_B = omegacw_mix_b[k];

		double ZETAQ_A = zetaQ_mix_a[jj];
		double ZETAQ_B = zetaQ_mix_b[jj];

		double ALPHA = alpha_mix[m];
		double XI = xi_mix[kk];
		double KAPPA = kappa_mix[mm];
		double MU = mu_mix[ll];
		double RAD = rad_mix[hh];

		write_posfile_mix_perc(ncells,xcf,ycf,zcf,zcoor,count);

		write_simCard(config,gamma,omegacc,omegacw,alpha,nsteps,ninfo,Lx,Ly,Lz,nsubsteps, bc, margin,relax_time,ncells,mu,lambda,kappa,rad,xi,wall_thickness,wall_kappa,SPOL,DPOL,JPOL,KPOL,zetaS,zetaQ,SNEM,KNEM,JNEM,WNEM,count);

		double ratio_a = OMEGACC_A/OMEGACC_B;
		double ratio_b = OMEGACW_A/OMEGACW_B;
		// write_summary(count,GAMMA_A,GAMMA_B,ZETAS_A,ZETAS_B,OMEGACC_A,OMEGACC_B,OMEGACW_A,OMEGACW_B,ALPHA,ZETAQ_A,ZETAQ_B);
		count++;

	}
	}
	}
	}
	}
	}
	}
	}
	}
	}

cout << "Computation done." << endl;

delete [] xcf;
delete [] ycf;
delete [] zcf;
return 0;

}







