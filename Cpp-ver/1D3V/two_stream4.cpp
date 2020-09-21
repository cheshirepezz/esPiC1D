#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <random>

int N_t       = 40000; // Number of time steps
int N_g       =   128; // Number of grid cells
int nppc      =    50; // Number of particles per cell
int rho0      =     1; // Background charge density
int ns        =     2; // Number of species

double v0     =   0.01; // Drift velocity
double v_th   =   0.1*v0; // Thermal velocity
double omegap =   1.0; // Plasma frequency
double qm[2]  = {- 1.0, 1.0}; // Charge to mass ratio electrons and positrons
double L      = 10.0;//3.0*sqrt(2.0)*M_PI*v0; // Domain lenght

using namespace std;

//****************************************************************************80
//
int main()
//
//****************************************************************************80
//
 /*
  * Last update:   21.12.2018
  * Code name:     two_stream4.cpp
  * Licensing:     GNU LGPL license.
  * Author:        Luca Pezzini
  * e-mail:        luca.pezzini@kuleuven.be
  * 
  * Physics:
  * From the theoretical point of view, we want a code that solves the
  * Vlasov-Maxwell system of equations, where we assume a one-dimensional
  * system and one degree of freedom for the particle motion.
  * Let us further assume that no magnetic fields are present.
  * The particles move according to the nonrelativistic Newton's equations.
  * A PiC method solves the equations above by employing a finite-difference
  * discretization in time and space.
  * 
  */
{
	//****************************************************************************80
	/*                                DATA STORAGE                                */
	//****************************************************************************80
	
  cout << setiosflags (ios::scientific);
	
	// Create a file .dat :
	/*
	string fname00   = "param.dat";
	ofstream fdata00 (fname00.c_str(), ios::out);
	fdata00 << setiosflags(ios::scientific);
	*/
	
	string fname0   = "energy3.dat";
	ofstream fdata0 (fname0.c_str(), ios::out);
	fdata0 << setiosflags(ios::scientific);
	
	//****************************************************************************80
	/*                                   START                                    */
	//****************************************************************************80
	cout << "\n";
	cout << " ------------------------------------------------------ \n";
  cout << "+                   TWO_STREAM.CPP                     +\n";
	cout << "+                                                      +\n";
  cout << "+                    C++ VERSION                       +\n";
  cout << "+      1-D TWO STREAM INSTABILITY SIMULATION           +\n";
  cout << "+                    PiC APPROACH                      +\n";
	cout << "+                                                      +\n";
	cout << " ------------------------------------------------------ \n";
	cout << " \n";
	
	//****************************************************************************80
	/*                            OUTPUT PARAMETER                                */
	//****************************************************************************80
	
	cout << "Domain lenght:                     " <<      L      << endl;
	cout << "Number of time steps:              " <<      N_t    << endl;
	cout << "Number of grid cells:              " <<      N_g    << endl;
	cout << "Number of particles per cell:      " <<      nppc   << endl;
	cout << "Background charge density:         " <<      rho0   << endl;
	cout << "Number of species:                 " <<      ns     << endl;
	cout << "Drift velocity:                    " <<      v0     << endl;
	cout << "Thermal velocity:                  " <<      v_th   << endl;
	cout << "Plasma frequency:                  " <<      omegap << endl;
  
	//****************************************************************************80
	/*                             GRID PARAMETER                                 */
	//****************************************************************************80

	double xn[N_g + 1]; // Nodes position array
	double xc[N_g]; // Centres position array
	double deltax = L/(double)N_g; // Cell size
	// Electric field array (nodes):
	double Enx[N_g + 1]; 
	double Eny[N_g + 1];
	double Enz[N_g + 1];
	// Current components array (nodes):
	double Jx[N_g + 1];
	double Jy[N_g + 1];
	double Jz[N_g + 1]; 
	// B field components (nodes):
	double Bny[N_g + 1];
	double Bnz[N_g + 1];
	// B field components (centres):
	double Bcy[N_g];
	double Bcz[N_g];

	// Initialise the node position array:
	for (int i = 0; i < N_g + 1; i++) 
	{
	  xn[i] = i*deltax;
	}
	
	// Initialise the node position array:
	for (int i = 1; i < N_g; i++) 
	{
	  xc[i] = i*0.5*deltax;
	}
	
	// Initialise E and J grid components array as 0:
	for (int i = 0; i < N_g + 1; i++)
	{
		Enx[i] = 0.0;
		Eny[i] = 0.0;
		Enz[i] = 0.0;
		Jx[i] = 0.0;
		Jy[i] = 0.0;
		Jz[i] = 0.0;
	}
	
	// Initialise B field components array (centres) as 0:
	for (int i = 0; i < N_g; i++)
	{
		Bcy[i] = 0.0;
		Bcz[i] = 0.0;
	}
	
	//****************************************************************************80
	/*                           PARTICLES PARAMETER                              */
	//****************************************************************************80
	
	int N_p = N_g*nppc; // Particles position array size
	double deltax_p = deltax/(double)nppc; // Particle spatial step
	// ns = 0 electrons; ns = 1 positrons
	double x_p[ns][N_p]; // Particles position array
	double qmac[ns]; //Charge of macroparticles
	double m_p[ns]; // Masses of the particles (abs. value)
	// Particles velocity components array:
	double vx[ns][N_p];
	double vy[ns][N_p];
	double vz[ns][N_p];
	// Particles E field components:
	double Epx = 0.0;
	double Epy = 0.0;
	double Epz = 0.0;
	// Particles B field components:
	double Bpx = 0.0;
	double Bpy = 0.0;
	double Bpz = 0.0;
	
	for (int i = 0; i < ns; i++)
	{
	 qmac[i] = ((rho0*deltax)/(double(nppc)))*(fabs(qm[i])/qm[i]);
	 m_p[i] = fabs(qmac[i]/qm[i]);
	}
	
	// Initialise the particle position array of each species:
	// the two species are overlap
	for (int i = 0; i < ns; i++) // Loop on the species of particles
	{
	  for (int j = 0; j < N_p; j++) // Loop on the particles position
	  {
	    x_p[i][j] = deltax_p*(j + 0.5);
	  }
	}
	
	// Construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> distribution (0.0,1.0);
	
	// Initialize velocities of each species:
	for (int i = 0; i < ns; i++)
	{
		for (int j = 0; j < N_p; j++)
		{
			vx[i][j] = v_th*distribution(generator);
			vz[i][j] = v_th*distribution(generator);
			
			// Filamentation instability:
			if (j % 2 == 0) 
			{
				vy[i][j] = v_th*distribution(generator) + v0;
			}
			else 
			{
				vy[i][j] = v_th*distribution(generator) - v0;
			}
		}
	}
	
  //****************************************************************************80
	/*                                ENERGY & TIME                               */
	//****************************************************************************80
	
	double deltat = L/(4.0*64.0); // Usually -> (deltat > v_th*deltax) we don't want particles skip cells
	double t[N_t + 1]; // Time array
	// Energy:
	double E_e[N_t + 1]; // Electric 
	double E_k[N_t + 1]; // Kinetic 
	double E_b[N_t + 1]; // Magnetic 
	// Maximum velocity by components:
	double vxmax[N_t + 1]; 
	double vymax[N_t + 1];
	double vzmax[N_t + 1];
	
	// Memorization symulation parameter:
	//fdata00 << L << " " << N_t  << " " << N_g << " " << nppc
	             //<< " " << rho0 << " " << ns  << " " << v0
							 //<< " " << v_th << " " << omegap << " " << deltat << endl;
	//fdata00 << endl << endl;
	
	// Initialise the time step array:
	for (int i = 0; i < N_t + 1; i++) 
	{
		t[i] = i*deltat;
	}
	
	// Initialise the energy and velocity component anytime as 0:
	for (int i = 0; i < N_t + 1; i++)
	{
		E_e[i] = 0.0;
		E_k[i] = 0.0;
		E_b[i] = 0.0;
		vxmax[i] = 0.0;
		vymax[i] = 0.0;
		vzmax[i] = 0.0;
	}
	
	// Calculate energy at t = 0:
	for (int l = 0; l < ns; l++) // Species loop
	{
		for (int j = 0; j < N_p; j++) // Particles loop  
		{
			E_k[0] += m_p[l]*(vx[l][j]*vx[l][j] + vy[l][j]*vy[l][j] + vz[l][j]*vz[l][j])/2.0; 
			
			// Select the maximum velocity:
			if (fabs(vx[l][j]) > vxmax[0]) 
			{
				vxmax[0] = vx[l][j];
			}
			if (fabs(vy[l][j]) > vymax[0]) 
			{
				vymax[0] = vy[l][j];
			}
			if (fabs(vz[l][j]) > vzmax[0]) 
			{
				vzmax[0] = vz[l][j];
			}
		}
	}
		
	fdata0 << t[0] << " " << E_e[0] << " "<< E_k[0] << " " << E_b[0] << " " << vxmax[0] << " " << vymax[0] << " " << vzmax[0] << endl;
 /* 
	// Fase space at t = 0:
	fname1[0] = "fase" + to_string(1) + ".dat"; // Change the name each time
	ofstream fdata1 (fname1[0].c_str(), ios::out);
	fdata1 << setiosflags(ios::scientific);
	
	
	for (int l = 0; l < ns; l++) // Species loop
	{
			for (int j = 0; j < N_p; j++) // Particles loop  
		  {
				fdata1 << x_p[l][j] << " " << vx[l][j] << " " << vy[l][j] << " " << vz[l][j] << endl;
			}
	}
	
	fdata1 << endl << endl;
	*/
	
	
	//****************************************************************************80
	/*                            FOUR-STEP PiC CYCLE                             */
	//****************************************************************************80

	//int f = 0;
	
	for (int i = 1; i < N_t + 1; i++) // Time loop
	{
		/*
		if (i % n == 0)
		{
			f += 1;
		  fname1[f] = "fase" + to_string(f + 1) + ".dat"; // Change the name each time
		  ofstream fdata1 (fname1[f].c_str(), ios::out);
		  fdata1 << setiosflags(ios::scientific);
		}
		*/
		/* 
		fname1[i] = "fase" + to_string(i + 1) + ".dat"; // Change the name each time
		ofstream fdata1 (fname1[i].c_str(), ios::out);
		fdata1 << setiosflags(ios::scientific);
		*/
		
		// Calculate B field on the nodes
		for (int a = 1; a < N_g; a++)
	  {
		  Bny[a] = (Bcy[a] + Bcy[a - 1])/2.0 ;
		  Bnz[a] = (Bcz[a] + Bcz[a - 1])/2.0;
		}
				
		// B field boundary conditions:
		Bny[0] = (Bcy[0] + Bcy[N_g - 1])/2.0 ;
		Bnz[0] = (Bcz[0] + Bcz[N_g - 1])/2.0;
		Bny[N_g] = Bny[0];
		Bnz[N_g] = Bnz[0];
		
		for (int l = 0; l < ns; l++) // Species loop
		{
			for (int j = 0; j < N_p; j++) // Particles loop  
		  {
				// For each particle, identify the two surrounding nodes k and k + 1:
			  int k = floor(x_p[l][j]/deltax); // Cut the first decimal number after the comma (positive number)
			  Epx = Enx[k]*(1.0 - (x_p[l][j] - xn[k])/deltax) + Enx[k + 1]*((x_p[l][j] - xn[k])/deltax);
				Epy = Eny[k]*(1.0 - (x_p[l][j] - xn[k])/deltax) + Eny[k + 1]*((x_p[l][j] - xn[k])/deltax);
				Epz = Enz[k]*(1.0 - (x_p[l][j] - xn[k])/deltax) + Enz[k + 1]*((x_p[l][j] - xn[k])/deltax); 
				Bpy = Bny[k]*(1.0 - (x_p[l][j] - xn[k])/deltax) + Bny[k + 1]*((x_p[l][j] - xn[k])/deltax);
				Bpz = Bnz[k]*(1.0 - (x_p[l][j] - xn[k])/deltax) + Bnz[k + 1]*((x_p[l][j] - xn[k])/deltax);
				    
				// Evolve the particle position and velocity:
				x_p[l][j] += vx[l][j]*deltat; 
				vx[l][j] += qm[l]*deltat*(Epx + vy[l][j]*Bpz - vz[l][j]*Bpy);
				vy[l][j] += qm[l]*deltat*(Epy - vx[l][j]*Bpz - vz[l][j]*Bpx);
				vz[l][j] += qm[l]*deltat*(Epz + vx[l][j]*Bpy - vy[l][j]*Bpx);
				
				//  Calculate kinetic energy:
				E_k[i] += m_p[l]*(vx[l][j]*vx[l][j] + vy[l][j]*vy[l][j] + vz[l][j]*vz[l][j])/2.0; 
				
				// Select the maximum velocity:
				if (fabs(vx[l][j]) >= vxmax[i]) 
				{
					vxmax[i] = vx[l][j];
				}
				if (fabs(vy[l][j]) >= vymax[i]) 
				{
					vymax[i] = vy[l][j];
				}
				if (fabs(vz[l][j]) >= vzmax[i]) 
				{
					vzmax[i] = vz[l][j];
				}
				
			  // Reflection boundary conditions for particles:
			  if (x_p[l][j] >= L)
			  {
			    x_p[l][j] -= L;
			  }
			  if (x_p[l][j] < 0.0)
			  {
			    x_p[l][j] += L;
			  }
				/*
				// Fase space memorization:
				//if (i % n == 0)
				//{
				  fdata1 << x_p[l][j] << " " << vx[l][j] << " " << vy[l][j] << " " << vz[l][j] << endl;
				//}
				*/
		  }
		}
		
		//if (i % n == 0)
		//{
		  //fdata1 << endl << endl;
		//}
		
		// Initialise J_g everywhere:
		for (int j = 0; j < N_g + 1; j++)
		{
			Jx[j] = 0.0;
			Jy[j] = 0.0;
			Jz[j] = 0.0;
		}
		
		// Each node will receive information from an unknown number of particles not farther
		// than a distance deltax:
		for (int l = 0; l < ns; l++)
		{
		  for (int j = 0; j < N_p; j++)
		  {
			  int k = floor(x_p[l][j]/deltax);
			  Jx[k] += qmac[l]*vx[l][j]*(1.0 - (x_p[l][j] - xn[k])/deltax)/deltax;
			  Jx[k + 1] += (qmac[l]*vx[l][j]*(x_p[l][j] - xn[k]))/(deltax*deltax);
				Jy[k] += qmac[l]*vy[l][j]*(1.0 - (x_p[l][j] - xn[k])/deltax)/deltax;
			  Jy[k + 1] += (qmac[l]*vy[l][j]*(x_p[l][j] - xn[k]))/(deltax*deltax);
				Jz[k] += qmac[l]*vz[l][j]*(1.0 - (x_p[l][j] - xn[k])/deltax)/deltax;
			  Jz[k + 1] += (qmac[l]*vz[l][j]*(x_p[l][j] - xn[k]))/(deltax*deltax);
		  }
		}
		
		// Current boundary conditions:
		Jx[0] += Jx[N_g]; // Adding the contribution of the last node
		Jx[N_g] = Jx[0]; // Imposing the continuity (ring)
		Jy[0] += Jy[N_g]; 
		Jy[N_g] = Jy[0];
		Jz[0] += Jz[N_g]; 
		Jz[N_g] = Jz[0];
		
		// Evolve the E field:
		for (int j = 1; j < N_g + 1; j++)
		{
		  Enx[j] -= deltat*Jx[j];
			Eny[j] -= deltat*((Bcz[j] - Bcz[j - 1])/deltax + Jy[j]);
			Enz[j] += deltat*((Bcy[j] - Bcy[j - 1])/deltax - Jz[j]);
		}
		
		// E (y,z) component boundary conditions:
		Enx[0] -= deltat*Jx[0];
		Eny[0] -= deltat*((Bcz[0] - Bcz[N_g - 1])/deltax + Jy[0]);
		Enz[0] += deltat*((Bcy[0] - Bcy[N_g - 1])/deltax - Jz[0]);
		Enx[N_g] = Enx[0];
		Eny[N_g] = Eny[0];
		Enz[N_g] = Enz[0];
		
		// Evolving the B field in the centres:
		for (int j = 0; j < N_g; j++)
		{
		  Bcy[j] += deltat*(Enz[j + 1] - Enz[j])/deltax;
			Bcz[j] -= deltat*(Eny[j + 1] - Eny[j])/deltax;
		}
		
		// Summation over the grid of the Electric energy for each time:
		for (int j = 0; j < N_g; j++) 
		{
			E_e[i] += deltax*(Enx[j]*Enx[j] + Eny[j]*Eny[j] + Enz[j]*Enz[j])/2.0;
			E_b[i] += deltax*(Bny[j]*Bny[j] + Bnz[j]*Bnz[j])/2.0;
		}
		//Energy memorization:
		fdata0 << t[i] << " " << E_e[i] << " "<< E_k[i] << " " << E_b[i] << " " << vxmax[i] << " " << vymax[i] << " " << vzmax[i] << endl;
	}
	
	fdata0 << endl << endl;
	
	//****************************************************************************80
	/*                                   END !                                    */
	//****************************************************************************80
	
	//fdata00.close();
	fdata0.close();
	//fdata1.close();

  cout << "\n";
	cout << " ----------------------------------------------------- \n";
	cout << "+                                                     +\n";
  cout << "+               NORMAL END OF EXECUTION !             +\n";
  cout << "+                                                     +\n";
	cout << " ----------------------------------------------------- \n";
	cout << " \n";

  return 0;
}