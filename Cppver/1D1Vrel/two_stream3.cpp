#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <random>

int N_t       = 3000; // Number of time steps
int N_g       =   64; // Number of grid cells
int nppc      =   200; // Number of particles per cell
int rho0      =    1; // Background charge density
int ns        =    2; // Number of species
int gamma_0   =    6; // Lorentz factor

double v0     =   sqrt(1.0 - 1.0/double(gamma_0*gamma_0)); // Drift velocity
double v_th   =   0.1; // Thermal velocity
double omegap =   1.0; // Plasma frequency
double qm[2]  = {- 1.0, 1.0}; // Charge to mass ratio electrons and positrons
double L      = 200.0; //sqrt(2)*M_PI*v_0*(gamma_0)**(3/2)/omegap; // Domain lenght

using namespace std;

//****************************************************************************80
//
int main()
//
//****************************************************************************80
//
 /*
  * Last update:   09.12.2018
  * Code name:     two_stream3.cpp
  * Licensing:     GNU LGPL license.
  * Author:        Luca Pezzini
  * e-mail:        luca.pezzini@kuleuven.be
  * 
  * Physics:
  * From the theoretical point of view, we want a code that solves the
  * relativistic Vlasov-Maxwell system of equations, where we assume a
  * one-dimensional system and one degree of freedom for the particle motion.
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
	
	string fname0   = "energyc2.5.dat";
	ofstream fdata0 (fname0.c_str(), ios::out);
	fdata0 << setiosflags(ios::scientific);
	
	//****************************************************************************80
	/*                                   START                                    */
	//****************************************************************************80
	
	cout << "\n";
	cout << " ------------------------------------------------------ \n";
  cout << "+                   TWO_STREAM3.CPP                    +\n";
	cout << "+                                                      +\n";
  cout << "+                    C++ VERSION                       +\n";
  cout << "+      1-D TWO STREAM INSTABILITY SIMULATION           +\n";
	cout << "+                 RELATIVISTIC VERSION                 +\n";
  cout << "+                    PiC APPROACH                      +\n";
	cout << "+                                                      +\n";
	cout << " ------------------------------------------------------ \n";
	cout << " \n";
	
	//****************************************************************************80
	/*                            OUTPUT PARAMETER                                */
	//****************************************************************************80
	
	cout << "Domain lenght:                     " << L << endl;
	cout << "Number of time steps:              " << N_t << endl;
	cout << "Number of grid cells:              " << N_g << endl;
	cout << "Number of particles per cell:      " << nppc << endl;
	cout << "Background charge density:         " << rho0 << endl;
	cout << "Number of species:                 " << ns << endl;
	cout << "Drift velocity:                    " << v0 << endl;
	cout << "Thermal velocity:                  " << v_th << endl;
	cout << "Lorentz factor:                    " << gamma_0 << endl;
	cout << "Plasma frequency:                  " << omegap << endl;
  
	//****************************************************************************80
	/*                             GRID PARAMETER                                 */
	//****************************************************************************80

	double x_g[N_g + 1]; // Nodes position array 
	double deltax = L/64.0; // Cell size
	double E_g[N_g + 1]; // Electric field array
	double J_g[N_g + 1]; // Current array

	for (int i = 0; i < N_g + 1; i ++) // Initialise the node position array
	{
	  x_g[i] = i*deltax;
	}
	
	for (int i = 0; i < N_g + 1; i ++) // Initialise E and J grid array as 0
	{
		E_g[i] = 0.0;
		J_g[i] = 0.0;
	}
	
	//****************************************************************************80
	/*                           PARTICLES PARAMETER                              */
	//****************************************************************************80
	
	int N_p = N_g*nppc; // Particles position array size
	double deltax_p = deltax/(double)nppc; // Particle spatial step
	// ns = 0 electrons; ns = 1 positrons
	double x_p[ns][N_p]; // Particles position array
	double v_p[ns][N_p]; // Particles velocity array
	double qmac[ns]; //Charge of macroparticles
	double m_p[ns]; // Masses of the particles (abs. value)
	double step; // Temporary variable for initialize part. pos. array
	double u_p[ns][N_p]; // Particles momentum
	double u0 = gamma_0*v0; // Relativistic drift momentum 
	
	for (int i= 0; i< ns; i++)
	{
	 qmac[i] = ((rho0*deltax)/(double(nppc)))*(fabs(qm[i])/qm[i]);
	 m_p[i] = fabs(qmac[i]/qm[i]);
	}
	
	// Initialise the particle position array of each species:
	// the two species are overlap
	for (int i = 0; i < ns; i++) // Loop on the species of particles
	{
		step = 0.0; 
	  for (int j = 0; j < N_p; j++) // Loop on the particles position
	  {
	    x_p[i][j] = deltax_p/2.0 + step;
		  step += deltax_p;
	  }
	}
	
	// Construct a trivial random generator engine from a time-based seed:
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> distribution (0.0,1.0);
		
	// Initialize particles momenta:
	for (int i = 0; i < ns; i++)
	{
		for (int j = 0; j < N_p; j++)
		{
			v_p[i][j] = v_th*distribution(generator);
			if (j % 2 == 0) 
			{
				u_p[i][j] = 1.0/(sqrt(1 - v_p[i][j]*v_p[i][j]))*v_p[i][j] + u0;
			}
			else 
			{
				u_p[i][j] = 1.0/(sqrt(1 - v_p[i][j]*v_p[i][j]))*v_p[i][j] - u0;
			}
		}
	}
	
  //****************************************************************************80
	/*                                ENERGY & TIME                               */
	//****************************************************************************80
	
	double deltat = deltax/4.0; // Usually -> (deltat > v_th*deltax) we don'twant particles skip cells
	double t[N_t]; // Time array
	double E_p; // Electric field at particles position
	double E_e[N_t +1]; // Electric enecrgy
	double E_k[N_t + 1]; // Kinetic energy
	double gamma_max[N_t + 1]; // Gamma max vector
	
	// Memorization symulation parameter:
	//fdata00 << L << " " << N_t << " "<< N_g << " " << nppc
	             //<< " " << rho0 << " " << ns << " "<< v0
							 //<< " " << v_th << " " << gamma_0
							 //<< " " << omegap << " " << deltat << endl;
	//fdata00 << endl << endl;
	
	// Initialise the time step array:
	for (int i = 0; i < N_t + 1; i++) 
	{
		t[i] = i*deltat;
	}
	
	// Initialise the energy component anywhere as 0:
	for (int i = 0; i < N_t + 1; i++)
	{
		E_e[i] = 0.0;
		E_k[i] = 0.0;
		gamma_max[i] = 0.0;
	}
	
	// Calculate energy at t = 0:
	for (int l = 0; l < ns; l++) // Species loop
	{
		for (int j = 0; j < N_p; j++) // Particles loop  
		{
			E_k[0] += m_p[l]*(sqrt(1.0 + u_p[l][j]*u_p[l][j]) - 1.0); // Kinetic energy
			
			// Select the maximum velocity:
			if (sqrt(1.0 + u_p[l][j]*u_p[l][j]) > gamma_max[0]) 
			{
				gamma_max[0] = sqrt(1.0 + u_p[l][j]*u_p[l][j]);
			}
		}
	}
			
	fdata0 << t[0] << " " << E_e[0] << " "<< E_k[0] << " " << gamma_max[0] << endl;

	/*
	// Fase space at t = 0:
	fname1[0] = "fase" + to_string(1) + ".dat"; // Change the name each time
	ofstream fdata1 (fname1[0].c_str(), ios::out);
	fdata1 << setiosflags(ios::scientific);
	

	for (int l = 0; l < ns; l++) // Species loop
	{
			for (int j = 0; j < N_p; j++) // Particles loop  
		  {
				fdata1 << x_p[l][j] << " " << u_p[l][j]/sqrt(1.0 + u_p[l][j]*u_p[l][j]) << endl;
			}
	}
	fdata1 << endl << endl;
	*/
	
	//****************************************************************************80
	/*                            FOUR-STEP PiC CYCLE                             */
	//****************************************************************************80
	
	for (int i = 1; i < N_t + 1; i++) // Time loop
	{
		 //fname1[i] = "fase" + to_string(i + 1) + ".dat"; // Change the name each time
		 //ofstream fdata1 (fname1[i].c_str(), ios::out);
		 //fdata1 << setiosflags(ios::scientific);
		 
		for (int l = 0; l < ns; l++) // Species loop
		{
			for (int j = 0; j < N_p; j++) // Particles loop  
		  {
				// For each particle, identify the two surrounding nodes k and k + 1:
			  int k = floor(x_p[l][j]/deltax); // Cut the first decimal number after the comma (positive number)
				// The E field on a particle located in a specific cell at position x_p is determined by the E field of the two surrounding nodes
		    E_p = E_g[k]*(1.0 - (x_p[l][j] - x_g[k])/deltax) + E_g[k + 1]*((x_p[l][j] - x_g[k])/deltax); 
			  x_p[l][j] += (u_p[l][j]*deltat)/sqrt(1.0 + u_p[l][j]*u_p[l][j]); // Evolve the particle position
				u_p[l][j] += qm[l]*E_p*deltat; // Evolve the particles velocity
				E_k[i] += m_p[l]*(sqrt(1.0 + u_p[l][j]*u_p[l][j]) - 1.0); // Calculate kinetic energy
				
				// Select the maximum velocity:
				if (sqrt(1.0 + u_p[l][j]*u_p[l][j]) > gamma_max[i]) 
				{
					gamma_max[i] = sqrt(1.0 + u_p[l][j]*u_p[l][j]);
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
				
				// Fase space memorization:
				//fdata1 << x_p[l][j] << " " << u_p[l][j]/sqrt(1.0 + u_p[l][j]*u_p[l][j]) << endl;
		  }
		}
		
		//fdata1 << endl << endl;
		
		// Initialise J_g everywhere:
		for (int j = 0; j < N_g + 1; j++)
		{
			J_g[j] = 0.0;
		}
		
		// Each node will receive information from an unknown number of particles not farther
		// than a distance deltax:
		for (int l = 0; l < ns; l++)
		{
		  for (int j = 0; j < N_p; j++)
		  {
			  int k = floor(x_p[l][j]/deltax);
			  J_g[k] += (qmac[l]*u_p[l][j]*(1.0 - (x_p[l][j] - x_g[k])/deltax))/(deltax*sqrt(1.0 + u_p[l][j]*u_p[l][j]));
			  J_g[k + 1] += (qmac[l]*u_p[l][j]*(x_p[l][j] - x_g[k]))/(deltax*deltax*sqrt(1.0 + u_p[l][j]*u_p[l][j]));
		  }
		}
		
		// The last node on the right of the domain actually corresponds to the 1st node on
		// the left and vice versa:
		J_g[0] += J_g[N_g]; // Adding the contribution of the last node
		J_g[N_g] = J_g[0]; // Imposing the continuity (ring)
		
		// Evolve the E field:
		for (int j = 0; j < N_g; j++) // N_g, instead of N_g + 1, because the 1st and the last are the same
		{
		  E_g[j] -= deltat*J_g[j];
		}
		
		// Summation over the grid of the Electric energy for each time:
		for (int j = 0; j < N_g; j++) 
		{
			E_e[i] += (E_g[j]*E_g[j]*deltax)/2.0;
		}
		
		//Energy memorization:
		fdata0 << t[i] << " " << E_e[i] << " "<< E_k[i] << " " << gamma_max[i] << " " << E_e[i] + E_k[i] <<  endl;	
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