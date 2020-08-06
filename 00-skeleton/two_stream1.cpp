#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <chrono>
#include <random>

int N_g    =  32; // number of grid cells
int L      =   1; // domain lenght
int nppc =  10; // number of particles per cell
int N_t    = 100; // symulation time step

using namespace std;

//****************************************************************************80
//
int main()
//
//****************************************************************************80
//
 /*
  * Last update:   04.11.2018
  * Code name:     two_stream.cpp
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

  cout << setiosflags (ios::scientific);
  string fname   = "fase_space.dat";
  ofstream fdata (fname.c_str(), ios::out);
  
	//****************************************************************************80
	/*                             GRID PARAMETER                                 */
	//****************************************************************************80

	double x_g[N_g + 1]; // nodes position array 
	double deltax = L/(double)N_g; // cell size
	double E_g[N_g + 1]; // electric field array
	double J_g[N_g + 1]; // current array

	for (int i = 0; i < N_g + 1; i += 1) // initialise the node position array
	{
	  x_g[i] = i*deltax;
	}
	
	for (int i = 0; i < N_g + 1; i += 1) // initialise E and J grid array as 0
	{
		E_g[i] = 0.0;
		J_g[i] = 0.0;
	}
	
	//****************************************************************************80
	/*                           PARTICLES PARAMETER                              */
	//****************************************************************************80
	
	int N_p = N_g*nppc; // particles position array size
	double deltax_p = deltax/(double)nppc; // particle spatial step
	double step = 0.0; // temporary variable for initialize part. pos. array
	double x_p[N_p]; // particles position array
	double v_p[N_p]; // particles velocity array
	double v_th = 0.001; // thermal velocity
	double charge_to_mass = 1.0; // charge to mass ratio: q_p/m_p = 1
	
	for (int i = 0; i < N_p; i += 1) // initialise the particle position array
	{
	  x_p[i] = deltax_p/2 + step;
		step += deltax_p;
	}
	
	// construct a trivial random generator engine from a time-based seed:
	
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  default_random_engine generator (seed);
  normal_distribution<double> distribution (0.0,1.0);
	
	// initialise the particle velocity array -> Maxwellian (random):
	
	for (int i = 0; i < N_p; i += 1) 
	{
		v_p[i] = v_th*distribution(generator); // random number normalized by v_th
		fdata << x_p[i] << " " << v_p[i] << endl;
	}
	
	/*
	for (int i = 0; i < N_p; i += 1) 
	{
		fdata << x_p[i] << " " << v_p[i] << endl;
	}
	*/
	
  //****************************************************************************80
	/*                              OTHER PARAMETER                               */
	//****************************************************************************80
	
	double deltat = deltax/2; // usually -> (deltat > v_th*deltax) we don'twant particles skip cells
	double t[N_t]; // time array
	int rho0 = 1; // background charge density
	double q_p = rho0*deltax/nppc; // constant number for all particles
	
	for (int i = 0; i < N_t; i += 1) // initialise the time step array
	{
		t[i] = i*deltat;
	}
	
	//****************************************************************************80
	/*                            FOUR-STEP PiC  CYCLE                            */
	//****************************************************************************80
	
	double E_p; // electric field at particles position 
	
	for (int i = 0; i < N_t; i += 1) // time loop
	{
	 	for (int j = 0; j < N_p; j += 1) // particles loop
		{
			int k = floor(x_p[j]/deltax); // cut the first decimal number after the comma (positive number)
		    E_p = E_g[k]*(1 - fabs(x_p[j] - x_g[k])/deltax) + E_g[k + 1]*(1 - fabs(x_p[j] - x_g[k])/deltax); // the E field on a particle located in a specific cell at position x_p is determined by the E field of the two surrounding nodes
			x_p[j] += deltat*v_p[j]; // Evolve the particle position
			v_p[j] += charge_to_mass*deltat*E_p; // Evolve the particles velocity
			
			// Reflection boundary conditions for particles:
			if(x_p[j] > L)
			{
			  x_p[j] -= L;
			}
			
			if (x_p[j] < 0)
			{
			  x_p[j] += L;
			}
		}
		
		// Initialise J_g everywhere
		for (int d = 0; d < N_g + 1; d += 1)
		{
			J_g[d] = 0;
		}
		// Each node will receive information from an unknown number of particles not farther than a distance deltax
		for (int j = 0; j < N_p;  j += 1)
		{
			int k = floor(x_p[j]/deltax);
			J_g[k] += q_p*v_p[j]*(1 - (x_p[j] - x_g[k])/deltax)/deltax;
			J_g[k + 1] += q_p*v_p[j]*(x_p[j] - x_g[k])/(deltax*deltax);
		}
		
		// The last node on the right of the domain actually corresponds to the 1st node on the left and vice versa
		J_g[0] += J_g[N_g];
		J_g[N_g] += J_g[0];
		
		// Evolve the E field
		for (int f; f < N_g; f += 1)
		{
		  E_g[f] -= deltat*J_g[f];
		}
	}
	
	
	for (int i = 0; i < N_p; i += 1) 
	{
		fdata << x_p[i] << " " << v_p[i] << endl;
	}
  
	
	

  fdata << endl << endl;
  fdata.close();

  cout << "\n";
  cout << " ----------------------------------------------------- \n";
  cout << "+                                                     +\n";
  cout << "+               NORMAL END OF EXECUTION !             +\n";
  cout << "+                                                     +\n";
  cout << " ----------------------------------------------------- \n";
  cout << " \n";
	//
  return 0;
}