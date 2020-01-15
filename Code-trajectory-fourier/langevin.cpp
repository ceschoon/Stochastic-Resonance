////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This program generates trajectories obeying a langevin equation.    //
//    It is used to describe the phenomenon of stochastic resonance.      //
//                                                                        //
//    This program has been designed for a project in a meteorology       //
//    course at the Université Libre de Bruxelles (PHYS-F450).            //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

/*  Execute this program with the following bash script:

#! /bin/bash

g++ -O3 -fopenmp -o langevin langevin.cpp utilities.cpp -DPLOT -DSIMULATION
time ./langevin > results.txt

 */

#include <random>
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <omp.h>
#include <stdlib.h>
#include "utilities.h"

using namespace std;

double const PI = 3.14159265358979;

////////////////////////////////////////////////////////////////////////////

// norm of fourier coefficients
void power_spectrum(vector<double>  x, double  time_step, 
                    vector<double> &y, double &freq_step)
{
	int N = x.size();
	y.resize(N);
	
	vector<complex<double>> X(N,0);
	vector<complex<double>> Y(N,0);
	
	freq_step = 2*PI/(N*time_step);
	
	for (int j=0; j<N; j++) X[j] = complex<double>(x[j],0);
	
	for (int k=0; k<N; k++) 
	{
		complex<double> I = complex<double>(0,1);
		for (int j=0; j<N; j++) Y[k] += time_step*X[j]*exp(-2*PI/N*j*k*I);
		Y[k] /= sqrt(2*PI);
		y[k] = norm(Y[k])/(2*PI)*pow(freq_step,2) ;
	}
}

////////////////////////////////////////////////////////////////////////////

// Potential
double U_double_well(double x, double lambda)
{
	return -lambda/2*x*x + 1.0/4*x*x*x*x;
}

double dU_double_well(double x, double lambda)
{
	return -lambda*x + x*x*x;
}

double d2U_double_well(double x, double lambda)
{
	return -lambda + 3*x*x;
}


////////////////////////////////////////////////////////////////////////////


int main()
{
	/////////////////////// Random number generator ////////////////////////
	
	random_device true_gen;
	int seed = true_gen();
	cout << endl;
	cout << "seed = " << seed << endl;
	
	// We create the rng engine inside of the parallelised loop
	// otherwise the same random number sequence is used by all threads
	
	uniform_real_distribution<double> dist01(0,1);
	
	////////////////////////////// Parameters //////////////////////////////
	
	// simulation parameters
	double T = 1e7;             // simulation time
	double dt = 1e-2;           // time step
	double nsteps = 50000;     // number of time steps between each recording
	int M = 1;                  // number of trajectories
	
	// noise parameters
	double q2 = 5e-2;          // amplitude squared
	double q = sqrt(q2);
	
	// potential parameters
	double lambda = 1;
	
	// forcing parameters
	double epsilon = 1e-1;      // amplitude
	double period = 1e5;        // periodicity
	double phi = 0;             // initial phase
	
	// Report
	cout << endl;
	cout << "T = " << T << endl;
	cout << "dt = " << dt << endl;
	cout << "M = " << M << endl;
	cout << "q = " << q << endl;
	cout << "lambda = " << lambda << endl;
	cout << "epsilon = " << epsilon << endl;
	cout << "period = " << period << endl;
	cout << "phi = " << phi << endl;
	
	int sys_result = system("mkdir -p data");
	
	//////////////////////////////// Theory ////////////////////////////////
	
	// Potential details (symmetric double well)
	double U_barrier = U_double_well(0,lambda)
	                 - U_double_well(sqrt(lambda),lambda);
	double d2U_sqrt_factor = sqrt(-d2U_double_well(0,lambda)*
	                               d2U_double_well(sqrt(lambda),lambda));
	
	// Transition rate
	double r = 1/(2*PI)*d2U_sqrt_factor*exp(-2*U_barrier/(q*q));
	
	cout << endl;
	cout << "q^2 = " << q*q << endl;
	cout << "resonant frequency: 2*r = " << 2*r << endl;
	cout << "signal frequency: omega = " << 2*PI/period << endl;
	cout << "Cramer's time = " << 1/r << endl;
	cout << "signal period = " << period  << endl;
	
	////////////////////////////// Simulation //////////////////////////////
	// Used to plot trajectory and do the fourier transform
	
	#ifdef SIMULATION
	//#pragma omp parallel for 
	for (int i=0; i<M; i++)
	{
		double t = 0;
		double x = 0;
		
		default_random_engine gen(seed+i);
		
		ofstream dumpFile("data/trajectory_"+to_string(i)+".dat");
		dumpFile << "# t" << " " << "x" << endl;
		
		while(t<T)
		{
			// generate a random number with normal distribution
			
			double ran1 = dist01(gen);
			double ran2 = dist01(gen);
			double ran3 = sqrt(-2*std::log(ran1)) * sin(ran2*2*PI);
			
			// advance by time step
			
			x += -dU_double_well(x,lambda)*dt         // potential
			  + epsilon*cos(2*PI/period*t+phi)*dt     // forcing
			  + sqrt(q*q*dt)*ran3;                    // noise
			t += dt;
			
			// print in file (every nsteps)
			
			if (int(t/(dt*nsteps)) != int((t-dt)/(dt*nsteps)))
				dumpFile << t << " " << x << endl;
		}
	}
	#endif
	
	// forcing
	
	ofstream dataFile("data/forcing.dat");
	
	for (double t=0; t<T; t+=dt*nsteps)
		// amplitude epsilon=0.5 for visibility on graph
		dataFile << t << " " << 0.5*cos(2*PI/period*t+phi) << endl;
	
	dataFile.close();
	
	
	// interesting time range to plot
	int n_time_steps_to_plot = 10*period/(dt*nsteps);
	n_time_steps_to_plot -= 1; // prettier graph
	
	//////////////////// gnuplot script for trajectories ///////////////////
	
	#ifdef PLOT
	for (int i=0; i<M; i++)
	{
		ofstream plotFile("data/plot_trajectory");
		
		plotFile << "set term svg enhanced mouse #size 600,500" << endl;
		plotFile << "set output 'trajectory_"+to_string(i)+".svg'" << endl;
		plotFile << endl;
		plotFile << "set title \"Integretion of the Langevin equation\" font \",20\"" << endl;
		plotFile << "set label \"epsilon = " << epsilon << " q^2 = " << q*q
			     //<< "\" at graph 0.05,0.90 font \",16\"" << endl; // top left
			     << "\" at graph 0.31,1.05 font \",16\"" << endl; // under title
		plotFile << "set xlabel \"t\" font \",20\"" << endl;
		plotFile << "set ylabel \"x\" font \",20\"" << endl;
		plotFile << endl;
		plotFile << "set format x '%.1t*10^{%T}'" << endl;
		plotFile << endl;
		plotFile << "set key off" << endl;
		plotFile << endl;
		
		plotFile << "plot \"forcing.dat\" "
		         << "every ::0::"+to_string(n_time_steps_to_plot)+" "
		         << "with lines, \\" << endl;
		plotFile << "     \"trajectory_"+to_string(i)+".dat\" "
		         << "every ::0::"+to_string(n_time_steps_to_plot)+" "
		         << "with lines" << endl;
		
		plotFile.close();
		
		sys_result = system("cd data; gnuplot plot_trajectory");
	}
	#endif
	
	//////////////////////// Plot of the potential /////////////////////////
	
	double dx = 0.01;
	double x0 = -1.6;
	double x1 =  1.6;
	
	int N = int((x1-x0)/dx) +1;
	vector<double> x(N,0);
	for (int i=0; i<N; i++) x[i] = x0+i*dx;
	vector<double> y(N,0);
	for (int i=0; i<N; i++) y[i] = U_double_well(x[i],lambda);
	
	// save in file
	
	dataFile.open("data/potential.dat");
	
	for (int i=0; i<N; i++)
		dataFile << x[i] << " " << y[i] << endl;
	
	dataFile.close();
	
	// plot using gnuplot
	
	ofstream plotFile;
	plotFile.open("data/plot_potential");
	
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'potential.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Double well potential\" font \",20\"" << endl;
	plotFile << "set label \"lambda = " << lambda
	         //<< "\" at graph 0.5,0.9 font \",16\"" << endl; // top left
	         << "\" at graph 0.42,1.05 font \",16\"" << endl; // under title
	plotFile << "set xlabel \"x\" font \",20\"" << endl;
	plotFile << "set ylabel \"U(x)\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set key off" << endl;
	plotFile << endl;
	
	plotFile << "plot \"potential.dat\" with lines" << endl;
	
	plotFile.close();
	
	sys_result = system("cd data; gnuplot plot_potential");
	
	
	/////////////////////////// Fourier Transform //////////////////////////
	
	#ifdef PLOT
	//#pragma omp parallel for 
	for (int i=0; i<M; i++)
	{
		// input time series x assumed to be a time serie evenly spaced in time
		vector<double> t;
		vector<double> x_raw;
		vector<double> x;
		vector<double> y;
		
		// TEST is success if ||^2 goes to about 0.06
		// as 0.5*cos() in forcing and 0.5 because cos is exp+exp-
		// factor (1/4)^2 = 1/16 = 0.0625
		
		//ifstream dataFile("data/forcing.dat"); //test
		ifstream dataFile("data/trajectory_"+to_string(i)+".dat");
		if (!dataFile) continue;
		readColumnVectorFromFile(dataFile, 0, t);
		readColumnVectorFromFile(dataFile, 1, x);
		
		double time_step = t[1]-t[0];
		double freq_step;
		power_spectrum(x,time_step,y,freq_step);
		
		ofstream dumpFile("data/spectrum_"+to_string(i)+".dat");
		dumpFile << "#freq abs_amplitude" << endl;
		
		int N = x.size();
		for (int k=0; k<N; k++)
		{
			double freq = k*freq_step;
			double amplitude = sqrt(y[k]);
			if (k<N) 
			{
				dumpFile << freq << " " << amplitude << endl;
			} 
		}
		
		// interesting range of frequencies to plot
		double signal_freq = 2*PI/period;
		double resonant_freq = 2*r;
		double max_freq_show = 3 * ((signal_freq>resonant_freq)?signal_freq:resonant_freq);
		double max_freq_DFT = 2*PI/time_step;
		int n_freq_show = N*max_freq_show/max_freq_DFT;
		n_freq_show -= 1; // prettier graph
		
		// plot using gnuplot
		
		ofstream plotFile;
		plotFile.open("data/plot_spectrum");
		
		plotFile << "set term svg enhanced mouse #size 600,500" << endl;
		plotFile << "set output 'spectrum_"+to_string(i)+".svg'" << endl;
		plotFile << endl;
		plotFile << "set title \"Fourier Transform of Trajectory\" font \",20\"" << endl;
		plotFile << "set xlabel \"omega\" font \",20\"" << endl;
		plotFile << "set ylabel \"amplitude of the response |A|\" font \",20\"" << endl;
		plotFile << endl;
		plotFile << "set format x '%.1t*10^{%T}'" << endl;
		plotFile << endl;
		plotFile << "set key off" << endl;
		plotFile << endl;
		
		plotFile << "plot \"spectrum_"+to_string(i)+".dat\" "
		         << "every ::0::"+to_string(n_freq_show)+" "
		         << "with lines" << endl;
		
		plotFile.close();
		
		sys_result = system("cd data; gnuplot plot_spectrum");
	}
	#endif
	
	////////////////////////////////////////////////////////////////////////
	
	
	return 0;
}