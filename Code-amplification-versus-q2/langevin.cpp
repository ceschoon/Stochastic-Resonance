////////////////////////////////////////////////////////////////////////////
//                                                                        //
//    This program computes the amplification of a signal for several     //
//    noise strenghs by stochastic resonance.                             //
//                                                                        //
//    This program has been designed for a project in a meteorology       //
//    course at the Université Libre de Bruxelles (PHYS-F450).            //
//                                                                        //
//    Author: Cedric Schoonen                                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

/*  Execute this program with the following bash script:

#! /bin/bash

g++ -O3 -fopenmp -o langevin langevin.cpp utilities.cpp -DSIMULATION
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

// returns norm of the fourier coefficient at freq
double power_at_freq(vector<double> x, double time_step, double freq)
{
	int N = x.size();
	complex<double> coef_freq = 0;
	
	vector<complex<double>> X(N,0);
	vector<complex<double>> Y(N,0);
	
	complex<double> I = complex<double>(0,1);
	for (int j=0; j<N; j++) X[j] = complex<double>(x[j],0);
	
	for (int j=0; j<N; j++) 
		coef_freq += time_step*X[j]*exp(-freq*j*time_step*I);
	coef_freq /= sqrt(2*PI);
	
	double freq_step = 2*PI/(N*time_step);
	
	return norm(coef_freq)/(2*PI)*pow(freq_step,2);
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
	double T = 1e8;             // simulation time
	double dt = 1e-2;           // time step
	double nsteps = 100000;     // number of time steps between each recording
	
	// noise parameters
	int M = 20;                 // number of noise strenghs to try
	double q2_min = 0;          // min noise strengh
	double q2_max = 0.2;        // max noise strengh
	
	vector<double> q2_vec(M,0);
	for (int i=0; i<M; i++) q2_vec[i] = q2_min+(q2_max-q2_min)/M*i;
	
	// potential parameters
	double lambda = 1;
	
	// forcing parameters
	double epsilon = 5e-3;      // amplitude
	double period = 1e5;        // periodicity
	double phi = 0;             // initial phase
	
	// Report
	cout << endl;
	cout << "T = " << T << endl;
	cout << "dt = " << dt << endl;
	cout << "M = " << M << endl;
	cout << "q2_min = " << q2_min << endl;
	cout << "q2_max = " << q2_max << endl;
	cout << "lambda = " << lambda << endl;
	cout << "epsilon = " << epsilon << endl;
	cout << "period = " << period << endl;
	cout << "phi = " << phi << endl;
	
	int sys_result = system("mkdir -p data");
	
	//////////////////////////////// Theory ////////////////////////////////
	
	double q = sqrt((q2_min+q2_max)/2);
	
	// Potential details (symmetric double well)
	double U_barrier = U_double_well(0,lambda)
	                 - U_double_well(sqrt(lambda),lambda);
	double d2U_sqrt_factor = sqrt(-d2U_double_well(0,lambda)*
	                               d2U_double_well(sqrt(lambda),lambda));
	
	// Transition rate
	double r = 1/(2*PI)*d2U_sqrt_factor*exp(-2*U_barrier/(q*q));
	
	cout << endl;
	cout << "Properties of mid-range noise, q^2 = (q2_min+q2_max)/2 :" << endl;
	cout << "q^2 = " << q*q << endl;
	cout << "resonant frequency: 2*r = " << 2*r << endl;
	cout << "signal frequency: omega = " << 2*PI/period << endl;
	cout << "Cramer's time = " << 1/r << endl;
	cout << "signal period = " << period  << endl;
	
	// Amplitude versus noise (symmetric double well)
	
	ofstream dataFile("data/amplification_versus_noise_theory.dat");
	dataFile << "#q2 abs_amplitude" << endl;
	
	int M_theory = 1000;
	vector<double> amplitudes_theory(M_theory,0);
	for (int i=0; i<M_theory; i++)
	{
		double q2 = q2_min + (q2_max-q2_min)/M_theory*i;
		double ri = 1/(2*PI)*d2U_sqrt_factor*exp(-2*U_barrier/q2);
		amplitudes_theory[i] = epsilon*(-sqrt(lambda)/q2)
		                     / sqrt(1+pow(2*PI/period,2)/pow(2*ri,2));
		dataFile << q2 << " " << abs( amplitudes_theory[i] ) << endl;
	}
	
	////////////////////////////// Simulation //////////////////////////////
	// Used to plot trajectory and do the fourier transform
	
	#ifdef SIMULATION
	#pragma omp parallel for 
	for (int i=0; i<M; i++)
	{
		//noise strengh sqrt
		double q = sqrt(q2_vec[i]);
		
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
	
	///////////////////// Amplitudes of signal frequency ///////////////////
	
	ofstream dumpFile("data/amplification_versus_noise.dat");
	dumpFile << "#q2 abs_amplitude" << endl;
	
	for (int i=0; i<M; i++)
	{
		// input time series x assumed to be a time serie evenly spaced in time
		vector<double> t;
		vector<double> x;
		vector<double> y;
		
		ifstream dataFile("data/trajectory_"+to_string(i)+".dat");
		if (!dataFile) continue;
		readColumnVectorFromFile(dataFile, 0, t);
		readColumnVectorFromFile(dataFile, 1, x);
		
		double signal_freq = 2*PI/period;
		double time_step = t[1]-t[0];
		
		double amplitude2 = power_at_freq(x,time_step,signal_freq);
		dumpFile << q2_vec[i] << " " << sqrt(amplitude2) << endl;
	}
	
	
	/////////////////// Plot of amplification versus noise /////////////////
	
	ofstream plotFile("data/plot_amplification_versus_noise");
	
	plotFile << "set term svg enhanced mouse #size 600,500" << endl;
	plotFile << "set output 'plot_amplification_versus_noise.svg'" << endl;
	plotFile << endl;
	plotFile << "set title \"Amplification versus noise\" font \",20\"" << endl;
	plotFile << "set label \"epsilon = " << epsilon << ", omega = " << 2*PI/period
		     //<< "\" at graph 0.05,0.90 font \",16\"" << endl; // top left
		     << "\" at graph 0.22,1.05 font \",16\"" << endl; // under title
	plotFile << "set xlabel \"q^2\" font \",20\"" << endl;
	plotFile << "set ylabel \"amplitude of the response |A|\" font \",20\"" << endl;
	plotFile << endl;
	plotFile << "set format x '%.1t*10^{%T}'" << endl;
	plotFile << endl;
	plotFile << "set key top right" << endl;
	plotFile << endl;
	
	plotFile << "plot \"amplification_versus_noise.dat\" with linespoints pointtype 7 pointsize 0.5 title 'simulation' ,\\" << endl;
	plotFile << "     \"amplification_versus_noise_theory.dat\" with lines title 'theory'" << endl;
	
	plotFile.close();
	
	sys_result = system("cd data; gnuplot plot_amplification_versus_noise");
	
	
	////////////////////////////////////////////////////////////////////////
	
	
	return 0;
}