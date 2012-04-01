#pragma once

#pragma unmanaged

#define _USE_MATH_DEFINES
#include<stdio.h>
#include <stdlib.h>
#include <time.h>
#include<math.h>
#include<vector>
#include"Body.h"

#include"Cloud.h"

#define isnan(x) ((x) != (x))

/* This program is a Picard-Maclaurin solver for Newton's N Body Problem. 
The Maclaurin polynomial are determined   from Picard iteration using Cauchy 
products. The algorithm is at www.math.jmu.edu/~jim/picard.htm.  You are 
asked to give the number of bodies, N,the mass of each body, the initial 
position (alpha[i,0],beta[i,0],gamma[i,0]) and the initial velocity 
(delta[i,0],rho[i,0],lambda[i,0]) of each body, the time step,  the time to 
run the simulation and the degree of the Maclaurin polynomial for the 
numerical solution.

Give the number of bodies (masses). */

#define N 200

/* Give the degree of the Maclaurin polynomials.*/

#define K1 1
#define DEF_RESOLUTION 600.
#define DEF_DAYS 900
#define DEF_OUTPUT false

#define gaussian 0.01720209895
#define gaussianInvert 58.132440867048959743368991607853

#define AUtoKM 149598000.
#define AUDtoKMS 1731.45833


class NBSolver
{


	//Constants
//	const double gaussian = 0.01720209895;
//	const double gaussianInvert = 58.132440867048959743368991607853;

public:
	/*double alpha[N+1][K1+1], beta[N+1][K1+1], gamma[N+1][K1+1];
	double delta[N+1][K1+1], rho[N+1][K1+1], lambda[N+1][K1+1];
	double mass[N+1], sigma[N+1][N+1][K1+1];*/

	Body bodies[N+1][K1+1];
	double mass[N+1];
	double sigma[N+1][N+1][K1+1];
	
private:
	double F[N+1][N+1][K1+1], G[N+1][N+1][K1+1];
	double H[N+1][N+1][K1+1];

public:
	//Computation vars:
	double h;
	int numTimesteps;
	int currentTimestep;

	NBSolver(void);
	~NBSolver(void);

	void Init();
	void RunLoop(int numLoops);
	void SumPolynomial(double t);
	void SetResolution(double resolutionInSec, double num_days);

	double GetElapsedDays();

private:
	
	void InitBodies();
	void InitBodies2();
	void InitCalculations();
	void LoopPolynomial(int k);

};

#pragma managed
