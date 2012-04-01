#include "NBSolver.h"

using namespace std;

#pragma unmanaged


NBSolver::NBSolver(void)
{
	numTimesteps = 0;
	h = 1;
	currentTimestep = 0;
}

NBSolver::~NBSolver(void)
{
}

void NBSolver::Init()
{
	/*Give the time step and time.*/

	SetResolution(DEF_RESOLUTION, DEF_DAYS);
/*
	double secInYear = 31556926.0;
	
	//seconds per timestep
	double resolutionInSec = DEF_RESOLUTION; 

	// Determine the number of time steps for this simulation.

	int secsInDay = 86400;
	double num_days = 203.;
	num_days = DEF_DAYS;

	numTimesteps = num_days * secsInDay / resolutionInSec;

	// h is in units of 1 earth-year circumference (radians)
	h = resolutionInSec / (86400. / gaussian);
*/
	currentTimestep = 0;

	//InitBodies2();
}

double NBSolver::GetElapsedDays()
{
	return currentTimestep * (h / gaussian);
}

void NBSolver::SetResolution(double resolutionInSec, double num_days)
{
	double secsInDay = 86400.;

	double elapsedDays = GetElapsedDays();

	if (resolutionInSec == 0)
		resolutionInSec = 1;

	numTimesteps = num_days * secsInDay / resolutionInSec;

	h = resolutionInSec / (secsInDay / gaussian);

	currentTimestep = elapsedDays * secsInDay / resolutionInSec;
	
	//if (currentTimestep > numTimesteps)
	//	  currentTimestep = numTimesteps;
}

void NBSolver::InitBodies()
{
	int i;
	int j;

	/*  Give the masses. */
	
	/*  mass of body Sun        =  1.
	mass of body Mercury    =  1.66013679527193035E-7
	mass of body Venus      =  2.44783959796682464E-6
	mass of body Earth      =  3.04043273871083524E-6
	mass of body Mars       =  3.22714936215392876E-7
	mass of body Jupiter    =  9.54790662147324233E-4
	mass of body Saturn     =  2.85877644368210402E-4
	mass of body Uranus     =  4.35540069686411149E-5
	mass of body Neptune    =  5.17759138448793649E-5
	mass of body Pluto      =  7.6923076923076926E-9 */
	
	//Mass of combined planet and moon systems
	double massSun = 1.989100E+30;
	double massMercury = 3.302000E+23;
	double massVenus = 4.868500E+24;
	double massEarth = 6.047090E+24;
	double massMars = 6.418500E+23;
	double massJupiter = 1.898993E+27;
	double massSaturn = 5.686028513220E+26;
	double massUranus = 8.684113190000E+25;
	double massNeptune = 1.024514700000E+26;
	double massPluto = 1.504000000000E+22;

	mass[1]  =  1.;
	mass[2]  =  massMercury / massSun;
	mass[3]  =  massVenus / massSun;
	mass[4]  =  massEarth / massSun;
	mass[5]  =  massMars / massSun;
	mass[6]  =  massJupiter / massSun;
	mass[7]  =  massSaturn / massSun;
	mass[8]  =  massUranus / massSun;
	mass[9]  =  massNeptune / massSun;
	mass[10] =  massPluto / massSun;
	mass[11] = massJupiter / massSun;
	/*mass[12] = 1000 / massSun;
	mass[13] = 1000 / massSun;
	mass[14] = 1000 / massSun;
	mass[15] = 1000 / massSun;*/
/*
	mass[1]       =  1.;
	mass[2]    =  1.66013679527193035E-7;
	mass[3]      =  2.44783959796682464E-6;
	mass[4]     =  3.04043273871083524E-6;
	mass[5]      =  3.22714936215392876E-7;
	mass[6]    =  9.54790662147324233E-4;
	mass[7]     =  2.85877644368210402E-4;
	mass[8]     =  4.35540069686411149E-5;
	mass[9]    =  5.17759138448793649E-5;
	mass[10]      =  7.6923076923076926E-9;
	
	*/
	/*    Set up the initial conditions  */
	//2000-Jan-01 00:00:00.0000

	//Uses planet/moon barycenters relative to solar system barycenter (J2000)
	/* SUN */
	bodies[1][0].x = -7.139147120601119E-03;
	bodies[1][0].y = -2.792019830318904E-03;
	bodies[1][0].z = 2.061825704683573E-04;
	bodies[1][0].vx = 5.374261885473961E-06 * gaussianInvert;
	bodies[1][0].vy = -7.410966640098345E-06 * gaussianInvert;
	bodies[1][0].vz = -9.422892899203406E-08 * gaussianInvert;
	
	/* MERCURY */
	
	bodies[2][0].x =  -1.478672244638274E-01;
	bodies[2][0].y =  -4.466929789292389E-01;
	bodies[2][0].z =  -2.313937957609384E-02;
	bodies[2][0].vx = 2.117424560897133E-02 * gaussianInvert;
	bodies[2][0].vy =  -7.105386129334717E-03 * gaussianInvert;
	bodies[2][0].vz = -2.522925839071368E-03 * gaussianInvert;
	
	/* VENUS */
	
	bodies[3][0].x =  -7.257693636228210E-01;
	bodies[3][0].y = -2.529582240743707E-02;
	bodies[3][0].z =  4.137802926392006E-02;
	bodies[3][0].vx =  5.189070601827867E-04 * gaussianInvert;
	bodies[3][0].vy =  -2.031355259129968E-02 * gaussianInvert;
	bodies[3][0].vz = -3.072686213835033E-04 * gaussianInvert;
	
	/* EARTH */
	
	bodies[4][0].x =  -1.756895992853660E-01; 
	bodies[4][0].y =    9.659716382992982E-01;
	bodies[4][0].z =   2.050240283368983E-04;
	bodies[4][0].vx =  -1.722463620800800E-02 * gaussianInvert;
	bodies[4][0].vy =  -3.020684845073816E-03 * gaussianInvert;
	bodies[4][0].vz =  -7.003423872759354E-08 * gaussianInvert;

	/* MARS */ 
	
	bodies[5][0].x =  1.383221919028136E+00;
	bodies[5][0].y =  -2.380173799124939E-02;
	bodies[5][0].z = -3.441183423421527E-02;
	bodies[5][0].vx = 7.533013358251247E-04 * gaussianInvert;
	bodies[5][0].vy =  1.517888772436052E-02 * gaussianInvert;
	bodies[5][0].vz =  2.996588796105037E-04 * gaussianInvert;
	
	/* JUPITER */
	
	bodies[6][0].x =  3.996320681110832E+00; 
	bodies[6][0].y =  2.932561823012092E+00;
	bodies[6][0].z =  -1.016168451332413E-01;
	bodies[6][0].vx =  -4.558099510764580E-03 * gaussianInvert;
	bodies[6][0].vy =  6.439346715906696E-03 * gaussianInvert;
	bodies[6][0].vz =  7.536243379766181E-05 * gaussianInvert;
	
	/* SATURN */
	
	bodies[7][0].x =  6.401418058908803E+00;
	bodies[7][0].y =  6.565252439589429E+00;
	bodies[7][0].z = -3.689199086691116E-01;
	bodies[7][0].vx = -4.285743775521677E-03 * gaussianInvert;
	bodies[7][0].vy =  3.884169867203816E-03 * gaussianInvert;
	bodies[7][0].vz = 1.027826777848790E-04 * gaussianInvert;
	
	/* URANUS */
	
	bodies[8][0].x = 1.442338133008372E+01;
	bodies[8][0].y = -1.373844069614078E+01;
	bodies[8][0].z = -2.379185362018481E-01;
	bodies[8][0].vx =  2.683753457289004E-03 * gaussianInvert;
	bodies[8][0].vy =  2.665032941399753E-03 * gaussianInvert;
	bodies[8][0].vz = -2.487077052909905E-05 * gaussianInvert;
	
	
	/* NEPTUNE */
	
	bodies[9][0].x =  1.680362717435433E+01;
	bodies[9][0].y = -2.499545061246260E+01;
	bodies[9][0].z = 1.274832850715894E-01;
	bodies[9][0].vx =  2.584745162555315E-03 * gaussianInvert;
	bodies[9][0].vy =  1.769356439106874E-03 * gaussianInvert;
	bodies[9][0].vz = -9.600253692009498E-05 * gaussianInvert;
	
	/* PLUTO */
	
	bodies[10][0].x = -9.884002870123060E+00;
	bodies[10][0].y = -2.796080827445471E+01;
	bodies[10][0].z = 5.851004954719751E+00;
	bodies[10][0].vx =  3.034076901648591E-03 * gaussianInvert;
	bodies[10][0].vy =  -1.545334774680257E-03 * gaussianInvert;
	bodies[10][0].vz = -7.122692805748236E-04 * gaussianInvert;

	/* ROCKET */
	
	bodies[11][0].x =  bodies[4][0].x * -1; 
	bodies[11][0].y =  bodies[4][0].y * -1; 
	bodies[11][0].z =  bodies[4][0].z * -1;
	bodies[11][0].vx =  bodies[4][0].vx * -1.2;
	bodies[11][0].vy =  bodies[4][0].vy;
	bodies[11][0].vz =  bodies[4][0].vz;
	
}
void NBSolver::InitBodies2()
{
	int i;
	int j;
	srand(7);

	//Mass of combined planet and moon systems
	double massSun = 1.989100E+30;
	double massMercury = 3.302000E+23;
	double massVenus = 4.868500E+24;
	double massEarth = 6.047090E+24;
	double massMars = 6.418500E+23;
	double massJupiter = 1.898993E+27;
	double massSaturn = 5.686028513220E+26;
	double massUranus = 8.684113190000E+25;
	double massNeptune = 1.024514700000E+26;
	double massPluto = 1.504000000000E+22;

	double avgvx = 0, avgvy = 0, avgvz = 0;
	for (int i = 1; i <= N; i++)
	{
		mass[i] = .1;
		double ratio = (i / (double)N);
		double angle = 2 * M_PI * ratio; //(rand() / (RAND_MAX+1.0)) * M_PI * 2;
		double radius = (rand() / (RAND_MAX+1.0)) * 14.5 + .5;
		double vel = ((rand() / (RAND_MAX+1.0)) * .4 + .05 * radius);
		double velAngle = ((rand() / (RAND_MAX+1.0)) - 0.5) * M_PI / 8;
		bodies[i][0].x =  cos(angle) * radius;
		bodies[i][0].y =  sin(angle) * radius;
		bodies[i][0].z =  ((rand() / (RAND_MAX+1.0)) - 0.5) * 10;
		double vx = cos(angle + M_PI / 2.0 + velAngle) * vel;
		double vy = sin(angle + M_PI / 2.0 + velAngle) * vel;
		double vz = 0;
		avgvx += vx;
		avgvy += vy;
		avgvz += vz;
		bodies[i][0].vx = vx;
		bodies[i][0].vy = vy;
		bodies[i][0].vz = vz;
		bodies[i][0].isDisabled = false;
	}
	avgvx /= N;
	avgvy /= N;
	avgvz /= N;

	for (int i = 1; i <= N; i++)
	{
		bodies[i][0].vx -= avgvx;
		bodies[i][0].vy -= avgvy;
		bodies[i][0].vz -= avgvz;
	}

}

void NBSolver::InitCalculations()
{
	int i, j;
	/*  Set up sigma from the initial conditions */

	double restitution = .5;
	double minDistance = 0.5;
	bool doCollision = false;
	for ( i = 1; i <=N; i++ )
	{
		Body left = bodies[i][0];
		if (left.isDisabled)
				continue;
		double dist = sqrt(left.x * left.x + left.y * left.y + left.z * left.z);
		/*double factor = pow(.9999, dist);
		bodies[i][0].vx *= factor;
		bodies[i][0].vy *= factor;
		bodies[i][0].vz *= factor;*/

        for ( j = i+1; j <= N; j++ )
		{
			if (bodies[j][0].isDisabled)
				continue;
			
			Body right = bodies[j][0];
			double a = left.x - right.x;
			double b = left.y - right.y;
			double c = left.z - right.z;

			double distance = sqrt(a*a + b*b + c*c);
			if (distance < minDistance)
				distance = minDistance;
			distance = 1./distance;
			sigma[i][j][0] = distance;

			sigma[j][i][0] = distance;
			
			if (doCollision && sigma[i][j][0] > 1/minDistance)
			{
				double totalMass = mass[i] + mass[j];
				bodies[i][0].vx = (restitution * mass[j] * (right.vx - left.vx) + mass[i] * left.vx + mass[j] * right.vx) / totalMass;
				bodies[i][0].vy = (restitution * mass[j] * (right.vy - left.vy) + mass[i] * left.vy + mass[j] * right.vy) / totalMass;
				bodies[i][0].vz = (restitution * mass[j] * (right.vz - left.vz) + mass[i] * left.vz + mass[j] * right.vz) / totalMass;

				bodies[j][0].vx = (restitution * mass[i] * (left.vx - right.vx) + mass[i] * left.vx + mass[j] * right.vx) / totalMass;
				bodies[j][0].vy = (restitution * mass[i] * (left.vy - right.vy) + mass[i] * left.vy + mass[j] * right.vy) / totalMass;
				bodies[j][0].vz = (restitution * mass[i] * (left.vz - right.vz) + mass[i] * left.vz + mass[j] * right.vz) / totalMass;

				if (restitution == 0)
				{
					mass[i] = totalMass;
					//mass[j] = 0.0;
					bodies[j][0].isDisabled = true;
				}
				else
				{
					double ax = (left.x + right.x)/2;
					double ay = (left.y + right.y)/2;
					double az = (left.z + right.z)/2;
					
					double vx = bodies[j][0].vx;
					double vy = bodies[j][0].vy;
					double vz = bodies[j][0].vz;
					double len = sqrt(vx*vx + vy*vy + vz*vz);

					bodies[j][0].x = ax + (vx / len) * minDistance * 0.51;
					bodies[j][0].y = ay + (vy / len) * minDistance * 0.51;
					bodies[j][0].z = az + (vz / len) * minDistance * 0.51;

					vx = bodies[i][0].vx;
					vy = bodies[i][0].vy;
					vz = bodies[i][0].vz;
					len = sqrt(vx*vx + vy*vy + vz*vz);

					bodies[i][0].x = ax + (vx / len) * minDistance * 0.51;
					bodies[i][0].y = ay + (vy / len) * minDistance * 0.51;
					bodies[i][0].z = az + (vz / len) * minDistance * 0.51;
				}
			}
		}
	}
}

void NBSolver::LoopPolynomial(int k)
{
//loop the bodies
	for (int i = 1; i <= N; i++ )
	{
		if (bodies[i][0].isDisabled)
			continue;
		//XYZ are easy

		bodies[i][k].x = bodies[i][k-1].vx / k;
		bodies[i][k].y = bodies[i][k-1].vy / k;
		bodies[i][k].z = bodies[i][k-1].vz / k;

		for (int j = i+1; j <= N; j++ )
		{
			if (bodies[j][0].isDisabled)
				continue;
			
			//Sigma is u
			//F is u2
			F[i][j][k-1] = 0;
			
			//G is u3
			G[i][j][k-1] = 0;
			for (int m = 0; m <= k-1; m++ )
			{
				F[i][j][k-1] += sigma[i][j][m]*sigma[i][j][k-1-m];
			}
			F[j][i][k-1] = F[i][j][k-1];
			
			//F & G must be separate loop otherwise compiler optimizations break the calculation
			for (int m = 0; m <= k-1; m++ )
			{
				G[i][j][k-1] += sigma[i][j][k-1-m]*F[i][j][m];
			}
			G[j][i][k-1] = G[i][j][k-1];
		}
		
		//Calculate vx, vy, vz
		bodies[i][k].vx = 0;
		bodies[i][k].vy = 0;
		bodies[i][k].vz = 0;

		for (int j = 1; j <= N; j++ )
		{
			if (bodies[j][0].isDisabled)
				continue;
			
			double Qx = 0;
			double Qy = 0;
			double Qz = 0;
			for (int m = 0; m <= k-1; m++ )
			{
				Body left = bodies[i][m];
				Body right = bodies[j][m];
				double g = G[i][j][k-1-m];
				Qx += (left.x - right.x)*g;
				Qy += (left.y - right.y)*g;
				Qz += (left.z - right.z)*g;
			}
			bodies[i][k].vx -= mass[j]*Qx;
			bodies[i][k].vy -= mass[j]*Qy;
			bodies[i][k].vz -= mass[j]*Qz;
		}

		bodies[i][k].vx /= k;
		bodies[i][k].vy /= k;
		bodies[i][k].vz /= k;

		//H is "Ajkm"
		for (int j = i+1; j <= N; j++ )
		{
			if (bodies[j][0].isDisabled)
				continue;
			
			H[i][j][k-1] = 0;
			//T is calculating u/sigma
			double T = 0;

			for (int m = 0; m <= k-1; m++ )
			{
				H[i][j][k-1] += 
					(bodies[i][k-1-m].x-bodies[j][k-1-m].x) * (bodies[i][m].vx-bodies[j][m].vx) + 
					(bodies[i][k-1-m].y-bodies[j][k-1-m].y) * (bodies[i][m].vy-bodies[j][m].vy) + 
					(bodies[i][k-1-m].z-bodies[j][k-1-m].z) * (bodies[i][m].vz-bodies[j][m].vz);
				
				T += H[i][j][m]*G[i][j][k-1-m];
			}

			sigma[i][j][k] = -T/k;
			sigma[j][i][k] = sigma[i][j][k];
		}
	}
}

void NBSolver::SumPolynomial(double t)
{
	int i, j, k;
	double x, y, z, vx, vy, vz;

	for ( i = 1; i <= N; i++ )
	{
		if (bodies[i][0].isDisabled)
				continue;
			
		x = bodies[i][K1-1].x + bodies[i][K1].x * t;
		y = bodies[i][K1-1].y + bodies[i][K1].y * t;
		z = bodies[i][K1-1].z + bodies[i][K1].z * t;
		vx = bodies[i][K1-1].vx + bodies[i][K1].vx * t;
		vy = bodies[i][K1-1].vy + bodies[i][K1].vy * t;
		vz = bodies[i][K1-1].vz + bodies[i][K1].vz * t;
		
		for ( k = 1; k <= K1-1; k++ )
		{
			x = x*t + bodies[i][K1-k-1].x;
			y = y*t + bodies[i][K1-k-1].y;
			z = z*t + bodies[i][K1-k-1].z;
			vx = vx*t + bodies[i][K1-k-1].vx;
			vy = vy*t + bodies[i][K1-k-1].vy;
			vz = vz*t + bodies[i][K1-k-1].vz;
		}

		/*for ( j = i+1; j <= N; j++ )
		{
			if (sigma[i][j][0] < .02)
			{
				double totalMass = mass[i] + mass[j];
				bodies[i][j].vx = (bodies[i][0].vx * mass[i] + bodies[j][0].vx * mass[j]) / totalMass;
				bodies[i][j].vy = (bodies[i][0].vy * mass[i] + bodies[j][0].vy * mass[j]) / totalMass;
				bodies[i][j].vz = (bodies[i][0].vz * mass[i] + bodies[j][0].vz * mass[j]) / totalMass;

				mass[i] = totalMass;
				mass[j] = 0.0;
				bodies[j][0].isDisabled = true;
			}
		}*/

		bodies[i][0].x = x;
		bodies[i][0].y = y;
		bodies[i][0].z = z;
		bodies[i][0].vx = vx;
		bodies[i][0].vy = vy;
		bodies[i][0].vz = vz;
	}
}

void NBSolver::RunLoop(int numLoops)
{
	//loop the timesteps
	for (int i=0; i<numLoops; i++)
	{
		/*if (currentTimestep > numTimesteps)
			break;*/

		InitCalculations();

		//loop the polynomials
		for (int k = 1; k <= K1; k++ )
		{
			LoopPolynomial(k);
		}
		
		SumPolynomial(h);

		currentTimestep++;
	}
}

#pragma managed
