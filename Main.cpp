/* Computational Fluid Dynamics
Spatial discretion:		Unstructured 2D grid
Physics model:			Euler method with Runga Kutta time stepping
Analysis:				Aerofoil coefficient of lift/drag for the NACA0012	*/

// Include interface files for functions (headers .h) 
#include "INITARRAYS.h"
#include "physicsmodel.h"
#include "Instrumentation and DEBUG.h"

// Copy and paste namespaces to this file (#Include) 
#include <iostream>		// For maths operations and output. 
#include <iomanip>		// For setprecision of output
#include <thread>

/*Keep in mind: 
	- Keep in mind that things may be bound by computation or data acquisistion
	- Profile the system and see what part is taking the longest, improve that.
	- CPU pipeline and branch predictor 
*/

/* Optimizing to do:																										Done? 			
	- Ensure loops are all ROW major (loops ALONG row first of 2D arrays)											    
	- Keep data close together in memory																					(DONE) [BASE]                                  
	- For arrays, prefer i++ to ++i as i++ gets element, then incrememnts which can be done simulatanously					(DONE) [BASE]
	- Prefer position dependent code																						(DONE) [BASE]
	- Prefer regular memory access patterns																			(kind of done) [BASE]
	- Use stack for everything you can as it is faster than heap															(DONE) [BASE]
	- Do not use globals as the compuer must assume it is visible to every .cpp, .h file (reloads it constantly)            (DONE) [BASE]
	- Prefer immutable data where you can (e.g. marked const which you promise the compiler not to change the variable)		(DONE) [BASE]

Code and data structure:
	- Keep data small																										(DONE) 
	- Minimize datatype sizes except for ints in math operations (e.g. always use 32 bit ints							   	
	  and use floats instead of doubles)	
	- Prefer 32-bit ints to all other sizes. It is sweet spot for integers as ALU (Arithmetic logic unit)					
	  inside the CPU can do one 64-bit operations at a time and two 32-bit operations at time.  
		o	64 bit may make some code 20x slower 
		o	8, 16-bit computations use conversion to 32 bits and back as the CPU works at a minimum of 32 bits! 
		o	Use small ints in arrays  (32 bits) |(DONE)
	- Avoid datatype conversions (e.g double to float and never do float to int)											(DONE)

	- Temporal cache coherency, ensure calculations are in an order such that a variable it not accessed, 
	  discard, then immediatly accessed again. 
	- Feed in data in 64 kb blocks, this is the size of one cacheline in the L1, L2, and L3 cache memory 
	  -> Avoid loading less than 64KB's of data from memory as there will be unused data in the cache. 	
	- Minimize flow control and avoid hard to predict branches, e.g. branches which are very random and gives BP 
	  a hard time predicting.
	  ->  Mark branches with LIKELY and UNLIKELY.
	- Avoid data dependencies as this doesnt allow for in-core (instruction level) parallelisation 
	  -> Reorginize equation structures to reduce data-dependency (get rid of required parethetsis) 
	- Prefer inline functions than function calls. Stops CPU requiring to jump to a new stack frame. 
	- Prefer unsigned to signed as unsigned is faster 
	- Strength reduction 
	  -> Reorginize code to allow for usage of fastest operations 
	- Minimize indirect writes by writing 64kb of data at a time. 
	- Loop unrolling, unrol a loop into its iterative calculations/multiply calculations in loop. 

*/

/* CPU i7 "bla bla", clockspeed 2.6 GHz, no multi-threading*/
/* Matlab: 814 seconds */
/* Optimizing done:	2n = 20				       Minimum       t_a,b        t_2 a,b         t_a+b          Differntial      % improval 
V0- Structured data (one big) [BASELINE]         6.1761       68.4345       128.007
  - Structured data (smaller data sets)          5.98996      60.11329      120.708		  128.54779           0.885          11.49 % 
V1- Reducing array size (double to float)	     4.04982	  40.8549       81.8442	      109.2894	          0.599          40.01 %
    to fit more data on cachelines and 
    ensuring no type conversions
- 
*/

/* Log: 
- Tried breaking down equations in FLUX into smaller equations -> Same/Slower
- Tried taking out function call for Q_array_reset -> Same/Slower... 
*/
// Start of code
int main()
{	
	/* Time code; start the clock. */
	Timer time_mainframe("Mainframe cycle total");

	/* Start of time loop */
	for (int t_i = 0; t_i < 10; t_i++)
	{
		/* Time code; start the clock. */
		Timer time_mainframecycle("Mainframe cycle");
		
		/* Establish data structures for simulation
		   Note, split into small segments for optimization...?? */
		INITARRAYS V;
		INITARRAYS::DS_FLUX dsQ;
		INITARRAYS::DS_DISS dsD;
		INITARRAYS::DS_TIME dsT;
		INITARRAYS::DS_BCON dsB;

		/* Extract grid data from files */
		V.IF1_extract();
		V.IF2_extract();

		/* Initialize data with initial conditions all cells and boundaries */
		V.init_conds(dsB);

		/* Start of simulation */
		for (int N = 0; N < V.NCYC; N++)
		{
			// Call boundary condition and time step functions
			BCON(dsB, V);
			TIME(dsT, dsB, V);

			// Apply artificial dissapation to solution but ONLY at first order.
			DISS(dsD, V);

			// Iterating through orders of Runga Kutta method 
			for (int NS = 0; NS < V.NRKS; NS++)
			{
				// Calculate flux through function for conservative variables
				FLUX(dsQ, dsB, V);

				// Updating solution (applying boundary conditions, fluxes etc to conservative variables)
				for (int i = 0; i < V.NC; i++)
				{
					V.COEFF = V.ST[NS] * (dsT.DT[i]);

					V.W[i][0] = std::abs(V.WS[i][0] - V.COEFF * (dsQ.Q[i][0] - dsD.D[i][0])); // abs as you cannoy have a -ve density
					V.W[i][1] = V.WS[i][1] - V.COEFF * (dsQ.Q[i][1] - dsD.D[i][1]);
					V.W[i][2] = V.WS[i][2] - V.COEFF * (dsQ.Q[i][2] - dsD.D[i][2]);
					V.W[i][3] = std::abs(V.WS[i][3] - V.COEFF * (dsQ.Q[i][3] - dsD.D[i][3])); // abs as you cannoy have a -ve energy

					// Calculate flow properties 
					V.U = V.W[i][1] / V.W[i][0];
					V.Vv = V.W[i][2] / V.W[i][0];
					V.Pr[i] = V.GMm1 * (V.W[i][3] - (0.5f * V.W[i][0]) * (V.U * V.U + V.Vv * V.Vv));

				} // End of "Updating Soln" 

				//Re-apply boundary conditions as solution has changed
				BCON(dsB, V);

			} // End of "Runga Kutta steps" 

			// Store old solution (Could be a problem if its not copied over!) 
			// Outside RK loop
			for (int i = 0; i < V.NC; i++)
			{
				V.WS[i][0] = V.W[i][0];
				V.WS[i][1] = V.W[i][1];
				V.WS[i][2] = V.W[i][2];
				V.WS[i][3] = V.W[i][3];
			}

			// Calculate aerodynamics of aerofoil from flow properties 
			// Reset force coefficients
			V.CX = 0.0f;
			V.CY = 0.0f;
			V.CM = 0.0f;

			for (int i = 0; i < V.ICU - 1; i++)
			{
				V.CP[i] = (V.Pr[V.K[i]] - 1) / V.CPD;
				V.XCM = 0.5f * (V.XP[V.IB[i]] + V.XP[V.IA[i]]) - 0.25f; // 0.25 is quater chord (normalized)
				V.YCM = 0.5f * (V.YP[V.IB[i]] + V.YP[V.IA[i]]);
				V.DCX = V.CP[i] * (V.YP[V.IB[i]] - V.YP[V.IA[i]]);
				V.DCY = -V.CP[i] * (V.XP[V.IB[i]] - V.XP[V.IA[i]]);
				V.CM = V.CM - (V.DCY * V.XCM) + (V.DCX * V.YCM);

				V.CX = V.CX + V.DCX;
				V.CY = V.CY + V.DCY;
				V.XF[i] = V.XCM + 0.25f;
			}

			// Calculate lift and drag coefficients
			V.CL[N] = (V.CY * V.CA) - (V.CX * V.SA);
			V.CD[N] = (V.CY * V.SA) + (V.CX * V.CA);

			// Print iterations
			//std::cout << "Iteration: " << N << "  CL = " << V.CL[N] << "  CD = " << V.CD[N] << "\n";

		} // End of iteration loop
 
	} // End of time loop for mainframe cycle

	return 0;
}



