#pragma once
/* Initialize variables, prepare spatial discretization data and pre-calculate any calculations for time*/

// External modules
#include <iostream>

// Internal modules
#include "Instrumentation and DEBUG.h"

// Constant for maths. 
// const double pi = (std::acos(-1));

const float pi = 3.14159265358979323846f;


// Initiate (V)ariable class to group variables together in memory 
// Making it public as I dont want a million getter functions. 
class INITARRAYS
{
private:
	// Initialize variables data type for extracting data. 
	std::string DC0, DC1, DC2, DC3, DC4, line;

public:

	// Grid properties (Need to check of i need to +-1 off ICU and NBD due to difference in indexing between matlab and every other language..) 
	// static to allow it to set size of array. (promises compiler not to change the stati variable for any instance of this class) 
	static const int ICL = 0;			// Aerofoil inner boundary (1st point)
	static const int ICU = 101;			// Aerofoil inner boundary (Last point)
	static const int NBD = 280;			// Outer boundary around aerofoil (Number of points)
	static const int NPTI = 2940;		// Total number of grid points
	static const int NC = 5600;			// Total number of cells
	static const int NEDGE = 8540;		// Total number of edges
	static const int NCYC = 2000;		// Number of iterations


	// Initialize arrays to hold grid data with data type. 
	int     K[NEDGE], P[NEDGE], IA[NEDGE], IB[NEDGE];
	float   XP[NPTI], YP[NPTI];
	float	W[NC][4], WS[NC][4], Pr[NC]; 
	float   CP[100], XF[100], CL[NCYC], CD[NCYC];

	// Main loop variables 
	float COEFF, U, Vv, CX, CY, CM, XCM, YCM, DCX, DCY;

	// Constructor and destructor of the initialization object. Using an initializer list! 
	INITARRAYS() : K{}, P{}, IA{}, IB{}, XP{}, YP{}, W{}, WS{}, Pr{}, CP{}, XF{}, CL{}, CD{}, COEFF(0.0f), U(0.0f), Vv(0.0f), CX(0.0f), 
				   CY(0.0f), CM(0.0f), XCM(0.0f), YCM(0.0f), DCX(0.0f), DCY(0.0f)
	{ d_print("Initialization object created\n") }
	~INITARRAYS() 
	{ d_print("Initialization object destroyed\n") }


	// Simulation properties 
	const float FMACH = 0.5f;			// Freestream Mach number
	const float AL = 2 * (pi / 180);	// Angle of incidence (rad)
	const float CFL = 5.4f;				// Courant number 

	// Runga kutta properties for time discretization
	const int NRKS	  = 4;							// Order of rungakutta method	
	const float ST[4] = {	(0.25f * CFL),			// Array of coefficients
							(0.33333f * CFL),
							(0.5f * CFL),
							(CFL)				};

	// Declare aerodynamic variables and initial conditions for simulation
	const float GM = 1.4f;								// Gamma of air
	const float GMm1 = GM - 1.0f;						// Gamma - 1 (Precalculated as it needs to be calcualted alot, could use macro..) 
	const float SA = std::sin(AL);						// Note sin is in radians. Precalcuted for speed. 
	const float CA = std::cos(AL);
	const float rho0 = 1.0f;							// Free stream normalized density
	const float P0 = 1.0f;								// Freestream normalized pressure 
	const float C0 = std::pow((GM * rho0 / P0), 0.5f);	// Freestream speed of sound.
	const float U0 = FMACH * C0 * CA;					//Free stream horizontal velocity
	const float V0 = FMACH * C0 * SA;					// Freestream vertical velocity 
	const float E0 = ((P0) / (rho0 * GMm1)) + (0.5f * ((U0 * U0) + (V0 * V0)));			// Freestream energy 
	const float CPD = GM * FMACH * FMACH * 0.5f;		// Freestream dynamic head. 

	/* Member functions for class to open file and extact data */
	void	IF1_extract();
	void	IF2_extract();

	/* Data structures for physics models.*/

	// FLUX model structure
	struct DS_FLUX
	{
	public:
		float Q[NC][4];				  // Declare array data 
		float DX, DY, U, Vv, QK, DPr; // Declare variables 

		/* Struct member functions */ 
		DS_FLUX() : Q{}, DX(0.0f), DY(0.0f), U(0.0f), Vv(0.0f), QK(0.0f), DPr(0.0f)
					{ d_print("Flux datastructure created")}
		~DS_FLUX()	{ d_print("Flux datastructured destroyed") }
		
		// initaits FLUX array at every function call. 
		void Qarray_reset(); 
	};

	struct DS_DISS
	{
	public:
		float D[NC][4], DW2[NC][4], DELW[4], T[4];		// Declare array data (Changes in simulation)
		float C2, U, Vv, DX, DY, DL2, A, EPS2, EPS4;	// Declare variables
		const float RK2 = 0.5f;							// 1st dissapation parameter, k2. // Const values 
		const float RK4 = 0.008f;						// 2nd dissapation parameter, k4. 

		/* Struct member functions */
		DS_DISS() : D{}, DW2{}, DELW{}, T{}, C2(0.0f), U(0.0f), Vv(0.0f), DX(0.0f), DY(0.0f), DL2(0.0f), A(0.0f), EPS2(0.0f), EPS4(0.0f)
				   { d_print("Diss datastructure created")}
		~DS_DISS() { d_print("Diss datastructured destroyed") }

		// Initiates dissapation array at every call 
		void Darray_reset();
	};

	struct DS_BCON
	{
		// Declare array data 
		float BW[NBD][5];

		// Declare variables 
		float	DX, DY, U, Vv, DH, STH, CTH, RE, CE, PE, XMACH, // BCON..
				S0, QN0, QT0, RI0, SE, QNE, QTE,
				RIE, C, SB, UB, VB, RB, PB, REB,
				QN, QTT;

		/* Struct member functions */
		DS_BCON() : BW{}, DX(0.0f), DY(0.0f), U(0.0f), Vv(0.0f), DH(0.0f), STH(0.0f), CTH(0.0f), RE(0.0f), CE(0.0f), PE(0.0f),
				    XMACH(0.0f), S0(0.0f), QN0(0.0f), QT0(0.0f), RI0(0.0f), SE(0.0f), QNE(0.0f), QTE(0.0f),
					RIE(0.0f), C(0.0f), SB(0.0f), UB(0.0f), VB(0.0f), RB(0.0f), PB(0.0f), REB(0.0f),
					QN(0.0f), QTT(0.0f)
					{ d_print("Boundary datastructure created")	}

		~DS_BCON()	{ d_print("Boundary datastructure destroyed") }
	};
	
	struct DS_TIME
	{
		float	 DT[NC], DT2[NC];
		float	 C2, DX, DY, DL2, U, Vv;

		/* Struct member functions */
		DS_TIME() : DT{}, DT2{}, C2(0.0f), DX(0.0f), DY(0.0f), DL2(0.0f), U(0.0f), Vv(0.0f)
					{ d_print("Time datastructure created")	}

		~DS_TIME()  { d_print("Time datastructure destroyed") }

		// Initate time array at every call
		void DTarray_reset();
	};


	// Establish initial conditions for simulation 
	void	init_conds(DS_BCON& B);

};

