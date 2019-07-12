#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "INITARRAYS.h"


void INITARRAYS::IF1_extract()
{
	// Get data from resource files. 
	std::ifstream inputfile("grid.txt");

	int count_i = 0; 

	// Check if file is open
	if (!inputfile.is_open()) { std::printf("\n\nNo file found...\n\n");	exit(-1); } // -1 is debug exit!	

	while (std::getline(inputfile, line) && true == !inputfile.eof())
	{
		std::stringstream ss(line);		// Marches down line? 
		std::getline(ss, DC0, '	');		// ?? 
		std::getline(ss, DC1, '	');		// K cell
		std::getline(ss, DC2, '	');		// IA 
		std::getline(ss, DC3, '	');		// IB 
		std::getline(ss, DC4, '	');		// P cell

		// Convert string into integer values (-1 to convert from matlab index system to C++ [start from 0]) 
		K[count_i] = std::stoi(DC1);
		K[count_i] -= 1;
		IA[count_i] = std::stoi(DC2);
		IA[count_i] -= 1;
		IB[count_i] = std::stoi(DC3);
		IB[count_i] -= 1;
		P[count_i] = std::stoi(DC4);

		// Cannot get negative index
		if (P[count_i] != 0)
		{
			P[count_i] -= 1;
		}
		
		// Incremenent counter for while loop
		count_i++;
	}
	inputfile.close();
}

void INITARRAYS::IF2_extract()
{
	// Get grid data
	std::ifstream inputfile2("grid2.txt");

	int count_j = 0;

	// Check if file is open
	if (!inputfile2.is_open()) { std::printf("\n\nNo file found...\n\n");	exit(-1); } // -1 is debug exit!	

	while (std::getline(inputfile2, line) && true == !inputfile2.eof())
	{

		std::stringstream ss(line);
		std::getline(ss, DC0, '	');		// XP
		std::getline(ss, DC1, '	');		// YP

		XP[count_j] = std::stof(DC0);
		YP[count_j] = std::stof(DC1);

		count_j++;
	}
	inputfile2.close();
}

void INITARRAYS::init_conds(DS_BCON& B)
{
	// Generate initial condition for working variable 
	for (int i = 0; i < NC; i++)
	{
		// For each sell, setting govern equation variables:
		W[i][0] = rho0;				// Density
		W[i][1] = rho0 * U0;		// Horizontal momentum
		W[i][2] = rho0 * V0;		// Vertical momentum
		W[i][3] = rho0 * E0;		// Potential of something
		Pr[i] = P0;				// Pressure 

		// Initialize store variables as initial solution.
		// This is previous iteration store for marchine solution
		WS[i][0] = rho0;			// Density
		WS[i][1] = rho0 * U0;		// Horizontal momentum
		WS[i][2] = rho0 * V0;		// Vertical momentum
		WS[i][3] = rho0 * E0;		// Potential of something

		// Initialization of boundary conditions at grid boundary
		// for each variable (freestream conditions)
		if (i < NBD)
		{
			B.BW[i][0] = rho0;
			B.BW[i][1] = rho0 * U0;
			B.BW[i][2] = rho0 * V0;
			B.BW[i][3] = rho0 * E0;
			B.BW[i][4] = P0;
		}
	}
}

void INITARRAYS::DS_FLUX::Qarray_reset()
{
	for (int i = 0; i < NC; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Q[i][j] = 0.0f;		// Assigning all values of Q as zero
		}
	}
}

void INITARRAYS::DS_TIME::DTarray_reset()
{
	// Pre-allocate time array with zeros
	for (int i = 0; i < NC; i++)
	{
		DT2[i] = 0.0f;
		DT[i] = 0.0f;
	}
}

void INITARRAYS::DS_DISS::Darray_reset()
{
	for (int i = 0; i < NC; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			// Preallocate dissapation arrays and reset them. 
			D[i][j] = 0.0f;
			DW2[i][j] = 0.0f;
		}
	}
}
