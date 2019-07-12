/* Implementation of physics model TIME function*/

#include	<iostream>
#include "physicsmodel.h"	// "Link" to interface 
#include "INITARRAYS.h"

void TIME(INITARRAYS::DS_TIME& TI, INITARRAYS::DS_BCON& B, INITARRAYS& V)
{
	// Reset TIME arrays to zero 
	TI.DTarray_reset();

	// Calculate time step around aerofoil boundary 
	for (int i = 0; i < (V.ICU - 1); i++)
	{
		TI.C2	= V.GM * V.Pr[V.K[i]] / V.W[V.K[i]][0];	// Speed of sound
		TI.DX	= V.XP[V.IB[i]] - V.XP[V.IA[i]];
		TI.DY	= V.YP[V.IB[i]] - V.YP[V.IA[i]];
		TI.DL2	= TI.DX * TI.DX + TI.DY * TI.DY;
		TI.U	= V.W[V.K[i]][1] / V.W[V.K[i]][0];		// Horizontal velocity
		TI.Vv	= V.W[V.K[i]][2] / V.W[V.K[i]][0];		// Vertical velocity

		// Timestep calculation
		TI.DT2[V.K[i]] = TI.DT2[V.K[i]] + (std::abs(TI.DY * TI.U - TI.DX * TI.Vv)) + (std::pow(std::abs(TI.C2 * TI.DL2), 0.5f));
	}

	// Calculate time step around outer boundary 
	for (int i = V.ICU - 1; i < V.NBD; i++)
	{
		TI.C2	= V.GM * B.BW[i][4] / B.BW[i][0];	// Speed of sound
		TI.DX	= V.XP[V.IB[i]] - V.XP[V.IA[i]];
		TI.DY	= V.YP[V.IB[i]] - V.YP[V.IA[i]];
		TI.DL2	= TI.DX * TI.DX + TI.DY * TI.DY;
		TI.U	= B.BW[i][1] / B.BW[i][0];			// Horizontal velocity
		TI.Vv	= B.BW[i][2] / B.BW[i][0];			// Vertical velocity

		// Timestep calculation
		TI.DT2[V.K[i]] = TI.DT2[V.K[i]] + (std::abs(TI.DY * TI.U - TI.DX * TI.Vv)) + (std::pow(std::abs(TI.C2 * TI.DL2), 0.5f));
	}

	// Calculating time step inside the domain 
	// Visits every edge and calculates velocity corresponding to the largest eigenvalue
	for (int i = V.NBD; i < V.NEDGE; i++)
	{
		TI.C2	= V.GM * (V.Pr[V.K[i]] + V.Pr[V.P[i]]) / (V.W[V.K[i]][0] + V.W[V.P[i]][0]);
		TI.U	= 0.5f * ((V.W[V.K[i]][1] / V.W[V.K[i]][0]) + (V.W[V.P[i]][1] / V.W[V.P[i]][0]));
		TI.Vv	= 0.5f * ((V.W[V.K[i]][2] / V.W[V.K[i]][0]) + (V.W[V.P[i]][2] / V.W[V.P[i]][0]));
		TI.DX	= V.XP[V.IB[i]] - V.XP[V.IA[i]];
		TI.DY	= V.YP[V.IB[i]] - V.YP[V.IA[i]];
		TI.DL2	= TI.DX * TI.DX + TI.DY * TI.DY;

		// Calculate timestep for both cells K and P
		TI.DT2[V.K[i]] = TI.DT2[V.K[i]] + (std::abs(TI.DY * TI.U - TI.DX * TI.Vv)) + (std::pow(std::abs(TI.C2 * TI.DL2), 0.5f));
		TI.DT2[V.P[i]] = TI.DT2[V.P[i]] + (std::abs(TI.DY * TI.U - TI.DX * TI.Vv)) + (std::pow(std::abs(TI.C2 * TI.DL2), 0.5f));
	}

	// Invert all timesteps
	for (int i = 0; i < V.NC; i++)
	{
		if (TI.DT2[i] > 0.0f)
		{
			TI.DT[i] = 1.0f / TI.DT2[i];
		}
		else
		{
			std::cout << "Negative time!!" << std::endl;
			std::cin.get();
			exit(-1);	// Breakpoint exit. 
		}

	}
}