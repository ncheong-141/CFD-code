/* Implementation of physics model FLUX function*/

#include	<iostream>
#include "physicsmodel.h"	// "Link" to interface 
#include "INITARRAYS.h"


void FLUX(INITARRAYS::DS_FLUX& FL, INITARRAYS::DS_BCON& B, INITARRAYS& V)
{
	// Resets FLUX array to zero
	FL.Qarray_reset();

	// Calculating flux contribution of aerofoil surface (e.g. pressure boundary cond)
	for (int i = 0; i < (V.ICU - 1); i++)
	{
		FL.DX = V.XP[V.IB[i]] - V.XP[V.IA[i]];
		FL.DY = V.YP[V.IB[i]] - V.YP[V.IA[i]];

		// Update fluxes for cells. 
		// No mass flow through physical boundary so no Q[i][0] or [3]
		FL.Q[V.K[i]][1] = FL.DY * V.Pr[V.K[i]];
		FL.Q[V.K[i]][2] = -FL.DX * V.Pr[V.K[i]];
	}

	// Calculating flux contribution of OUTER BOUNDARY
	// Using BW as they contain the boundary conditions. 
	// Note this might not be -1.. 
	for (int i = V.ICU - 1; i < V.NBD; i++)
	{
		FL.U  = B.BW[i][1] / B.BW[i][0];		// Horizontal velocity
		FL.Vv = B.BW[i][2] / B.BW[i][0];		// Vertical velocity
		FL.DX = V.XP[V.IB[i]] - V.XP[V.IA[i]];
		FL.DY = V.YP[V.IB[i]] - V.YP[V.IA[i]];
		FL.QK = FL.DY * FL.U - FL.DX * FL.Vv;

		// Updata fluxes for cells. Note, no cell P as P doesnt exist at boundaries.
		FL.Q[V.K[i]][0] = FL.Q[V.K[i]][0] + (FL.QK * B.BW[i][0]); // Mass flux (density)
		FL.Q[V.K[i]][1] = FL.Q[V.K[i]][1] + (FL.QK * B.BW[i][1]) + (B.BW[i][4] * FL.DY);
		FL.Q[V.K[i]][2] = FL.Q[V.K[i]][2] + (FL.QK * B.BW[i][2]) - (B.BW[i][4] * FL.DX);
		FL.Q[V.K[i]][3] = FL.Q[V.K[i]][3] + (FL.QK * (B.BW[i][3] + B.BW[i][4]));			
	}

	// Calculating flux contribution inside the domain 
	for (int i = V.NBD; i < V.NEDGE; i++)
	{
		// Note averaging velocities between cells K and P.
		FL.U = 0.5f * ((V.W[V.K[i]][1] / V.W[V.K[i]][0])
			+ (V.W[V.P[i]][1] / V.W[V.P[i]][0]));
		FL.Vv = 0.5f * ((V.W[V.K[i]][2] / V.W[V.K[i]][0])
			+ (V.W[V.P[i]][2] / V.W[V.P[i]][0]));
		FL.DX = V.XP[V.IB[i]] - V.XP[V.IA[i]];
		FL.DY = V.YP[V.IB[i]] - V.YP[V.IA[i]];
		FL.QK = FL.DY * FL.U - FL.DX * FL.Vv;
		FL.DPr = V.Pr[V.K[i]] + V.Pr[V.P[i]];

		// Applying fluxes to both cells K and P 
		FL.Q[V.K[i]][0] = FL.Q[V.K[i]][0] + (0.5f * FL.QK * (V.W[V.K[i]][0] + V.W[V.P[i]][0])); 	// Mass flux in K
		FL.Q[V.P[i]][0] = FL.Q[V.P[i]][0] - (0.5f * FL.QK * (V.W[V.K[i]][0] + V.W[V.P[i]][0])); 	// Mass flux in P

		FL.Q[V.K[i]][1] = FL.Q[V.K[i]][1] + (0.5f * (FL.QK * (V.W[V.K[i]][1] + V.W[V.P[i]][1]) + (FL.DPr * FL.DY))); // Mass flux in horizontal momentum 
		FL.Q[V.P[i]][1] = FL.Q[V.P[i]][1] - (0.5f * (FL.QK * (V.W[V.K[i]][1] + V.W[V.P[i]][1]) + (FL.DPr * FL.DY)));

		FL.Q[V.K[i]][2] = FL.Q[V.K[i]][2] + (0.5f * (FL.QK * (V.W[V.K[i]][2] + V.W[V.P[i]][2]) - (FL.DPr * FL.DX))); // Mass flux in vertical momentum
		FL.Q[V.P[i]][2] = FL.Q[V.P[i]][2] - (0.5f * (FL.QK * (V.W[V.K[i]][2] + V.W[V.P[i]][2]) - (FL.DPr * FL.DX)));

		FL.Q[V.K[i]][3] = FL.Q[V.K[i]][3] + (0.5f * FL.QK * (V.W[V.K[i]][3] + V.W[V.P[i]][3] + FL.DPr)); 	// Mass flux in energy
		FL.Q[V.P[i]][3] = FL.Q[V.P[i]][3] - (0.5f * FL.QK * (V.W[V.K[i]][3] + V.W[V.P[i]][3] + FL.DPr)); 	// 
	}
}