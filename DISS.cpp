/* Implementation of physics model DISS function*/

#include	<iostream>
#include "physicsmodel.h"	// "Link" to interface 
#include "INITARRAYS.h"

void DISS(INITARRAYS::DS_DISS& DI, INITARRAYS& V)
{

	// Reset DISS array to zeros. 
	DI.Darray_reset();

	// Not considering boundary conditions so no first two loops. 

	// Inner fluid dissapation calculations. 
	for (int i = V.NBD; i < V.NEDGE; i++)
	{
		// Calculate difference in governing equaiton variables 
		DI.DELW[0] = V.W[V.P[i]][0] - V.W[V.K[i]][0];
		DI.DELW[1] = V.W[V.P[i]][1] - V.W[V.K[i]][1];
		DI.DELW[2] = V.W[V.P[i]][2] - V.W[V.K[i]][2];
		DI.DELW[3] = V.W[V.P[i]][3] - V.W[V.K[i]][3] + V.Pr[V.P[i]] - V.Pr[V.K[i]];

		// Cell K calculations
		DI.DW2[V.K[i]][0] = DI.DW2[V.K[i]][0] + DI.DELW[0];
		DI.DW2[V.K[i]][1] = DI.DW2[V.K[i]][1] + DI.DELW[1];
		DI.DW2[V.K[i]][2] = DI.DW2[V.K[i]][2] + DI.DELW[2];
		DI.DW2[V.K[i]][3] = DI.DW2[V.K[i]][3] + DI.DELW[3];

		// Cell P calculations
		DI.DW2[V.P[i]][0] = DI.DW2[V.P[i]][0] - DI.DELW[0];
		DI.DW2[V.P[i]][1] = DI.DW2[V.P[i]][1] - DI.DELW[1];
		DI.DW2[V.P[i]][2] = DI.DW2[V.P[i]][2] - DI.DELW[2];
		DI.DW2[V.P[i]][3] = DI.DW2[V.P[i]][3] - DI.DELW[3];

		// Flow properties 
		DI.C2 = V.GM * (V.Pr[V.K[i]] + V.Pr[V.P[i]]) / (V.W[V.K[i]][0] + V.W[V.P[i]][0]);
		DI.U = 0.5f * ((V.W[V.K[i]][1] / V.W[V.K[i]][0]) + (V.W[V.P[i]][1] / V.W[V.P[i]][0]));
		DI.Vv = 0.5f * ((V.W[V.K[i]][2] / V.W[V.K[i]][0]) + (V.W[V.P[i]][2] / V.W[V.P[i]][0]));
		DI.DX = V.XP[V.IB[i]] - V.XP[V.IA[i]];
		DI.DY = V.YP[V.IB[i]] - V.YP[V.IA[i]];
		DI.DL2 = DI.DX * DI.DX + DI.DY * DI.DY;
		DI.A = (std::abs(DI.DY * DI.U - DI.DX * DI.Vv)) + (std::pow(std::abs(DI.C2 * DI.DL2), 0.5f));		// Estimation of largest eigenvalue

		// Dissapation activation controllers 
		DI.EPS2 = DI.RK2 * std::abs((V.Pr[V.P[i]] - V.Pr[V.K[i]]) / (V.Pr[V.P[i]] + V.Pr[V.K[i]]));
		DI.EPS4 = std::fmax(0.0f, (DI.RK4 - DI.EPS2));

		// Calculate corresponding dissapation 
		DI.T[0] = DI.A * (DI.EPS2 * DI.DELW[0] - (DI.EPS4 * (DI.DW2[V.P[i]][0] - DI.DW2[V.K[i]][0])));
		DI.T[1] = DI.A * (DI.EPS2 * DI.DELW[1] - (DI.EPS4 * (DI.DW2[V.P[i]][1] - DI.DW2[V.K[i]][1])));
		DI.T[2] = DI.A * (DI.EPS2 * DI.DELW[2] - (DI.EPS4 * (DI.DW2[V.P[i]][2] - DI.DW2[V.K[i]][2])));
		DI.T[3] = DI.A * (DI.EPS2 * DI.DELW[3] - (DI.EPS4 * (DI.DW2[V.P[i]][3] - DI.DW2[V.K[i]][3])));

		// For K cells: 
		DI.D[V.K[i]][0] = DI.D[V.K[i]][0] + DI.T[0];
		DI.D[V.K[i]][1] = DI.D[V.K[i]][1] + DI.T[1];
		DI.D[V.K[i]][2] = DI.D[V.K[i]][2] + DI.T[2];
		DI.D[V.K[i]][3] = DI.D[V.K[i]][3] + DI.T[3];

		// For P cells:
		DI.D[V.P[i]][0] = DI.D[V.P[i]][0] - DI.T[0];
		DI.D[V.P[i]][1] = DI.D[V.P[i]][1] - DI.T[1];
		DI.D[V.P[i]][2] = DI.D[V.P[i]][2] - DI.T[2];
		DI.D[V.P[i]][3] = DI.D[V.P[i]][3] - DI.T[3];
	}
}