/* Implementation of physics model BCON function*/

#include	<iostream>
#include "physicsmodel.h"	// "Link" to interface 
#include "INITARRAYS.h"

// Calculates boundary conditioons for each cell at the boundarys

void BCON(INITARRAYS::DS_BCON& B, INITARRAYS& V)
{
	// Loop over outer boundary points. 
	for (int i = V.ICU - 1; i < V.NBD; i++)
	{
		B.DX  = V.XP[V.IB[i]] - V.XP[V.IA[i]];
		B.DY  = V.YP[V.IB[i]] - V.YP[V.IA[i]];
		B.DH  = std::pow(((B.DX * B.DX) + (B.DY * B.DY)), 0.5f);
		B.STH = B.DY / B.DH;
		B.CTH = B.DX / B.DH;
		B.RE  = V.W[V.K[i]][0];
		B.U   = V.W[V.K[i]][1] / B.RE;
		B.Vv  = V.W[V.K[i]][2] / B.RE;
		B.PE = V.Pr[V.K[i]];
		B.CE = std::pow((std::abs(V.GM * B.PE / B.RE)), 0.5f);
		B.XMACH = std::pow((B.U * B.U + B.Vv * B.Vv), 0.5f) / B.CE;

		// Freestream variables, 0 is freestream properties. 
		B.S0 = std::pow((V.P0 / V.rho0), V.GM);		// Entropy of freestream, near boundary must stay constant
		B.QN0 = (V.U0 * B.STH) - (V.V0 * B.CTH);		// Resolving velocities 
		B.QT0 = (V.U0 * B.CTH) + (V.V0 * B.STH);
		B.RI0 = B.QN0 - (2.0f * V.C0 / V.GMm1);		// Internal energy

		// Inside domain, E = local properties
		B.SE = std::pow((B.PE / B.RE), V.GM);
		B.QNE = (B.U * B.STH) - (B.Vv * B.CTH);
		B.QTE = (B.U * B.CTH) + (B.Vv * B.STH);
		B.RIE = B.QNE + (2.0f * B.CE / V.GMm1);		// Internal energy

		// CHeck if flow is subsonic or supersonic using mach number
		if (B.XMACH < 1.0f)
		{
			// Averaging Freestream data from last iteration
			B.QN = 0.5f * (B.RIE + B.RI0);
			B.C = V.GMm1 * 0.25f * (B.RIE - B.RI0);
			// if so, thren calculates from outside freestream ( something to do with projections of usable data..)
			if (B.QN < 0.0f)
			{
				B.QTT = B.QT0;
				B.SB = B.S0;
			}
			else //Calculate it from the local properties. 
			{
				B.QTT = B.QTE;
				B.SB = B.SE;
			}
		}
		else // Supersonic flow
		{
			if (B.QN0 < 0.0f)
			{
				B.QN  = B.QN0;
				B.C   = V.C0;
				B.QTT = B.QT0;
				B.SB  = B.S0;
			}
			else
			{
				B.QN  = B.QNE;
				B.C   = B.CE;
				B.QTT = B.QTE;
				B.SB  = B.SE;
			}
		}

		// Calculate primary variables at boundaries 
		B.UB = (B.QN * B.STH) + (B.QTT * B.CTH);
		B.VB = (B.QTT * B.STH) - (B.QN * B.CTH);
		B.RB = std::pow((B.C * B.C / V.GM / B.SB), (1.0f / V.GMm1));
		B.PB = B.C * B.C * B.RB / V.GM;
		B.REB = (B.PB / V.GMm1) + (0.5f * B.RB * ((B.UB * B.UB) + (B.VB * B.VB)));

		// Calculate conservative variables for boundaries! 
		B.BW[i][0] = B.RB;
		B.BW[i][1] = B.RB * B.UB;
		B.BW[i][2] = B.RB * B.VB;
		B.BW[i][3] = B.REB;
		B.BW[i][4] = B.PB;
	}
}