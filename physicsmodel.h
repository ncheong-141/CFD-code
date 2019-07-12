#pragma once
/* Physics model functions for simulation*/

#include <iostream>
#include "INITARRAYS.h"

// Inputting data from INITARRAYS.
/* The & operator, this syntax means that any use of object in the function body
   refers to the actual object which was passed into the function and not a copy.
   All modifications will modify this object and be visible once the function has completed. */
void FLUX(INITARRAYS::DS_FLUX& FL, INITARRAYS::DS_BCON& B, INITARRAYS& V);
void TIME(INITARRAYS::DS_TIME& TI, INITARRAYS::DS_BCON& B, INITARRAYS& V);
void BCON(INITARRAYS::DS_BCON& B,  INITARRAYS& V);
void DISS(INITARRAYS::DS_DISS& DI, INITARRAYS& V);