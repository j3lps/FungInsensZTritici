#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "standardFunctions.h"
#include "main.h"

//_____________________________________________________________________________________________

//RUNGE-KUTTA CATEGORY 4

void rk4(double pPrev[], double pNext[], double dpdt[])
/* Reads in the previous variable values in p[], the differential equations with p[] values 
as 'initial conditions' and the timestep h, and the end reads out the new variable values */
{
	unsigned int		i;					/* counter */
	double				*k1, *k2, *k3, *k4;

	k1 = (double *)calloc(nDiff, sizeof(double));
	k2 = (double *)calloc(nDiff, sizeof(double));
	k3 = (double *)calloc(nDiff, sizeof(double));
	k4 = (double *)calloc(nDiff, sizeof(double));

	for (i = 0; i <= nDiff - 1; ++ i)
	{
		k1[i] = dpdt[i];
		k2[i] = dpdt[i] + 0.5 * STEPSIZE * k1[i];
		k3[i] = dpdt[i] + 0.5 * STEPSIZE * k2[i];
		k4[i] = dpdt[i] + STEPSIZE * k3[i];
		pNext[i] = pPrev[i] + STEPSIZE * ((k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6);
	}
	free(k1);
	free(k2);
	free(k3);
	free(k4);
}

int integerPower(int base, int power)
{
	int i, result = 1;

	for(i = 0; i < power; i++)
	{
		result *= base;
	}
	return(result);
}
