#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "main.h"
#include "canopyDevelopment.h"

/* Function specifies canopy growth */
double canopyGrowth(double t)
{
	double	growth = 0.0, LAI = 0.0, SAI = 0.0, diffLAI = 0.0;
	double	a = 8, b = 165, c = 1300, m = 2300, n = 300;
	/* NOTE: this function has parameters in common with Canopy_senescence(t)
	See H:\Femke\Note books\Work notes - hand written\Fungicide resistance
	\Syngenta - wheat integrated control\ Wrok progress - lab book\ ...
	... June 2015; p.3 for further details. */

	LAI = (a / (1 + exp(-(t - c) / b)));
	diffLAI = (a * exp(-(t - c) / b)) / (b * (1 + exp(-(t - c) / b)) *  (1 + exp(-(t - c) / b)));
	SAI = (a / (1 + exp(-(t - m) / n)));

	growth = diffLAI / (LAI - SAI);

	return growth;
}

/* Function specifies canopy senescence */
double canopySenescence(double t)
{
	double	senescence = 0.0, LAI = 0.0, SAI = 0.0, diffSAI = 0.0;
	double	a = 8, b = 165, c = 1300, m = 2300, n = 300;
	/* NOTE: this function has parameters in common with Canopy_senescence(t)
	See H:\Femke\Note books\Work notes - hand written\Fungicide resistance
	\Syngenta - wheat integrated control\ Wrok progress - lab book\ ...
	... June 2015; p.3 for further details. */

	LAI = (a / (1 + exp(-(t - c) / b)));
	diffSAI = (a * exp(-(t - m) / n)) / (n * (1 + exp(-(t - m) / n)) *  (1 + exp(-(t - m) / n)));
	SAI = (a / (1 + exp(-(t - m) / n)));

	senescence = diffSAI / (LAI - SAI);

	return senescence;
}

void	initialiseCanopy(double *pNew)
{
	double initialLAI = 0.000046;
	// Set the initial crop area
	pNew[0] = initialLAI;
	// Set the total LAI to have the same crop area
	pNew[3 * NGENOTYPES + 2] = initialLAI;
}

