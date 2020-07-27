#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <tchar.h>

#include "fungicideDynamics.h"
#include "main.h"

void sprayDosages()
{                 // these are doses for at risk fungicide. 
	doseHR.T1 = 1.0; //1.0 //1.37 for full dose ? 
	doseHR.T2 = 1.0;
	
	//doseHR.T1 = 0.0;
	//doseHR.T2 = 0.0;
}

void sprayTiming()
{
	timing.T1 = T_EMERGENCE_L3;			// Emergence of leaf 3
//	timing.T1 = T_EMERGENCE_L2;			// Emergence of leaf 3
	timing.T2 = T_FLAG_LEAF_EMERGENCE;	// Emergence of leaf 1
}

void initialiseLeafAreasAtTimeOfSpray()
{
	leafArea.T1 = 0.0;	
	leafArea.T2 = 0.0;	
}

void applySprayProgram()
{
	initialiseLeafAreasAtTimeOfSpray();
	sprayTiming();
	sprayDosages();
}

double sprayDoseHR(double t, double *pVar)
{
	double	latentAreaAsc = 0.0, latentAreaPyc = 0.0, infectiousArea = 0.0;
	double	el0, el1;
	double	DECAY_HR = 0.011;

	// Determine the leaf area present to intercept fungicides at time of spray
	if(t >= timing.T1 && t < timing.T1 + STEPSIZE)
	{
		leafArea.T1 = pVar[3 * NGENOTYPES + 2];
	}
	if(t >= timing.T2 && t < timing.T2 + STEPSIZE)
	{
		leafArea.T2 = pVar[3 * NGENOTYPES + 2];
	}	
	
	el0 = - DECAY_HR * (t - timing.T1);
	el1 = - DECAY_HR * (t - timing.T2);
	
	if(t >= timing.T2)
	{
		effectiveDoseHR = (doseHR.T1 / leafArea.T1) * exp(el0) + (doseHR.T2 / leafArea.T2) * exp(el1);
	}
	else if(t >= timing.T1)
	{
		effectiveDoseHR = (doseHR.T1 / leafArea.T1) * exp(el0);
	}
	else effectiveDoseHR = 0.0;

	return effectiveDoseHR;
}

double sprayEffectHR(double t, double *pVar, double sensitivityEffect)
{
	// Calculates the reduction of the life history parameter affected by the fungicide
	double DOSE_RESPONSE_HR = 300, MAX_REDUCTION_HR = 0.65;
//  double DOSE_RESPONSE_HR = 300, MAX_REDUCTION_HR = 0.65;
// The default value of MAX_REDUCTION_HR for a systemic fungicide is 0.65. (IE+LP+SP affected)
// The default value of MAX_REDUCTION_HR for a protectant fungicide is 0.999, which gives the same level of disease control in the first year. (only IE affected)
	double Effect, temp;

	temp = - DOSE_RESPONSE_HR * sprayDoseHR(t, pVar);
	Effect = 1 - MAX_REDUCTION_HR * (1 - sensitivityEffect) * (1 - exp(temp));

	return Effect;
}

double sprayEffectHR2nd(double t, double *pVar, double sensitivityEffect) // This is copied from sprayEffectHR. To make a second fungicide be present, as a mixing partner which is not evolving, I add a dose of a second fungicide. I keep sprayDoseHR the same and assume the same dose of each fungicdie is used (eg 0.5 and 0.5)
{
	// Calculates the reduction of the life history parameter affected by the fungicide
	double DOSE_RESPONSE_HR2 = 300, MAX_REDUCTION_HR2 = 0.78;                  // these determine the efficacy of the mixing partner. I assume is is identical to the at risk fungicide. 
	//  double DOSE_RESPONSE_HR = 300, MAX_REDUCTION_HR = 0.65;
	// The default value of MAX_REDUCTION_HR for a systemic fungicide is 0.65. (IE+LP+SP affected)
	// The default value of MAX_REDUCTION_HR for a protectant fungicide is 0.999, which gives the same level of disease control in the first year. (only IE affected)
	double Effect2, temp;

	// Dose of mixing partner (as proportion of the main fungicide)
	double dose2 = 0.5;

	temp = -DOSE_RESPONSE_HR2 * sprayDoseHR(t, pVar);
	Effect2 = 1 - MAX_REDUCTION_HR2 * (1 - sensitivityEffect) * (1 - exp(dose2*temp));

	return Effect2;
}

double sprayEffectPartRes(double t, double *pVar, double sensitivityEffect)
{
	// Effect of fungicide on the partially resistant strain
	// Calculates the reduction of the life history parameter affected by the fungicide
	double DOSE_RESPONSE_HR = 300, MAX_REDUCTION_HR = 0.65;
	double Effect, temp;

	temp = - DOSE_RESPONSE_HR * sprayDoseHR(t, pVar);
	Effect = 1 - MAX_REDUCTION_HR * (1 - sensitivityEffect) * (1 - exp(temp));
	
	return Effect;
}
