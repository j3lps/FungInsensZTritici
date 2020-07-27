#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "differentialEquations.h"
#include "main.h"
#include "canopyDevelopment.h"
#include "pycnidiospores.h"
#include "ascospores.h"

/*------------------- DEFINING THE DIFFERENTIAL EQUATIONS ----------------------*/
void differentialEquationsHealthyCanopy(double *pVar, double *dpdt, double t)
{
	// Rename the variables for ease of interpretation
	double	H, A;			// Healthy leaf area (GAI)
	
	A = pVar[3 * NGENOTYPES + 2];
	H = pVar[0];
	
	//********************** Healthy area index development **************************************************************
	dpdt[3 * NGENOTYPES + 2] = canopyGrowth(t) * A;					// LAI with just growth, not senescence
	dpdt[0] = canopyGrowth(t) * A - canopySenescence(t) * H;		// Healthy area index; growth and senescence may occur
}

void differentialEquationsFullModel(double *pVar, double *dpdt, double t, unsigned int **geneSummary, double *freqAscospores, double *freqLesions, double *IEAscospores, double *IELesions, double *LPRateAscospores, double *LPRateLesions, double *SporeProdLesions, double *SporeProdAscospores)
{
	/* Naming conventions
	pVar[0]											healthy area index
	pVar[1]-pVar[NGENOTYPES]						latent area index
	pVar[NGENOTYPES + 1] - pVar[2*NGENOTYPES]		infectious area index
	pVar[0]											Total leaf area index

	See 'Genotype order and syntax.txt' for genotype order
	*/
	double	*freqPseudothecia, infectiousSum = 0.0;
	double	probLandOnHealthyTissue = 1.0, probLandOnHost = 1.0, probLesionsMeet = 0.0;
	double	latentAscosporeTissues = 0.0, latentPycnidiosporeTissues = 0.0, infectiousTissues = 0.0;
	double	*mutationFrequencies, totalAscosporeInfection = 0.0, totalPycnidiosporeInfection = 0.0;
	double	totalSenescence = 0.0, totalLesionDeath = 0.0;
	unsigned int	var, nGeno;

	freqPseudothecia = (double *)calloc(NGENOTYPES, sizeof(double));

	// Calculate the genotype specific mutation frequencies
	mutationFrequencies = (double *)calloc(NGENOTYPES, sizeof(double));
	mutation(SporeProdLesions, geneSummary, mutationFrequencies, pVar);

	// Calc probability of landing on a host plant (rather than the ground)
	probLandOnHost = 1 - exp(-pVar[3 * NGENOTYPES + 2]);

	// Calc probability of landing on a healthy leaf
	probLandOnHealthyTissue = pVar[0] / pVar[3 * NGENOTYPES + 2];

	// Calc total reduction in healthy area due to infection by ascospores
	for (var = 1; var <= NGENOTYPES; ++var){
		totalAscosporeInfection += ascosporeDensity(t) * freqAscospores[var - 1] * probLandOnHost * probLandOnHealthyTissue * IEAscospores[var - 1];
	}

	// Calc total reduction in healthy area due to infection by pycnidiospores
	for (var = NGENOTYPES + 1; var <= 2 * NGENOTYPES; ++var){
		totalPycnidiosporeInfection += probLandOnHost * probLandOnHealthyTissue * IELesions[var - 1 - NGENOTYPES] * mutationFrequencies[var - 1 - NGENOTYPES];
	}

	//********************** Healthy area index development **********************************************************
	dpdt[0] = canopyGrowth(t) * pVar[3 * NGENOTYPES + 2] - canopySenescence(t) * pVar[0] - totalAscosporeInfection  * pVar[0] - totalPycnidiosporeInfection  * pVar[0]; 		// Healthy area index; growth and senescence may occur

	//******* Latent tissues arising from infection by ascospores ****************************************************
	for (var = 1; var <= NGENOTYPES; ++var){
		dpdt[var] = ascosporeDensity(t) * freqAscospores[var - 1] * probLandOnHost * probLandOnHealthyTissue * IEAscospores[var - 1] * pVar[0] - LPRateAscospores[var - 1] * pVar[var] - canopySenescence(t) * pVar[var];
	}

	//******* Latent tissues arising from infection by pycnidiospores ************************************************
	for (var = NGENOTYPES + 1; var <= 2 * NGENOTYPES; ++var){
		dpdt[var] = probLandOnHost * probLandOnHealthyTissue * IELesions[var - 1 - NGENOTYPES] * mutationFrequencies[var - 1 - NGENOTYPES] * pVar[0] - LPRateLesions[var - 1 - NGENOTYPES] * pVar[var] - canopySenescence(t) * pVar[var];
	}

	//******* Infectious tissues *************************************************************************************
	for (var = 2 * NGENOTYPES + 1; var <= 3 * NGENOTYPES; ++var){
		dpdt[var] = LPRateAscospores[var - 1 - 2 * NGENOTYPES] * pVar[var - 2 * NGENOTYPES] + LPRateLesions[var - 1 - 2 * NGENOTYPES] * pVar[var - NGENOTYPES] - RATE_INFECTIOUS_PERIOD * pVar[var];
	}

	//******** Dead infectious tissues *******************************************************************************
	for (var = 1; var <= 2 * NGENOTYPES; ++var){
		totalSenescence += canopySenescence(t) * pVar[var];
	}
	for (var = 2 * NGENOTYPES + 1; var <= 3 * NGENOTYPES; ++var){
		totalLesionDeath += RATE_INFECTIOUS_PERIOD  * pVar[var];
	}
	dpdt[3 * NGENOTYPES + 1] = canopySenescence(t) * pVar[0] + totalSenescence + totalLesionDeath;
	
	//******** Total leaf area index *********************************************************************************
	dpdt[3 * NGENOTYPES + 2] = canopyGrowth(t) * pVar[3 * NGENOTYPES + 2];


	//******** BETWEEN SEASON DYNAMICS *******************************************************************************
	//******** ASEXUAL REPRODUCTION **********************************************************************************
	//******** Pycnidiospore production for between season transfer **************************************************
	// V(t) * I_k(t) * rho_pyc_k
	for (var = 3 * NGENOTYPES + 3; var <= 4 * NGENOTYPES + 2; ++var){
		dpdt[var] = pseudotheciaViability(t) * mutationFrequencies[var - 3 * NGENOTYPES - 3];
	}

	//******** SEXUAL REPRODUCTION ***********************************************************************************
	//******** Calculate lesion frequencies **************************************************************************
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
		infectiousSum += pVar[2 * NGENOTYPES + 1 + nGeno];
	}
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
		if (infectiousSum == 0){
			freqLesions[nGeno] = 0;
		}
		else{
			freqLesions[nGeno] = pVar[2 * NGENOTYPES + 1 + nGeno] / infectiousSum; // area under V(t) I_k(t) rho_k(t) divided by area for all genotypes combined
		}
	}
	
	//******** Calculate pseudothecia frequencies *********************************************************************
	// Apply mating
	matingEvent(pVar, freqLesions, freqPseudothecia, geneSummary);

	//******** Ascospore production for between season transfer ******************************************************
	// V(t) * q_k,x(t) * rho_asc_k
	probLesionsMeet = calcProbLesionsMeet(pVar);
	for (var = 4 * NGENOTYPES + 3; var <= 5 * NGENOTYPES + 2; ++var){
		dpdt[var] = pseudotheciaViability(t) * probLesionsMeet * freqPseudothecia[var - 4 * NGENOTYPES - 3] * SporeProdAscospores[var - 3 - 4 * NGENOTYPES];
	}

	free(freqPseudothecia);
	free(mutationFrequencies);
}

