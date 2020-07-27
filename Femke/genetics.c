#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <tchar.h>

#include "genetics.h"
#include "main.h"
#include "standardFunctions.h"
#include "fungicideDynamics.h"

void createGenotypes(unsigned int** geneSummary)
{
	/* In all arrays referring to genotypes the genotypes are ordered in the same systematic order, i.e. elements 0
	of any genotype array always refers to the pathogen genotype which is fully sensitive to fungicide treatment and
	which is recognised by a resistant host. In a 4 gene case this is represented by ABCD.
	0 = fungicide sensitive; avirulent (elicitors produced and recognised by the resistant host)
	1 = fungicide resistant; virulent
	Note that within each genotype, the first N_FUNG_RES_GENES are associated with fungicide resistance, whereas the
	remainder of the NGENES, N_VIR_GENES, are associated with virulence

	The complete order is given by :
	0		ABCD		0000	fully sensitive to fungicide treatment and recognised by a resistant host
	1		ABCd		0001
	2		ABcD		0010
	3		ABcd		0011
	4		AbCD		0100
	5		AbCd		0101
	6		AbcD		0110
	7		Abcd		0111
	8		aBCD		1000
	9		aBCd		1001
	10		aBcD		1010
	11		aBcd		1011
	12		abCD		1100
	13		abCd		1101
	14		abcD		1110
	15		abcd		1111	fully resistant to fungicide and fully virulent	*/

	unsigned int nGene, nGeno, repeats = 1, repeatsCount, counter1, counter2, counter3;

	for (nGene = 0; nGene < NGENES; ++nGene) {
		counter1 = 0;
		counter2 = NGENOTYPES / (pow(2.0, (double)(nGene + 1)));
		counter3 = counter2;
		repeatsCount = 0;
		for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno) {
			if (repeatsCount < repeats) {
				if (nGeno >= counter1 && nGeno < counter2) {
					geneSummary[nGeno][nGene] = 0;
				}
				else if (nGeno >= counter2 && nGeno < counter2 + counter3) {
					geneSummary[nGeno][nGene] = 1;
				}
				else {
					counter1 = counter2 + counter3;
					counter2 = counter1 + counter3;
					++repeatsCount;
					--nGeno;
				}
			}
		}
		repeats *= 2;
	}
}

void getGenotypicInfectionParameters(double* IEAscospores, double* IELesions, double* LPRateAscospores, double* LPRateLesions, double* SporeProdLesions, double* SporeProdAscospores, unsigned int** geneSummary, double* pVar, double t)
{
	/* This function determines the infection efficiency, spore production, and latent period of each genotype,
	given the cultivar and fungicides present in the simulation.
	First, loop through genotypes and determine the number of genes conferring fungicide sensitivity and the number of
	genes conferring avirulence as they reduce the base infection efficiency of the pathogen.
	Note that the first	N_FUNG_RES_GENES are associated with fungicide resistance whereas the rest N_VIR_GENES are associated with
	virulence. */

	// These base infection parameters are on a sensitive cultivar with no fungicide
	double	BASE_IE_ASC = 0.02, BASE_LP_RATE_ASC = 1.0 / 536.0;
	double	BASE_IE_LES = 0.02, BASE_LP_RATE_LES = 1.0 / 473, BASE_SPORE_PROD_LESIONS = 0.2, BASE_SPORE_PROD_ASCOSPORES = 0.001 * 0.2;
	
	// AVIRULENCE_EFFECT determines cultivar resistance, SENSITIVITY_EFFECT determines effect of partial fungicide resistance, COST_FUNG_RES is the cost of fungicide resistance
	// NB: SENSITIVITY_EFFECT does not change the efficacy of the fungicide, that is changed by MAX_REDUCTION_HR in fungicideDynamics.c
	// but changing AVIRULENCE_EFFECT *does* change the efficacy of the cultivar.
	double	AVIRULENCE_EFFECT = 1.0, SENSITIVITY_EFFECT = 0.4, COST_FUNG_RES = 0.015;

	// CHECK: Can remove N_VIR_GENES, as we're not doing selection. Just need to update parameters based on the RR of the cultivar.

	// Relating Tau (AVIRULENCE_EFFECT) to resistance rating, for 1 QTL. 
	// RR = 3, 4, 5, 6, 7
	// Tau = 1, 0.96, 0.94, 0.93, 0.91 
	unsigned int RR;
	if (processNumber == 0) {
		RR = 0;
	}
	else RR = (processNumber - 1) % 5;
	switch (RR) {
	case 0:
		AVIRULENCE_EFFECT = 1.0;
		break;
	case 1:
		AVIRULENCE_EFFECT = 0.96;
		break;
	case 2:
		AVIRULENCE_EFFECT = 0.94;
		break;
	case 3:
		AVIRULENCE_EFFECT = 0.93;
		break;
	case 4:
		AVIRULENCE_EFFECT = 0.91;
		break;
	}

	unsigned int nGeno, nGene, fungicideSensitiveGenes, avirulenceGenes, fungicideResistanceGenes, virulenceGenes;
	double	sensitivityEffect, sensitivityEffect2nd, TotalFungicideEffect, Total2ndFungicideEffect, fitnessCostFungRes, fitnessCostVirulence;
	double	baseLPRateLes, baseLPRateAsc, baseSporeProdLes, baseSporeProdAsc;
	double	max_reduction;

	// Initialise external parameters
	RATE_INFECTIOUS_PERIOD = 1.0 / 456;

	// Set the base latent period and spore production for each genotype
	baseLPRateAsc = BASE_LP_RATE_ASC;
	baseLPRateLes = BASE_LP_RATE_LES;
	baseSporeProdAsc = BASE_SPORE_PROD_ASCOSPORES;
	baseSporeProdLes = BASE_SPORE_PROD_LESIONS;
	// Loop over each genotype and update the spore production and latent period as necessary
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno) {
		
		// Work out whether this genotype is sensitive or insensitive to the fungicide
		fungicideSensitiveGenes = 0;
		avirulenceGenes = 0;
		for (nGene = 0; nGene < N_FUNG_RES_GENES; ++nGene) {
			if (geneSummary[nGeno][nGene] == 0) {
				++fungicideSensitiveGenes; // Determine the number genes conferring fungicide sensitivity
			}
		}
		for (nGene = N_FUNG_RES_GENES; nGene < NGENES; ++nGene) {
			if (geneSummary[nGeno][nGene] == 0) {
				++avirulenceGenes; // Determine the number genes conferring avirulence
			}
		}

		// If fungicides are applied, work out the effect of each fungicide
		if (processNumber >= 6) {
			// Work out the effect of the fungicide(s)
			sensitivityEffect = SENSITIVITY_EFFECT;
			// This pathogen genotype has a fungicide resistance gene
			if (fungicideSensitiveGenes == 0) {
				// If processNumber < 16 then the fungicide insensitive strain is absolutely resistant; the fungicide has no effect
				if (processNumber < 16) {
					// 1.0 implies the fungicide has no effect.
					TotalFungicideEffect = 1.0;
				}
				else {
					// For partial resistance the following function works out the effect of the fungicide
					TotalFungicideEffect = sprayEffectPartRes(t, pVar, sensitivityEffect);
				}
			}
			// This pathogen strain is fully sensitive to the fungicide
			else if (fungicideSensitiveGenes == 1) {
				// Total reduction in life cycle para due to fungicide treatment - HR = high-risk
				TotalFungicideEffect = sprayEffectHR(t, pVar, 0.0);
			}
			else {
				printf("You cannot be here!");
				exit(1);
			}
		}
		else {
			// Processes 0-5 have no fungicide control applied
			TotalFungicideEffect = 1.0;
		}

		// Work out whether this strain has a fitness cost
		fungicideResistanceGenes = N_FUNG_RES_GENES - fungicideSensitiveGenes;
		fitnessCostFungRes = 1 - COST_FUNG_RES * fungicideResistanceGenes;

		// Determine whether a mixing partner is being added. 
		if ((processNumber >= 11 && processNumber <= 15) || (processNumber >= 21 && processNumber <= 25)) {
			// Calculate the effect of the 2nd fungicide at time t
			Total2ndFungicideEffect = sprayEffectHR2nd(t, pVar, 0.0);
		}
		else {
			// For simulations without a mixing partner, set Total2ndFungicideEffect = 1.0; no effect
			Total2ndFungicideEffect = 1.0;
		}

		// QoI fungicide affects IE, LP and SP, while protectant (mixing partner) only affects IE
		IEAscospores[nGeno] = BASE_IE_ASC * AVIRULENCE_EFFECT * TotalFungicideEffect * Total2ndFungicideEffect * fitnessCostFungRes;
		LPRateAscospores[nGeno] = baseLPRateAsc * AVIRULENCE_EFFECT * TotalFungicideEffect * fitnessCostFungRes;
		SporeProdAscospores[nGeno] = baseSporeProdAsc * AVIRULENCE_EFFECT * TotalFungicideEffect * fitnessCostFungRes;
		IELesions[nGeno] = BASE_IE_LES * AVIRULENCE_EFFECT * TotalFungicideEffect * Total2ndFungicideEffect * fitnessCostFungRes;
		LPRateLesions[nGeno] = baseLPRateLes * AVIRULENCE_EFFECT * TotalFungicideEffect * fitnessCostFungRes;
		SporeProdLesions[nGeno] = baseSporeProdLes * AVIRULENCE_EFFECT * TotalFungicideEffect * fitnessCostFungRes;

	}
}
