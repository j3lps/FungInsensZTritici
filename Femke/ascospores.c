#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ascospores.h"
#include "main.h"
#include "simulations.h"

//================ ASCOSPORE DENSITY FUNCTION =========================================================
double ascosporeDensity(double t)
{
	// Calculates the ascospores released per m-2
	double	SHAPE_PARAMETER1 = 0.0007, SHAPE_PARAMETER2 = 0.01;
	double	ascospores;

	ascospores = SHAPE_PARAMETER1 * t * t * exp(- SHAPE_PARAMETER2 * t);
	
	// This should turn off the epidemic:
	//ascospores = 0.0;

	return ascospores;
}

//================ SEXUAL REPRODUCTION ================================================================

void createOneGenePunnettSquare()
{
	/* Here we define the probability of a cross between two haploid parents leading to offspring 
	of a specific genotype for a one gene case.
	The resultant information will be stored in the global 3D array 'OneGenePunnettSquare'.
	The resultant punnett square should be
	{
		A offspring
		{
				A		a
			A	{1.0	0.5},
			a	{0.5	0.0}
		},
		a offspring
		{
				A		a
			A	{0.0	0.5},
			a	{0.5	1.0}
		}
	} */

	unsigned int	outGenotype, firstGenotype, secondGenotype;
	
	// Loop through each genotype, whereby A = 0 and a = 1
	for(outGenotype = 0; outGenotype < 2; ++outGenotype){
		// Loop through the first parent's genotype
		for(firstGenotype = 0; firstGenotype < 2; ++firstGenotype){
			//Loop through the second parent's genotype
			for(secondGenotype = 0; secondGenotype < 2; ++secondGenotype){
				if(outGenotype == firstGenotype){
					oneGenePunnettSquare[outGenotype][firstGenotype][secondGenotype] += 0.5;
				}
				if(outGenotype == secondGenotype){			
					oneGenePunnettSquare[outGenotype][firstGenotype][secondGenotype] += 0.5;
				}
			}
		}
	}
}

void matingEvent(double *pVar, double *freqLesions, double *freqPseudothecia, unsigned int **geneSummary)
{
	unsigned int nGene, offspringGeno, parent1Geno, parent2Geno;
	unsigned int offspringAllele, parent1Allele, parent2Allele;
	double	propSingleMating, propAllMatings;

	for (offspringGeno = 0; offspringGeno < NGENOTYPES; ++offspringGeno){ 
		// Stores the frequency of a specific genotype within the ascospores
		propAllMatings = 0.0;

		/* Loop through all genotype crosses and work out the proportion of each that gives offspringGeno. 
		Add to propAllMatings. */
		for (parent1Geno = 0; parent1Geno < NGENOTYPES; ++parent1Geno){
			for (parent2Geno = 0; parent2Geno < NGENOTYPES; ++parent2Geno){
				// Stores the proportion of the offspring genotype resulting from this single mating event
				propSingleMating = 1.0;
				for (nGene = 0; nGene < NGENES; ++nGene){
					/* Determine whether the gene confers for fungicide sensitivity / avirulence (0) or
					fungicide resistance / virulence (1) within both the offspring and the parents */
					offspringAllele = geneSummary[offspringGeno][nGene];
					parent1Allele = geneSummary[parent1Geno][nGene];
					parent2Allele = geneSummary[parent2Geno][nGene];
					/* For each individual gene work out the probability of the Allele of the offspring occuring 
					after a mating between parent1 and parent2. Multiply the probabilities for all genes, which
					assumes that all genes are sufficiently distant from each other. The probabilities can be 
					read out from the one gene punnett square. */
					propSingleMating *= oneGenePunnettSquare[offspringAllele][parent1Allele][parent2Allele];
				}
				propAllMatings += (propSingleMating * freqLesions[parent1Geno] * freqLesions[parent2Geno]);
			}
		}
		freqPseudothecia[offspringGeno] = propAllMatings;
	}
}

double pseudotheciaViability(double t)
{
	/* This function calculates the probability that pseudothecia produced at time t release
	viable ascospores during the next host growing season.
	The function assumes that pseudothecia are more likely to remain viable until the next 
	host growing season if they have been produced later during the host growing season. 
	It is assumed that pseudothecia produced before anthesis are very unlikely to produce 
	viable ascospores at the start of the next growing season, whereas it is highly likely 
	that pseudothecia produced just before the end of the growing season will produce 
	viable ascospores at the start of the next growing season.*/
	
	double viability;
	double VIABILITY_MIDPOINT = 2600, VIABILITY_RATE = 100.0;

	viability = 1 / (1 + exp( -(t - VIABILITY_MIDPOINT) / VIABILITY_RATE ) );
//	viability = 1;

	return viability;
}

void calcAscosporeSums(double time, double *pCur, double *pNew, double *ascosporeSum)
{
	// Calculates the area under Int (V(t) I_k(t) dt) between time 0 and end of simulation see numerator of eqn 7
	// infectiousSum[NGENOTYPES] stores denominator of eqn 7
	unsigned int nGeno;
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno)
	{
		ascosporeSum[nGeno] += 0.5 * STEPSIZE * (pCur[4 * NGENOTYPES + 3 + nGeno] + pNew[4 * NGENOTYPES + 3 + nGeno]);
	}
}

void sexualReproductionBetweenSeasons(double *freqAscospores, double *ascosporeSum)
{
	/* Here we apply sexual reproduction during the season to determine the genotype frequencies of the primary inoculum in
	year x+1 based on the infectious tissue densities and the genotype specific lesion spore production in year x
	q_k_x+1 = int(V(t)*q_k_x(t)*rho_asc_i) / (sum_j=1_N(int(V(t)*q_j_x(t)*rho_asc_j))). */
	unsigned int nGeno;

	// Determine genotype frequencies at the start of the next season
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
		ascosporeSum[NGENOTYPES] += ascosporeSum[nGeno];
	}
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
		freqAscospores[nGeno] = ascosporeSum[nGeno] / ascosporeSum[NGENOTYPES]; // area under V(t) I_k(t) rho_k(t) divided by area for all genotypes combined
	}
}

double calcProbLesionsMeet(double *pVar)
{
	/* Calculates the probability that lesions meet and undergo a mating event. The probability increases if the overall disease 
	severity increases. */
	double probLesionsMeet, x = 0.03;

	probLesionsMeet = 1 - exp(- x * calculateSeverity(pVar));
	
	return probLesionsMeet;
} 
