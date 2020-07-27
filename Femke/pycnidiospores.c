#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <tchar.h>

#include "pycnidiospores.h"
#include "main.h"

void	mutation(double * SporeProdLesions, unsigned int **geneSummary, double *mutationFrequencies, double *pVar)
{
//	double	noMutationProb, MUTATION_PROB_RESISTANCE = 0.0000001, MUTATION_PROB_VIRULENCE = 0.001;
	double	noMutationProb, MUTATION_PROB_RESISTANCE = 0.0000001, MUTATION_PROB_VIRULENCE = 0.0;
	unsigned int offspringGeno, parentGeno, nResMutations, nVirMutations, marker, nGene, offspringAllele, parentAllele;
//	unsigned int pascalElement, coeff;
		
	/* First calculate the probability of not undergoing any mutations, denoted by S in the model description. For 3
	genes S=1-3*omega-3*omega^2-omega^3. Using the Pascal triangle you can work out the terms making up S. 
	
									# mutations from 0 to NGENES
	#genes		1							1		1
				2						1		2		1
				3					1		3		3		1
				4				1		4		6		4		1		==>		S = 1 -4*o -6*o^2 -4*o^3 -o^4
				etc
	The kth coefficient of the nth row of the triangle of Pascal can be calculated as n over k, which is given by
	n! / (k! * (n-k)!). The k elements are numbered from 0 to the total number of genes
	
	If the virulence genes and fungicide resistance genes do not have the same mutation probability this become 
	slighlty more complex. Assume omega is the probability for fungicide resistance mutations and xi is the 
	probability for virulence mutations, then
	S = 1 - (omega - 2 * xi) - (xi^2 - 2 * omega * xi) - omega * xi^2
	*/

	for (offspringGeno = 0; offspringGeno < NGENOTYPES; ++offspringGeno){
		noMutationProb = 1.0;
		for (parentGeno = 0; parentGeno < NGENOTYPES; ++parentGeno){
			// Stores the number of mutations required to get from the parent genotype to the offspring genotype
			nResMutations = 0;
			nVirMutations = 0;
			for (nGene = 0; nGene < N_FUNG_RES_GENES; ++nGene){
				/* Determine whether the gene has undergone mutation and if so add to the mutation count */
				offspringAllele = geneSummary[offspringGeno][nGene];
				parentAllele = geneSummary[parentGeno][nGene];
				if (offspringAllele != parentAllele){
					++nResMutations;
				}
			}
			for (nGene = N_FUNG_RES_GENES; nGene < NGENES; ++nGene){
				/* Determine whether the gene has undergone mutation and if so add to the mutation count */
				offspringAllele = geneSummary[offspringGeno][nGene];
				parentAllele = geneSummary[parentGeno][nGene];
				if (offspringAllele != parentAllele){
					++nVirMutations;
				}
			}
			if ((nResMutations + nVirMutations) != 0){
				mutationFrequencies[offspringGeno] += pow(MUTATION_PROB_RESISTANCE, (double)nResMutations) * pow(MUTATION_PROB_VIRULENCE, (double)nVirMutations) * SporeProdLesions[parentGeno] * pVar[2 * NGENOTYPES + 1 + parentGeno];
				noMutationProb -= pow(MUTATION_PROB_RESISTANCE, (double)nResMutations) * pow(MUTATION_PROB_VIRULENCE, (double)nVirMutations);
			}
			else{
				marker = parentGeno;
			}
		}
		mutationFrequencies[offspringGeno] += noMutationProb * SporeProdLesions[marker] * pVar[2 * NGENOTYPES + 1 + marker];
	}
}
