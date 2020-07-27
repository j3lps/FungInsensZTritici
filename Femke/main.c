#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <tchar.h>

#include "main.h"
#include "canopyDevelopment.h"
#include "differentialEquations.h"
#include "ascospores.h"
#include "genetics.h"
#include "fungicideDynamics.h"
#include "simulations.h"

/*--------------------------------- PARAMETER DECLARATION -----------------------------------------------------------*/
double	*pCur, *pNew, *diffs;			// Stores old and new variable values and differential equations
double	*freqLesions, *freqAscospores;	// Stores genotype frequencies of lesions and ascospores
unsigned int	**geneSummary;			// 2D array pointer storing for each genotype which genes are virulent/ resistant vs. avirulent/ sensitive

unsigned int	NGENOTYPES, nDiff, i, j;	// Number of genotypes and diff equations present in the population
double	oneGenePunnettSquare[2][2][2] = { { {0.0, 0.0} , {0.0, 0.0} }, { {0.0, 0.0} , {0.0, 0.0} } };
double	RATE_INFECTIOUS_PERIOD;
char			cultivar = 'r';				// The cultivar is either fully susceptible (s) or carries some resistance (r)

struct	Dosages	doseHR;				// Stores the dose of each high risk fungicide spray
struct	Dosages	doseHR2nd;         // stores the dose of the second fungicde, the mixing partner. 
struct	Sprays	timing;
struct	Sprays	leafArea;			// Stores the leaf area at the time of a given spray application
double	effectiveDoseHR;			// Effective daily dose concentration

double	HADDiseaseAbsent, HADDiseasePresent, HADLossPercentage; // HAD determined according to healthy and latent area index
double	TotalDiseasedGAI, TotalHealthyGAI, TotalGAILoss;
double	TotalFungEffect;
double  Total2ndFungEffect; // create a second fungicide 
double	fungResStartSeason, fungResEndSeason;
double	virulenceStartSeason, virulenceEndSeason;
double  selectionRatioFungicideResistance, selectionRatioVirulence;
double  timeSpent;

unsigned int	effectiveLife, nGene, processNumber;

// This program runs simulations of fungicide resistance developing in Zymoseptoria tritici
// 26 simulations are coded, dependent on the command line argument
// 0: No pathogen
// 1:5 : R3 : R7, no fungicide
// 6:10 : R3 : R7, with fungicide, absolute resistance
// 11:15 : R3 : R7, with 2 fungicides, absolute resistance
// 16:20 : R3 : R7, with fungicide, partial resistance
// 21:25 : R3 : R7, with 2 fungicides, partial resistance

int main(int argc, char* argv[])
{
	// This is used to track the time of a simulation
	clock_t begin, end;
	begin = clock();
	
	if (argc == 2) {
		processNumber = atoi(argv[1]);
		printf("Process number is %d \n", processNumber);
	}
	else processNumber = 0;

	//-------- DETERMINE NUMBER OF GENOTYPES PRESENT IN THE POPULATION ------
	NGENOTYPES = integerPower(2, NGENES);	// The number of genotypes is 2^NGENES in a haploid pathogen
	nDiff = 5 * NGENOTYPES + 3;				// The number of differential equations: healthy tissues; NGENOTYPE classes for latent lesions from ascospores, latent lesions from pycnidiospores, infectious tissues and senesced infectious tissues; total leaf area index.
	
	//------ ASSIGN 2D ARRAY LENGTH FOR GENOTYPES -----------------------------
	// Allocate a one-dimensional array of pointers to double
	geneSummary = (unsigned int **) calloc(NGENOTYPES,sizeof(unsigned int*));
	// For each of the pointers to double, allocate an array of doubles
	for(i = 0; i < NGENOTYPES; ++i){
		geneSummary[i] = (unsigned int *)calloc(NGENES, sizeof(unsigned int));
	}
	// For every gene, and genotype, state whether genotype x for gene y is wild-type (0) or resistant (1)
	createGenotypes(geneSummary);
	// You can now access gene 3 of genotype 30 using geneSummary[30][3]
	
	// Create a punnet square for a cross between genotypes from one gene
	createOneGenePunnettSquare();

	// Run the simulation
	runMultipleSeasonsResistanceOnly(pNew, pCur, diffs, geneSummary, freqAscospores);

	free(geneSummary);

	end = clock();
	timeSpent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("run time in minutes = %.2f\n", timeSpent/60.0);

	return(0); /* Add a breakpoint so that the console screen stays up until you continue debugging with F5. 
	This means that you can read whatever you decide to	print to the screen */
}
