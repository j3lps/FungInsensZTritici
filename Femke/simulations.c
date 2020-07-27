#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <tchar.h>

#include "simulations.h"
#include "main.h"
#include "differentialEquations.h"
#include "standardFunctions.h"
#include "fungicideDynamics.h"
#include "canopyDevelopment.h"
#include "ascospores.h"
#include "pycnidiospores.h"
#include "genetics.h"

void canopyDynamics(double *pPrev, double *pNext, double *dpdt)
{
	double *pTemp, time = 0.0, outputTime = 0.0, outputInterval = 1.0;

	FILE *Hfp = NULL;
	fopen_s(&Hfp, "healthyCanopyDynamics.txt", "w"); // Delete old material stored in the file
	fclose(Hfp);
	fopen_s(&Hfp, "healthyCanopyDynamics.txt", "a"); // Open file to start appending

	HADDiseaseAbsent = 0.0;
	TotalHealthyGAI = 0.0;

	//------ ASSIGN ARRAY LENGTHS FOR VARIABLES AND DIFFERENTIAL EQUATIONS ---
	/*Each array needs to store info for a healthy class, 2 * NGENOTYPE latent classes (lesions resulting from
	ascospore infection have a different LP than those resulting from infection by and pycnidispores) and
	NGENOTYPE infectious classes */
	pPrev = (double *)calloc(nDiff, sizeof(double));
	pNext = (double *)calloc(nDiff, sizeof(double));
	dpdt = (double *)calloc(nDiff, sizeof(double));

	//Initialise canopy areas
	initialiseCanopy(pNext);

	outputTime = time + outputInterval;

	fprintf(Hfp, "Time (DD)\t H\n");
	fprintf(Hfp, "%.0f\t %6f\n", time, pNext[0]);

	while (time < T_END)
	{
		time += STEPSIZE;

		pTemp = pPrev; pPrev = pNext; pNext = pTemp;
		// Calculate the derivates of the healthy canopy
		differentialEquationsHealthyCanopy(pPrev, dpdt, time);
		// Calculate the HAI at the next time point
		rk4(pPrev, pNext, dpdt);
		
		// Update the healthy area duration, if in the grain filling period
		calcHADDiseaseAbsent(time, pPrev, pNext);
		// Update the healthy area duratin over the whole time from crop emergence
		calcTotalHealthyGAI(time, pPrev, pNext);
		
		if (time >= outputTime){
			fprintf(Hfp, "%.0f\t %6f\n", time, pNext[0]);
			outputTime = floor(outputTime) + outputInterval;
		}
	}

	fclose(Hfp);
	free(pPrev);
	free(pNext);
	free(dpdt);
}

void epidemicParameterFit(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores)
{
	// SET MUTATION_PROBABILITY TO 0.0 IN PYCNIDIOSPORES.C BEFORE RUNNING!!!!!!!!!!!!!!!!! 
	// SET CULTIVAR TO 's' (SUSCEPTIBLE)
	// SET FUNGICIDE DOSAGES TO 0
	
	double time = 0.0, *pTemp, outputInterval = 1.0, timeOutput = 0.0, severity = 0.0;
	double *freqLesions, *IEAscospores, *IELesions, *LPRateAscospores, *LPRateLesions, *SporeProdLesions, *SporeProdAscospores;
	unsigned int eqn;
	
	FILE *fpDynamics = NULL;
	fopen_s(&fpDynamics, "EpidemicProgressFitting.txt", "w"); // Delete old material stored in the file
	fclose(fpDynamics);
	fopen_s(&fpDynamics, "EpidemicProgressFitting.txt", "a"); // Open file to start appending

	FILE *fpHAD = NULL;
	fopen_s(&fpHAD, "HADAndSeverityFit.txt", "w"); // Delete old material stored in the file
	fclose(fpHAD);
	fopen_s(&fpHAD, "HADAndSeverityFit.txt", "a"); // Open file to start appending

	// Determine healthy area dynamics in absence of disease
	canopyDynamics(pCur, pNew, diffs);

	HADDiseasePresent = 0.0;
	HADLossPercentage = 0.0;

	//------ ASSIGN ARRAY LENGTHS FOR VARIABLES AND DIFFERENTIAL EQUATIONS ---
	/*Each array needs to store info for a healthy class, 2 * NGENOTYPE latent classes (lesions resulting from
	ascospore infection have a different LP than those resulting from infection by and pycnidispores) and
	NGENOTYPE infectious classes */
	pCur = (double *)calloc(nDiff, sizeof(double));
	pNew = (double *)calloc(nDiff, sizeof(double));
	diffs = (double *)calloc(nDiff, sizeof(double));

	freqLesions = (double *)calloc(NGENOTYPES, sizeof(double));

	//----- CALCULATE GENOTYPE SPECIFIC INFECTION PARAMETERS -----------------
	IEAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	IELesions = (double *)calloc(NGENOTYPES, sizeof(double));
	LPRateAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	LPRateLesions = (double *)calloc(NGENOTYPES, sizeof(double));
	SporeProdLesions = (double *)calloc(NGENOTYPES, sizeof(double));
	SporeProdAscospores = (double *)calloc(NGENOTYPES, sizeof(double));

	//------------------- DEFINE FUNGICIDE TREATMENT PROGRAM ---------------
	applySprayProgram();
	
	//-------------------- SET INITIAL GENOTYPE FREQUENCIES -----------------
	/* The only pathogen strain present is fully susceptibe to fungicide treatment and fully avirulent */
	freqAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
//	freqAscospores[(NGENOTYPES / 2) - 1] = 1.0;
	freqAscospores[0] = 1.0;
	
	//Initialise canopy areas
	initialiseCanopy(pNew);
	
	timeOutput = time + outputInterval;
	
	// Print column labels
	fprintf(fpDynamics, "Time (DD)\tH\t");
	for (eqn = 1; eqn <= NGENOTYPES; ++eqn){
		fprintf(fpDynamics, "latentAsc%i\t", eqn);
	}
	for (eqn = 1; eqn <= NGENOTYPES; ++eqn){
		fprintf(fpDynamics, "latentPyc%i\t", eqn);
	}
	for (eqn = 1; eqn <= NGENOTYPES; ++eqn){
		fprintf(fpDynamics, "infectious%i\t", eqn);
	}
	fprintf(fpDynamics, "Severity\n");

	// Print first value
	fprintf(fpDynamics, "%.0f\t");
	for (eqn = 0; eqn <= 3 * NGENOTYPES; ++eqn){
		fprintf(fpDynamics, "%12f\t", pNew[eqn]);
	}
	if(calculateSeverity(pNew) >= 0.1){
		severity = calculateSeverity(pNew);
	}
	else{
		severity = 0.0;
	}
	fprintf(fpDynamics, "%6f\n", severity);

	while (time < T_END){
		
		time += STEPSIZE;

		pTemp = pCur; pCur = pNew; pNew = pTemp;
		getGenotypicInfectionParameters(IEAscospores, IELesions, LPRateAscospores, LPRateLesions, SporeProdLesions, SporeProdAscospores, geneSummary, pCur, time);
		differentialEquationsFullModel(pCur, diffs, time, geneSummary, freqAscospores, freqLesions, IEAscospores, IELesions, LPRateAscospores, LPRateLesions, SporeProdLesions, SporeProdAscospores);
		rk4(pCur, pNew, diffs);

		calcHADDiseasePresent(time, pCur, pNew);

//		Used for dose-reponse curve parameter fit
		if (time >= T_EMERGENCE_L2 + 532 && time < T_EMERGENCE_L2 + 532 + STEPSIZE)
		{
			if (calculateSeverity(pNew) >= 0.1){
				severity = calculateSeverity(pNew);
			}
			else{
				severity = 0.0;
			}
			fprintf(fpHAD, "Severity for dose-response fit = %6f\n", severity);
		}
		
		if (time >= timeOutput){
			fprintf(fpDynamics, "%.0f\t", time);
			for (eqn = 0; eqn < 3 * NGENOTYPES + 1; ++eqn){
				fprintf(fpDynamics, "%12f\t", pNew[eqn]);
			}
			
			if (calculateSeverity(pNew) >= 0.1){
				severity = calculateSeverity(pNew);
			}
			else{ 
				severity = 0.0;  
			}
			fprintf(fpDynamics, "%6f\n", severity);
			timeOutput = floor(timeOutput) + outputInterval;
		}
	}

	HADLossPercentage = (1 - (HADDiseasePresent / HADDiseaseAbsent)) * 100;
	fprintf(fpHAD, "HAD in absence of disease = %6f\nHAD in presence of disease = %6f\nHAD loss = %.2f%%\n", HADDiseaseAbsent, HADDiseasePresent, HADLossPercentage);

	free(pCur);
	free(pNew);
	free(diffs);
	free(freqAscospores);
	fclose(fpDynamics);
	fclose(fpHAD);
}

void runOneSeason(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores, unsigned int season)
{
	double time = 0.0, *pTemp, outputInterval = 1.0, timeOutput = 0.0, severity = 0.0;
	double *IEAscospores, *IELesions, *LPRateAscospores, *LPRateLesions, *SporeProdLesions, *SporeProdAscospores;
	double *freqLesions, *ascosporeSum, *pycnidiosporeSum;
	unsigned int nGeno, nGene, avirulenceGenes, eqn;
	double *pFungEffect;
	FILE *fpSingleSeason = NULL;

	// Specify whether to write the first season to file every day (1) or not (0)
	unsigned int printDaily = 0;			
	if (processNumber <= 10) printDaily = 1;

	if (printDaily && season == 1) {
		char myFile[FILE_NAME_BUFFER];
		int n = sprintf_s(myFile, FILE_NAME_BUFFER, "Dynamics%d.csv", processNumber);
		fopen_s(&fpSingleSeason, myFile, "w"); // Delete old material stored in the file
		fclose(fpSingleSeason);
		fopen_s(&fpSingleSeason, myFile, "a"); // Open file to start appending
	}
	HADDiseasePresent = 0.0;
	HADLossPercentage = 0.0;
	TotalDiseasedGAI = 0.0;
	TotalGAILoss = 0.0;
	TotalFungEffect = 0.0;
	
	// What is this doing??
	pFungEffect = (double *)calloc(3, sizeof(double));
	pFungEffect[0] = 1; pFungEffect[1] = 1; pFungEffect[2] = 1;

	//------ ASSIGN ARRAY LENGTHS FOR VARIABLES AND DIFFERENTIAL EQUATIONS ---
	/*Each array needs to store info for a healthy class, 2 * NGENOTYPE latent classes (lesions resulting from
	ascospore infection have a different LP than those resulting from infection by and pycnidispores) and
	NGENOTYPE infectious classes */
	pCur = (double *)calloc(nDiff, sizeof(double));
	pNew = (double *)calloc(nDiff, sizeof(double));
	diffs = (double *)calloc(nDiff, sizeof(double));

	//----- ASSIGN ARRAY LENGHTS FOR GENOTYPE FREQUENCIES OF LESIONS ---------- // CHECK, not sure what these variables are doing. Why are they only NGENOTYPES long? Shouldn't it be ngenes?
	freqLesions = (double *)calloc(NGENOTYPES, sizeof(double));
	pycnidiosporeSum = (double *)calloc(NGENOTYPES + 1, sizeof(double)); // Stores area under the genotype specific pycnidiospore production curve
	ascosporeSum = (double *)calloc(NGENOTYPES + 1, sizeof(double)); // Stores area under the genotype specific infectious tissue density curve

	//----- CALCULATE GENOTYPE SPECIFIC INFECTION PARAMETERS -----------------
	IEAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	IELesions = (double *)calloc(NGENOTYPES, sizeof(double));
	LPRateAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	LPRateLesions = (double *)calloc(NGENOTYPES, sizeof(double));
	SporeProdLesions = (double *)calloc(NGENOTYPES, sizeof(double));
	SporeProdAscospores = (double *)calloc(NGENOTYPES, sizeof(double));

	//------------------- DEFINE FUNGICIDE TREATMENT PROGRAM ------------------
	applySprayProgram();

	//Initialise canopy areas
	initialiseCanopy(pNew);

	timeOutput = time + outputInterval;

	// Record the fungicide resistance frequencies at the start of the season
	fungResStartSeason = 0.0;
	for (nGeno = NGENOTYPES / 2; nGeno < NGENOTYPES; ++nGeno){
		fungResStartSeason += freqAscospores[nGeno];
	}

	if (N_VIR_GENES == 0){	// 1 resistance QTL
		virulenceStartSeason = 1;
	}
	else if (N_VIR_GENES == 1){	// 1 resistance QTL
		virulenceStartSeason = freqAscospores[1] + freqAscospores[3];
	}
	else if (N_VIR_GENES == 2){
		virulenceStartSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[5] + freqAscospores[6];
	}
	else if (N_VIR_GENES == 3){
		virulenceStartSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[9] + freqAscospores[10] + freqAscospores[12];
	}
	else if (N_VIR_GENES == 4){
		virulenceStartSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[8] + freqAscospores[17] + freqAscospores[18] + freqAscospores[20] + freqAscospores[24];
	}
	else if (N_VIR_GENES == 5){
		virulenceStartSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[8] + freqAscospores[16] + freqAscospores[33] + freqAscospores[34] + freqAscospores[36] + freqAscospores[40] + freqAscospores[48];
	}
	else{
		printf("Error: the number of cultivar resistance QTLs does not exist");
		exit(10);
	}

	if (printDaily && season == 1) {
		// Print column labels
		fprintf(fpSingleSeason, "Time,H");
		for (eqn = 1; eqn <= NGENOTYPES; ++eqn) {
			fprintf(fpSingleSeason, ", latentAsc%i", eqn);
		}
		for (eqn = 1; eqn <= NGENOTYPES; ++eqn) {
			fprintf(fpSingleSeason, ", latentPyc%i", eqn);
		}
		for (eqn = 1; eqn <= NGENOTYPES; ++eqn) {
			fprintf(fpSingleSeason, ", infectious%i", eqn);
		}
		fprintf(fpSingleSeason, ", Severity\n");

		// Print initial values
		fprintf(fpSingleSeason, "%.0f");
		for (eqn = 0; eqn <= 3 * NGENOTYPES; ++eqn) {
			fprintf(fpSingleSeason, ", %12f", pNew[eqn]);
		}
		if (calculateSeverity(pNew) >= 0.1) {
			severity = calculateSeverity(pNew);
		}
		else {
			severity = 0.0;
		}
		fprintf(fpSingleSeason, ", %6f\n", severity);
	}
	
	// Run a single season
	while (time < T_END){

		time += STEPSIZE;

		// CHECK: What is this?? Shuffling the fungicide effect? This sequence loses the original pFungEffect[2]: c(0,1,2) -> c(1,0,0)
		pFungEffect[2] = pFungEffect[0]; pFungEffect[0] = pFungEffect[1]; pFungEffect[1] = pFungEffect[2];
		pTemp = pCur; pCur = pNew; pNew = pTemp;
		// Calculate the current infection parameters
		getGenotypicInfectionParameters(IEAscospores, IELesions, LPRateAscospores, LPRateLesions, SporeProdLesions, SporeProdAscospores, geneSummary, pCur, time);
		differentialEquationsFullModel(pCur, diffs, time, geneSummary, freqAscospores, freqLesions, IEAscospores, IELesions, LPRateAscospores, LPRateLesions, SporeProdLesions, SporeProdAscospores);
		rk4(pCur, pNew, diffs);
		pFungEffect[1] = sprayEffectHR(time, pNew, 0.0);

		calcHADDiseasePresent(time, pCur, pNew);
		calcTotalDiseasedGAI(time, pCur, pNew);
		calcTotalFungEffect(time, pFungEffect);
		
		calcAscosporeSums(time, pCur, pNew, ascosporeSum);

		if (printDaily && season == 1) {
			if (time >= timeOutput) {
				fprintf(fpSingleSeason, "%.0f", time);
				for (eqn = 0; eqn < 3 * NGENOTYPES + 1; ++eqn) {
					fprintf(fpSingleSeason, ", %12f", pNew[eqn]);
				}

				if (calculateSeverity(pNew) >= 0.1) {
					severity = calculateSeverity(pNew);
				}
				else {
					severity = 0.0;
				}
				fprintf(fpSingleSeason, ", %6f\n", severity);
				timeOutput = floor(timeOutput) + outputInterval;
			}
		}
	}

	printf("%6f\t%6f\n", HADDiseaseAbsent, HADDiseasePresent);
	HADLossPercentage = (1 - (HADDiseasePresent / HADDiseaseAbsent)) * 100;	
	TotalGAILoss = (1 - (TotalDiseasedGAI / TotalHealthyGAI)) * 100;

	// Calculate genotype frequencies at start of the next growing season
	sexualReproductionBetweenSeasons(freqAscospores, ascosporeSum);

	fungResEndSeason = 0.0;
	for (nGeno = NGENOTYPES / 2; nGeno < NGENOTYPES; ++nGeno){
		fungResEndSeason += freqAscospores[nGeno];
	}
	selectionRatioFungicideResistance = fungResEndSeason / fungResStartSeason;

	if (N_VIR_GENES == 0){	// 1 resistance QTL
		virulenceEndSeason = 1;
	}
	else if (N_VIR_GENES == 1){	// 1 resistance QTL
		virulenceEndSeason = freqAscospores[1] + freqAscospores[3];
	}
	else if (N_VIR_GENES == 2){
		virulenceEndSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[5] + freqAscospores[6];
	}
	else if (N_VIR_GENES == 3){
		virulenceEndSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[9] + freqAscospores[10] + freqAscospores[12];
	}
	else if (N_VIR_GENES == 4){
		virulenceEndSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[8] + freqAscospores[17] + freqAscospores[18] + freqAscospores[20] + freqAscospores[24];
	}
	else if (N_VIR_GENES == 5){
		virulenceEndSeason = freqAscospores[1] + freqAscospores[2] + freqAscospores[4] + freqAscospores[8] + freqAscospores[16] + freqAscospores[33] + freqAscospores[34] + freqAscospores[36] + freqAscospores[40] + freqAscospores[48];
	}
	else{
		printf("Error: the number of cultivar resistance QTLs does not exist");
		exit(10);
	}

	free(pCur);
	free(pNew);
	free(diffs);
	free(freqLesions);
	free(ascosporeSum);
	free(pycnidiosporeSum);
	free(IEAscospores);
	free(IELesions);
	free(LPRateAscospores);
	free(LPRateLesions);
	free(SporeProdLesions);
	free(SporeProdAscospores);
	if (printDaily && season == 1) fclose(fpSingleSeason);
}

double calculateSeverity(double *pNew)
{
	double severity = 0.0, totalInfectious = 0.0;
	unsigned int eqn;

	if (pNew[3 * NGENOTYPES + 2] < 0.1){
		severity = 0.0;
	}
	else{
		for (eqn = 2 * NGENOTYPES + 1; eqn <= 3 * NGENOTYPES; ++eqn){
			totalInfectious += pNew[eqn];						// Total density of infectious tissues
		}
		severity = (totalInfectious / pNew[3 * NGENOTYPES + 2]) * 100.0;
	}

	return severity;
}

void calcHADDiseaseAbsent(double time, double *pPrev, double *pNext)
{
	// Calculate Healthy area duration (HAD) over the grain-fill period
	if (time >= T_ANTHESIS){
		HADDiseaseAbsent += 0.5 * STEPSIZE * (pPrev[0] + pNext[0]);
	}
}

void calcHADDiseasePresent(double time, double *pCur, double *pNew)
{
	unsigned int eqn;
	double sumCur, sumNew;

	// Calculate Healthy area duration (HAD) over the grain-fill period
	if (time >= T_ANTHESIS){
		sumCur = 0.0;
		sumNew = 0.0;
		for (eqn = 0; eqn <= 2 * NGENOTYPES; ++eqn){		// healthy and latent tissues contribute to HAD
			sumCur += pCur[eqn];
			sumNew += pNew[eqn];
		}
		HADDiseasePresent += 0.5 * STEPSIZE * (sumCur + sumNew);
	}
}

void calcTotalHealthyGAI(double time, double *pPrev, double *pNext)
{
	// Calculate Healthy area duration (HAD) over the grain-fill period
	if (time >= T_CROP_EMERGENCE){
		TotalHealthyGAI += 0.5 * STEPSIZE * (pPrev[0] + pNext[0]);
	}
}

void calcTotalDiseasedGAI(double time, double *pCur, double *pNew)
{
	unsigned int eqn;
	double sumCur, sumNew;

	// Calculate Healthy area duration (HAD) over the grain-fill period
	if (time >= T_CROP_EMERGENCE){
		sumCur = 0.0;
		sumNew = 0.0;
		for (eqn = 0; eqn <= 2 * NGENOTYPES; ++eqn){		// healthy and latent tissues contribute to HAD
			sumCur += pCur[eqn];
			sumNew += pNew[eqn];
		}
		TotalDiseasedGAI += 0.5 * STEPSIZE * (sumCur + sumNew);
	}
}

void runMultipleSeasonsResistanceOnly(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores)
{
	double time = 0.0, outputInterval = 1.0, timeOutput = 0.0, severity = 0.0;
	unsigned int season = 1, nGeno;
	char buff[1024]; // File buffer

	// Create a file to store the results at the end of each season
	FILE *fpMultipleSeasons = NULL;
	char first[] = "Sim";
	char middle[10];
	sprintf_s(middle, sizeof(middle), "%d", processNumber);
	char end[] = ".csv";
	char fileName[50];
	snprintf(fileName, sizeof(fileName), "%s%s%s", first, middle, end);
	fopen_s(&fpMultipleSeasons, fileName, "w"); // Delete old material stored in the file
	fclose(fpMultipleSeasons);
	fopen_s(&fpMultipleSeasons, fileName, "a"); // Open file to start appending
	setvbuf(stdout, buff, _IOLBF, 1024);

	// Determine HAD in the absence of disease, so that we can calculate HAD loss
	canopyDynamics(pCur, pNew, diffs);

	//-------------------- SET INITIAL GENOTYPE FREQUENCIES -----------------
	// The initial pathogen strain present is fully sensitive to fungicide treatment and 	
	// fully avirulent, representing a susceptible culitvar controlled by fungicides 
	freqAscospores = (double *)calloc(NGENOTYPES, sizeof(double)); 
	freqAscospores[0] = 1.0;
	// If running no pathogen, then can turn off the pathogen here by setting the frequency of all genotypes to zero.
	if (processNumber == 0) freqAscospores[0] = 0.0;

	HADLossPercentage = 0.0;

	// Write the header for the output file
	fprintf(fpMultipleSeasons, "Season, HAD_loss");
	for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno) {
		fprintf(fpMultipleSeasons, ", Freq%d",nGeno);
	}
	fprintf(fpMultipleSeasons, "\n");

	// Run for 100 seasons
	while (season < 101)
	{
		runOneSeason(pNew, pCur, diffs, geneSummary, freqAscospores, season);
		fprintf(fpMultipleSeasons, "%i, %6f", season, HADLossPercentage);
		for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
			fprintf(fpMultipleSeasons, ", %.12f ", freqAscospores[nGeno]);
		}
		fprintf(fpMultipleSeasons, "\n");
		++season;
	}
	
	free(freqAscospores);
	fclose(fpMultipleSeasons);
}

void runMultipleSeasonsVirulenceOnly(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores)
{
	double time = 0.0, outputInterval = 1.0, timeOutput = 0.0, severity = 0.0;
	unsigned int season = 1, nGeno;
	char buff[1024]; // File buffer

	FILE *fpMultipleSeasons = NULL;
	char first[] = "Sim";
	char middle[10];
	sprintf_s(middle,sizeof(middle),"%d",processNumber);
	char end[] = ".txt";
	char fileName[50]; 
	snprintf(fileName,sizeof(fileName),"%s%s%s",first,middle,end);
	fopen_s(&fpMultipleSeasons, fileName, "w"); // Delete old material stored in the file
	fclose(fpMultipleSeasons);
	fopen_s(&fpMultipleSeasons, fileName, "a"); // Open file to start appending
	setvbuf(stdout, buff, _IOLBF, 1024);

	// Determine HAD in the absence of disease
	canopyDynamics(pCur, pNew, diffs);

	//-------------------- SET INITIAL GENOTYPE FREQUENCIES -----------------
	// The initial pathogen strain present is fully sensitive to fungicide treatment and 	
	// but no fungicide is applied and the pathogen is fully avirulent
	freqAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	freqAscospores[0] = 1.0; // 1.0

	HADLossPercentage = 0.0;

	fprintf(fpMultipleSeasons, "Season\t%% HAD loss\tSelection ratio\n");

	//	while (HADLossPercentage < HAD_LOSS_THRESHOLD)
	while (season < 101)
	{
		runOneSeason(pNew, pCur, diffs, geneSummary, freqAscospores, season);
		fprintf(fpMultipleSeasons, "%i\t%6f\t%6f\t", season, HADLossPercentage, selectionRatioVirulence);
		for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
			fprintf(fpMultipleSeasons, "%.12f\t", freqAscospores[nGeno]);
		}
		fprintf(fpMultipleSeasons, "\n");
		++season;
	}

	effectiveLife = season - 2;
	fprintf(fpMultipleSeasons, "Effective life = %i\n", effectiveLife);

	free(freqAscospores);
	fclose(fpMultipleSeasons);
}

void runMultipleSeasonsIntegratedControl(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores)
{
	double time = 0.0, outputInterval = 1.0, timeOutput = 0.0, severity = 0.0, total_sens_freq = 0.0, total_avir_freq = 0.0;
	unsigned int season = 1, nGeno, fung_replacement_number = 0, cv_replacement_number = 0, i, j;
	char buff[1024]; // File buffer

	FILE *fpMultipleSeasons = NULL;
	fopen_s(&fpMultipleSeasons, "1QTL_1_dose_tau_0_7Case3.txt", "w"); // Delete old material stored in the file
	fclose(fpMultipleSeasons);
	fopen_s(&fpMultipleSeasons, "1QTL_1_dose_tau_0_7Case3.txt", "a"); // Open file to start appending
	setvbuf(stdout, buff, _IOLBF, 1024);

	// Determine HAD in the absence of disease
	canopyDynamics(pCur, pNew, diffs);

	//-------------------- SET INITIAL GENOTYPE FREQUENCIES -----------------
	// The initial pathogen strain present is  fully sensitive to fungicide treatment and 	
	// fully avirulent
	freqAscospores = (double *)calloc(NGENOTYPES, sizeof(double));
	freqAscospores[0] = 1.0; 

	HADLossPercentage = 0.0;

	fprintf(fpMultipleSeasons, "season\t%% HAD loss\tSelection ratio fung res\tSelection ratio virulence\tTotal %% GAI loss\tTotal Fung Effect\n");

//	while (HADLossPercentage < HAD_LOSS_THRESHOLD)
	while (season < 51)
	{
		runOneSeason(pNew, pCur, diffs, geneSummary, freqAscospores, season);
		fprintf(fpMultipleSeasons, "%i\t%6f\t%6f\t%6f\t%6f\t%6f\t", season, HADLossPercentage, selectionRatioFungicideResistance, selectionRatioVirulence, TotalGAILoss, TotalFungEffect);
		for (nGeno = 0; nGeno < NGENOTYPES; ++nGeno){
			fprintf(fpMultipleSeasons, "%.12f\t", freqAscospores[nGeno]);
		}
		fprintf(fpMultipleSeasons, "\n");
		++season;

/*		if (fung_replacement_number < 2){
			for (i = 0; i < NGENOTYPES / 2; ++i){
				total_sens_freq += freqAscospores[i];
			}
			if (total_sens_freq < 0.95){
				for (j = 0; j < NGENOTYPES / 2; ++j){
					freqAscospores[j] += freqAscospores[(NGENOTYPES / 2) + j];
					freqAscospores[(NGENOTYPES / 2) + j] = 0;
				}
				++fung_replacement_number;
			}
			total_sens_freq = 0.0;
		}

		if (cv_replacement_number < 2){
			total_avir_freq  = freqAscospores[0] + freqAscospores[NGENOTYPES / 2];
			if (total_avir_freq < 0.95){
				for (j = 1; j < NGENOTYPES / 2; ++j){
					freqAscospores[0] += freqAscospores[j];
					freqAscospores[NGENOTYPES / 2] += freqAscospores[(NGENOTYPES / 2) + j];
					freqAscospores[j] = 0;
					freqAscospores[(NGENOTYPES / 2) + j] = 0;
				}
				++cv_replacement_number;
			}
			total_avir_freq = 0.0;
		}
*/
	}

	effectiveLife = season - 2;
	fprintf(fpMultipleSeasons, "Effective life = %i\n", effectiveLife);

	free(freqAscospores);
	fclose(fpMultipleSeasons);
}

void calcTotalFungEffect(double time, double *pFungEffect)
{
	// Calculate Healthy area duration (HAD) over the grain-fill period
	if (time >= T_CROP_EMERGENCE) {
		TotalFungEffect += 0.5 * STEPSIZE * ((1 - pFungEffect[0]) + (1 - pFungEffect[1]));
	}
}
