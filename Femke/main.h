#ifndef MAIN_H //if MAIN_H is not defined (BLAH_H is just a label for the following code) 
#define MAIN_H //if MAIN_H is not previously defined then define it :)    

#define		STEPSIZE				0.1		// For integrator
#define		T_CROP_EMERGENCE		515.0	// Crop emergence date since sowing at 1 Oct
#define		T_EMERGENCE_L3			1456.0	// GS32; Senescence onset DD for first leaf
#define		T_EMERGENCE_L2			1588.0	// Emergence date of L2; spray time used for dose-response curve parameter fitting
#define		T_FLAG_LEAF_EMERGENCE	1700.0	// DD of full emergence of flag leaf
#define		T_ANTHESIS				2066.0
#define		T_END					3100.0	// DD of end of simulation representing harvest or time at which GAI=0
#define		NGENES					2		// Number of genes associated with virulence and/or fungicide resistance
#define		N_FUNG_RES_GENES		1		// Total number of genes associated with fungicide resistance
#define		N_VIR_GENES				NGENES - N_FUNG_RES_GENES // Total number of genes associated with virulence
#define		NSPRAYS					2		// Number of sprays applied per season
#define		FILE_NAME_BUFFER		500
#define		HAD_LOSS_THRESHOLD		5.0	

// ======================= GENETICS ====================================================================
extern	unsigned int	NGENOTYPES, nDiff;					// Number of genotypes present in the population
extern	double			oneGenePunnettSquare[2][2][2];
extern	double			RATE_INFECTIOUS_PERIOD;
extern	char			fungicideResistance;
extern	char			cultivar;

// ======================== FUNGICIDE TREATMENT ========================================================
extern struct Dosages		// This structure is used to define the dose of individual sprays
{
	double		T1;
	double		T2;
} doseHR;			// A high and a low rsik fungicide are considered

extern struct Sprays		// This structure is used to define whether the sprays are active or not and at what times the sprays will be applied if active
{
	double		T1;
	double		T2;
} timing, leafArea;

extern	double		effectiveDoseHR;		// Effective daily dose concentration
extern	double		HADDiseaseAbsent, HADDiseasePresent, HADLossPercentage; // HAD determined according to healthy and latent area index
extern	double		TotalDiseasedGAI, TotalHealthyGAI, TotalGAILoss;
extern	double		TotalFungEffect;
extern	double		selectionRatioFungicideResistance, selectionRatioVirulence;
extern	double		fungResStartSeason, fungResEndSeason;
extern	double		virulenceStartSeason, virulenceEndSeason;

extern unsigned int effectiveLife, processNumber;

#endif 
