#ifndef SIMULATIONS_H 
#define SIMULATIONS_H 

void	canopyDynamics(double *pPrev, double *pNext, double *dpdt);
void	epidemicParameterFit(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores);
void	runOneSeason(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores, unsigned int season);
double	calculateSeverity(double *pNew);
void	calcHADDiseaseAbsent(double time, double *pPrev, double *pNext);
void	calcHADDiseasePresent(double time, double *pCur, double *pNew);
void	calcTotalHealthyGAI(double time, double *pPrev, double *pNext);
void	calcTotalDiseasedGAI(double time, double *pPrev, double *pNext);
void	runMultipleSeasonsResistanceOnly(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores);
void	runMultipleSeasonsVirulenceOnly(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores);
void	runMultipleSeasonsIntegratedControl(double *pNew, double *pCur, double *diffs, unsigned int **geneSummary, double *freqAscospores);
void	calcTotalFungEffect(double time, double *pFungEffect);

#endif