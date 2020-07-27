#ifndef PYCNIDIOSPORES_H 
#define PYCNIDIOSPORES_H 

void	mutation(double *SporeProdLesions, unsigned int **geneSummary, double *mutationFrequencies, double *pVar);
void	asexualReproductionBetweenSeasons(double *freqAscospores, double *pycnidiosporeSum);
void	calcPycnidiosporeSums(double time, double *pCur, double *pNew, double *pycnidiosporeSum);

#endif 