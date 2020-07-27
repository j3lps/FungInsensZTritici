#ifndef GENETICS_H 
#define GENETICS_H 

void createGenotypes(unsigned int **geneSummary);
void getGenotypicInfectionParameters(double *IEAscospores, double *IELesions, double *LPRateAscospores, double *LPRateLesions, double *SporeProdLesions, double *SporeProdAscospores, unsigned int **geneSummary, double *pVar, double t);

#endif 