#ifndef		DIFFERENTIAL_EQUATIONS_H
#define		DIFFERENTIAL_EQUATIONS_H

void	differentialEquationsHealthyCanopy(double *pVar, double *dpdt, double t);
void	differentialEquationsFullModel(double *pVar, double *dpdt, double t, unsigned int **geneSummary, double *freqAscospores, double *freqLesions, double *IEAscospores, double *IELesions, double *LPRateAscospores, double *LPRateLesions, double *SporeProdLesions, double *SporeProdAscospores);

#endif	/* Header file is included */
