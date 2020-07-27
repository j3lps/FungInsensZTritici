#ifndef		SEXUAL_REPRODUCTION_H
#define		SEXUAL_REPRODUCTION_H

double	ascosporeDensity(double t);
void	createOneGenePunnettSquare();
void	matingEvent(double *pVar, double *freqLesions, double *freqPseudothecia, unsigned int **geneSummary);
double	pseudotheciaViability(double t);
void	calcAscosporeSums(double time, double *pCur, double *pNew, double *ascosporeSum);
void	sexualReproductionBetweenSeasons(double *freqAscospores, double *ascosporeSum);
double	calcProbLesionsMeet(double *pVar);

#endif	/* Header file is included */
