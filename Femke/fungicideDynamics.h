#ifndef FUNGICIDE_DYNAMICS_H //if MAIN_H is not defined (BLAH_H is just a label for the following code) 
#define FUNGICIDE_DYNAMICS_H //if MAIN_H is not previously defined then define it :)    

void sprayDosages();
void sprayTiming();
void initialiseLeafAreasAtTimeOfSpray();
void applySprayProgram();
double sprayDoseHR(double t, double *pVar);
double sprayEffectHR(double t, double *pVar, double sensitivityEffect);
double sprayEffectHR2nd(double t, double *pVar, double sensitivityEffect);
double sprayEffectPartRes(double t, double *pVar, double sensitivityEffect);

#endif 