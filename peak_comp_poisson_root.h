#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <iomanip>

// ROOT libraries
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMath.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#define S32K        32768 //2^15
#define NPAR        20
#define NSPECT      32
#define NSIMDATA    16

//forward declarations
void   find_chisqMin(); // ROOT minimization class
void   plotSpectra();   // ROOT canvas
int    readMCA(FILE* inp,char* filename,float inpHist[NSPECT][S32K]);
int    readFMCA(FILE* inp,char* filename,float inpHist[NSPECT][S32K]);

//global variables
float  expHist[NSPECT][S32K];
float  simHist[NSIMDATA][NSPECT][S32K];
double chisq;
double scaledSimHist[NSIMDATA][NSPECT][S32K];
double bgHist[NSPECT][S32K];
float  resultsHist[NSPECT][S32K];

//global flags
int verbosity; //0=only print chisq, 1=print parameters, 2=print debug info

// new ROOT stuff
TApplication *theApp;
TH1D *data[NSPECT],*sim[NSPECT];
double lrchisq(const double *par); // likelihood ratio chisq for spectrum sp
double pchisq(const double *par); // pearson chisq
double nchisq(const double *par); // neyman chisq
double expCurrent[S32K],simCurrent[NSIMDATA][S32K];
int spCurrent;

// fit parameters
double aFinal[NPAR][NSPECT]; // final parameters
