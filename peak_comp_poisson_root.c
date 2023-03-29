#include "peak_comp_poisson_root.h"
#include "read_config.c"

int main(int argc, char *argv[]){

  int i,j,k,m;
  FILE *expData, *simData, *results, *scalingFactors;

  if(argc != 2){
    printf("\npeak_comp_poisson parameter_file\n");
    printf("Compares the .mca spectra designated in the parameter file "
           "specified and generates cool statistics.\n\n");
    exit(-1);
  }

  // initialize values
  addBackground = 0;
  for(i = 0; i < NSPECT; i++){
    for(j = 0; j < S32K; j++){
      expHist[i][j] = 0.;
      for(k = 0; k < NSIMDATA; k++){
        simHist[k][i][j] = 0.;
        scaledSimHist[k][i][j] = 0.;
      }
    }
  }

  for(i = 0; i < NPAR; i++){
    for(j = 0; j < NSPECT; j++){
      aFinal[i][j] = 0.;
    }
  }

  // initialize ROOT stuff
  char simName[132];
  char dataName[132];
  for(i = 0; i < NSPECT; i++){
    sprintf(simName, "sim_%2d", i);
    sim[i] = new TH1D(simName, ";;", S32K, 0, S32K - 1);
  }
  for(i = 0; i < NSPECT; i++){
    sprintf(dataName, "data_%2d", i);
    data[i] = new TH1D(dataName, ";;", S32K, 0, S32K - 1);
  }

  readConfigFile(argv[1]); // grab data from the config file

  // check that the number of spectra being compared is fine
  if(numSpectra >= NSPECT){
    printf("ERROR: A spectrum number specified in the parameter file is larger "
           "than the maximum value of %i.  Reduce it or increase NSPECT in "
           "peak_comp.h and recompile.\n",
           NSPECT);
    exit(-1);
  }

  // read in the .mca files
  // check file extension of exp data and copy into float histogram
  if((expData = fopen(expDataName, "r")) == NULL){
    printf("ERROR: Cannot open the experiment data file %s!\n", expDataName);
    exit(-1);
  }
  const char *dot = strrchr(expDataName, '.'); // get the file extension
  if(strcmp(dot + 1, "mca") == 0){
    readMCA(expData, expDataName, expHist);
  }else if(strcmp(dot + 1, "fmca") == 0){
    readFMCA(expData, expDataName, expHist);
  }else{
    printf("ERROR: Improper type of input file: %s\n", expDataName);
    printf(
        "Integer array (.mca) and float array (.fmca) files are supported.\n");
    exit(-1);
  }
  fclose(expData);

  for(i = 0; i < NSPECT; i++){
    for(j = 0; j < S32K; j++){
      for(m=0;m<rebinFactor;m++){
        data[i]->Fill(j*rebinFactor + m, expHist[i][j]);
      }
    }
  }
  
  for(i = 0; i < numSimData; i++){
    if((simData = fopen(simDataName[i], "r")) == NULL){
      printf("ERROR: Cannot open the simulated data file %s!\n", simDataName[i]);
      exit(-1);
    }
    const char *dots = strrchr(simDataName[i], '.'); // get the file extension
    if(strcmp(dots + 1, "mca") == 0){
      readMCA(simData, simDataName[i], simHist[i]);
    }else if(strcmp(dots + 1, "fmca") == 0){
      readFMCA(simData, simDataName[i], simHist[i]);
    }else{
      printf("ERROR: Improper type of input file: %s\n", simDataName[i]);
      printf("Integer array (.mca) and float array (.fmca) files are supported.\n");
      exit(-1);
    }
    fclose(simData);
  }

  //set limits, taking re-binning into account
  for(i = 0; i < numSpectra; i++){
    fitStartCh[i] = startCh[i]/rebinFactor;
    fitEndCh[i] = endCh[i]/rebinFactor;
  }

  /*for(i = 0; i < numSpectra; i++)
    for(j = 0; j < S32K; j++)
      for(k=0;k<numSimData;k++){
        printf("sim hist %i sp %i ch %i: %f\n",k,spectrum[i],j,simHist[k][spectrum[i]][j]);
      }*/

  /*// read into ROOT
  for(i = 0; i < numSpectra; i++)
    for(j = 0; j < S32K; j++)
      sim[i]->Fill(j, simHist[0][i][j]);*/

  if(verbosity>0){
    printf("Spectra read in...\n");
  }
  
  // do the fit
  double chisq = find_chisqMin();

  // scale simulated data
  for(i = 0; i < numSpectra; i++){
    for(j = 0; j < S32K; j++){
      for(m=0;m<rebinFactor;m++){
        for(k=0;k<numSimData;k++){
          scaledSimHist[k][spectrum[i]][j*rebinFactor + m]=0.;
          scaledSimHist[k][spectrum[i]][j*rebinFactor + m] += aFinal[k+3][i] * simHist[k][spectrum[i]][j];
        }
      }
    }
  }

  // add background to simulated data
  for(i = 0; i < numSpectra; i++){
    for(j = 0; j < S32K; j++){
      for(m=0;m<rebinFactor;m++){
        if(addBackground == 1){
          bgHist[spectrum[i]][j*rebinFactor + m] = aFinal[0][i] + aFinal[1][i]*(double)j + aFinal[2][i]*(double)j*(double)j;
        }else if(addBackground == 2){
          bgHist[spectrum[i]][j*rebinFactor + m] = aFinal[0][i] + aFinal[1][i]*(double)j;
        }else if(addBackground == 3){
          bgHist[spectrum[i]][j*rebinFactor + m] = aFinal[0][i];
        }else{
          bgHist[spectrum[i]][j*rebinFactor + m] = 0.;
        }
      }
    }
  }

  // fit result histogram
  for(i = 0; i < numSpectra; i++){
    for(j = 0; j < S32K; j++){
      //no need to take re-binning into account here as
      //this was already done in the earlier steps (bgHist,scaledSimHist)
      resultsHist[spectrum[i]][j] = bgHist[spectrum[i]][j];
      for(k=0;k<numSimData;k++){
        resultsHist[spectrum[i]][j] += (float)scaledSimHist[k][spectrum[i]][j];
      }
    }
  }

  // print output
  if(verbosity>0){
    printf("Fit chisq: %.15f\n", chisq);
  }else{
    printf("%.15f\n", chisq);
  }

  // write results
  if(saveResults == 1){
    //write spectra
    sprintf(str, "fit_poisson.fmca");
    if((results = fopen(str, "w")) == NULL){
      printf("ERROR: Cannot open the output file %s!\n", str);
      exit(-1);
    }
    for(i = 0; i < numSpectra; i++){
      fwrite(resultsHist[spectrum[i]], S32K * sizeof(float), 1, results);
    }
    fclose(results);

    //save scaling factors
    sprintf(str, "scalingFactors.dat");
    if((scalingFactors = fopen(str, "w")) == NULL){
      printf("ERROR: Cannot open the output file %s!\n", str);
      exit(-1);
    }
    for(i = 0; i < numSpectra; i++){
      fprintf(scalingFactors, "%d %.9f\n", i + 1, aFinal[0][i]);
    }
    fclose(scalingFactors);
  }

  /*printf("Scaling factors:\n");
  for(i = 0; i < numSpectra; i++){
    printf("%d %.9f\n", i + 1, aFinal[0][i]);
  }*/

  // plot results
  if(plotMode >= 0){
    theApp = new TApplication("App", &argc, argv);
    plotSpectra();
  }

  return 0; // great success
}

double lrchisq(const double *par){
  // chisq from likelihood ratio test
  // for more information see Baker and Cousins pg. 439 and Appendix
  double yi = 0.;      // model
  double ni = 0.;      // experiment
  double lrchisq = 0.; // likelihood ratio chisq
  int i = 0;
  int j = 0;

  for(i = fitStartCh[spCurrent]; i <= fitEndCh[spCurrent]; i++){
    ni = expCurrent[i]; // events in ith bin

    // calculate model in the ith bin
    yi=0;
    for(j=0; j<numSimData; j++){
      
      if((j>=1)&&(useRelIntensities)&&(relIntensityAvailable[j])){
        yi += par[3]*par[j+3] * simCurrent[j][i]; //amplitude relative to the 1st ampliutude
      }else{
        yi += par[j+3] * simCurrent[j][i];
      }
      
    }
    // add background if neccesary
    if(addBackground == 1){
      yi += par[0] + par[1]*(double)i + par[2]*(double)i*(double)i;
    }else if(addBackground == 2){
      yi += par[0] + par[1]*(double)i;
    }else if(addBackground == 3){
      yi += par[0];
    }

    //model cannot give less than 0 counts
    if(yi<0.0){
      return 1E10;
    }

    // evaluate chisq given input parameters
    if((ni > 0.)&&(yi != 0.)){
      lrchisq += (yi - ni + ni * log(ni / yi));
    }else{
      lrchisq += yi; // the log(0) case
    }
  }
  lrchisq *= 2.;

/*   printf("Parameters: %f %f %f %f, computed lrchisq: %f\n",par[0],par[1],par[2],par[3],lrchisq);
  getc(stdin); */
  return lrchisq;
}

double find_chisqMin(){
  int i = 0;
  int j = 0;
  int k = 0;
  char str[256];
  double intExp,intSim;
  double chisq = 0.;

  if(verbosity>0)
    printf("Fitting data...\n");

  for(i = 0; i < numSpectra; i++){
    if(verbosity>0){
      printf("-Spectrum %i-\n",i+1);
    }
    // for more information see minimizer class documentation
    // https://root.cern.ch/root/html/ROOT__Math__Minimizer.html
    char minName[132] = "Minuit";
    char algoName[132] = ""; // default (Migard for Minuit)
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    // set tolerance , etc...
    min->SetMaxFunctionCalls(10000000); // for Minuit
    min->SetMaxIterations(100000);
    min->SetTolerance(0.001);
    if(verbosity>1){
      min->SetPrintLevel(1);
    }else{
      min->SetPrintLevel(0); // set to 1 for more info
    }

    // communication to the likelihood ratio function
    // (via global parameters)
    for(j = 0; j < S32K; j++){
      expCurrent[j] = (double)expHist[spectrum[i]][j];
    }
    for(j = 0; j < S32K; j++){
      for(k=0;k<numSimData;k++)
        simCurrent[k][j] = (double)simHist[k][spectrum[i]][j];
    }
    spCurrent = i; //needed for fit

    //calculate integrals
    intSim=0.;
    intExp=0.;
    for(j = fitStartCh[i]; j <= fitEndCh[i]; j++){
      for(k=0;k<numSimData;k++){
        intSim += simCurrent[k][j];
        //printf("ch %i sim val: %f\n",j,simCurrent[k][j]);
      }
      intExp += (double)expHist[spectrum[i]][j]; 
    }
    //printf("intExp: %f, intSim: %f\n",intExp,intSim);

    // create function wrapper for minmizer
    // a IMultiGenFunction type
    ROOT::Math::Functor lr(&lrchisq, 3+numSimData); // likelihood ratio chisq

    // step size and starting variables
    // may need to change for best performance
    // under different running conditions
    double ratio = intExp/intSim;
    double variable[3+NSIMDATA];
    double step[3+NSIMDATA];

    for(j=0;j<3+numSimData;j++){
      variable[j] = ratio/2.;
      step[j] = ratio/100.;
      //printf("Variable %i: val %f, step %f\n",j,variable[j],step[j]);
    }
    variable[1]/=100000.; //slope for linear background is usually very small
    variable[2]/=100000.;

    min->SetFunction(lr);

    // Set pars for minimization
    for(j=0;j<3+numSimData;j++){
      sprintf(str, "a%i", j);
      if((j>=4)&&(useRelIntensities)&&(ril[j-3]!=rih[j-3])){
        min->SetLimitedVariable(j, str, variable[j], step[j],ril[j-3],rih[j-3]); //set relative intensity
      }else if((j>=3)&&(forcePosAmp==1)){
        min->SetLimitedVariable(j, str, variable[j], step[j],0.0,5.0*ratio); //j>=2 for amplitudes
      }else if((j==1)&&(forceNegSlopeBG==1)){
        min->SetLimitedVariable(j, str, -1.0*variable[j], step[j]/100., -1000.0, 0.0); //slope not allowed to be positive
      }else{
        min->SetVariable(j, str, variable[j], step[j]);
      }
      //printf("Variable %i initialized to %f, step size %f",j,variable[j], step[j]);
      //printf(".\n");
    }

    // do the minimization
    min->Minimize();

    // grab parameters from minimum
    const double *xs = min->X();

    //if(verbosity>0)
    //  printf("Spectrum %i fit parameters: %f %f %f\n",i,xs[0],xs[1],xs[2]);

    // assuming 3 parameters
    // save pars
    for(j = 0; j < 3+numSimData; j++){
      if((j>=3)&&(useRelIntensities)&&(relIntensityAvailable[j-3])){
        aFinal[j][i] = xs[j]*xs[3]; //intensity relative to first intensity
      }else{
        aFinal[j][i] = xs[j];
      }
    }

    //print amplitudes
    if(verbosity>0){
      for(j = 3; j < 3+numSimData; j++){
        printf("Spectrum %i, amplitude %i: %f\n",i+1,j-1,aFinal[j][i]);
      }
    }

    chisq += min->MinValue();

  }

  return chisq;
}

void plotSpectra(){
  
  int i, j, k;
  TH1D *results[NSPECT];
  TH1D *resultsBGData[NSPECT];
  TH1D *resultsSimData[NSIMDATA][NSPECT];
  TCanvas *c = new TCanvas("c1", "", 1618, 1000);

  // initialize and fill results histo
  char resultsName[132],resultsBGName[132],resultsSimDataName[132];
  for(i = 0; i < NSPECT; i++){
    sprintf(resultsName, "results_%2d", i);
    results[i] = new TH1D(resultsName, ";;", S32K, 0, S32K - 1);
    if(plotMode>=1){
      sprintf(resultsBGName, "resultsBG_%2d", i);
      resultsBGData[i] = new TH1D(resultsBGName, ";;", S32K, 0, S32K - 1);
      for(k=0;k<numSimData;k++){
        sprintf(resultsSimDataName, "resultsSimData_%2d_%2d", i, k);
        resultsSimData[k][i] = new TH1D(resultsSimDataName, ";;", S32K, 0, S32K - 1);
      }
    }
    
  }
  // be careful with indicies here
  for(i = 0; i < numSpectra; i++){
    for(j = 0; j < S32K; j++){
      results[spectrum[i]]->Fill(j, resultsHist[spectrum[i]][j]);
      if(plotMode>=1){
        resultsBGData[spectrum[i]]->Fill(j, bgHist[spectrum[i]][j]);
        for(k=0;k<numSimData;k++)
          resultsSimData[k][spectrum[i]]->Fill(j, scaledSimHist[k][spectrum[i]][j] + bgHist[spectrum[i]][j]);
      }
    }
  }

  // display limits
  Double_t low[NSPECT], high[NSPECT];

  // divide canvas to accomodate all spectra
  // assumes number to plot <= 3*2=6
  c->Divide(3, (int)(ceil(numSpectra/3.0)), 1E-6, 1E-6, 0);

  // formatting the plot area
  // be careful with indicies
  for(i = 1; i <= numSpectra; i++){
    c->GetPad(i)->SetRightMargin(0.01);
    c->GetPad(i)->SetTopMargin(0.025);
    Int_t ind = i - 1;
    low[i] = startCh[ind];
    high[i] = endCh[ind];
  }

  // plot
  for(i = 1; i <= numSpectra; i++){
    c->cd(i);
    Int_t ind = i - 1;
    // data in black
    data[spectrum[ind]]->SetLineStyle(1);
    data[spectrum[ind]]->SetLineWidth(2);
    data[spectrum[ind]]->SetLineColor(12);
    data[spectrum[ind]]->GetXaxis()->SetRangeUser(low[i], high[i]);
    data[spectrum[ind]]->SetStats(0);
    data[spectrum[ind]]->Draw("HIST");

    // simulation in red
    results[spectrum[ind]]->SetLineStyle(1);
    results[spectrum[ind]]->SetLineWidth(2);
    results[spectrum[ind]]->SetLineColor(46);
    results[spectrum[ind]]->Draw("HIST SAME");
    if(plotMode>=1){
      //plot individual simulated data
      for(k=0;k<numSimData;k++){
        resultsSimData[k][spectrum[ind]]->SetLineStyle(1);
        resultsSimData[k][spectrum[ind]]->SetLineWidth(2);
        resultsSimData[k][spectrum[ind]]->SetLineColor(799+k*20);
        resultsSimData[k][spectrum[ind]]->Draw("HIST SAME");
      }
      if(plotMode!=2){
        //plot background
        resultsBGData[spectrum[ind]]->SetLineStyle(1);
        resultsBGData[spectrum[ind]]->SetLineWidth(1);
        resultsBGData[spectrum[ind]]->SetLineColor(920);
        resultsBGData[spectrum[ind]]->Draw("HIST SAME");
      }
      
    }
    
  }

  c->cd(1);
  // x1,y1,x2,y2
  TLegend *leg;
  if(plotMode>=1){
    leg = new TLegend(0.70, 0.85 - numSimData*0.05, 0.90, 0.95);
  }else{
    leg = new TLegend(0.70, 0.85, 0.90, 0.95);
  }
  leg->AddEntry(data[0], "Experiment", "l");
  leg->AddEntry(results[spectrum[0]], "Simulation", "l");
  if(plotMode>=1){
    for(k=0;k<numSimData;k++){
      sprintf(resultsName, "Branch %i", k+1);
      leg->AddEntry(resultsSimData[k][spectrum[k+1]], resultsName, "l");
    }
  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->Draw();

  c->SetBorderMode(0);
  c->Update();

  theApp->Run(kTRUE);
}

int readMCA(FILE *inp, char *filename, float inpHist[NSPECT][S32K]){
  
  if(verbosity>0){
    printf("Reading %i spectra in MCA file: %s\n",endSpectrum,filename);
  }
  
  int mcaHist[NSPECT][S32K];
  for(int i = 0; i <= endSpectrum; i++){
    if(fread(mcaHist[i], S32K * sizeof(int), 1, inp) != 1){
      printf("ERROR: Error reading file %s!\n", filename);
      exit(-1);
    }
  }
  
  for(int i = 0; i < NSPECT; i++){
    int bin = -1;
    for(int j = 0; j < S32K; j++){
      if(j%rebinFactor==0){
        bin++;
      }
      inpHist[i][bin] += (float)mcaHist[i][j];
    }
  }

  return 1;
}

int readFMCA(FILE *inp, char *filename, float inpHist[NSPECT][S32K]){

  if(verbosity>0){
    printf("Reading %i spectra in FMCA file: %s\n",endSpectrum,filename);
  }

  float mcaHist[NSPECT][S32K];
  for(int i = 0; i <= endSpectrum; i++){
    if(fread(mcaHist[i], S32K * sizeof(float), 1, inp) != 1){
      printf("ERROR: Error reading file %s!\n", filename);
      exit(-1);
    }
  }

  for(int i = 0; i < NSPECT; i++){
    int bin = -1;
    for(int j = 0; j < S32K; j++){
      if(j%rebinFactor==0){
        bin++;
      }
      inpHist[i][bin] += mcaHist[i][j];
    }
  }

  return 1;
}
