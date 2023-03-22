FILE *config;
int spectrum[NSPECT], startCh[NSPECT], endCh[NSPECT];
int numSpectra, endSpectrum, maxNumCh, numSimData;
double ril[NSIMDATA], rih[NSIMDATA];
int addBackground; // 0=no, 1=constant background
int plotOutput;    // 0=no, 1=yes, 2=detailed
int forcePosAmp;   // 0=no, 1=yes
int forceNegSlopeBG; //0=no, 1=yes
int useRelIntensities; //0=no, 1=yes
int relIntensityAvailable[NSIMDATA]; //0=no, 1=yes
int plotMode;      // 0 = normal, 1=detailed
int saveResults;        // 0=no, 1=yes
char expDataName[256], simDataName[NSIMDATA][256]; // file names
char str[256], str1[256], str2[256];

void readConfigFile(const char *fileName) {
  int index = 0;
  numSpectra = 0;
  endSpectrum = 0;
  maxNumCh = 0;
  numSimData = 0;
  saveResults = 0;
  if((config = fopen(fileName, "r")) == NULL){
    printf("ERROR: Cannot open the config file %s!\n", fileName);
    exit(-1);
  }
  while (!(feof(config))) // go until the end of file is reached
  {
    if(fgets(str, 256, config) != NULL){
      if(index < NSPECT)
        // spectrum, channel range and step function parameter data
        if(sscanf(str, "%i %i %i", &spectrum[index], &startCh[index], &endCh[index]) == 3){
          if(spectrum[index] > endSpectrum)
            endSpectrum = spectrum[index];
          if((endCh[index] - startCh[index] + 1) > maxNumCh)
            maxNumCh = endCh[index] - startCh[index] + 1;
          index++;
          numSpectra++;
        }
      
      //simulated data parameters
      if(sscanf(str, "%s %s %lf %lf", str1, str2, &ril[numSimData], &rih[numSimData]) == 4){
        if(strcmp(str1, "SIMULATED_DATA") == 0){
          if(numSimData == 0){
            printf("ERROR: cannot set relative intensities for the first simulated dataset (%s).\n",str2);
            printf("Intensities should only be specified for later datasets, relative to the first dataset.\n");
            exit(-1);
          }
          if(numSimData<NSIMDATA){
            if((verbosity>0)&&(numSimData>0))
              printf("Relative intensities provided for simulated data %i.\n",numSimData+1);
            strcpy(simDataName[numSimData], str2);
            relIntensityAvailable[numSimData]=1;
            numSimData++;
          }
        }
      }else if (sscanf(str, "%s %s", str1, str2) == 2){
        if (strcmp(str1, "SIMULATED_DATA") == 0){
          if(numSimData<NSIMDATA){
            if((verbosity>0)&&(numSimData>0))
              printf("Relative intensities not provided for simulated data %i.\n",numSimData+1);
            strcpy(simDataName[numSimData], str2);
            relIntensityAvailable[numSimData]=0;
            numSimData++;
          }
        }
      }

      if(sscanf(str, "%s %s", str1, str2) == 2){ // single parameter data
        if(strcmp(str1, "EXPERIMENT_DATA") == 0)
          strcpy(expDataName, str2);
        
        if(strcmp(str1, "ADD_BACKGROUND") == 0){
          if(strcmp(str2, "quad") == 0)
            addBackground = 1; //quadratic
          else if(strcmp(str2, "lin") == 0)
            addBackground = 2; //linear
          else if(strcmp(str2, "const") == 0)
            addBackground = 3; //constant
          else
            addBackground = 0;
        }
        if(strcmp(str1, "FORCE_POSITIVE_AMPLITUDE") == 0){
          if(strcmp(str2, "yes") == 0)
            forcePosAmp = 1;
          else
            forcePosAmp = 0;
        }
        if(strcmp(str1, "FORCE_NEGATIVE_SLOPE_BACKGROUND") == 0){
          if(strcmp(str2, "yes") == 0)
            forceNegSlopeBG = 1;
          else
            forceNegSlopeBG = 0;
        }
        if(strcmp(str1, "USE_RELATIVE_INTENSITIES") == 0){
          if(strcmp(str2, "yes") == 0)
            useRelIntensities = 1;
          else
            useRelIntensities = 0;
        }
        if((strcmp(str1, "PLOT_MODE") == 0)||(strcmp(str1, "PLOT_OUTPUT") == 0)){
          if (strcmp(str2, "yes") == 0)
            plotMode = 0;
          else if (strcmp(str2, "detailed") == 0)
            plotMode = 1;
          else if (strcmp(str2, "detailed_nobg") == 0)
            plotMode = 2;
          else
            plotMode = -1;
        }
        if(strcmp(str1, "VERBOSE") == 0){
          if (strcmp(str2, "yes") == 0)
            verbosity = 1;
          else if (strcmp(str2, "debug") == 0)
            verbosity = 2;
          else
            verbosity = 0;
        }
        if(strcmp(str1, "SAVE_RESULTS") == 0){
          if (strcmp(str2, "yes") == 0)
            saveResults = 1;
          else
            saveResults = 0;
        }
      }
      if(sscanf(str, "%s %s", str1, str2) == 1){ // listing of simulated data
        if(strcmp(str1, "<---END_OF_PARAMETERS--->") == 0)
          break;
      }
    }
  }
  fclose(config);

  if(verbosity>0){
    printf("Taking experiment data from file: %s\n", expDataName);
    if(numSimData > 0){
      printf("Taking simulated data from file(s): %s", simDataName[0]);
      for(int i=1;i<numSimData;i++)
        {
          printf(", %s", simDataName[i]);
          if(useRelIntensities){
            printf(" with relative intensity ranging from %.2lf to %.2lf",ril[i],rih[i]);
          }
        }
      printf("\n");
    }else{
      printf("ERROR: no simulated data specified!\n");
      exit(-1);
    }
    for (index = 0; index < numSpectra; index++)
      printf("Will compare spectrum %i from channels %i to %i.\n",
            spectrum[index], startCh[index], endCh[index]);
    if(addBackground == 0)
      printf("Will not add background to simulated data.\n");
    if(addBackground == 2)
      printf("Will add a linear background to simulated data.\n");
    if(plotOutput == 0)
      printf("Will not plot output data.\n");
    if(plotOutput == 1)
      printf("Will plot output data.\n");
    if(plotOutput == 2)
      printf("Will plot detailed output data.\n");
    if(saveResults == 1)
      printf("Saving scaled fit results.\n");
    printf("Finished reading parameter file...\n");
  }
  
}
