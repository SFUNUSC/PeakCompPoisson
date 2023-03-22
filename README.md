# PeakCompPoisson

A peak comparison program using the Poisson maximum likelihood ratio chi-square method.

## Description

Uses Poisson maximum likelihood ratio chi-square for the comparison of experimental to GEANT4-simulated RDM/DSAM lineshapes for TIP experiment analysis. Utilizes the [ROOT](https://root.cern/) Minuit function minimization libraries for numerical minimization of chi-square statistic w.r.t. model parameters. Particularly useful for low statistics data sets.

### Useful References: 

* A. Chester et al. Nucl. Inst. and Meth. A 882 (2017) 69-83.
* A. Chester et al. Phys. Rev. C. 96 (2017) 011302.
* S. Baker and R. D. Cousins. Nucl. Inst. and Meth. 221 (1984) 437-442.
* G. Cowan. Statistical Data Analysis. Clarendon Press, 1998.
* F. James. Minuit Function Minimization and Error Analysis Reference Manual. CERN 3 1994. Version 94.1. CERN Program Library entry D506.

## Installation

Requires [ROOT](https://root.cern/) (tested with v6.x) to be installed with environment variables set up properly.  Environment variable setup can be done by adding to your `.bashrc` (and then reloading the terminal):

```
#ROOT configuration in .bashrc
export ROOTSYS=/path/to/root
export ROOTINC=$ROOTSYS/include
export ROOTLIB=$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTLIB
```

Once this is done, use `make` to compile.  Tested using g++ and GNU make on Ubuntu 16.04 and CentOS 7.

## Use

The `peak_comp_poisson` code takes a parameter file as its only argument.

See the [sample parameter file](sample_parameters.dat) for an example.

Valid parameters that can be used are:

|**Parameter**|**Description**|
|:-----------:|:-------------:|
| EXPERIMENT_DATA | Contains the path to the .mca or .fmca file containing experiment data (**required**). |
| ADD_BACKGROUND | no: Don't add any background (default).<br>const: Add a constant background to the simulated data.<br>lin: Add a linear background to the simulated data.<br>quad: Add a quadratic background to the simulated data. |
| FORCE_NEGATIVE_SLOPE_BACKGROUND | yes: Force the background term to decrease (or stay constant) with increasing energy.<br>no: Leave background parameters free (default). |
| PLOT_OUTPUT | no: Show chisq stats only.<br>yes: Show a plot of the simulated and experimental data alongside chisq stats.<br>detailed: Same as 'yes', except plot all simulated datasets and background as well.<br>detailed_nobg: Same as 'detailed', except don't show background. |
| VERBOSE | no: Only print the fit chisq value (useful for interfacing with bash scripts).<br>yes: Print the results of reading in the parameter file and fitting.<br>debug: Print the above plus debug info. |
| SAVE_RESULTS | no: Do not save fit results (default).<br>yes: Save fit results to a file (`fit_poisson.fmca`) alongside scaling factors (`scalingFactors.dat`). |
| SIMULATED_DATA | Contains the path(s) to the .mca or .fmca file(s) containing simulated data (**required**).<br>Multiple datasets may be specified, each on a separate line with the same format.<br>Two numbers can be specified after the filename, corresponding to low and high limits for the relative intensity of that data, relative to any other datasets specified.  **Note**: intensity limits are only valid for the second and later files, the first file will be assumed to have relative intensity 1. |
| SPECTRUM  START_CHANNEL END_CHANNEL | Contains a list of the spectra to analyze in the experiment and simulated data, along with the channel range to analyze (**required**).  Multiple spectra may be specified, each on a separate line with the same format. |


## Credits

Original author: Aaron Chester

Current maintainer: Jonathan Williams
