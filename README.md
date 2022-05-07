# constituents-extraction

Curve fitting to extract Aggregate and Monomer constituents
from frame vs intensity data for biological samples.

I am not a biologist, so I am just going to call them `peak1` and `peak2`.


## Contents

- [`fit.py`]: Python script that does the curve fitting
- [`data/`]: frame vs intensity raw data
- [`output/`]: results of the curve fitting
  - PNG/PDF plots
  - TSV files containing the parameter values

[`fit.py`]: fit.py
[`data/`]: data/
[`output/`]: output/


## Technical details

1. Data is loaded and normalised to between 0 and 1.
2. Unsophisticated heuristics are used to make initial guesses
   for the fitting parameters.
3. Curves are fitted using `scipy.optimize`.
4. Results are exported to [`data/`].

- A straight line `y = m x + b` is used to model the background intensity.
- An [exponentially modified Gaussian] is used to model each peak in intensity.
  This is a convolution of the density functions of
  a normal distribution with an exponential distribution.
  We use the chromatography version (see <<https://doi.org/10.1002/cem.1343>>):
  ````
  f(x) =
          h sigma/tau sqrt(pi)/2
          . exp[1/2 (sigma/tau)^2 - (x - mu)/tau]
          . erfc[1/sqrt(2) (sigma/tau - (x - mu)/sigma)]
  ````
  where
  - `h` is Gaussian amplitude
  - `mu` is Gaussian mean
  - `sigma` is Gaussian standard deviation
  - `tau` is exponential relaxation time.

![Plot of 2-peak fit for Sample 3 data.][plot]

[exponentially modified Gaussian]:
  https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
[plot]: output/fit-Sample_3.txt.png
