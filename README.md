# constituents-extraction

Curve fitting to extract Aggregate and Monomer constituents
from frame vs intensity data for a biological sample.

I am not a biologist, so I am just going to call them `peak1` and `peak2`.


## Contents

- [`fit.py`]: Python scipt that does the curve fitting.
- [`data/`]: frame vs intensity raw data.
- [`output/`]: results of the curve fitting.

[`fit.py`]: fit.py
[`data/`]: data/
[`output/`]: output/


## Technical details

1. Data is loaded and normalised to between 0 and 1.
2. Unsophisticated heuristics are used to make initial guesses
   for the fitting parameters.
3. Curves are fitted using `scipy.optimize`.
4. Results are exported to [`data/`].

![Plot of 2-peak fit for Sample 3 data.][plot]

[plot]: output/fit-Sample_3.txt.png
