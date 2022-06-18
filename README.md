# constituents-extraction

Curve fitting to extract aggregate and monomer constituents
from frame vs intensity data for biological samples.

I ([@yawnoc]) am not a biologist,
so I am just going to call them `peak1` and `peak2`.

[@yawnoc]: https://github.com/yawnoc


## License

Copyright 2022 Conway

This repository is licensed under 'MIT No Attribution' (MIT-0), see [LICENSE].
In short, this means you can do whatever you want with it.

[LICENSE]: LICENSE


## Contents

- [`fit.py`]: Python script that does the curve fitting
- [`data/`]: frame vs intensity raw data
- [`output/`]: results of the curve fitting
  - PNG/PDF plots
  - TSV files containing the parameter values

PDF output is not version-controlled, and can be found in [Releases].

[`fit.py`]: fit.py
[`data/`]: data/
[`output/`]: output/
[Releases]: https://github.com/yawnoc/constituents-extraction/releases


## Technical details

### Strategy

1. Data is loaded and normalised to between 0 and 1.
2. Unsophisticated heuristics are used to make initial guesses
   for the fitting parameters. <br>
   ![Plot of 2-peak parameter guesses for Sample 3 data.][guess-plot]
3. Curves are fitted using `scipy.optimize`. <br>
   ![Plot of 2-peak fit for Sample 3 data.][fit-plot]
4. Results are exported to [`output/`].

### Model

A straight line $y = m x + b$ is used to model the background intensity.

An [exponentially modified Gaussian] is used to model each peak in intensity.
This is a convolution of the density functions of
a normal distribution and an exponential distribution.
We use the chromatography version:

$$
\DeclareMathOperator{\erfc}{erfc}
f(x) =
        \frac{h \sigma}{\tau} \sqrt{\frac{\pi}{2}}
        \exp \left[ \frac{1}{2} \left( \frac{\sigma}{\tau} \right)^2 - \frac{x - \mu}{\tau} \right]
        \erfc \left[ \frac{1}{\sqrt{2}} \left( \frac{\sigma}{\tau} - \frac{x - \mu}{\sigma} \right) \right]
$$

where
- $h$ is Gaussian amplitude
- $\mu$ is Gaussian mean
- $\sigma$ is Gaussian standard deviation
- $\tau$ is exponential relaxation time.

[exponentially modified Gaussian]:
  https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution
[guess-plot]: output/guess-Sample_3.txt.png
[fit-plot]: output/fit-Sample_3.txt.png
