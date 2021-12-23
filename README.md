GMM-NEA 
======

## Overview

`R` codes implementing the methodologies reported in the paper "*Reconstruction of Nuclear Ensemble Approach Electronic Spectra using Probabilistic Machine Learning*" (L. Cerdán and D. Roca-Sanjuán, 2022). It contains a fully functional version of the algorithm to reconstruct the NEA spectra using both auto-delta and GMM-NEA.

## Usage

It requires the installation of [R](https://cran.r-project.org/) and/or the [Rstudio IDE](https://www.rstudio.com/products/rstudio/).


```
$ Rscript GMM_NEA.R input_file molecule outliers ci
```

The argument `input_file` specifies the name of the `.csv` file containing the values for the vertical excitation energies (VEE) and oscillator strengths (f) for all transitions (columns) and geometries (rows). An example of this kind of file can be found in the folder `Data`.

The argument `molecule` specifies the molecule identifier to add to the output files

The argument `outliers` specifies the program to look for outliers (optional argument). By default, it is set to `FALSE`.

The argument `ci` tells the program to compute the confidence intervals (CI) of the reconstruction. By default, it is set to `FALSE`.

Example to run the script with default options: 

```
$ Rscript GMM_NEA.R target_uracile.csv U6OH_with_outliers
```

Example to run the script looking for outliers and computing the CI:

```
$ Rscript GMM_NEA.R target_uracile.csv U6OH outliers ci
```


## Reference

L. Cerdán and D. Roca-Sanjuán, Reconstruction of Nuclear Ensemble Approach Electronic Spectra using Probabilistic
Machine Learning. *J. Chem. Theory Comput.*,  2022, XXX, XXX−XXX
