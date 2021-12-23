GMM-NEA 
======

## Overview

`R` code implementing the methodologies reported in the paper "*Reconstruction of Nuclear Ensemble Approach Electronic Spectra using Probabilistic Machine Learning*" (L. Cerdán and D. Roca-Sanjuán, 2022). 
It contains a fully functional version of the algorithm to reconstruct the NEA spectra using both auto-delta and GMM-NEA.

## Requirements and Dependencies

It requires the installation of [R](https://cran.r-project.org/) and/or the [Rstudio IDE](https://www.rstudio.com/products/rstudio/).

Before running this script, make sure that the following packages are installed: `'stats'`, `'scales'`, `'robustbase'`, `'matrixStats'`, `'mclust'`, `'foreach'`, `'doParallel'`
, `'parallel'`

## Usage

1) Save the source files `GMM_NEA.R` and `GMM_NEA_helper_funcs.R` into a folder. 

2) Into that same folder save the `.csv` file containing the values of the vertical excitation energies (*VEE*) and oscillator strengths (*f*) 
for all transitions (columns) and geometries/configurations (rows) to reconstruct the NEA spectrum. The first half of the columns contain the *VEE*, the second half the *f*s. 
An example of this kind of file can be found in the folder `Data`.

3) Access a terminal (OS or RStudio), go to the folder containing the source files, and type:

```
$ Rscript GMM_NEA.R input_file molecule [outliers] [ci]
```

The argument `input_file` specifies the name of the `.csv` file containing the *VEE* and *f* values.

The argument `molecule` specifies the molecule name/identifier to add to the output files

The *optional* argument `outliers` specifies the program to look for outliers (optional argument).

The *optional* argument `ci` tells the program to compute the confidence intervals (CI) of the reconstruction.

Examples: 

```
$ Rscript GMM_NEA.R vee_f_uracil_radical.csv U6OH outliers ci
```

<img src="figures/Spectra_eV_auto_d_vs_GMM_NEA_100_geoms_U6OH.png" width="400"/>

```
$ Rscript GMM_NEA.R vee_f_benzene.csv benzene
```

<img src="figures/Spectra_eV_auto_d_vs_GMM_NEA_250_geoms_benzene.png" width="400"/>


## Output

Once the computation is over, the following files are saved into the folder `./molecule`:

- auto_d_spectra_eV_*n_geoms*_*molecule*.csv

- GMM_NEA_spectra_eV_*n_geoms*_*molecule*.csv

- Spectra_eV_auto_d_vs_GMM_NEA_*n_geoms*_*molecule*.png


## Reference

L. Cerdán and D. Roca-Sanjuán, Reconstruction of Nuclear Ensemble Approach Electronic Spectra using Probabilistic
Machine Learning. *J. Chem. Theory Comput.*,  2022, XXX, XXX−XXX
