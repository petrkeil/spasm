# spasm - **Sp**ecies **as**sociation **m**etrics for various data types

*Petr Keil*

This package contains code for several species association metrics and analyses.

## Beware, this is not a regular R package!

**Almost all functions in this package have much better and more mature counterpars
in other existing packages!** 'spasm' is just an assortment of code and wrappers that I needed
for a review of species association metrics, and putting the code to the form of 
a package is my way of organizing and documenting the code.

The existing packages that do a better job than `spasm` are:

- `EcoSimR` for cooccurrence analysis and null models.
- `spatstat` for point pattern analysis.
- `vegan` for dissimilarities, and all sorts of association indices.
- `vario` for community variograms.

**Before using a `spasm` function, please have a look at these packages - they may
do what you need in a more standard and better way.**

## Installing spasm

The easiest way to install spasm is using the ```install_github``` function from
the ```devtools``` package:

```{r}
library(devtools)
install_github("petrkeil/spasm")
```

## Overview

The functions provided here are classified according to the type of ecological
data that they use. 

### Spatially implicit site x species incidence matrices

These are functions whose names begin with ```C_```. A typical example is Stone & Roberts'
C-score.

### Spatially implicit site x species abundance matrices

Names of these functions begin with ```CA_```.

### Spatially explicit site x species matrices

So far these are not provided here, but in Dan McGlinn's package ```vario``` [here](https://github.com/dmcglinn/vario).

### Point patterns

The following methods are implemented, togehter with their respective null models:

- The P-M classification of point pattern overlap
- Bivariate pair correlation function 


## What needs to be cited

The package uses some example data, whose source should be credited, if these data are re-used. Specifically:

- Data from **Sonoran desert long-term plots** come from Rodriguez-Buritica S., Raichle H., Webb R.H., Turner R.M., and Venable D.L. (2013) One hundred and six years of population and community dynamics of Sonoran Desert Laboratory perennials. *Ecology*, 94:976.
- A **collection of site by species incidence (binary) matrices** from Atmar, W. and B.D. Patterson. 1995. Nestedness temperature calculator. [An Internet gopher publication in VisualBasic by AICS Research Inc, University Park, New Mexico, and The Field Museum, Chicago, posted on fmppr.fmnh.org 70]
- A **collection of site by species abundance matricess** from Ulrich W. and Gotelli N.J. (2010) Null model analysis of species associations using abundance data. **Ecology**, 91: 3384-3397.
- **Harvard forest plot** from Orwig D, Foster D, Ellison A. 2015. Harvard Forest CTFS-ForestGEO Mapped Forest Plot since 2014. Harvard Forest Data Archive: HF253.
