# Irregular2dWavelets

This is an R package to implement the two-dimensional wavelet decomposition approach for irregular shapes and grids proposed in "Spatial Multiresolution Analysis of the Effect of PM2.5 on Birth Weights" by Antonelli et. al (2016). The paper can be found at the following link

http://arxiv.org/pdf/1609.05186v1.pdf

The package can be read in using the following lines of code:

```{r, echo=TRUE, message=FALSE}
library(devtools)
install_github(repo = "jantonelli111/Irregular2dWavelets")
library(Irregular2dWavelets)
```

If you wish to see the vignette associated with the package that helps to illustrate its usage through a data example, use the following lines of code. This will take a minute or two longer to build.

```{r, echo=TRUE, message=FALSE}
library(devtools)
install_github(repo = "jantonelli111/Irregular2dWavelets", build_vignettes = TRUE)
library(Irregular2dWavelets)
```
Then, to view the vignette simply type into R

```{r, echo=TRUE, message=FALSE}
vignette("Irregular2dWavelets", package="Irregular2dWavelets")
```
The PDF of the vignette can also be found in the vignettes folder of this github repository.
