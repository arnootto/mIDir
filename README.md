# The mIDir package
An R package to implement the methodology described in *A refreshing take on the inverted Dirichlet via a modeparameterization with some statistical illustrations* (2025).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/mIDir")
```
## Example
Code to reproduce Example 5.1: Cantaloupe spectra data (V3 vs V4)
```{r}
# Fit the mIDir and cmIDir models to V3 and V4 of the cantaloupe spectra data
library(mIDir)
library(rrcov)
data("fruit")
X <- fruit[,c("V3","V4")]*10 #scale the data
a=ML.mID(X)
b=ML.cMID(X)
```
Code to reproduce Example 5.2: Australian Institute of Sport data
```{r}
# Fit the mixtures of mIDir to the Australian Institute of Sport data
library(mIDir)
library(sn)
data("ais")
c=ML.mMID(X = as.data.frame(cbind(ais$LBM,ais$Wt,ais$BMI,ais$WCC,ais$Bfat)),k=2, method = "BFGS",initialization = "random.soft",max.iter = 1000)
```
