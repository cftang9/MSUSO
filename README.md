# Multiple Ordinal Dominance Curves and Uniform Stochastic Ordering

This repository contains R programs for the article, “Multiple Ordinal Dominance Curves and Uniform Stochastic Ordering”. 
<!-- This article has been submitted for publication. -->

Prior to using R programs on this repository, please download the main R program [EGJ_USO_Library.R](https://raw.githubusercontent.com/cftang9/MSUSO/master/EGJ_USO_Library.r). 

## Part 1. Reproducing simulation results in the manuscript

The following simulations generate random samples both from the mixture normal distribution and ODCs. Since both distingusing distribution methods and GOF tests depends on the ODC between consecutive distributions. it also suffices to generating random samples from ODCs. In the manuscripts, we consider G_1 from [Tang et. al. (2017)](../master/README.md#reference) on the left of Figure 1 below for star-shaped ODC and G_1^{-1} for non-star-shaped ODC.
The sequence of the ODCs from [Wang and Tang (2020)](../master/README.md#reference) on the right of Figure 1 are for power curves comparison.
The computation times in the followings are based on on a computer with a 3.0GHz processor and 64GB of memory. 

```R``` codes for ODC plots in Figure 1: [Figure_1_ODCs.R](https://github.com/cftang9/MSUSO/blob/master/Figure_1_ODCs_Plot.r).
![Figure 1](../master/Figure_1_ODC_Plot.png)

### 1.1. Table 1: Equal Distribution Test under USO with k=3

The ```R``` codes to reproduce Table 1 is attached: [Testing_Equality_k3.R](https://github.com/cftang9/MSUSO/blob/master/Testing_Equality_k3.R). The calculation took approximately 10 minutes. 

<!--1.1 Table 1 in Section 3 of the manuscript 
To reproduce Table 1, which involves four classic copulas: Clayton, Frank, Gumbel, and Gaussian, please run this R program:
[Clayton_Frank_Gumbel_and_Gaussian_n=100.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D100.R).
But be aware of that, because the number of replications is 10,000, this program might take a long time to finish. As stated in our manuscript, our calculation of Table 1 took approximately 73 minutes on a computer with a 3.1GHz processor and 16GB of memory. 
-->

### 1.2. Table 2: Size and Power of Goodness-of-fit Tests for USO with k=3
[Testing_GOF_k3.R](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k3.R) provides the size and power studies with k=3 sample with sample size n=200 which requires 6 hours approximately. 

### 1.3. Figure 2: Power Curves Goodness-of-fit Tests for USO with k=3
The power curves comparison in Figure 2, (```R``` codes [Testing_GOF_k3_PC.R](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k3_PC.R))
with k=3 samples and equal sample sizes n=200, requires 6 hours totally on a computer with a 3.0GHz processor and 64GB of memory. 

### 1.4. Table 3: Distinguishing Distributions under USO with k=3

The ```R``` codes to reproduce Table 3 is attached: [Testing_Jump_k3.R](https://github.com/cftang9/MSUSO/blob/master/Testing_Jump_k3.R). The calculation took approximately 10 minutes. 


## Part 2. Reproducing simulation results in Web Appendix
In addition to the simulation results in the manuscript, more simulations results are provided in the supplementary materials with ```R``` codes attached in the followings. 

### 2.1. Tables: Equal Distribution Test under USO with k=4,5
Other than k=3 samples, we applied the distinghishing distribution methods to samples
[k=4](https://github.com/cftang9/MSUSO/blob/master/Testing_Equality_k4.R),
and [k=5](https://github.com/cftang9/MSUSO/blob/master/Testing_Equality_k5.R)
with sample sizes n=60,100,and 200. All the calculations took less than 10 minuntes. 

### 2.2. Tables: Goodness-of-fit Tests for USO with k=4,5
In addition to k=3 in [Table 2](../master/README.md#22-table-2-size-and-power-of-goodness-of-fit-tests-for-uso-with-k3), 
we provide the size and power study for 
[k=4](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k4.R),
and [k=5](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k5.R)
samples with sample sizes n=60, 100, and 200. All the calculations took less than 10 minuntes. 

### 2.3. Figures: Power Curves for Goodness-of-fit Tests for USO with k=4, 5
We also consider more settings for power curves comparison for GOF tests with k=4,5 samples and sample sizes n= 200 and 400 with ```R``` codes attached.
* For k=4 with sample sizes [n=200](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k4_PC.R) which took approximately 8 hours on a computer with a 3.0GHz processor and 64GB of memory, respectively. 
* For k=5 with sample sizes [n=200](https://github.com/cftang9/MSUSO/blob/master/Testing_GOF_k5_PC.R) which took approximately 10 hours, respectively. 

### 2.4. Tables: Distinguishing Distributions under USO with k=4,5
Other than k=3 samples, we applied the distinghishing distribution methods to samples
[k=4](https://github.com/cftang9/MSUSO/blob/master/Testing_Jump_k4.R),
and [k=5](https://github.com/cftang9/MSUSO/blob/master/Testing_Jump_k5.R)
with sample sizes n=60,100,and 200. All the calculations took less than 10 minuntes. 


## Part 3. MFAP4 data analysis
We applied both distinguishing distribution methods and GOF tests to microfibrillar-associated protein 4 (MFAP4) data with clinical cohort characteristics and MFAP4 serum levels from [Bracht et. al. (2016)](../master/README.md#reference) in [MFAP4.xlsx](https://static-content.springer.com/esm/art%3A10.1186%2Fs12967-016-0952-3/MediaObjects/12967_2016_952_MOESM1_ESM.xlsx). We grouped the MFAP4 levels in fibrosis stages and saved in ```R``` data form ```data_MFAP4``` in [MFAP4.Rdata](../master/MFAP4.Rdata). 

Here we provide the empirical estimators and estimators under USO for ODCs between consecutive fibrosis stages with 
```R``` codes attached: [Figure_3](../master/Figure_3_MFAP4.R).
![Figure 3](../master/Figure_3_MFAP4.png)

### 3.1. Table 3 part 1: Equal test for USO for MFAP4 levels

### 3.2. Table 3 part 2: Goodness-of-fit test for USO for MFAP4 levels

The second part of Table 3 () provides the departure ```M1s```, ```M2s```, and ```Mss``` of consecutive distributions from USO in L_p norm with p=1,2, and supremum norms, respectively. The critical values ```boot.cv.Skps``` and ```boot.cv.Wkps``` for the cumulated test statistics ```Skp``` and ```Wkp```, respectively, are provided. The Bonferroni-corrected critical values are saved in function ```Bon.cvs``` with L_p norms ```Bon.cv.p1```, ```$Bon.cv.p2```, and ```Bon.cv.ps```. 

### 3.3. Table 3 part 3: Distinghishing Fibrosis stages from MFAP4 levels
The first part (first 5 columns) of Table 3 provides the differences ```D1s```, ```D2s```, and ```Dss``` of distributions from equality in L_p norm with p=1,2, and supremum norms, respectively. The thresholds for each L_p differences are provided to determine if the consecutive distributions are distinct. 

## Reference: 
1. Thilo Bracht, Christian Mölleken, Maike Ahrens, Gereon Poschmann, Anders Schlosser, Martin Eisenacher, Kai Stühler, Helmut E. Meyer, Wolff H. Schmiegel, Uffe Holmskov, Grith L. Sorensen and Barbara Sitek (2016). [Evaluation of the biomarker candidate MFAP4 for non-invasive assessment of hepatic fibrosis in hepatitis C patients.](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-016-0952-3) *Journal of Translational Medicine*. 14:201.
2. Dewei Wang, Chuan-Fa Tang, and Joshua M. Tebbs (2020). [More powerful goodness-of-fit tests for uniform stochastic ordering.](http://www.sciencedirect.com/science/article/pii/S0167947319302531) *Computational Statistics & Data Analysis*. 144:106898.
3. Chuan-Fa Tang, Dewei Wang, and Joshua M. Tebbs (2017). [Nonparametric goodness-of-fit tests for uniform stochastic ordering](https://projecteuclid.org/euclid.aos/1513328583) *The Annals of Statistics*. 45:2565-2589.


