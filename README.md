# Multiple Ordinal Dominance Curves and Uniform Stochastic Ordering

This repository contains R programs for the article, “Multiple Ordinal Dominance Curves and Uniform Stochastic Ordering”. 
<!-- This article has been submitted for publication. -->

Prior to using R programs on this repository, please download the main R program [MUSOLibrary.R](https://raw.githubusercontent.com/cftang9/MSUSO/master/MUSOLibrary.R?token=AK5HQA6Z4FIJ4GDV5CVOAYLAVBZ6S). 

## Part 1. Reproducing simulation results

The following simulations generate random samples both from the mixture normal distribution and ODCs. Since both distingusing distribution methods and GOF tests depends on the ODC between consecutive distributions. it also suffices to generating random samples from ODCs. In the manuscripts, we consider R_1 from [Tang et. al. (2017)](../master/README.md#reference) on the right of Figure 1 below for unequal but USO ordered distributions. On the other hand, the ODC R_1^{-1} [Wang and Tang (2020)](../master/README.md#reference) is used to generate random samples which do not satisfy USO. The sequence of the ODCs on the right of Figure 1 is for power curves comparison. The computation times in the followings are based on on a computer with a 3.0GHz processor and 64GB of memory. 

```R``` codes for ODC plots in Figure 1: [Figure_1_ODCs.R](../master/Figure_1_ODCs.R).
![Figure 1](../master/Figure_1_ODCs.png)

### 1.1. Table 1: Distinguishing Distributions under USO with k=3

The ```R``` codes to reproduce Table 1 is attached: [Table_1_DD_k3.R](https://github.com/cftang9/MSUSO/blob/master/Table_1_DD_k3.R). The calculation took approximately 10 minutes. 

<!--2.1 Table 1 in Section 3 of the manuscript 
To reproduce Table 1, which involves four classic copulas: Clayton, Frank, Gumbel, and Gaussian, please run this R program:
[Clayton_Frank_Gumbel_and_Gaussian_n=100.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D100.R).
But be aware of that, because the number of replications is 10,000, this program might take a long time to finish. As stated in our manuscript, our calculation of Table 1 took approximately 73 minutes on a computer with a 3.1GHz processor and 16GB of memory. 
-->

### 1.2. Table 2: Size and Power of Goodness-of-fit Tests for USO with k=3
[Table_2](https://github.com/cftang9/MSUSO/blob/master/Table_2_GOF_k3.R) provides the size and power studies with k=3 sample with sample size n=200 which requires 6 hours approximately. 

### 1.3. Figure 2: Power Curves Goodness-of-fit Tests for USO with k=3
The power curves comparison in [Figure_2](https://github.com/cftang9/MSUSO/blob/master/Figure_2_GOF_PowerCurves_k3_200.R)
with k=3 sample with sample size n=200 which requires 6 hours. 

## Part 2. Reproducing simulation results in Web Appendix
In addition to the simulation results in the manuscript in [Part 2](../master/README.md#part-2-reproducing-simulation-results), more simulations results are provided in the supplementary materials with ```R``` codes attached in the followings. 

### 2.1. Tables: Distinguishing Distributions under USO with k=2,4,5
Other than k=3 samples, we applied the distinghishing distribution methods to samples
[k=2](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k2.R),
[k=4](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k4.R),
and [k=5](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k5.R)
with sample sizes n=200 and 400. All the calculations took less than 10 minuntes. 

### 2.2. Tables: Goodness-of-fit Tests for USO with k=2, 4, 5
In addition to k=3 in [Table 2](../master/README.md#22-table-2-size-and-power-of-goodness-of-fit-tests-for-uso-with-k3), 
we provide the size and power study for 
[k=2](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k2.R),
[k=4](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k4.R),
and [k=5](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k5.R)
samples with sample sizes n=200 and 400 with ```R``` codes attached. The calculation took approximately 8 hours. 

### 2.3. Figures: Power Curves for Goodness-of-fit Tests for USO with k=2, 3, 4, 5
We also consider more settings for power curves comparison for GOF tests with k=2,3,4,5 samples and sample sizes n= 200 and 400 with ```R``` codes attached.
* For k=3 with sample size [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k3_400.R) which took approximately 8 hours on a computer with a 3.0GHz processor and 64GB of memory. 
* For k=2 with sample sizes [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k2_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k2_400.R), they took approximately 4 and 6 hours, respectively.  
* For k=4 with sample sizes [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k4_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k4_400.R) which took approximately 8 and 10 hours, respectively. 
* For k=5 with sample sizes [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k5_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k5_400.R) which took approximately 10 and 12 hours, respectively. 

## Part 3. MFAP4 data analysis
We applied both distinguishing distribution methods and GOF tests to microfibrillar-associated protein 4 (MFAP4) data with clinical cohort characteristics and MFAP4 serum levels from [Bracht et. al. (2016)](../master/README.md#reference) in [MFAP4.xlsx](https://static-content.springer.com/esm/art%3A10.1186%2Fs12967-016-0952-3/MediaObjects/12967_2016_952_MOESM1_ESM.xlsx). We grouped the MFAP4 levels in fibrosis stages and saved in ```R``` data form ```data_MFAP4``` in [MFAP4.Rdata](../master/MFAP4.Rdata). 

Here we provide the empirical estimators and estimators under USO for ODCs between consecutive fibrosis stages with 
```R``` codes attached: [Figure_3](../master/Figure_3_MFAP4.R).
![Figure 3](../master/Figure_3_MFAP4.png)

### 3.1. Table 3 part 1: Equal test for USO for MFAP4 levels


### 3.2. Table 3 part 2: Goodness-of-fit test for USO for MFAP4 levels

The second part of Table 3 () provides the departure ```M1s```, ```M2s```, and ```Mss``` of consecutive distributions from USO in L_p norm with p=1,2, and supremum norms, respectively. The critical values ```boot.cv.Skps``` and ```boot.cv.Wkps``` for the cumulated test statistics ```Skp``` and ```Wkp```, respectively, are provided. The Bonferroni-corrected critical values are saved in function ```Bon.cvs``` with L_p norms ```Bon.cv.p1```, ```$Bon.cv.p2```, and ```Bon.cv.ps```. 
```R
set.seed(05222021)
MGOFUSO(data_MFAP4)
```
```R
$M1s
[1] 0.11992756 0.12407684 0.04938947 0.12970074

$M2s
[1] 0.14243105 0.14830742 0.06981125 0.15658425

$Mss
[1] 0.3488498 0.3322567 0.2897778 0.3263503

$Skps
[1] 0.4230946 0.5171340 1.2972346

$Wkps
[1] 0.1297007 0.1565843 0.3488498

$boot.cv.Skps
[1] 1.304033 1.588324 3.441980

$boot.cv.Wkps
[1] 0.5882955 0.6911207 1.3172881

$decision.Skps
[1] FALSE FALSE FALSE

$decision.Wkps
[1] FALSE FALSE FALSE

$decision.Bon
[1] FALSE FALSE FALSE
```
```R
Bon.cvs(k=5)
```
```R
$Bon.cv.p1
[1] 0.7283533

$Bon.cv.p2
[1] 0.8414987

$Bon.cv.ps
[1] 1.579907
```



### 3.1. Table 3 part 1: Distinghishing Fibrosis stages from MFAP4 levels
The first part (first 5 columns) of Table 3 provides the differences ```D1s```, ```D2s```, and ```Dss``` of distributions from equality in L_p norm with p=1,2, and supremum norms, respectively. The thresholds for each L_p differences are provided to determine if the consecutive distributions are distinct. 
```R
set.seed(05222021)
MDDUSO(data_MFAP4)
```
```R
$D1s
[1] 0.6544701 0.9131620 1.4065704 0.7376320

$D2s
[1] 0.7043039 0.9764447 1.5331956 0.7869467

$Dss
[1] 1.092246 1.555362 2.311166 1.036642

$thresholds
      95%       95%       95% 
0.8187272 0.9042299 1.4685763 

$distinction.p1
[1] 2 3

$distinction.p2
[1] 2 3

$distinction.ps
[1] 2 3
```


## Reference: 
1. Thilo Bracht, Christian Mölleken, Maike Ahrens, Gereon Poschmann, Anders Schlosser, Martin Eisenacher, Kai Stühler, Helmut E. Meyer, Wolff H. Schmiegel, Uffe Holmskov, Grith L. Sorensen and Barbara Sitek (2016). [Evaluation of the biomarker candidate MFAP4 for non-invasive assessment of hepatic fibrosis in hepatitis C patients.](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-016-0952-3) *Journal of Translational Medicine*. 14:201.
2. Dewei Wang, Chuan-Fa Tang, and Joshua M. Tebbs (2020). [More powerful goodness-of-fit tests for uniform stochastic ordering.](http://www.sciencedirect.com/science/article/pii/S0167947319302531) *Computational Statistics & Data Analysis*. 144:106898.
3. Chuan-Fa Tang, Dewei Wang, and Joshua M. Tebbs (2017). [Nonparametric goodness-of-fit tests for uniform stochastic ordering](https://projecteuclid.org/euclid.aos/1513328583) *The Annals of Statistics*. 45:2565-2589.


