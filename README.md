# Multisample distinction and goodness-of-fit tests for uniform stochastic ordering via ordinal dominance curve

This repository contains R programs for the article, “Multisample distinction and goodness-of-fit tests for uniform stochastic ordering via ordinal dominance curve”. 
<!-- This article has been submitted for publication. -->

Prior to using R programs on this repository, please download the main R program [MUSOLibrary.R](https://raw.githubusercontent.com/cftang9/MSUSO/master/MUSOLibrary.R?token=AK5HQA6Z4FIJ4GDV5CVOAYLAVBZ6S). 

<!--
which requires installing `R` packages `Rcpp` and `copula`. We would like to point out that loading or executing functions in `Rcpp` packages may encounter some technical problems for Windows users if your `R` software was recently updated to the latest version. One may run these codes in `Rstudio` and follow what it suggests to solve the problem.  After successfully loading the main R program, the function `IndvsPQD` will automate critical value calculations for the practitioner. 
--> 

## Part 1. Illustration

Here we consider three mixture normal distributions. Each of them are mixed with the standard normal distribution with proportion p and a normal distribution with mean m and variance 1 with proportion (1-p). Here we generate three independent samples from mixture normal with p=(0.2,0.2,0.2), m=c(2,2.6,3.2) with sample sizes n=(150, 200, 250) using the ```MixNormal``` function which return a ```list``` of vectors of random samples:
```R
set.seed(20210521)
Data_demo = MixNormal(n=c(150,200,250),p=c(0.2,0.2,0.2),m=c(2,2.6,3.2))
```
These three distributions are simple ordered in the sense of uniform stochastic ordering (USO). 

### 1.1. Distinguishing Distribution Methods under USO

To perform the distinguishing distribution methods, we input the ```Data_demo``` in the function ``` MDDUSO```: 
```R
set.seed(20210521)
MDDUSO(Data_demo)
```

The function ```MDDUSO``` returns the scaled Lp strengh of USO, ```D1s, D2s, Dss``` for L1, L2, and sup norm, respectively for each ordinal dominance curve (ODC) from each consecutive distributions. The thresholds for each Lp strengh is reported by ```thresholds```. Then ```distinction.p1``` reports the locations where the inequality between ith consecutive distributions holds under USO for L1 (L2, sup norm) strength ```D1s``` (```D2s```, ```Dss```) larger than the first (second, third) variable of ```threshold```. 
<!--
For each collecion of Lp strength, the sum and maximum are reported by ```Skps``` and ```Wkps```. The bootstrapped critical values for ```Skps``` and ```Wkps``` are given by ```boot.cv.Skps``` and ```boot.cv.Wkps```. Lastly, ```decision.p1``` reports ```TRUE``` if USO is rejected by ```Skps``` > ```boot.cv.Skps``` and reports ```FALSE``` otherwise. Similar to ```decision.Skps```, the  ```decision.Skps``` reports ```TRUE``` if USO is rejected by ```Skps``` > ```boot.cv.Skps``` and reports ```FALSE``` otherwise.-->

```R
$D1s
[1] 1.081647 1.405280

$D2s
[1] 1.188057 1.698258

$Dss
[1] 1.85164 2.84605

$thresholds
      95%       95%       95% 
0.7186769 0.8021878 1.3424391 

$distinction.p1
[1] 1 2

$distinction.p2
[1] 1 2

$distinction.ps
[1] 1 2
```

### 1.2. Goodness-of-fit tests for USO

To perform the GOF tests, we input the ```Data_demo``` in the function ``` MDDUSO```: 
```R
set.seed(20210521)
MGOFUSO(Data_demo)
```
Input the ```Data_demo``` in the function ``` MGOFUSO```, it returns the scaled Lp departures, ```M1s, M2s, Mss``` for L1, L2, and sup norm, respectively, from USO for each ODC from each consecutive ODCs. For each collecion of Lp departures, the sum and maximum are reported by ```Skps``` and ```Wkps```. The bootstrapped critical values for ```Skps``` and ```Wkps``` are given by ```boot.cv.Skps``` and ```boot.cv.Wkps```. Lastly, ```decision.Skps``` reports ```TRUE``` if USO is rejected by ```Skps``` > ```boot.cv.Skps```, and reports ```FALSE``` otherwise. Similar to ```decision.Skps```, the  ```decision.Skps``` reports ```TRUE``` if USO is rejected by ```Skps``` > ```boot.cv.Skps```, and reports ```FALSE``` otherwise. The Bonferroni corrected method ```decnsion.Bon``` reports ```TRUE``` if USO is rejected when ```Wkps``` is larger than the Bonferroni corrected critical values, and reports ```FALSE``` otherwise. 

```R
$M1s
[1] 0.08392172 0.08892899

$M2s
[1] 0.1108007 0.1164118

$Mss
[1] 0.3454021 0.3187576

$Skps
[1] 0.1728507 0.2272125 0.6641597

$Wkps
[1] 0.08892899 0.11641176 0.34540211

$boot.cv.Skps
[1] 0.7735035 0.9545352 2.0628668

$boot.cv.Wkps
[1] 0.5441434 0.6587912 1.3409020

$decision.Skps
[1] FALSE FALSE FALSE

$decision.Wkps
[1] FALSE FALSE FALSE

$decision.Bon
[1] FALSE FALSE FALSE
```

### 1.3 For your own data
Please use the following ```R``` commands after naming the data ```my_data```:
```R
source(https://raw.githubusercontent.com/cftang9/MSUSO/master/MUSOLibrary.R?token=AK5HQA6Z4FIJ4GDV5CVOAYLAVBZ6S)
# name your data by my_data
MDDUSO(my_data)
MDDUSO(my_data)
```

## Part 2. Reproducing simulation results

### 2.0. Figure 1: Ordinal dominance curves
```R``` codes for ODC plots in Figure 1: [Figure_1_ODCs.R](../master/Figure_1_ODCs.R).
![Figure 1](../master/Figure_1_ODCs.png)

### 2.1. Table 1: Distinguishing Distributions under USO with k=3

[Table_1_DD_k3.R](https://github.com/cftang9/MSUSO/blob/master/Table_1_DD_k3.R)

<!--2.1 Table 1 in Section 3 of the manuscript 
To reproduce Table 1, which involves four classic copulas: Clayton, Frank, Gumbel, and Gaussian, please run this R program:
[Clayton_Frank_Gumbel_and_Gaussian_n=100.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D100.R).
But be aware of that, because the number of replications is 10,000, this program might take a long time to finish. As stated in our manuscript, our calculation of Table 1 took approximately 73 minutes on a computer with a 3.1GHz processor and 16GB of memory. 
-->

### 2.2. Table 2: Size and Power of Goodness-of-fit Tests for USO with k=3

The 


[Table_2_GOF_k3.R](https://github.com/cftang9/MSUSO/blob/master/Table_2_GOF_k3.R)

### 2.3. Figure 2: Power Curves Goodness-of-fit Tests for USO with k=3

[Figure_2_GOF_PowerCurves_k3_200.R](https://github.com/cftang9/MSUSO/blob/master/Figure_2_GOF_PowerCurves_k3_200.R)

## Part 3. Reproducing simulation results in Web Appendix
<!--
### 3.0. Figure A.1: Ordinal dominance curves
```R``` codes for ODC plots in Figure 1: [Supp_Figure_A.1_ODCs.R](../master/Supp_Figure_1_ODCs.R).
![Supp Figure 1](../master/Supp_Figure_1_ODCs.png)
-->
### 3.1. Tables: Distinguishing Distributions under USO with k=2,4,5

[Supp_Table_DD_k2.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k2.R)
[Supp_Table_DD_k4.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k4.R)
[Supp_Table_DD_k5.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_DD_k5.R)

<!--
### 2.2 Tables C.1 and C.2 in Web Appendix C

Table 1 considers n=100. We also included the same table but with n=50 and 200 in Web Appendix C. To reproduce those two tables. Please run [Clayton_Frank_Gumbel_and_Gaussian_n=50.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D50.R)
and
[Clayton_Frank_Gumbel_and_Gaussian_n=200.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D200.R), respectively.
-->
<!--
### 2.3 Tables C.3-C.5 in Web Appendix C

In addition to the Clayton, Frank, Gaussian, and Gumbel copulas, we have also considered the FGM and CA copulas and a restricted bivariate t distribution family. The results are presented in Tables C.3-C.5 in Web Appendix C. To reproduce these tables, please run
[Restricted_t_FGM_and_CA n=50.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D100.R),
[Restricted_t_FGM_and_CA n=100.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D50.R),
and
[Restricted_t_FGM_and_CA_n=200.R](https://raw.githubusercontent.com/cftang9/PQD/master/Restricted_t_FGM_and_CA_n%3D200.R).

## Part 3: To reproduce the real data analysis results in Section 4 of the manuscript
We applied all tests in this manuscript to three data applications. To reproduce the results of our analysis (Table 2 and Figures 2-4), please run the R program for each. The data included in the CSV file will be automatically read by the corresponding R program.
-->

### 3.2. Tables: Goodness-of-fit Tests for USO with k=2, 4, 5

[Supp_Table_GOF_k2.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k2.R)
[Supp_Table_GOF_k4.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k4.R)
[Supp_Table_GOF_k5.R](https://github.com/cftang9/MSUSO/blob/master/Supp_Table_GOF_k5.R)


### 3.3. Figures: Power Curves for Goodness-of-fit Tests for USO with k=2, 4, 5

We put  similar settings with two samples with sample size with ```R``` codes for k=3 with [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k3_400.R) 

k=2 with [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k2_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k2_400.R)

k=4 with [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k4_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k4_400.R)

k=5 with [n=200](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k5_200.R) and [n=400](https://github.com/cftang9/MSUSO/blob/master/Supp_Figure_GOF_PowerCurves_k5_400.R)



## Part 4. MFAP4 data analysis
We applied both distinguishing distribution methods and GOF tests to microfibrillar-associated protein 4 (MFAP4) data with clinical cohort characteristics and MFAP4 serum levels from [Bracht et. al. (2016)](../master/README.md#reference) in [MFAP4.xlsx](https://static-content.springer.com/esm/art%3A10.1186%2Fs12967-016-0952-3/MediaObjects/12967_2016_952_MOESM1_ESM.xlsx)
We grouped the MFAP4 levels in fibrosis stages and saved as a data ```data_MFAP4``` in [MFAP4.Rdata](../master/MFAP4.Rdata). 

```R``` codes for ODC plots in Figure 3: [Figure_3_MFAP4.R](../master/Figure_3_MFAP4.R).
![Figure 3](../master/Figure_3_MFAP4.png)

### 4.1. Table 3 part 1: Distinghishing Fibrosis stages from MFAP4 levels
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
### 4.2. Table 3 part 2: Goodness-of-fit test for USO for MFAP4 levels
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


## Reference: 
1. Thilo Bracht, Christian Mölleken, Maike Ahrens, Gereon Poschmann, Anders Schlosser, Martin Eisenacher, Kai Stühler, Helmut E. Meyer, Wolff H. Schmiegel, Uffe Holmskov, Grith L. Sorensen and Barbara Sitek (2016). [Evaluation of the biomarker candidate MFAP4 for non-invasive assessment of hepatic fibrosis in hepatitis C patients.](https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-016-0952-3) *Journal of Translational Medicine*. 14:201.
2. Dewei Wang, Chuan-Fa Tang, and Joshua M. Tebbs (2020). [More powerful goodness-of-fit tests for uniform stochastic ordering.](http://www.sciencedirect.com/science/article/pii/S0167947319302531) *Computational Statistics & Data Analysis*. 144:106898.
3. Chuan-Fa Tang, Dewei Wang, and Joshua M. Tebbs (2017). [Nonparametric goodness-of-fit tests for uniform stochastic ordering](https://projecteuclid.org/euclid.aos/1513328583) *The Annals of Statistics*. 45:2565-2589.)


<!--
<img src="https://render.githubusercontent.com/render/math?math=e^{i \pi} = -1">
-->


<!--
## Part 3: To reproduce the real data analysis results in Section 4 of the manuscript
We applied all tests in this manuscript to three data applications. To reproduce the results of our analysis (Table 2 and Figures 2-4), please run the R program for each. The data included in the CSV file will be automatically read by the corresponding R program.


### 3.1 Twins Data

Data: [TwinsData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.csv) 
(R program: [TwinsData.R](https://raw.githubusercontent.com/cftang9/PQD/master/TwinsData.R))

### 3.2 Education data

Data: [EducationData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.csv)
(R program: [EducationData.R](https://raw.githubusercontent.com/cftang9/PQD/master/EducationData.R))


### 3.3 Stock Data

Data: [StockData.csv](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.csv) 
(R program: [StockData.R](https://raw.githubusercontent.com/cftang9/PQD/master/StockData.R))

-->




