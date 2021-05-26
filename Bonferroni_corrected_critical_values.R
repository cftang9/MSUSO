rm(list=ls(all=TRUE))
#source("MUSOLibrary.R"); 

B0 = 10000; n0 = 3000
Bon.p1 = array(,B0); Bon.p2 = array(,B0); Bon.ps = array(,B0); 
set.seed(05242021)
for(b0 in 1:B0){
  X10 = sort(runif(n0)); X20 = sort(runif(n0)); 
  R10=rep(0,n0+1)
  for(i in 2:(n0+1))
  {
    R10[i]=sum(X10<=X20[i-1])/n0
  }
  u10 = seq(0,1,by=1/n0); 
  r10=c((1 - R10[2:(n0+1)])/(1 - u10[1:n0]));
  Mr10 = cummin(r10); 
  MR10 = c(1-(1-u10[1:n0])*Mr10,1);
  
  ub10=1-Mr10-R10[2:(n0+1)]+Mr10*(1:n0)/n0;
  lb10=1-Mr10-R10[2:(n0+1)]+Mr10*(0:(n0-1))/n0;
  
  c1 = sqrt((n0*n0)/(n0+n0)); 
  
  Bon.p1[b0] = c1*(sum((ub10^(1+1)-lb10^(1+1))[Mr10>0]/(1+1)/Mr10[Mr10>0]))^(1/1);
  Bon.p2[b0] = c1*(sum((ub10^(2+1)-lb10^(2+1))[Mr10>0]/(2+1)/Mr10[Mr10>0]))^(1/2);
  Bon.ps[b0] = c1*max(ub10); #print(b0)
  
}
alpha = 0.05; 
k = 2:7; 
Bon.cv.p1 = quantile(Bon.p1,1-alpha/(k-1)); Bon.cv.p1; 
Bon.cv.p2 = quantile(Bon.p2,1-alpha/(k-1)); Bon.cv.p2; 
Bon.cv.ps = quantile(Bon.ps,1-alpha/(k-1)); Bon.cv.ps; 

# (n=3000); 
# > Bon.cv.p1 = quantile(Bon.p1,1-alpha/(k-1)); Bon.cv.p1; 
# 95%     97.5% 98.33333%    98.75%       99% 99.16667% 
#   0.5787743 0.6567475 0.7048794 0.7283533 0.7541543 0.7732230 
# > Bon.cv.p2 = quantile(Bon.p2,1-alpha/(k-1)); Bon.cv.p2; 
# 95%     97.5% 98.33333%    98.75%       99% 99.16667% 
#   0.6731821 0.7596477 0.8093516 0.8414987 0.8702084 0.8854050 
# > Bon.cv.ps = quantile(Bon.ps,1-alpha/(k-1)); Bon.cv.ps; 
# 95%     97.5% 98.33333%    98.75%       99% 99.16667% 
#   1.339493  1.454522  1.525212  1.579907  1.613301  1.641733 



