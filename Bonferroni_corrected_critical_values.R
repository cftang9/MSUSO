rm(list=ls(all=TRUE))
source("MUSOLibrary.R"); 

B0 = 10000; n0 = 3000
Bon.p1 = array(,B0); Bon.p2 = array(,B0); Bon.ps = array(,B0); 
set.seed(05242021)
for(b0 in 1:B0){
  X0 = runif(n0); Y0 = runif(n0); 
  uj0 = seq(0,1,by=1/n0); cj0 = sqrt((n0*n0)/(n0+n0)); 
  
  ### calcuating test statistics 
  n.uj0 = length(uj0); n.Y0 = length(Y0); 
  emp.Rj0 = Rmna(X0,Y0) 
  r.emp.Rj0 = (1-emp.Rj0[1:n.Y0])/(1-uj0[1:n.Y0]); 
  Mrj0 = cummin(r.emp.Rj0); 
  MRj0 = 1-Mrj0[1:n.Y0]*(1-uj0[1:n.Y0]); MRj0[n.uj0] = 1; 
  
  ubj0 = (MRj0[1:n.Y0] + Mrj0/n.Y0 - emp.Rj0[1:n.Y0]); 
  lbj0 = (MRj0[1:n.Y0] - emp.Rj0[1:n.Y0]);
  Bon.p1[b0] = cj0*(sum((ubj0^(1+1)-lbj0^(1+1))[Mrj0>0]/(1+1)/Mrj0[Mrj0>0]))^(1/1);
  Bon.p2[b0] = cj0*(sum((ubj0^(2+1)-lbj0^(2+1))[Mrj0>0]/(2+1)/Mrj0[Mrj0>0]))^(1/2);
  Bon.ps[b0] = cj0*max(ubj0,lbj0); print(b0)
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



