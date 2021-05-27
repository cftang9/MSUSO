rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

####################################################################################################
### n=200
n = 200; 
nv = c(n,n,n,n); k = length(nv)

B0 = 10000; nv0 = nv
Wk10 = array(NA,B0); Wk20 = array(NA,B0); Wks0 = array(NA,B0); 
set.seed(02132021)
for(b0 in 1:B0){
  X0.data = list(X10 = runif(nv0[1]), 
                 X20 = runif(nv0[2]),
                 X30 = runif(nv0[3]), 
                 X40 = runif(nv0[4])
  )
  temp = DDDST(X0.data)
  Wk10[b0] = max(temp$D1s); 
  Wk20[b0] = max(temp$D2s); 
  Wks0[b0] = max(temp$Dss); 
}

jump.cv.p1.200 = quantile(Wk10,0.95) 
jump.cv.p2.200 = quantile(Wk20,0.95) 
jump.cv.ps.200 = quantile(Wks0,0.95) 

#######################################################

B = 1000; 

#######################################################
### RR

Table.DD.K4.RR.200 = array(,c(4,9)); 
for(j in 1:4){
  Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
  set.seed(02142021)
  for(b in 1:B){
    X.data = Gen.RR(nv,case=j)
    temp = DDDST(X.data)
    Jump.p1[b,] = c(temp$D1s>jump.cv.p1.200);
    Jump.p2[b,] = c(temp$D2s>jump.cv.p2.200);
    Jump.ps[b,] = c(temp$Dss>jump.cv.ps.200); print(b)
  }
  
  exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
  under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
  over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
  TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
  
  exact.jump=c( (j-1)%%2, (j-1)%/%2, (j-1)%/%4); 
  
  for(b in 1:B){
    over.rate.p1 = (all(Jump.p1[b,] >= exact.jump) & any(Jump.p1[b,] > exact.jump))/B + over.rate.p1; 
    exact.rate.p1 = all(Jump.p1[b,] == exact.jump)/B + exact.rate.p1; 
    under.rate.p1 = (all(Jump.p1[b,] <= exact.jump) & any(Jump.p1[b,] < exact.jump))/B + under.rate.p1; 
    over.rate.p2 = (all(Jump.p2[b,] >= exact.jump) & any(Jump.p2[b,] > exact.jump))/B + over.rate.p2; 
    exact.rate.p2 = all(Jump.p2[b,] == exact.jump)/B + exact.rate.p2; 
    under.rate.p2 = (all(Jump.p2[b,] <= exact.jump) & any(Jump.p2[b,] < exact.jump))/B + under.rate.p2; 
    over.rate.ps = (all(Jump.ps[b,] >= exact.jump) & any(Jump.ps[b,] > exact.jump))/B + over.rate.ps; 
    exact.rate.ps = all(Jump.ps[b,] == exact.jump)/B + exact.rate.ps; 
    under.rate.ps = (all(Jump.ps[b,] <= exact.jump) & any(Jump.ps[b,] < exact.jump))/B + under.rate.ps;
  }
  TP.p1 = mean(Jump.p1%*%exact.jump); TP.p2 = mean(Jump.p2%*%exact.jump); TP.ps = mean(Jump.ps%*%exact.jump)
  FP.p1 = mean(Jump.p1%*%(1-exact.jump)); FP.p2 = mean(Jump.p2%*%(1-exact.jump)); FP.ps = mean(Jump.ps%*%(1-exact.jump))
  
  C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
  C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
  C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
  
  Table.DD.K4.RR.200[j,] = c(C1,C2,C3); 
}

Table.DD.K4.MN.200 = array(,c(8,9)); 
for(j in 1:8){
  Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
  set.seed(02142021)
  for(b in 1:B){
    X.data = Gen.DD(nv,case=j)
    temp = DDDST(X.data)
    Jump.p1[b,] = c(temp$D1s>jump.cv.p1.200);
    Jump.p2[b,] = c(temp$D2s>jump.cv.p2.200);
    Jump.ps[b,] = c(temp$Dss>jump.cv.ps.200); 
  }
  
  exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
  under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
  over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
  TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
  
  exact.jump=c( (j-1)%%2, (j-1)%/%2, (j-1)%/%4); 
  
  for(b in 1:B){
    over.rate.p1 = (all(Jump.p1[b,] >= exact.jump) & any(Jump.p1[b,] > exact.jump))/B + over.rate.p1; 
    exact.rate.p1 = all(Jump.p1[b,] == exact.jump)/B + exact.rate.p1; 
    under.rate.p1 = (all(Jump.p1[b,] <= exact.jump) & any(Jump.p1[b,] < exact.jump))/B + under.rate.p1; 
    over.rate.p2 = (all(Jump.p2[b,] >= exact.jump) & any(Jump.p2[b,] > exact.jump))/B + over.rate.p2; 
    exact.rate.p2 = all(Jump.p2[b,] == exact.jump)/B + exact.rate.p2; 
    under.rate.p2 = (all(Jump.p2[b,] <= exact.jump) & any(Jump.p2[b,] < exact.jump))/B + under.rate.p2; 
    over.rate.ps = (all(Jump.ps[b,] >= exact.jump) & any(Jump.ps[b,] > exact.jump))/B + over.rate.ps; 
    exact.rate.ps = all(Jump.ps[b,] == exact.jump)/B + exact.rate.ps; 
    under.rate.ps = (all(Jump.ps[b,] <= exact.jump) & any(Jump.ps[b,] < exact.jump))/B + under.rate.ps;
  }
  TP.p1 = mean(Jump.p1%*%exact.jump); TP.p2 = mean(Jump.p2%*%exact.jump); TP.ps = mean(Jump.ps%*%exact.jump)
  FP.p1 = mean(Jump.p1%*%(1-exact.jump)); FP.p2 = mean(Jump.p2%*%(1-exact.jump)); FP.ps = mean(Jump.ps%*%(1-exact.jump))
  
  C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
  C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
  C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
  
  Table.DD.K4.MN.200[j,] = c(C1,C2,C3); 
}


####################################################################################################
### n=400
n = 400; 
nv = c(n,n,n,n); k = length(nv)

B0 = 10000; nv0 = nv
Wk10 = array(NA,B0); Wk20 = array(NA,B0); Wks0 = array(NA,B0); 
set.seed(02132021)
for(b0 in 1:B0){
  X0.data = list(X10 = runif(nv0[1]), 
                 X20 = runif(nv0[2]),
                 X30 = runif(nv0[3]),
                 X40 = runif(nv0[4])
  )
  temp = DDDST(X0.data)
  Wk10[b0] = max(temp$D1s); 
  Wk20[b0] = max(temp$D2s); 
  Wks0[b0] = max(temp$Dss); print(b0)
}

jump.cv.p1.400 = quantile(Wk10,0.95) 
jump.cv.p2.400 = quantile(Wk20,0.95) 
jump.cv.ps.400 = quantile(Wks0,0.95) 

#######################################################

B = 1000; 

#######################################################
### RR

Table.DD.K4.RR.400 = array(,c(4,9)); 
for(j in 1:4){
  Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
  set.seed(02142021)
  for(b in 1:B){
    X.data = Gen.RR(nv,case=j)
    temp = DDDST(X.data)
    Jump.p1[b,] = c(temp$D1s>jump.cv.p1.400);
    Jump.p2[b,] = c(temp$D2s>jump.cv.p2.400);
    Jump.ps[b,] = c(temp$Dss>jump.cv.ps.400); print(b)
  }
  
  exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
  under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
  over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
  TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
  
  exact.jump=c( (j-1)%%2, (j-1)%/%2, (j-1)%/%4); 
  
  for(b in 1:B){
    over.rate.p1 = (all(Jump.p1[b,] >= exact.jump) & any(Jump.p1[b,] > exact.jump))/B + over.rate.p1; 
    exact.rate.p1 = all(Jump.p1[b,] == exact.jump)/B + exact.rate.p1; 
    under.rate.p1 = (all(Jump.p1[b,] <= exact.jump) & any(Jump.p1[b,] < exact.jump))/B + under.rate.p1; 
    over.rate.p2 = (all(Jump.p2[b,] >= exact.jump) & any(Jump.p2[b,] > exact.jump))/B + over.rate.p2; 
    exact.rate.p2 = all(Jump.p2[b,] == exact.jump)/B + exact.rate.p2; 
    under.rate.p2 = (all(Jump.p2[b,] <= exact.jump) & any(Jump.p2[b,] < exact.jump))/B + under.rate.p2; 
    over.rate.ps = (all(Jump.ps[b,] >= exact.jump) & any(Jump.ps[b,] > exact.jump))/B + over.rate.ps; 
    exact.rate.ps = all(Jump.ps[b,] == exact.jump)/B + exact.rate.ps; 
    under.rate.ps = (all(Jump.ps[b,] <= exact.jump) & any(Jump.ps[b,] < exact.jump))/B + under.rate.ps;
  }
  
  TP.p1 = mean(Jump.p1%*%exact.jump); TP.p2 = mean(Jump.p2%*%exact.jump); TP.ps = mean(Jump.ps%*%exact.jump)
  FP.p1 = mean(Jump.p1%*%(1-exact.jump)); FP.p2 = mean(Jump.p2%*%(1-exact.jump)); FP.ps = mean(Jump.ps%*%(1-exact.jump))
  
  C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
  C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
  C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
  
  Table.DD.K4.RR.400[j,] = c(C1,C2,C3); 
}


Table.DD.K4.MN.400 = array(,c(8,9)); 
for(j in 1:8){
  Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
  set.seed(02142021)
  for(b in 1:B){
    X.data = Gen.DD(nv,case=j)
    temp = DDDST(X.data)
    Jump.p1[b,] = c(temp$D1s>jump.cv.p1.200);
    Jump.p2[b,] = c(temp$D2s>jump.cv.p2.200);
    Jump.ps[b,] = c(temp$Dss>jump.cv.ps.200); print(b)
  }
  
  exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
  under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
  over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
  TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
  
  exact.jump=c( (j-1)%%2, (j-1)%/%2, (j-1)%/%4); 
  
  for(b in 1:B){
    over.rate.p1 = (all(Jump.p1[b,] >= exact.jump) & any(Jump.p1[b,] > exact.jump))/B + over.rate.p1; 
    exact.rate.p1 = all(Jump.p1[b,] == exact.jump)/B + exact.rate.p1; 
    under.rate.p1 = (all(Jump.p1[b,] <= exact.jump) & any(Jump.p1[b,] < exact.jump))/B + under.rate.p1; 
    over.rate.p2 = (all(Jump.p2[b,] >= exact.jump) & any(Jump.p2[b,] > exact.jump))/B + over.rate.p2; 
    exact.rate.p2 = all(Jump.p2[b,] == exact.jump)/B + exact.rate.p2; 
    under.rate.p2 = (all(Jump.p2[b,] <= exact.jump) & any(Jump.p2[b,] < exact.jump))/B + under.rate.p2; 
    over.rate.ps = (all(Jump.ps[b,] >= exact.jump) & any(Jump.ps[b,] > exact.jump))/B + over.rate.ps; 
    exact.rate.ps = all(Jump.ps[b,] == exact.jump)/B + exact.rate.ps; 
    under.rate.ps = (all(Jump.ps[b,] <= exact.jump) & any(Jump.ps[b,] < exact.jump))/B + under.rate.ps;
  }
  TP.p1 = mean(Jump.p1%*%exact.jump); TP.p2 = mean(Jump.p2%*%exact.jump); TP.ps = mean(Jump.ps%*%exact.jump)
  FP.p1 = mean(Jump.p1%*%(1-exact.jump)); FP.p2 = mean(Jump.p2%*%(1-exact.jump)); FP.ps = mean(Jump.ps%*%(1-exact.jump))
  
  C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
  C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
  C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
  
  Table.DD.K4.MN.400[j,] = c(C1,C2,C3); 
}

### Table 1. 
Table.DD.K4.RR.200
Table.DD.K4.RR.400
Table.DD.K4.MN.200
Table.DD.K4.MN.400
