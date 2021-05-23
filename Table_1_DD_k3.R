rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

### n=200
n = 200; 
nv = c(n,n,n); k = length(nv)


B0 = 10000; nv0 = nv
Wk10 = array(NA,B0); Wk20 = array(NA,B0); Wks0 = array(NA,B0); 
set.seed(02132021)
for(b0 in 1:B0){
  X0.data = list(X10 = runif(nv0[1]), 
                 X20 = runif(nv0[2]),
                 X30 = runif(nv0[3])
  )
  temp = DDDST(X0.data)
  Wk10[b0] = temp$MaxDps[1]; 
  Wk20[b0] = temp$MaxDps[2]; 
  Wks0[b0] = temp$MaxDps[3]; 
}

jump.cv.p1 = quantile(Wk10,0.95) 
jump.cv.p2 = quantile(Wk20,0.95) 
jump.cv.ps = quantile(Wks0,0.95) 


#######################################################

Table.K3.200 = array(,c(4,9)); 

B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 

Table.K3.200[1,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.200[2,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.200[3,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.200[4,] = c(C1,C2,C3); 

Table.K3.200




#######################################################

Table.K3.200 = array(,c(4,9)); 


B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 

Table.K3.200[1,] = c(C1,C2,C3); 


B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.200[3,] = c(C1,C2,C3); 

B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.200[4,] = c(C1,C2,C3); 

Table.K3.200

######################################################################################################

### n = 400

n = 400; 
nv = c(n,n,n); k = length(nv)
Table.K3.400 = array(,c(4,9)); 

B0 = 10000; nv0 = nv
Wk10 = array(NA,B0); Wk20 = array(NA,B0); Wks0 = array(NA,B0); 
set.seed(02132021)
for(b0 in 1:B0){
  X0.data = list(X10 = runif(nv0[1]), 
                 X20 = runif(nv0[2]),
                 X30 = runif(nv0[3])
  )
  temp = DDDST(X0.data)
  Wk10[b0] = temp$MaxDps[1]; 
  Wk20[b0] = temp$MaxDps[2]; 
  Wks0[b0] = temp$MaxDps[3]; 
}

jump.cv.p1 = quantile(Wk10,0.95) 
jump.cv.p2 = quantile(Wk20,0.95) 
jump.cv.ps = quantile(Wks0,0.95) 

###########################################################


B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 

Table.K3.400[1,] = c(C1,C2,C3); 


B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[2,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[3,] = c(C1,C2,C3); 


B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[4,] = c(C1,C2,C3); 

Table.K3.400



###########################################################

Table.K3.400 = array(,c(4,9)); 

B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 

Table.K3.400[1,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,0); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[2,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(0,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[3,] = c(C1,C2,C3); 



B = 1000; 
Jump.p1 = array(NA,c(B,(k-1))); Jump.p2 = array(NA,c(B,(k-1))); Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); M2s.data = array(NA,c(B,(k-1))); Mss.data = array(NA,c(B,(k-1))); 
set.seed(02142021)
for(b in 1:B){
  if(1==1){
    J = c(1,1); signal = 0.6; 
    exact.jump=J; 
    k = length(J) + 1; 
    p = rep(1,k) * 0.2; 
    n = nv
    m = signal*cumsum(c(1,J)) - signal + 2
    X.data = MixNormal(n=n,p=p,m=m)
  }
  temp = DDDST(X.data)
  Jump.p1[b,] = c(temp$D1s>jump.cv.p1);
  Jump.p2[b,] = c(temp$D2s>jump.cv.p2);
  Jump.ps[b,] = c(temp$Dss>jump.cv.ps); 
}

exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
under.rate.p1 = 0; under.rate.p2 = 0; under.rate.ps = 0; 
over.rate.p1 = 0; over.rate.p2 = 0; over.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
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
TP.p1 = mean(Jump.p1%*%J); TP.p2 = mean(Jump.p2%*%J); TP.ps = mean(Jump.ps%*%J)
FP.p1 = mean(Jump.p1%*%(1-J)); FP.p2 = mean(Jump.p2%*%(1-J)); FP.ps = mean(Jump.ps%*%(1-J))

C1 = c(exact.rate.p1, under.rate.p1, over.rate.p1); 
C2 = c(exact.rate.p2, under.rate.p2, over.rate.p2); 
C3 = c(exact.rate.ps, under.rate.ps, over.rate.ps); 
Table.K3.400[4,] = c(C1,C2,C3); 


Table.K3.400