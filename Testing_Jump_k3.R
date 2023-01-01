rm(list=ls(all=TRUE))
source("EGJ_USO_Library.r")

n = 100; nv = c(n,n,n); 
#nv = c(97, 176, 135,  67,  67); 
k = length(nv); 

set.seed(02132021)

B0 = 10000; 
TPs0 = array(,c(B0,3))
set.seed(101020220)
for(b0 in 1:B0){
  X10 = rMixUSO(n=nv[1],pl=0.0)
  X20 = rMixUSO(n=nv[2],pl=0.0)
  X30 = rMixUSO(n=nv[3],pl=0.0)
  
  R120 = Rmn(X10,X20); #LSMR120 = LSMRmn(R120); 
  R230 = Rmn(X20,X30); #LSMR230 = LSMRmn(R230); 
  
  D12s0 = DP(R120,m=nv[1],n=nv[2])
  D23s0 = DP(R230,m=nv[2],n=nv[3])
  
  Ds0 = rbind(D12s0,D23s0)

  TPs0[b0,] = apply(Ds0,2,max)
  print(b0)
}
alpha = 0.05
qB = apply(TPs0,2,function(x){quantile(x=x,1-alpha)})
B = 1000; #n = 200; n = 400;

Jump.p1 = array(NA,c(B,(k-1))); 
Jump.p2 = array(NA,c(B,(k-1))); 
Jump.ps = array(NA,c(B,(k-1))); 
# M1s.data = array(NA,c(B,(k-1))); 
# M2s.data = array(NA,c(B,(k-1))); 
# Mss.data = array(NA,c(B,(k-1))); 

Jump.p11 = array(NA,c(B,(k-1))); 
Jump.p21 = array(NA,c(B,(k-1))); 
Jump.ps1 = array(NA,c(B,(k-1))); 
# M1s.data1 = array(NA,c(B,(k-1))); 
# M2s.data1 = array(NA,c(B,(k-1))); 
# Mss.data1 = array(NA,c(B,(k-1))); 


jumps = c(0,0)
set.seed(02142021)
for(b in 1:B){
  X1 = rRLS_Jumps(n=nv[1],jumps=jumps,at=0);
  X2 = rRLS_Jumps(n=nv[2],jumps=jumps,at=1);
  X3 = rRLS_Jumps(n=nv[3],jumps=jumps,at=2);
  
  
  R12 = Rmn(X1,X2); #LSMR120 = LSMRmn(R120); 
  R23 = Rmn(X2,X3); #LSMR230 = LSMRmn(R230); 
  
  D12s = DP(R12,nv[1],nv[2])
  D23s = DP(R23,nv[2],nv[3])
  
  M12s = MP(R12,nv[1],nv[2])
  M23s = MP(R23,nv[2],nv[3])
  
  Ds = rbind(D12s,D23s)
  Ms = rbind(M12s,M23s)
  
  k = length(nv)
  cs = array(,k-1)
  for(j in 1:(k-1)){
    cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1]))
  }
  
  Jump.p1[b,] = c(Ds[,1]>qB[1]);
  Jump.p2[b,] = c(Ds[,2]>qB[2]);
  Jump.ps[b,] = c(Ds[,3]>qB[3]); 
  
  delta.p1 = c(0,sort(Ds[,1]));  #c(0,sort(Ds[,1]/log(cs))); 
  delta.p2 = c(0,sort(Ds[,2])); 
  delta.ps = c(0,sort(Ds[,3])); 
  
  JP1 = array(,c(k,k-1)); QP1 = array(,k); 
  JP2 = array(,c(k,k-1)); QP2 = array(,k); 
  JPs = array(,c(k,k-1)); QPs = array(,k); 
  
  constp1 = log(log(cs))*2/3#0.652; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 0.652/0.585#log(log(cs))/2; #log(log(cs))/2
  constp2 = log(log(cs))*3/4#0.725; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 0.725/0.676#log(log(cs))/2; #log(log(cs))/2
  constps = log(log(cs))*1; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 1.201/1.353#log(log(cs))/2; #log(log(cs))/2
  
  for(j in 1:k){
    JP1[j,] = c(Ds[,1]>delta.p1[j]); 
    JP2[j,] = c(Ds[,2]>delta.p2[j]);
    JPs[j,] = c(Ds[,3]>delta.ps[j]);
    
    QP1[j] = sum( (1-JP1[j,])*(Ds[,1]/cs) + (JP1[j,])*((Ms[,1]+log(cs)*constp1)/cs));
    QP2[j] = sum( (1-JP2[j,])*(Ds[,2]/cs) + (JP2[j,])*((Ms[,2]+log(cs)*constp2)/cs));
    QPs[j] = sum( (1-JPs[j,])*(Ds[,3]/cs) + (JPs[j,])*((Ms[,3]+log(cs)*constps)/cs));
  }
  
  delta.p1.star = delta.p1[which.min(QP1)]
  delta.p2.star = delta.p2[which.min(QP2)]
  delta.ps.star = delta.ps[which.min(QPs)]
  
  Jump.p11[b,] = c(Ds[,1]>delta.p1.star);
  Jump.p21[b,] = c(Ds[,2]>delta.p2.star);
  Jump.ps1[b,] = c(Ds[,3]>delta.ps.star); 
}

exact.jump = jumps; 
exact.rate.p1 = 0; exact.rate.p2 = 0; exact.rate.ps = 0; 
TP.p1 = 0; TP.p2 = 0; TP.ps = 0; 
for(b in 1:B){
  exact.rate.p1 = all(Jump.p1[b,] == exact.jump)/B + exact.rate.p1; 
  exact.rate.p2 = all(Jump.p2[b,] == exact.jump)/B + exact.rate.p2; 
  exact.rate.ps = all(Jump.ps[b,] == exact.jump)/B + exact.rate.ps; 
}
TP.p1 = mean(Jump.p1%*%jumps); TP.p2 = mean(Jump.p2%*%jumps); TP.ps = mean(Jump.ps%*%jumps)
FP.p1 = mean(Jump.p1%*%(1-jumps)); FP.p2 = mean(Jump.p2%*%(1-jumps)); FP.ps = mean(Jump.ps%*%(1-jumps))

exact.rate.p1; 
exact.rate.p2; 
exact.rate.ps; 

TP.p1; TP.p2; TP.ps; 
FP.p1; FP.p2; FP.ps; 

### delta

exact.rate.p11 = 0; exact.rate.p21 = 0; exact.rate.ps1 = 0; 
TP.p11 = 0; TP.p21 = 0; TP.ps1 = 0; 
for(b in 1:B){
  exact.rate.p11 = all(Jump.p11[b,] == exact.jump)/B + exact.rate.p11; 
  exact.rate.p21 = all(Jump.p21[b,] == exact.jump)/B + exact.rate.p21; 
  exact.rate.ps1 = all(Jump.ps1[b,] == exact.jump)/B + exact.rate.ps1; 

}
TP.p11 = mean(Jump.p11%*%jumps); TP.p21 = mean(Jump.p21%*%jumps); TP.ps1 = mean(Jump.ps1%*%jumps)
FP.p11 = mean(Jump.p11%*%(1-jumps)); FP.p21 = mean(Jump.p21%*%(1-jumps)); FP.ps1 = mean(Jump.ps1%*%(1-jumps))

exact.rate.p11; 
exact.rate.p21; 
exact.rate.ps1; 

TP.p11; TP.p21; TP.ps1; 
FP.p11; FP.p21; FP.ps1; 


TPcomp = rbind(c(TP.p1, TP.p11),
               c(TP.p2, TP.p21),
               c(TP.ps, TP.ps1)
)
colnames(TPcomp) = c("original TP", "TP with delta_star")
rownames(TPcomp) = c("p=1", "p=2", "p=inf")

FPcomp = rbind(c(FP.p1, FP.p11),
               c(FP.p2, FP.p21),
               c(FP.ps, FP.ps1)
)
colnames(FPcomp) = c("original FP", "FP with delta_star")
rownames(FPcomp) = c("p=1", "p=2", "p=inf")


EXcomp = rbind(c(exact.rate.p1, exact.rate.p11),
               c(exact.rate.p2, exact.rate.p21),
               c(exact.rate.ps, exact.rate.ps1)
)
colnames(EXcomp) = c("Exact", "Exact with delta_star")
rownames(EXcomp) = c("p=1", "p=2", "p=inf")

TPcomp = round(TPcomp,3); TPcomp; 
FPcomp = round(FPcomp,3); FPcomp; 
EXcomp = round(EXcomp,3); EXcomp; 

rbind(apply(Jump.p1,2,mean),apply(Jump.p11,2,mean))
rbind(apply(Jump.p2,2,mean),apply(Jump.p21,2,mean))
rbind(apply(Jump.ps,2,mean),apply(Jump.ps1,2,mean))

#c
paste(EXcomp[1,1], "&", TPcomp[1,1], "&", FPcomp[1,1], "&", 
      EXcomp[2,1], "&", TPcomp[2,1], "&", FPcomp[2,1], "&", 
      EXcomp[3,1], "&", TPcomp[3,1], "&", FPcomp[3,1])
#delta
paste(EXcomp[1,2], "&", TPcomp[1,2], "&", FPcomp[1,2], "&",
      EXcomp[2,2], "&", TPcomp[2,2], "&", FPcomp[2,2], "&",
      EXcomp[3,2], "&", TPcomp[3,2], "&", FPcomp[3,2])



