rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")

n = 60
nv = c(n,n,n,n,n); k=length(nv); 

B0 = 10000; 
UPs0 = array(,c(B0,3)); TPs0 = array(,c(B0,3)); MPs0 = array(,c(B0,3))
set.seed(101020220)
for(b0 in 1:B0){
  X10 = runif(nv[1]); 
  X20 = runif(nv[2]); 
  X30 = runif(nv[3]); 
  X40 = runif(nv[4]); 
  X50 = runif(nv[5]); 
  
  R120 = Rmn(X10,X20); #LSMR120 = LSMRmn(R120); 
  R230 = Rmn(X20,X30); #LSMR230 = LSMRmn(R230); 
  R340 = Rmn(X30,X40); #LSMR340 = LSMRmn(R340); 
  R450 = Rmn(X40,X50); #LSMR450 = LSMRmn(R450); 
  
  D12s0 = DP(R120,nv[1],nv[2]);
  D23s0 = DP(R230,nv[2],nv[3]);
  D34s0 = DP(R340,nv[3],nv[4]);
  D45s0 = DP(R450,nv[4],nv[5]);
  
  Ds0 = rbind(D12s0,D23s0,D34s0,D45s0);
  UPs0[b0,] = apply(Ds0,2,max);
  TPs0[b0,] = apply(Ds0,2,sum);
  MPs0[b0,] = D12s0; 
  print(b0);
}
alpha = 0.05
qT = apply(TPs0,2,function(x){quantile(x=x,1-alpha)})
qU = apply(UPs0,2,function(x){quantile(x=x,1-alpha)})
qB = apply(MPs0,2,function(x){quantile(x=x,1-alpha/(k-1))})


B = 1000
jumps = c(0,0,0,0)
UPs = array(,c(B,3)); TPs = array(,c(B,3))
set.seed(101020221)
for(b in 1:B){
  X1 = rRLS_Jumps(n=nv[1],jumps=jumps,at=0);
  X2 = rRLS_Jumps(n=nv[2],jumps=jumps,at=1);
  X3 = rRLS_Jumps(n=nv[3],jumps=jumps,at=2);
  X4 = rRLS_Jumps(n=nv[4],jumps=jumps,at=3);
  X5 = rRLS_Jumps(n=nv[5],jumps=jumps,at=4);
  
  R12 = Rmn(X1,X2); LSMR12 = LSMRmn(R12); 
  R23 = Rmn(X2,X3); LSMR23 = LSMRmn(R23); 
  R34 = Rmn(X3,X4); LSMR34 = LSMRmn(R34); 
  R45 = Rmn(X4,X5); LSMR45 = LSMRmn(R45); 
  
  D12s = DP(R12,nv[1],nv[2])
  D23s = DP(R23,nv[2],nv[3])
  D34s = DP(R34,nv[3],nv[4])
  D45s = DP(R45,nv[4],nv[5])
  
  Ds = rbind(D12s,D23s,D34s,D45s)
  UPs[b,] = apply(Ds,2,max)
  TPs[b,] = apply(Ds,2,sum)
  print(b)
}

paste(mean(TPs[,1] > qT[1]), "&", mean(UPs[,1] > qU[1]), "&",
      mean(TPs[,2] > qT[2]), "&", mean(UPs[,2] > qU[2]), "&",
      mean(TPs[,3] > qT[3]), "&", mean(UPs[,3] > qU[3]))

