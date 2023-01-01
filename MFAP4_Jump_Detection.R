rm(list=ls(all=TRUE))
working.dir = getwd(); 
setwd(".."); parent = getwd(); 
setwd(working.dir)
source(paste0(parent,"/MUSOLibrary.R"))
#source(paste0(parent,"/CorreCurves.R"))
library(readxl)
raw.data <- read_excel("MFAP4.xlsx")

data = cbind(raw.data$`Fibrosis Stage`, raw.data$`MFAP4 U/mL`)
X1 = sort(data[data[,1]==(1-1), 2]); 
X2 = sort(data[data[,1]==(2-1), 2]); 
X3 = sort(data[data[,1]==(3-1), 2]); 
X4 = sort(data[data[,1]==(4-1), 2]); 
X5 = sort(data[data[,1]==(5-1), 2]); 

X.data = list(X1 = X1,
              X2 = X2,
              X3 = X3,
              X4 = X4,
              X5 = X5
)
k = length(X.data); 
nv = array(NA,k);
for(j in 1:k){nv[j] = length(X.data[[j]])}
us = list(); cs = 0; 
for(j in 1:(k-1)){
  us[[j]] = seq(0,1,by=1/nv[j+1]); 
  cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1])); 
}

#nv = c(97, 176, 135,  67,  67)

nv0 = nv

set.seed(02132021)

B0 = 10000; 
M1s0 = array(NA,c(B0,(k-1))); M2s0 = array(NA,c(B0,(k-1))); Mss0 = array(NA,c(B0,(k-1))); 
for(b0 in 1:B0){
  X0.data = list(X10 = runif(nv0[1]), 
                 X20 = runif(nv0[2]),
                 X30 = runif(nv0[3]), 
                 X40 = runif(nv0[4]),
                 X50 = runif(nv0[5])
  )
  M1s = NA; M2s = NA; Mss = NA; 
  us = list(); fun.Mrs = list(); cs = 0; 
  for(j in 1:(k-1)){
    us[[j]] = seq(0,1,by=1/nv0[j+1]); 
    cs[j] = sqrt((nv0[j]*nv0[j+1])/(nv0[j]+nv0[j+1])); 
  }
  for(j in 1:(k-1)){
    ### loading data and parameters
    X = X0.data[[j]]; 
    Y = X0.data[[j+1]]; 
    uj = us[[j]]; 
    cj = cs[j]; 
    ### calcuating test statistics 
    n.uj = length(uj); n.Y = length(Y); 
    emp.Rj = ecdf(X)(quantile(Y,uj)); 
    r.emp.Rj = (1-emp.Rj[2:n.uj])/(1-uj[1:n.Y]); 
    Mrj = cummin(r.emp.Rj); 
    MRj = 1-Mrj[1:n.Y]*(1-uj[1:n.Y]); MRj[n.uj] = 1; 
    
    ubj = (1-(1:n.Y)/n.Y); 
    lbj = (1-(0:(n.Y-1))/n.Y); 
    M1j = cj*(-(sum( (1-Mrj)^1/(1+1)*(ubj^(1+1)-lbj^(1+1)))))^(1/1)
    M2j = cj*(-(sum( (1-Mrj)^2/(2+1)*(ubj^(2+1)-lbj^(2+1)))))^(1/2)
    Msj = cj * max((1-Mrj)*lbj)
    
    ### saving statsitics and slope function r
    M1s[j] = M1j; 
    M2s[j] = M2j; 
    Mss[j] = Msj; 
    fun.Mrs[[j]] = approxfun(uj,c(Mrj,Mrj[n.Y]),f=1); 
  }
  M1s0[b0,] = M1s; 
  M2s0[b0,] = M2s; 
  Mss0[b0,] = Mss; 
}

Wk10 = apply(M1s0, 1, max); Sk10 = apply(M1s0, 1, sum); 
Wk20 = apply(M2s0, 1, max); Sk20 = apply(M2s0, 1, sum); 
Wks0 = apply(Mss0, 1, max); Sks0 = apply(Mss0, 1, sum); 

jump.cv.p1 = quantile(Wk10,0.95); quantile(Sk10,0.95);  
jump.cv.p2 = quantile(Wk20,0.95); quantile(Sk10,0.95); 
jump.cv.ps = quantile(Wks0,0.95); quantile(Sk10,0.95); 
# jump.cv.p1 # 0.8297206 for n=50 # 0.8209637 for n=100 # 0.8057219 for n=200 # 0.8164940 for n=400 # 0.8126778 for n=1000 
# jump.cv.p2 # 0.9131916 for n=50 # 0.9055811 for n=100 # 0.8941796 for n=200 # 0.9013435 for n=400 # 0.8937240 for n=1000 
# jump.cv.ps # 1.5000000 for n=50 # 1.4849240 for n=100 # 1.5000000 for n=200 # 1.4849240 for n=400 # 1.4758050 for n=1000 

B = 1; #n = 200; n = 400;

Jump.p1 = array(NA,c(B,(k-1))); 
Jump.p2 = array(NA,c(B,(k-1))); 
Jump.ps = array(NA,c(B,(k-1))); 
M1s.data = array(NA,c(B,(k-1))); 
M2s.data = array(NA,c(B,(k-1))); 
Mss.data = array(NA,c(B,(k-1))); 

Jump.p11 = array(NA,c(B,(k-1))); 
Jump.p21 = array(NA,c(B,(k-1))); 
Jump.ps1 = array(NA,c(B,(k-1))); 
M1s.data1 = array(NA,c(B,(k-1))); 
M2s.data1 = array(NA,c(B,(k-1))); 
Mss.data1 = array(NA,c(B,(k-1))); 

M21p1s.data1 = array(NA,c(B,(k-1))); 
M21p2s.data1 = array(NA,c(B,(k-1))); 
M21pss.data1 = array(NA,c(B,(k-1))); 


set.seed(02142021)
for(b in 1:B){
  us = list(); fun.Mrs = list(); cs = 0; 
  for(j in 1:(k-1)){
    us[[j]] = seq(0,1,by=1/nv[j+1]); 
    cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1])); 
  }
  
  for(j in 1:(k-1)){
    ### loading data and parameters
    X = X.data[[j]]; 
    Y = X.data[[j+1]]; 
    uj = us[[j]]; 
    cj = cs[j]; 
    ### calcuating test statistics 
    n.uj = length(uj); n.Y = length(Y); 
    emp.Rj = ecdf(X)(quantile(Y,uj)); 
    r.emp.Rj = (1-emp.Rj[2:n.uj])/(1-uj[1:n.Y]); 
    Mrj = cummin(r.emp.Rj); 
    MRj = 1-Mrj[1:n.Y]*(1-uj[1:n.Y]); MRj[n.uj] = 1; 
    
    ubj = (1-(1:n.Y)/n.Y); 
    lbj = (1-(0:(n.Y-1))/n.Y); 
    M1j = cj*(-(sum( (1-Mrj)^1/(1+1)*(ubj^(1+1)-lbj^(1+1)))))^(1/1)
    M2j = cj*(-(sum( (1-Mrj)^2/(2+1)*(ubj^(2+1)-lbj^(2+1)))))^(1/2)
    Msj = cj * max((1-Mrj)*lbj)
    
    ### saving statsitics and slope function r
    M1s.data[b,j] = M1j; 
    M2s.data[b,j] = M2j; 
    Mss.data[b,j] = Msj; 
    
    rmn = 0; rmn = Rmn(X,Y); 
    LSRmn = 0; LSrmn = LSMRmn(rmn); 
    M21p1s.data1[b,j] = Mp21(rmn,LSrmn,m=length(X),n=length(Y),p=1)
    M21p2s.data1[b,j] = Mp21(rmn,LSrmn,m=length(X),n=length(Y),p=2)
    M21pss.data1[b,j] = Ms21(rmn,LSrmn,m=length(X),n=length(Y))
    
  }
  
  Jump.p1[b,] = c(M1s.data[b,]>jump.cv.p1);
  Jump.p2[b,] = c(M2s.data[b,]>jump.cv.p2);
  Jump.ps[b,] = c(Mss.data[b,]>jump.cv.ps); 
  
  
  
  delta.p1 = c(0,sort(M1s.data[b,])); 
  delta.p2 = c(0,sort(M2s.data[b,])); 
  delta.ps = c(0,sort(Mss.data[b,])); 
  
  JP1 = array(,c(k,k-1)); QP1 = array(,k); 
  JP2 = array(,c(k,k-1)); QP2 = array(,k); 
  JPs = array(,c(k,k-1)); QPs = array(,k); 
  
  
  constp1 = log(log(cs))*2/3 #0.652; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 0.652/0.585#log(log(cs))/2; #log(log(cs))/2
  constp2 = log(log(cs))*3/4#log(log(cs)) #0.725; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 0.725/0.676#log(log(cs))/2; #log(log(cs))/2
  constps = log(log(cs))*1#log(log(cs)) #1.201; #log(log(cs)); #log(log(cs)) #1/2#sqrt(log(log(cs))/2) # 1.201/1.353#log(log(cs))/2; #log(log(cs))/2
  
  for(j in 1:k){
    JP1[j,] = c(M1s.data[b,]>delta.p1[j]); 
    JP2[j,] = c(M2s.data[b,]>delta.p2[j]);
    JPs[j,] = c(Mss.data[b,]>delta.ps[j]);
    
    QP1[j] = sum( (1-JP1[j,])*(M1s.data[b,]/1) + (JP1[j,])*((M21p1s.data1+log(cs)*constp1)/1));
    QP2[j] = sum( (1-JP2[j,])*(M2s.data[b,]/1) + (JP2[j,])*((M21p2s.data1+log(cs)*constp2)/1));
    QPs[j] = sum( (1-JPs[j,])*(Mss.data[b,]/1) + (JPs[j,])*((M21pss.data1+log(cs)*constps)/1));
  }
  
  delta.p1.star = delta.p1[which.min(QP1)]
  delta.p2.star = delta.p2[which.min(QP2)]
  delta.ps.star = delta.ps[which.min(QPs)]
  
  Jump.p11[b,] = c(M1s.data[b,]>delta.p1.star);
  Jump.p21[b,] = c(M2s.data[b,]>delta.p2.star);
  Jump.ps1[b,] = c(Mss.data[b,]>delta.ps.star); 
}

JP1comp = rbind(Jump.p1, Jump.p11)
rownames(JP1comp) = c("original", "+delta")
colnames(JP1comp) = c("12", "23", "34", "45")

JP2comp = rbind(Jump.p2, Jump.p21)
rownames(JP2comp) = c("original", "+delta")
colnames(JP2comp) = c("12", "23", "34", "45")

JPscomp = rbind(Jump.ps, Jump.ps1)
rownames(JPscomp) = c("original", "+delta")
colnames(JPscomp) = c("12", "23", "34", "45")

JP1comp
JP2comp
JPscomp

par(mfrow=c(1,3))
plot(delta.p1,QP1,type="s")
plot(delta.p2,QP2,type="s")
plot(delta.ps,QPs,type="s")
par(mfrow=c(1,1))

