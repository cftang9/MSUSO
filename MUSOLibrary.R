#### Preparing ODCs ####
#### Linear part and inverse of the linear part. ####
Li = function(u,a,b,s){
  La = 1-(1-a)*s; Lb = 1-(1-b)*s; 
  Ind = c(a<=u&u<=b); 
  Li = (La + s*(u-a))*Ind
  return(Li); 
}
InvLi = function(v,a,b,s){
  La = 1-(1-a)*s; Lb = 1-(1-b)*s; 
  Ind = c(La<=v&v<=Lb); 
  InvLi = (a + (v-La)/s)*Ind
  return(InvLi); 
}
#### Quadratic Bezier function and it's inverse function. ####
tu = function(u,xv3){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  tu = array(,length(u))
  if(2*x2 == x1+x3){print("Warning!")}
  if(2*x2 != x1+x3){
    tu = ( x1 - x2 + 
             sqrt( abs(x2^2 - x1*x3 + (x1 - 2*x2 + x3)*u) )
    ) / (x1 - 2*x2 + x3) * 
      c(x1 <= u & u <= x3)
  }
  return(tu)
}
QB = function(u,xv3,yv3){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  w = tu(u,xv3); 
  QB = ((1-w)^2*y1 + 2*(1-w)*w*y2 + w^2*y3)*c(x1<=u&u<=x3); 
  return(QB)
}
InvQB = function(v,xv3,yv3){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  w = tu(v,yv3); 
  InvQB = ((1-w)^2*x1 + 2*(1-w)*w*x2 + w^2*x3)*c(y1<=v&v<=y3); 
  return(InvQB)
}
GlueQB = function(u,b1,a2,s1,s2,d=1/3){
  distance = a2-b1; 
  yb1 = 1 - (1-b1)*s1; ya2 = 1- (1-a2)*s2; 
  b1r = b1 + distance*d; yb1r = 1 - (1-b1r)*s1; 
  a2l = a2 - distance*d; ya2l = 1 - (1-a2l)*s2; 
  mab = (b1r+a2l)/2; ymab = (yb1r+ya2l)/2; 
  QB1 = QB(u,c(b1,b1r,mab),c(yb1,yb1r,ymab));
  QB2 = QB(u,c(mab,a2l,a2),c(ymab,ya2l,ya2));
  GlueQB = QB1 + QB2*c(u>mab); 
  return(GlueQB)
}
InvGlueQB = function(v,b1,a2,s1,s2,d=1/3){
  distance = a2-b1; 
  yb1 = 1 - (1-b1)*s1; ya2 = 1- (1-a2)*s2; 
  b1r = b1 + distance*d; yb1r = 1 - (1-b1r)*s1; 
  a2l = a2 - distance*d; ya2l = 1 - (1-a2l)*s2; 
  mab = (b1r+a2l)/2; ymab = (yb1r+ya2l)/2; 
  InvQB1 = InvQB(v,c(b1,b1r,mab),c(yb1,yb1r,ymab));
  InvQB2 = InvQB(v,c(mab,a2l,a2),c(ymab,ya2l,ya2));
  InvGlueQB = InvQB1 + InvQB2*c(v>ymab); 
  return(InvGlueQB)
}


#### Nonstrictly star-shaped ODC and it's inverse. ####
USOR = function(u,av,bv,sv){ #not include zero 
  nu = length(u); na = length(av); 
  USOR = array(0,nu)
  for(linear in 1:na){
    USOR = Li(u,av[linear],bv[linear],sv[linear]) + USOR; 
  }
  for(strictly in 1:(na-1)){
    USOR = GlueQB(u,bv[strictly],av[strictly+1],sv[strictly],sv[strictly+1])*c(u<av[strictly+1] & u>bv[strictly]) + USOR; 
  }
  return(USOR)
}
InvUSOR = function(v,av,bv,sv){ #not include zero
  nv = length(v); na = length(av); 
  yav = 1-(1-av)*sv; ybv = 1-(1-bv)*sv; 
  InvUSOR = array(0,nv)
  for(linear in 1:na){
    InvUSOR = InvLi(v,av[linear],bv[linear],sv[linear]) + InvUSOR; 
  }
  for(strictly in 1:(na-1)){
    Ind = c(v<yav[strictly+1] & v>ybv[strictly])
    InvUSOR[Ind] = InvGlueQB(v[Ind],bv[strictly],av[strictly+1],sv[strictly],sv[strictly+1]); 
  }
  return(InvUSOR)
}
rUSOR = function(n,av,bv,sv){
  rUSOR = InvUSOR(runif(n),av,bv,sv);
  return(rUSOR); 
}


#### ############################################

NSR = function(u,delta1=0.8,delta2=0.8){
  nu = length(u); #nv = length(v); 
  r = (1+delta1)/(1-delta1); r.inv = (1-delta1)/(1+delta1); 
  cc = 2*delta1/(1+delta1); 
  a = 4/(5*r+4); b = 1/5; 
  u1 = 7/8*(1-delta2); u2 = a + (1-a)*(1-delta2)
  v1 = 7/8*(1-delta2) + (1/8)*(a+(1-a)*(1-delta2)) ; v2 = b + (1-b)*(1-delta2);
  m1 = 1-delta2; m2 = (r-1 - m1*r*(1-r))/(r^2-1); 
  NSR = (USOR(u,c(0,v2), c(u1,1), c(1,0))-u)*(1-1/r)+u
  return(NSR)
}

QBP = function(u,xv3,yv3,prop){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  QBP = array(0,length(u));
  Ind = c(x1<=u&u<=x3); 
  QBP[Ind] = (QB(u[Ind],xv3,yv3)-u[Ind])*prop + u[Ind]
  #w = tu(u[Ind],xv3); 
  #QBP[Ind] = prop*(((1-w)^2*y1 + 2*(1-w)*w*y2 + w^2*y3) - u[Ind]) + u[Ind]; 
  return(QBP)
}

InvQBP = function(v,xv3,yv3,prop){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  
  c1 = prop/((x1+x3-2*x2)^2);
  c2 = y1*(x3-x2)^2 + y3*(x1-x2)^2 + 2*y2*(x1-x2)*(x3-x2) + (x2^2-x1*x3)*(y1-2*y2+y3)
  c3 = -y1*(x3-x2) + y3*(x1-x2) + y2*(x3-x1); 
  c4 = x2^2-x1*x3; 
  c5 = x1+x3-2*x2; 
  c6 =  ( (y1+y3-2*y2) + ((1-prop)/prop)*(x1+x3-2*x2) )*(x1+x3-2*x2) ; 
  
  c10 = 1/((x1+x3-2*x2)^2);
  c60 = (y1+y3-2*y2)*(x1+x3-2*x2); 
  
  IndP = c( abs( 1/c5*( ( (y1*c5/(c10*c60) - ( c2*c5/c60 - (c3*c5/c60)^2 -c4 ))^(1/2) - c3*c5/c60 )^2 - c4 ) - x1) < 10^(-10) & 
              abs( 1/c5*( ( (y3*c5/(c10*c60) - ( c2*c5/c60 - (c3*c5/c60)^2 -c4 ))^(1/2) - c3*c5/c60 )^2 - c4 ) - x3) < 10^(-10)
  )
  
  py1 = (y1 - x1)*prop + x1; 
  py3 = (y3 - x3)*prop + x3; 
  
  InvQBP = array(0,length(v)); 
  if(IndP == 1){
    Ind = c(py1<=v&v<=py3); 
    InvQBP[Ind] = 1/c5*( ( (1)*(v[Ind]*c5/(c1*c6) - ( c2*c5/c6 - (c3*c5/c6)^2 -c4 ))^(1/2) - c3*c5/c6 )^2 - c4 ); 
  }
  if(IndP != 1){
    Ind = c(py1<=v&v<=py3); 
    InvQBP[Ind] = 1/c5*( ( (-1)*(v[Ind]*c5/(c1*c6) - ( c2*c5/c6 - (c3*c5/c6)^2 -c4 ))^(1/2) - c3*c5/c6 )^2 - c4 ); 
  }
  return(InvQBP)
}

InvQBPP = function(v,xv3,yv3,prop){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  c1 = prop/((x1+x3-2*x2)^2);
  c2 = y1*(x3-x2)^2 + y3*(x1-x2)^2 + 2*y2*(x1-x2)*(x3-x2) + (x2^2-x1*x3)*(y1-2*y2+y3)
  c3 = -y1*(x3-x2) + y3*(x1-x2) + y2*(x3-x1); 
  c4 = x2^2-x1*x3; 
  c5 = x1+x3-2*x2; 
  c6 = (y1+y3-2*y2)*(x1+x3-2*x2) + ((1-prop)/prop)*(x1+x3-2*x2)^2; 
  
  c10 = 1/((x1+x3-2*x2)^2);
  c60 = (y1+y3-2*y2)*(x1+x3-2*x2); 
  
  py1 = (y1 - x1)*prop + x1; 
  py2 = (y2 - x2)*prop + x2; 
  py3 = (y3 - x3)*prop + x3; 
  
  InvQBPP = array(0,length(v)); 
  Ind = c(py1<=v&v<=py3); 
  InvQBPP[Ind] = 1/c5*( ( (1)*(v[Ind]*c5/(c1*c6) - ( c2*c5/c6 - (c3*c5/c6)^2 -c4 ))^(1/2) - c3*c5/c6 )^2 - c4 ); 
  return(InvQBPP)
}
InvQBPN = function(v,xv3,yv3,prop){
  x1 = xv3[1]; x2 = xv3[2]; x3 = xv3[3]; 
  y1 = yv3[1]; y2 = yv3[2]; y3 = yv3[3]; 
  c1 = prop/((x1+x3-2*x2)^2);
  c2 = y1*(x3-x2)^2 + y3*(x1-x2)^2 + 2*y2*(x1-x2)*(x3-x2) + (x2^2-x1*x3)*(y1-2*y2+y3)
  c3 = -y1*(x3-x2) + y3*(x1-x2) + y2*(x3-x1); 
  c4 = x2^2-x1*x3; 
  c5 = x1+x3-2*x2; 
  c6 = (y1+y3-2*y2)*(x1+x3-2*x2) + ((1-prop)/prop)*(x1+x3-2*x2)^2; 
  
  c10 = 1/((x1+x3-2*x2)^2);
  c60 = (y1+y3-2*y2)*(x1+x3-2*x2); 
  
  py1 = (y1 - x1)*prop + x1; 
  py2 = (y2 - x2)*prop + x2; 
  py3 = (y3 - x3)*prop + x3; 
  
  InvQBPN = array(0,length(v)); 
  Ind = c(py1<=v&v<=py3); 
  InvQBPN[Ind] = 1/c5*( ( (-1)*(v[Ind]*c5/(c1*c6) - ( c2*c5/c6 - (c3*c5/c6)^2 -c4 ))^(1/2) - c3*c5/c6 )^2 - c4 ); 
  return(InvQBPN)
}

GlueQBP = function(u,b1,a2,s1,s2,prop,d=1/3){
  #s1 = 1; s2 = 0;
  distance = a2 - b1; 
  yb1 = 1 - (1-b1)*s1; ya2 = 1- (1-a2)*s2; 
  b1r = b1 + distance*d; yb1r = 1 - (1-b1r)*s1; 
  a2l = a2 - distance*d; ya2l = 1 - (1-a2l)*s2; 
  mab = (b1r+a2l)/2; ymab = (yb1r+ya2l)/2; 
  
  QB1 = QBP(u,c(b1,b1r,mab),c(yb1,yb1r,ymab),prop);
  QB2 = QBP(u,c(mab,a2l,a2),c(ymab,ya2l,ya2),prop);
  GlueQBP = QB1 + QB2*c(u>mab); 
  return(GlueQBP)
}

InvGlueQBP = function(v,b1,a2,s1,s2,prop,d=1/3){
  #s1 = 1; s2 = prop; 
  distance = a2 - b1; 
  yb1 = 1 - (1-b1)*s1; ya2 = 1 - (1-a2)*s2; 
  b1r = b1 + distance*d; yb1r = 1 - (1-b1r)*s1; 
  a2l = a2 - distance*d; ya2l = 1 - (1-a2l)*s2; 
  mab = (b1r+a2l)/2; ymab = (yb1r+ya2l)/2; 
  
  pyb1 = prop*(yb1 - b1) + b1; 
  pya2 = prop*(ya2 - a2) + a2;
  pyb1r = prop*(yb1r - b1r) + b1r; 
  pya2l = prop*(ya2l - a2l) + a2l; 
  pymab = prop*(ymab - mab) + mab; 
  
  Ind1 = c(pyb1<=v & v<=pymab); 
  Ind2 = c(pymab<v & v<=pya2); 
  
  InvGlueQBP2 = array(0,length(v))
  InvGlueQBP2[Ind1] = InvQBPN(v[Ind1],c(b1,b1r,mab),c(yb1,yb1r,ymab),prop) + InvGlueQBP2[Ind1]; 
  InvGlueQBP2[Ind2] = InvQBPN(v[Ind2],c(mab,a2l,a2),c(ymab,ya2l,ya2),prop) + InvGlueQBP2[Ind2]; 
  
  return(InvGlueQBP2)
}

InvNSR = function(v,delta1=0.8,delta2=0.8){
  nv = length(v); 
  r = (1+delta1)/(1-delta1); r.inv = (1-delta1)/(1+delta1); 
  cc = 2*delta1/(1+delta1); 
  a = 4/(5*r+4); b = 1/5; 
  u1 = 7/8*(1-delta2); u2 = a + (1-a)*(1-delta2)
  v1 = 7/8*(1-delta2) + (1/8)*(a+(1-a)*(1-delta2)) ; v2 = b + (1-b)*(1-delta2);
  m1 = 1-delta2; m2 = (r-1 - m1*r*(1-r))/(r^2-1); 
  
  a2 = v2; b1 = u1; s1 = 1; s2 = 1/r;
  
  distance = a2 - b1; 
  yb1 = 1 - (1-b1)*s1; ya2 = 1- (1-a2)*s2; 
  #b1r = b1 + distance*d; yb1r = 1 - (1-b1r)*s1; 
  #a2l = a2 - distance*d; ya2l = 1 - (1-a2l)*s2; 
  #mab = (b1r+a2l)/2; ymab = (yb1r+ya2l)/2; 
  
  prop = 1 - 1/r; 
  
  #pyb1 = prop*(yb1 - b1) + b1; 
  #pya2 = prop*(ya2 - a2) + a2; 
  
  #pyb1r = prop*(yb1r - b1r) + b1r; 
  #pya2l = prop*(ya2l - a2l) + a2l; 
  #pymab = prop*(ymab - mab) + mab; 
  
  Ind = c(yb1<v & v<ya2); 
  
  InvNSR = array(0,nv); 
  InvNSR = Li(v,0,b1,1) + InvNSR; 
  InvNSR = InvLi(v,a2,1,(1-prop)) + InvNSR; 
  
  InvNSR[Ind] = InvGlueQBP(v[Ind],b1,a2,s1=1,s2=0,prop); 
  
  return(InvNSR)
}

Rmn = function(X,Y){
  m = length(X)
  n = length(Y)
  Rmn = c(array(0,(n)),1)
  for (i in 2:(n+1)){
    Rmn[i] = sum(X<=Y[i-1])/m
  }
  return(Rmn)
}
LSMRmn = function(Rmn_data){
  m = length(Rmn_data)-1
  u = seq(0,1,by=1/m)
  u_slope = array(0,m+1)
  LSMRmn = array(1,m+1)
  alpha = array(0,m+1)
  u_slope[1:m] = (1-Rmn_data[1:m])/(1-u[1:m])
  alpha = cummin(u_slope)
  LSMRmn = 1-(1-u)*alpha
  return(LSMRmn)
}

Mp21 = function(Rmn_data, LSMRmn_data,m,n,p=1){
  cmn = sqrt(m*n/(m+n))
  #Mp21 = cmn*TrapArea(LSMRmn_data-Rmn_data,p)
  u = seq(0,1,by=1/n)
  slope = (((1-LSMRmn_data[1:n])/(1-u[1:n])))
  nss = sum(slope>0)
  d = (LSMRmn_data[1:nss] - Rmn_data[1:nss])
  mid = ((slope[1:nss]/n+d[1:nss])^(1+p) - (d[1:nss])^(1+p))/(slope[1:nss]*(p+1))
  Mp21 = cmn*(sum(mid))^(1/p)
  return(Mp21)
}
Ms21 = function(Rmn_data, LSMRmn_data,m,n){
  cmn = sqrt(m*n/(m+n))
  u = seq(0,1,by=1/n)
  slope = (1-LSMRmn_data[1:n])/(1-u[1:n])
  Ms21 = cmn*(max(LSMRmn_data[1:n]+slope/n-Rmn_data[1:n]))
  return(Ms21)
}


# Park's test statistics and critical values. 
ParkCV = function(N=2,nk=2){
  TABLE=c(2.706, 4.231, 5.435, 6.498, 7.480, 8.407, 9.295, 10.152, 10.985,
          4.577, 7.059, 9.177, 11.127, 12.974, 14.750, 16.475, 18.159, 19.809,
          6.175, 9.669, 12.752, 15.637, 18.400, 21.080, 23.904, 26.671, 29.393,
          7.665, 12.196, 16.268, 20.115, 23.925, 27.827, 31.647, 35.403, 39.108,
          9.094, 14.675, 19.751, 24.585, 29.627, 34.546, 39.374, 44.131, 48.832
  );
  TABLE = array(TABLE,c(9,5))
  ParkCV = TABLE[nk-1,N-1];
  return(ParkCV); 
}
ParkLRT = function(X,tk){
  N = length(X[,1]); n = length(X[1,]); k = length(tk)-1; 
  mij = array(,c(N,k)); nij = array(,c(N,k)); 
  #X = array(runif(N*n),c(N,n))
  for(i in 1:N){
    for(j in 1:k){
      nij[i,j] = sum(X[i,]>tk[j] & X[i,]<=tk[j+1])
    }
    mij[i,] = cumsum(nij[i,])
  }
  theta = (mij-nij)/mij ; 
  Isotheta = array(,c(N,k))
  for(j in 1:k){
    Isotheta[,j] = pava(theta[,j])
  }
  Isotheta = Isotheta[,-1]; 
  nij = nij[,-1]; mij = mij[,-1]; theta = theta[,-1];  
  LogNonLikelihood = sum((mij-nij)*log(theta) + nij*log(1-theta))
  LogIsoLikelihood = sum((mij-nij)*log(Isotheta) + nij*log(1-Isotheta))
  ParkLRT = -2*(LogIsoLikelihood-LogNonLikelihood)
  return(ParkLRT); 
}

#### Generating random variables from USO distributions ####

#### Test statistic functions ####
TestBon = function(Ms, cv){
  Ms.max = max(Ms); 
  TestBon = as.numeric(Ms.max>cv)
  return(TestBon)
}
TestMax = function(Ms, cv){
  Ms.max = max(Ms); 
  TestMax = as.numeric(Ms.max>cv)
  return(TestMax)
}
TestSum = function(Ms, cv){
  Ms.sum = sum(Ms); 
  TestSum = as.numeric(Ms.sum>cv)
  return(TestSum)
}
TestCom = function(Ms, cv, bv){
  k = length(Ms)+1; #number of samples
  TestCom = c( sum(Ms) > cv | max(Ms) > bv
  ); 
  return(TestCom)
}

TestSubBon = function(test.st, Bon.cvs, threshold){
  k = length(test.st)+1; #number of samples
  ActSample = c(test.st > threshold); na = sum(ActSample); 
  TestSubBon = F; 
  if(na > 0) TestSubBon = c(max(test.st[c(1:(k-1))*ActSample]) > Bon.cvs[na] );
  return(TestSubBon)
}
TestSubSum = function(test.st, cvs, threshold){
  k = length(test.st)+1; #number of samples
  ActSample = c(test.st > threshold); na = sum(ActSample);
  TestSubSum = F; 
  if(na > 0) TestSubSum = c(sum(test.st[c(1:(k-1))*ActSample]) > cvs[sum(2^((k-2):0)*ActSample)+1]); 
  return(TestSubSum)
}
TestSubMax = function(test.st, cvs, threshold){
  k = length(test.st)+1; #number of samples
  ActSample = c(test.st > threshold); na = sum(ActSample); 
  TestSubMax = F; 
  if(na > 0) TestSubMax = c(max(test.st[c(1:(k-1))*ActSample]) > cvs[sum(2^((k-2):0)*ActSample)+1]); 
  return(TestSubMax)
}
TestSubCom = function(test.st, sum.cvs, max.cvs, threshold){
  k = length(test.st)+1; #number of samples
  ActSample = c(test.st > threshold); na = sum(ActSample);
  TestSubCom = F; 
  if(na > 0){ 
    TestSubCom = c( sum(test.st[c(1:(k-1))*ActSample]) > sum.cvs[sum(2^((k-2):0)*ActSample)+1] | 
                      max(test.st[c(1:(k-1))*ActSample]) > max.cvs[sum(2^((k-2):0)*ActSample)+1]
    ); 
  }
  return(TestSubCom)
}

DLp12 = function(x,s,p){
  n = length(x); 
  x1 = x[-n]; x2 = x[-1]; s = s[-n]; 
  ind.positive.s = c(s<0); 
  x1 = x1[ind.positive.s]; x2 = x2[ind.positive.s]; s = s[ind.positive.s]
  DLp12 = (sum((-s)^p*((1-x1)^(p+1) - (1-x2)^(p+1))/(p+1)))^(1/p); 
  return(DLp12)
}

DLs12 = function(x,s,p){
  n = length(x); 
  x1 = x[-n]; x2 = x[-1]; s = s[-n]; 
  ind.positive.s = c(s<0); 
  x1 = x1[ind.positive.s]; x2 = x2[ind.positive.s]; s = s[ind.positive.s]
  DLs12 = max( 0-s*(1-x1), 0-s*(1-x2)); 
  return(DLs12)
}

DLps12 = function(x,s){
  n = length(x); 
  x1 = x[-n]; x2 = x[-1]; s = s[-n]; 
  ind.positive.s = c(s<0); 
  x1 = x1[ind.positive.s]; x2 = x2[ind.positive.s]; s = s[ind.positive.s]
  DLps12 = 0; 
  p = 1; DLps12[1] = (sum((-s)^p*((1-x1)^(p+1) - (1-x2)^(p+1))/(p+1)))^(1/p); 
  p = 2; DLps12[2] = (sum((-s)^p*((1-x1)^(p+1) - (1-x2)^(p+1))/(p+1)))^(1/p); 
  DLps12[3] = max( 0-s*(1-x1), 0-s*(1-x2)); 
  return(DLps12)
}

invf = function(u,a,m,b,Ra,Rm,Rb){
  n = length(u)
  Ram  = Ra-Rm
  Ramb = Ra-2*Rm+Rb
  invf = array(,n)
  for (i in 1:n){
    cs = (Ram+sqrt(Ram^2-Ramb*(Ra-u[i])))/(Ramb)
    am = a-m
    amb = a-2*m+b
    invf[i] = ((cs*amb - am)^2-(m^2-a*b))/amb
  }
  return(invf)
}

RLS = function(u,delta=0){
  n = length(u)
  RLS = array(,n)
  a1 = 7/8*delta
  m1 = delta
  b1 = delta + 1/(8*18)
  a2 = delta + 7/(8*18)
  m2 = delta + 1/18
  b2 = delta + 1/9
  a3 = delta + 4/9
  m3 = delta + 1/2
  b3 = 7/8*delta + 9/16
  Ra1 = 7/8*delta
  Rm1 = delta
  Rb1 = delta + 1/18
  Ra2 = delta + 7/18
  Rm2 = delta + 8/18
  Rb2 = delta + (65)/(18*8)
  Ra3 = delta + (71)/(18*8)
  Rm3 = delta + 1/2
  Rb3 = 7/8*delta + 9/16
  for (i in 1:n){
    if( u[i]>=0 & u[i]<=a1){
      RLS[i] = u[i]
    }
    if( u[i]>a1 & u[i]<=b1){
      w = ftf(u=u[i],a=a1,m=m1,b=b1)
      RLS[i] = (1-w)^2*Ra1 + 2*(1-w)*w*Rm1 + w^2*Rb1
    }
    if( u[i]>b1 & u[i]<=a2){
      RLS[i] = delta + (u[i]-delta)*8
    }
    if( u[i]>a2 & u[i]<=b2){
      w = ftf(u=u[i],a=a2,m=m2,b=b2)
      RLS[i] = (1-w)^2*Ra2 + 2*(1-w)*w*Rm2 + w^2*Rb2
    }
    if( u[i]>b2 & u[i]<=a3){
      RLS[i] = (8/18+delta) + (u[i] - (1/18+delta))/8
    }
    if( u[i]>a3 & u[i]<=b3){
      w = ftf(u=u[i],a=a3,m=m3,b=b3)
      RLS[i] = (1-w)^2*Ra3 + 2*(1-w)*w*Rm3 + w^2*Rb3
    }
    if( u[i]>b3 & u[i]<=1){
      RLS[i] = u[i]
    }
    
  }
  return(RLS)
}

RIS = function(u,delta=0){
  n = length(u)
  RIS = array(,n)
  Ra1 = 7/8*delta
  Rm1 = delta
  Rb1 = delta + 1/(8*18)
  Ra2 = delta + 7/(8*18)
  Rm2 = delta + 1/18
  Rb2 = delta + 1/9
  Ra3 = delta + 4/9
  Rm3 = delta + 1/2
  Rb3 = 7/8*delta + 9/16
  a1 = 7/8*delta
  m1 = delta
  b1 = delta + 1/18
  a2 = delta + 7/18
  m2 = delta + 8/18
  b2 = delta + (65)/(18*8)
  a3 = delta + (71)/(18*8)
  m3 = delta + 1/2
  b3 = 7/8*delta + 9/16
  for (i in 1:n){
    if( u[i]>=0 & u[i]<=a1){
      RIS[i] = u[i]
    }
    if( u[i]>a1 & u[i]<=b1){
      RIS[i] = invf(u=u[i],a=Ra1,m=Rm1,b=Rb1,Ra=a1,Rm=m1,Rb=b1)
    }
    if( u[i]>b1 & u[i]<=a2){
      RIS[i] = delta + (u[i]-delta)/8
    }
    if( u[i]>a2 & u[i]<=b2){
      RIS[i] = invf(u=u[i],a=Ra2,m=Rm2,b=Rb2,Ra=a2,Rm=m2,Rb=b2)
    }
    if( u[i]>b2 & u[i]<=a3){
      RIS[i] = (1/18+delta) + (u[i] - (8/18 + delta))*8
    }
    if( u[i]>a3 & u[i]<=b3){
      RIS[i] = invf(u=u[i],a=Ra3,m=Rm3,b=Rb3,Ra=a3,Rm=m3,Rb=b3)
    }
    if( u[i]>b3 & u[i]<=1){
      RIS[i] = u[i]
    }
  }
  return(RIS)
}

#rm(list=ls(all=TRUE))

gen.sample.sim.table = function(nv=c(200,200), case=1){
  k = length(nv); 
  ### k=2 under H0 from AoS paper
  if(k == 2 & case == 1){
    X.data = list(
      X1 = runif(nv[1]),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 2){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 3){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.8),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 4){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 5){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      X2 = runif(nv[2])
    )
  }
  # k=2 under H1 from AoS paper
  if(k == 2 & case == 6){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta=0.25),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 7){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      X2 = runif(nv[2])
    )
  }
  if(k == 2 & case == 8){
    X.data = list(
      X1 = rnorm(nv[1],0,2),
      X2 = rnorm(nv[2])
    )
  }
  if(k == 2 & case == 9){
    X.data = list(
      X1 = rnorm(nv[1],0,1/2),
      X2 = rnorm(nv[2])
    )
  }
  
  
  ### k=3 under H0
  if(k == 3 & case == 1){
    X.data = list(
      X1 = runif(nv[1]),
      X2 = runif(nv[2]), 
      X3 = runif(nv[3])
    )
  }
  if(k == 3 & case == 2){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3])
    )
  }
  if(k == 3 & case == 3){
    X.data = list(
      # X1 = InvNSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      # X2 = InvNSR(runif(nv[2]),delta1=0.4,delta2=0.4),
      # X3 = runif(nv[3])
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = NSR(runif(nv[3]),delta1=0.4,delta2=0.4)
    )
  }
  
  ### k=3 under H1
  if(k == 3 & case == 4){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = sort(runif(nv[2])),
      X3 = NSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4)
    )
  }
  if(k == 3 & case == 5){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3])
    )
  }
  if(k == 3 & case == 6){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = InvNSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4)
    )
  }
  ### k=4 under H0  
  if(k == 4 & case == 1){
    X.data = list(
      X1 = runif(nv[1]),
      X2 = runif(nv[2]), 
      X3 = runif(nv[3]),
      X4 = runif(nv[4])
    )
  }
  if(k == 4 & case == 2){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4])
    )
  }
  if(k == 4 & case == 3){
    X.data = list(
      # X1 = InvNSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      # X2 = InvNSR(runif(nv[2]),delta1=0.4,delta2=0.4),
      # X3 = runif(nv[3]),
      # X4 = runif(nv[4])
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = NSR(runif(nv[3]),delta1=0.4,delta2=0.4),
      X4 = NSR(runif(nv[4]),delta1=0.4,delta2=0.4)
    )
  }
  
  ### k=4 under H1
  if(k == 4 & case == 4){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = sort(runif(nv[2])),
      X3 = NSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4),
      X4 = NSR(sort(runif(nv[4])),delta1=0.4,delta2=0.4)
    )
  }
  if(k == 4 & case == 5){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4])
    )
  }
  if(k == 4 & case == 6){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = InvNSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4),
      X4 = InvNSR(sort(runif(nv[4])),delta1=0.4,delta2=0.4)
    )
  }
  
  
  ### k=5 under H0  
  if(k == 5 & case == 1){
    X.data = list(
      X1 = runif(nv[1]),
      X2 = runif(nv[2]), 
      X3 = runif(nv[3]),
      X4 = runif(nv[4]),
      X5 = runif(nv[5])
    )
  }
  if(k == 5 & case == 2){
    X.data = list(
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4]),
      X5 = runif(nv[5])
    )
  }
  if(k == 5 & case == 3){
    X.data = list(
      # X1 = InvNSR(runif(nv[1]),delta1=0.8,delta2=0.4),
      # X2 = InvNSR(runif(nv[2]),delta1=0.4,delta2=0.4),
      # X3 = runif(nv[3]),
      # X4 = runif(nv[4]),
      # X5 = runif(nv[5])
      X1 = InvNSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = NSR(runif(nv[3]),delta1=0.4,delta2=0.4),
      X4 = NSR(runif(nv[4]),delta1=0.4,delta2=0.4),
      X5 = NSR(runif(nv[5]),delta1=0.4,delta2=0.4) 
    )
  }
  
  ### k=5 under H1
  if(k == 5 & case == 4){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = sort(runif(nv[2])),
      X3 = NSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4),
      X4 = NSR(sort(runif(nv[4])),delta1=0.4,delta2=0.4),
      X5 = NSR(sort(runif(nv[5])),delta1=0.4,delta2=0.4)
    )
  }
  if(k == 5 & case == 5){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4]),
      X5 = runif(nv[5])
    )
  }
  if(k == 5 & case == 6){
    X.data = list(
      X1 = NSR(runif(nv[1]),delta1=0.4,delta2=0.4),
      X2 = runif(nv[2]),
      X3 = InvNSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4),
      X4 = InvNSR(sort(runif(nv[4])),delta1=0.4,delta2=0.4),
      X5 = InvNSR(sort(runif(nv[5])),delta1=0.4,delta2=0.4)
    )
  }
  return(X.data)
}


gen.sample.power.curve = function(nv=c(200,200),delta=0,case=1){
  #delta is between 0 and 0.5
  k = length(nv); 
  if(k==2 & case==1){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2])
    )
  }
  
  if(k==3 & case==1){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = runif(nv[3])
    )
  }
  
  if(k==3 & case==2){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = RLS(runif(nv[3]),delta)
    )
  }
  
  if(k==4 & case==1){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4])
    )
  }
  
  if(k==4 & case==2){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = RLS(runif(nv[3]),delta),
      X4 = RLS(runif(nv[4]),delta)
    )
  }
  
  if(k==5 & case==1){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = runif(nv[3]),
      X4 = runif(nv[4]),
      X4 = runif(nv[5])
    )
  }
  
  if(k==5 & case==2){
    X.data = list(
      X1 = RIS(runif(nv[1]),delta),
      X2 = runif(nv[2]),
      X3 = RLS(runif(nv[3]),delta),
      X4 = RLS(runif(nv[4]),delta),
      X5 = RLS(runif(nv[5]),delta)
    )
  }
  
}

MGOFUSO = function(X.data){
  k = length(X.data); 
  nv = array(NA,k); M1s = NA; M2s = NA; Mss = NA; 
  for(j in 1:k){
    nv[j] = length(X.data[[j]])
  }
  us = list(); fun.Mrs = list(); cs = 0; 
  for(j in 1:(k-1)){
    us[[j]] = seq(0,1,by=1/nv[j+1]); 
    cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1])); 
  }
  #emp.Rs = list(); Maj.Rs = list(); 
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
    Mrj = cummin((1-emp.Rj[2:n.uj])/(1-uj[1:n.Y])); 
    MRj = 1-Mrj[1:n.Y]*(1-uj[1:n.Y]); MRj[n.uj] = 1; 
    ubj = 1-Mrj-emp.Rj[2:n.uj]+Mrj*(1:n.Y)/n.Y; 
    lbj = 1-Mrj-emp.Rj[2:n.uj]+Mrj*(0:(n.Y-1))/n.Y; 
    M1j = cj*(sum((ubj^(1+1)-lbj^(1+1))[Mrj>0]/(1+1)/Mrj[Mrj>0]))^(1/1);
    M2j = cj*(sum((ubj^(2+1)-lbj^(2+1))[Mrj>0]/(2+1)/Mrj[Mrj>0]))^(1/2);
    Msj = cj*max(ubj); 
    ### saving statsitics and slope function r
    M1s[j] = M1j; 
    M2s[j] = M2j; 
    Mss[j] = Msj; 
    fun.Mrs[[j]] = approxfun(uj,c(Mrj,Mrj[n.Y]),f=1); 
  }
  
  M1s; M2s; Mss; 
  
  Sk1 = sum(M1s); Sk1;  
  Sk2 = sum(M2s); Sk2; 
  Sks = sum(Mss); Sks; 
  
  Wk1 = max(M1s); Wk1; 
  Wk2 = max(M2s); Wk2; 
  Wks = max(Mss); Wks; 
  
  BB = 1000; 
  boot.M1s = array(NA,c(BB,(k-1))); 
  boot.M2s = array(NA,c(BB,(k-1))); 
  boot.Mss = array(NA,c(BB,(k-1))); 
  
  L = 1000; 
  s = seq(0,1,1/L); 
  ss = s[1:L];
  
  s.hat.rs = list(); s.MRs = list(); 
  
  for(j in 1:(k-1)){
    s.hat.rs[[j]] = fun.Mrs[[j]](s); 
    s.MRs[[j]] = 1-s.hat.rs[[j]]*(1-s)
  }
  boot.M1s = array(NA,c(BB,(k-1))); boot.M2s = array(NA,c(BB,(k-1))); boot.Mss = array(NA,c(BB,(k-1))); 
  for(bb in 1:BB){
    boot.Bs = list(); boot.emp.Bs = list(); 
    for(j in 1:k){
      ### simulating data from uniform distribution
      boot.Bs[[j]] = runif(nv[j]); 
      boot.emp.Bs[[j]] = ecdf(boot.Bs[[j]]); 
    }
    for(j in 1:(k-1)){
      ### loading data; 
      cj = cs[j]; 
      boot.emp.Bj1 = boot.emp.Bs[[j]]; 
      boot.emp.Bj2 = boot.emp.Bs[[j+1]]; 
      s.MRj1 = s.MRs[[j]]; 
      s.hat.rj = s.hat.rs[[j]];
      ### calculaing bounds from uniform distribution
      boot.Lj = cj*(boot.emp.Bj1(s.MRj1) - s.MRj1 - s.hat.rj*(boot.emp.Bj2(s) - s))
      R.boot.Lj = boot.Lj + s; 
      r.boot.Lj = (1-R.boot.Lj[1:L])/(1-ss); 
      hat.r.boot.Lj = cummin(r.boot.Lj); 
      MR.boot.Lj = c(1-(1-ss)*hat.r.boot.Lj, 1); 
      DR.boot.Lj = MR.boot.Lj - R.boot.Lj; 
      boot.M1j = mean(DR.boot.Lj); 
      boot.M2j = (mean((DR.boot.Lj)^2))^(1/2); 
      boot.Msj = max(DR.boot.Lj); 
      ### save simulated statistics
      boot.M1s[bb, j] = boot.M1j; 
      boot.M2s[bb, j] = boot.M2j; 
      boot.Mss[bb, j] = boot.Msj; 
    }
  }
  boot.cv.Sk1 = quantile(apply(boot.M1s,1,sum),0.95);
  boot.cv.Sk2 = quantile(apply(boot.M2s,1,sum),0.95);
  boot.cv.Sks = quantile(apply(boot.Mss,1,sum),0.95);
  boot.cv.Wk1 = quantile(apply(boot.M1s,1,max),0.95);
  boot.cv.Wk2 = quantile(apply(boot.M2s,1,max),0.95);
  boot.cv.Wks = quantile(apply(boot.Mss,1,max),0.95);
  
  MGOFUSO = list(
    M1s = M1s, M2s = M2s, Mss = Mss, 
    Skps = c(Sk1, Sk2, Sks), 
    Wkps = c(Wk1, Wk2, Wks), 
    boot.cv.Skps = as.numeric(c(boot.cv.Sk1, boot.cv.Sk2, boot.cv.Sks)),
    boot.cv.Wkps = as.numeric(c(boot.cv.Wk1, boot.cv.Wk2, boot.cv.Wks)),
    decision.Skps = as.logical(c(Sk1>boot.cv.Sk1, Sk2>boot.cv.Sk2, Sks>boot.cv.Sks)), 
    decision.Wkps = as.logical(c(Wk1>boot.cv.Wk1, Wk2>boot.cv.Wk2, Wks>boot.cv.Wks))
  )
  return(MGOFUSO)
}

Bon.cvs = function(k){
  if(k==2){
    Bon.cv.p1 = 0.580; 
    Bon.cv.p2 = 0.676; 
    Bon.cv.ps = 1.350; 
  }
  if(k==3){
    Bon.cv.p1 = 0.6585602; # (n0=3000)
    Bon.cv.p2 = 0.7628116; # (n0=3000)
    Bon.cv.ps = 1.4588240; # (n0=3000)
  }
  if(k==4){
    Bon.cv.p1 = 0.6993266; # (n0=3000)
    Bon.cv.p2 = 0.8052981; # (n0=3000)
    Bon.cv.ps = 1.5228830; # (n0=3000)
  }
  if(k==5){
    Bon.cv.p1 = 0.7287816; # (n0=3000)
    Bon.cv.p2 = 0.8350943; # (n0=3000)
    Bon.cv.ps = 1.5731570; # (n0=3000)
  }
  if(k==6){
    Bon.cv.p1 = 0.7721178; # (n0=3000)
    Bon.cv.p2 = 0.8778637; # (n0=3000)
    Bon.cv.ps = 1.6285020; # (n0=3000)
  }
  return(Bon.cvs = list(Bon.cv.p1 = Bon.cv.p1,
                        Bon.cv.p2 = Bon.cv.p2,
                        Bon.cv.ps = Bon.cv.ps))
}


MixNormal = function(n=c(100,100),p=c(0.25,0.25),m=c(1,2)){
  k = length(n); 
  MixNormal = list(); 
  for(j in 1:k){
    Ind = 0; 
    Ind = rbinom(n[j], 1, p[j]); 
    MixNormal[[j]] = sort(Ind*rnorm(n[j],0,1) + (1-Ind)*rnorm(n[j],m[j],1))
  }
  return(MixNormal)  
}


MixGamma = function(n=c(100,100),p=c(0.25,0.25),m=c(1,2),s=1,r=1){
  k = length(n); 
  MixGamma = list(); 
  for(j in 1:k){
    Ind = 0; 
    Ind = rbinom(n[j], 1, p[j]); 
    MixGamma[[j]] = sort(Ind * rgamma(n[j], shape = s, rate = r) + 
                           (1-Ind)*( m[j] + rgamma(n[j], shape = s, rate = r)) 
    )
  }
  return(MixGamma)  
}

MUSODecision = function(X.data,BB=1000,unif.n=1000,alpha=0.05){
  
  BB=1000; unif.n=1000; alpha=0.05
  
  k = length(X.data); 
  if(k==2){
    Bon.cv.p1 = 0.580; 
    Bon.cv.p2 = 0.676; 
    Bon.cv.ps = 1.350; 
  }
  if(k==3){
    Bon.cv.p1 = 0.6585602; # (n0=3000)
    Bon.cv.p2 = 0.7628116; # (n0=3000)
    Bon.cv.ps = 1.4588240; # (n0=3000)
  }
  if(k==4){
    Bon.cv.p1 = 0.6993266; # (n0=3000)
    Bon.cv.p2 = 0.8052981; # (n0=3000)
    Bon.cv.ps = 1.5228830; # (n0=3000)
  }
  if(k==5){
    Bon.cv.p1 = 0.7287816; # (n0=3000)
    Bon.cv.p2 = 0.8350943; # (n0=3000)
    Bon.cv.ps = 1.5731570; # (n0=3000)
  }
  if(k==6){
    Bon.cv.p1 = 0.7721178; # (n0=3000)
    Bon.cv.p2 = 0.8778637; # (n0=3000)
    Bon.cv.ps = 1.6285020; # (n0=3000)
  }
  
  nv = array(NA,k); M1s = NA; M2s = NA; Mss = NA; 
  for(j in 1:k){
    nv[j] = length(X.data[[j]])
  }
  
  us = list(); fun.Mrs = list(); cs = 0; 
  for(j in 1:(k-1)){
    us[[j]] = seq(0,1,by=1/nv[j+1]); 
    cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1])); 
  }
  L = 2000; 
  switch.more.uniform = TRUE; 
  switch.boot.exact.Lp = TRUE; 
  
  emp.Rs = list(); Maj.Rs = list(); 
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
    ubj = 1-Mrj-emp.Rj[2:n.uj]+Mrj*(1:n.Y)/n.Y; 
    lbj = 1-Mrj-emp.Rj[2:n.uj]+Mrj*(0:(n.Y-1))/n.Y; 
    M1j = cj*(sum((ubj^(1+1)-lbj^(1+1))[Mrj>0]/(1+1)/Mrj[Mrj>0]))^(1/1);
    M2j = cj*(sum((ubj^(2+1)-lbj^(2+1))[Mrj>0]/(2+1)/Mrj[Mrj>0]))^(1/2);
    Msj = cj*max(ubj); 
    ### saving statsitics and slope function r
    M1s[j] = M1j; 
    M2s[j] = M2j; 
    Mss[j] = Msj; 
    fun.Mrs[[j]] = approxfun(uj,c(Mrj,Mrj[n.Y]),f=1); 
  }
  
  M1s; M2s; Mss; 
  
  Sk1 = sum(M1s); Sk1;  
  Sk2 = sum(M2s); Sk2; 
  Sks = sum(Mss); Sks; 
  
  Wk1 = max(M1s); Wk1; 
  Wk2 = max(M2s); Wk2; 
  Wks = max(Mss); Wks; 
  
  
  boot.M1s = array(NA,c(BB,(k-1))); 
  boot.M2s = array(NA,c(BB,(k-1))); 
  boot.Mss = array(NA,c(BB,(k-1))); 
  
  
  s = seq(0,1,1/L); 
  ss = s[1:L];
  
  s.hat.rs = list(); s.MRs = list(); 
  
  for(j in 1:(k-1)){
    s.hat.rs[[j]] = fun.Mrs[[j]](s); 
    s.MRs[[j]] = 1-s.hat.rs[[j]]*(1-s)
  }
  boot.M1s = array(NA,c(BB,(k-1))); boot.M2s = array(NA,c(BB,(k-1))); boot.Mss = array(NA,c(BB,(k-1))); 
  for(bb in 1:BB){
    boot.Bs = list(); boot.emp.Bs = list(); 
    for(j in 1:k){
      ### simulating data from uniform distribution
      switch = 1; 
      if(switch.more.uniform==0){boot.Bs[[j]] = runif(nv[j]);}
      if(switch.more.uniform==1){boot.Bs[[j]] = runif(unif.n);}
      boot.emp.Bs[[j]] = ecdf(boot.Bs[[j]]); 
    }
    for(j in 1:(k-1)){
      ### loading data; 
      cj = cs[j]; lambdaj = nv[j+1]/(nv[j]+nv[j+1]); 
      boot.emp.Bj1 = boot.emp.Bs[[j]]; 
      boot.emp.Bj2 = boot.emp.Bs[[j+1]]; 
      s.MRj1 = s.MRs[[j]]; 
      s.hat.rj = s.hat.rs[[j]];
      ### calculaing bounds from uniform distribution
      if(switch.more.uniform==0){boot.Lj = cj*(boot.emp.Bj1(s.MRj1) - s.MRj1 - s.hat.rj*(boot.emp.Bj2(s) - s));}
      if(switch.more.uniform==1){boot.Lj = sqrt(lambdaj)*sqrt(unif.n)*(boot.emp.Bj1(s.MRj1) - s.MRj1) - s.hat.rj*sqrt(1-lambdaj)*sqrt(1000)*(boot.emp.Bj2(s) - s);}
      R.boot.Lj = boot.Lj + s; 
      r.boot.Lj = (1-R.boot.Lj[1:L])/(1-ss); 
      hat.r.boot.Lj = cummin(r.boot.Lj); 
      MR.boot.Lj = c(1-(1-ss)*hat.r.boot.Lj, 1); 
      DR.boot.Lj = MR.boot.Lj - R.boot.Lj; 
      DR.boot.rj = (0-DR.boot.Lj)/(1-s); DR.boot.rj[L+1] = 0; 
      
      if(switch.boot.exact.Lp==0){
        temp = DLps12(s,DR.boot.rj); 
        boot.M1j = temp[1]; 
        boot.M2j = temp[2]; 
        boot.Msj = temp[3]; 
      }
      if(switch.boot.exact.Lp==1){
        boot.M1j = mean(DR.boot.Lj); 
        boot.M2j = (mean((DR.boot.Lj)^2))^(1/2); 
        boot.Msj = max(DR.boot.Lj); 
      }
      
      ### save simulated statistics
      boot.M1s[bb, j] = boot.M1j; 
      boot.M2s[bb, j] = boot.M2j; 
      boot.Mss[bb, j] = boot.Msj; 
    }
  }
  boot.cv.Sk1 = quantile(apply(boot.M1s,1,sum),1-alpha);
  boot.cv.Sk2 = quantile(apply(boot.M2s,1,sum),1-alpha);
  boot.cv.Sks = quantile(apply(boot.Mss,1,sum),1-alpha);
  boot.cv.Wk1 = quantile(apply(boot.M1s,1,max),1-alpha);
  boot.cv.Wk2 = quantile(apply(boot.M2s,1,max),1-alpha);
  boot.cv.Wks = quantile(apply(boot.Mss,1,max),1-alpha);
  
  DSks = as.logical(c(Sk1>boot.cv.Sk1, Sk2>boot.cv.Sk2, Sks>boot.cv.Sks));
  DWks = as.logical(c(Wk1>boot.cv.Wk1, Wk2>boot.cv.Wk2, Wks>boot.cv.Wks));
  DBon = c(Wk1>Bon.cv.p1, Wk2>Bon.cv.p2, Wks>Bon.cv.ps);
  return(MUSODecision = list(Sks=DSks,Wks=DWks,Bon=DBon))
}

DDDST = function(X.data){
  k = length(X.data)
  nv = array(NA,k); 
  for(j in 1:k){
    nv[j] = length(X.data[[j]])
  }
  M1s = NA; M2s = NA; Mss = NA; 
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
    M1s[j] = M1j; 
    M2s[j] = M2j; 
    Mss[j] = Msj; 
    fun.Mrs[[j]] = approxfun(uj,c(Mrj,Mrj[n.Y]),f=1); 
  }
  return(list(D1s = M1s, D2s = M2s, Dss = Mss, MaxDps = c(max(M1s), max(M2s), max(Mss))))
}


