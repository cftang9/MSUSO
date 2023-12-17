p_DR = function(x,m,k, Asy=F){
  #x=10; m=7; k=3
  if(Asy==T){
    p_DR = pnorm(x,mean=((m-1)*sum(1/seq(2:m))),sd=sqrt((m-1)*sum(3/seq(2:m)-(1/seq(2:m))^2)),lower.tail=FALSE)
  }
  if(Asy==F){
    R1 = array(0, 1+k*(m-1)); R1[(m-1) + 1] = 1; #R1; 
    u = seq(0,(m-1)); #u; 
    wu = dbinom(u,size=(m-1),prob=(1/k)); #wu; sum(wu); 
    for(i in 1:(k-1)){
      R2 = array(0, 1+k*(m-1));
      for(l in ((m-1)+1):(k*(m-1)+1) ){
        R2[l] = sum( wu * R1[ l-u ])
      }
      R1 = R2; 
    }
    p_DR = sum(R2[((m-1)+1):(1+k*(m-1))]*pchisq(x,df=(0:((k-1)*(m-1))),lower.tail=FALSE))
  }
  return(p_DR)
}


DR_Eqd_USO = function(X_data, m){
  #X_data <- list(rnorm(100),rnorm(100),rnorm(100)); m = 7; 
  
  k = length(X_data);
  nv = array(NA,k);
  Emp_Fj = list();
  X_pool = c();
  for(j in 1:k){
    nv[j] = length(X_data[[j]]);
    Emp_Fj[[j]] = ecdf(X_data[[j]]);
    X_pool = c(X_pool,X_data[[j]]);
  }
  n_pool = sum(nv);
  Emp_X_pool = ecdf(X_pool);
  LogEmpLikRatio = 0;
  
  M = quantile(X_pool,(0:m)/m)
  M[1] = M[1]-0.1;
  # Range = range(X_pool)
  # M = seq(Range[1],Range[2],by=(Range[2]-Range[1])/m)
  # 
  for(i in 1:(m-1)){
    i1 = i; i2 = i+1;
    x = as.numeric(M[i1]); y = as.numeric(M[i2]); 
    Emp_theta_j = array(,k);
    njx = array(,k); njy = array(,k);
    for(j in 1:k){
      Emp_theta_j[j] = (1-Emp_Fj[[j]](y))/(1-Emp_Fj[[j]](x));
      njx[j] = nv[j]*(1-Emp_Fj[[j]](x));
      njy[j] = nv[j]*(1-Emp_Fj[[j]](y)); 
    }
    Rn1 = 0; Rn2 = 0;
    Rd1 = 0; Rd2 = 0;
    if( all(njx>0) ){
      Iso_theta_j = Iso::pava(y = Emp_theta_j, w = njx, decreasing=FALSE);
      Eqd_theta_j = rep((1-Emp_X_pool(y))/(1-Emp_X_pool(x)),k);
      for(j in 1:k){
        if(njy[j] > 0 & njx[j] > njy[j]){
          Rn1 = Rn1 + njy[j]*log(Eqd_theta_j[j])
          Rd1 = Rd1 + njy[j]*log(Iso_theta_j[j])
          Rn2 = Rn2 + (njx[j]-njy[j])*log(1-Eqd_theta_j[j])
          Rd2 = Rd2 + (njx[j]-njy[j])*log(1-Iso_theta_j[j])
        }
      }
    }
    LogEmpLikRatio = -2*(Rn1+Rn2-Rd1-Rd2) + LogEmpLikRatio #(3.1)
  }

  if(1==0){
    # m = 7; 
    if(k==3){cDR = 12.11926 }
    if(k==4){cDR = 14.16423 }
    if(k==5){cDR = 15.88740 }
    #return(list(LogEmpLikRatio = LogEmpLikRatio, dDR = as.numeric(LogEmpLikRatio>cDR )))
  }
    return(list(LogEmpLikRatio = LogEmpLikRatio))
}

if(1==0){
  B = 10000; test = array(,B)
  for(b in 1:B){
    X_data <- list(runif(60),runif(60),runif(60),runif(60),runif(60));
    test[b] = DR_Eqd_USO(X_data)$LogEmpLikRatio;
  }
  quantile(test,0.95)
}


DR_GOF_USO = function(X_data,m=4){
  #X_data <- list(rnorm(100),rnorm(100),rnorm(100)); m = 7; 
  
  k = length(X_data);
  nv = array(NA,k);
  Emp_Fj = list();
  X_pool = c();
  for(j in 1:k){
    nv[j] = length(X_data[[j]]);
    Emp_Fj[[j]] = ecdf(X_data[[j]]);
    X_pool = c(X_pool,X_data[[j]]);
  }
  n_pool = sum(nv);
  Emp_X_pool = ecdf(X_pool);
  LogEmpLikRatio = 0;
  
  M = quantile(X_pool,(0:m)/m)
  M[1] = M[1]-0.1; 
  
  
  for(i in 1:(m-1)){
    i1 = i; i2 = i+1;
    x = as.numeric(M[i1]); y = as.numeric(M[i2]); 
    Emp_theta_j = array(,k);
    njx = array(,k); njy = array(,k);
    for(j in 1:k){
      Emp_theta_j[j] = (1-Emp_Fj[[j]](y))/(1-Emp_Fj[[j]](x));
      njx[j] = nv[j]*(1-Emp_Fj[[j]](x));
      njy[j] = nv[j]*(1-Emp_Fj[[j]](y)); 
    }
    Rn1 = 0; Rn2 = 0;
    Rd1 = 0; Rd2 = 0;
    if( all(njx>0) ){
      Iso_theta_j = Iso::pava(y = Emp_theta_j, w = njx);
      for(j in 1:k){
        if(njy[j] > 0){
          Rn1 = Rn1 + njy[j]*(Iso_theta_j[j])
          Rd1 = Rd1 + njy[j]*(Emp_theta_j[j])
        }
        if(njx[j] > njy[j]){
          Rn2 = Rn2 + (njx[j]-njy[j])*(1-Iso_theta_j[j])
          Rd2 = Rd2 + (njx[j]-njy[j])*(1-Emp_theta_j[j])
        }
      }
    }
    LogEmpLikRatio = -2*(Rn1+Rn2-Rd1-Rd2) + LogEmpLikRatio #(3.1)
  }
  
  if(1==1){
    
    if(k==3 & m == 4 ){cDR =  9.177}
    if(k==3 & m == 7 ){cDR = 14.750}
    if(k==3 & m == 10){cDR = 19.809}
    
    if(k==4 & m == 4 ){cDR = 12.752}
    if(k==4 & m == 7 ){cDR = 21.080}
    if(k==4 & m == 10){cDR = 29.393}
    
    if(k==5 & m == 4 ){cDR = 16.268}
    if(k==5 & m == 7 ){cDR = 27.827}
    if(k==5 & m == 10){cDR = 39.108}
    
    if(k==10 & m == 4 ){cDR = 33.545}
    if(k==10 & m == 7 ){cDR = 61.179}
    if(k==10 & m == 10){cDR = 87.632}
    
  }
  return(list(LogEmpLikRatio = LogEmpLikRatio, dDR = as.numeric(LogEmpLikRatio>cDR )))
}



#X_data <- list(rnorm(100),rnorm(50),rnorm(50)); DR_GOF_USO(X_data);

