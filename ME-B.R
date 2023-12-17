# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# hello <- function() {
#   print("Hello, world!")
# }

MEB_Eqd_USO = function(X_data){
  k = length(X_data);
  nv = array(NA,k);
  Emp_Fj = list();
  X_pool = c();
  for(j in 1:k){
    nv[j] = length(X_data[[j]]);
    Emp_Fj[[j]] = ecdf(X_data[[j]]);
    X_pool = c(X_pool,X_data[[j]]);
  }
  n = sum(nv);
  X_pool = sort(X_pool);
  Emp_X_pool = ecdf(X_pool);
  LogEmpLikRatio = 0;

  for(i in 1:(n^2)){
    i1 = (i-1)%/%n + 1; i2 = (i-1)%%n + 1;
    if(i1 < i2){
      x = X_pool[i1]; y = X_pool[i2];
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
        Eqd_theta_j = rep((1-Emp_X_pool(y))/(1-Emp_X_pool(x)),k);
        for(j in 1:k){
          if(njy[j] > 0){
            Rn1 = Rn1 + njy[j]*(Eqd_theta_j[j])
            Rd1 = Rd1 + njy[j]*(Iso_theta_j[j])
          }
          if(njx[j] > njy[j]){
            Rn2 = Rn2 + (njx[j]-njy[j])*(1-Eqd_theta_j[j])
            Rd2 = Rd2 + (njx[j]-njy[j])*(1-Iso_theta_j[j])
          }
        }
      }
      LogEmpLikRatio = -2/n^2*(Rn1+Rn2-Rd1-Rd2) + LogEmpLikRatio
    }
  }
  #CV = c(0.777, 1.112, 1.373, 1.490)
  # if(k <= 5){
  #   return(list(EL = LogEmpLikRatio, CV = CV[(k-1)], Reject = c(LogEmpLikRatio>CV[k-1]) ))
  # }
  # if(k > 5) {return(LogEmpLikRatio);}
  
  if(k==3 & nv[1]== 60){CV=0.7945713}
  if(k==3 & nv[1]==100){CV=0.7795464}
  if(k==3 & nv[1]==200){CV=0.7743946}
  if(k==4 & nv[1]== 60){CV=0.9045977}
  if(k==4 & nv[1]==100){CV=0.9104321}
  if(k==4 & nv[1]==200){CV=0.9184299}
  if(k==5 & nv[1]== 60){CV=1.039331}
  if(k==5 & nv[1]==100){CV=1.042572}
  if(k==5 & nv[1]==200){CV=1.041397}
  
  return(list(EL = LogEmpLikRatio, CV = CV, Reject = as.numeric(LogEmpLikRatio>CV) ))
}

if(1==0){
  B = 10000; test = array(,B)
  for(b in 1:B){
    X_data <- list(runif(60),runif(60),runif(60));
    test[b] = MEB_Eqd_USO(X_data)$LogEmpLikRatio;
  }
  quantile(test,0.95)
}

#X_data <- list(rnorm(100),rnorm(50),rnorm(50)); MEB_Eqd_USO(X_data);

MEB_GOF_USO = function(X_data){
  k = length(X_data);
  nv = array(NA,k);
  Emp_Fj = list();
  X_pool = c();
  for(j in 1:k){
    nv[j] = length(X_data[[j]]);
    Emp_Fj[[j]] = ecdf(X_data[[j]]);
    X_pool = c(X_pool,X_data[[j]]);
  }
  n = sum(nv);
  X_pool = sort(X_pool);
  Emp_X_pool = ecdf(X_pool);
  LogEmpLikRatio = 0;

  for(i in 1:(n^2)){
    i1 = (i-1)%/%n + 1; i2 = (i-1)%%n + 1;
    if(i1 < i2){
      x = X_pool[i1]; y = X_pool[i2];
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
      LogEmpLikRatio = -2/n^2*(Rn1+Rn2-Rd1-Rd2) + LogEmpLikRatio
    }
  }
  
  if(k==3 & nv[1]== 60){CV=0.9186258}
  if(k==3 & nv[1]==100){CV=0.9482728}
  if(k==3 & nv[1]==200){CV=0.9159945}
  if(k==4 & nv[1]== 60){CV=1.288373}
  if(k==4 & nv[1]==100){CV=1.301387}
  if(k==4 & nv[1]==200){CV=1.319737}
  if(k==5 & nv[1]== 60){CV=1.662549}
  if(k==5 & nv[1]==100){CV=1.688474}
  if(k==5 & nv[1]==200){CV=1.655239}

  return(list(EL = LogEmpLikRatio, CV = CV, Reject = as.numeric(LogEmpLikRatio>CV) ))

}

#X_data <- list(rnorm(100),rnorm(50),rnorm(50)); MEB_GOF_USO(X_data);

rMW = function(n,a=1){
  U = runif(n); 
  rMW = sqrt(-2*log(U)) * (U>exp(-1/2)) + 
        sqrt(-2/a*(log(U)+(1-a)/2)) * (U<=exp(-1/2))
  return(rMW)
}

#plot(Rmn(rMW(1000,a=1.5),rMW(1000,a=1)))



