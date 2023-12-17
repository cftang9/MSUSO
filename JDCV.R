JD_CV = function(nv){
  B0 = 10000; k = length(nv); 
  UPs0 = array(,c(B0,3)); 
  set.seed(102220220)
  for(b0 in 1:B0){
    X0 = list()
    for(i in 1:k){
      X0[[i]] = runif(nv[i]); 
    }
    R0 = list(); D0 = array(,c((k-1),3)); 
    for(j in 1:(k-1)){
      R0[[j]] = Rmn(X0[[j]],X0[[j+1]]); 
      D0[j,] = DP(R0[[j]],nv[j],nv[j+1]); 
    }
    UPs0[b0,] = apply(D0,2,max);
  }
  alpha = 0.05;
  qU = apply(UPs0,2,function(x){quantile(x=x,1-alpha)});
  return(list(qU=qU))
}