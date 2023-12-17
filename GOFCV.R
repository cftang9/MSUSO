GOFTest_CV = function(nv){
  B0 = 10000; k = length(nv); 
  SPs0 = array(,c(B0,3)); WPs0 = array(,c(B0,3)); MPs0 = array(,c(B0,3))
  set.seed(102220220)
  for(b0 in 1:B0){
    X0 = list()
    for(i in 1:k){
      X0[[i]] = runif(nv[i]); 
    }
    R0 = list(); M0 = array(,c((k-1),3)); 
    for(j in 1:(k-1)){
      R0[[j]] = Rmn(X0[[j]],X0[[j+1]]); 
      M0[j,] = MP(R0[[j]],nv[j],nv[j+1]); 
    }
    SPs0[b0,] = apply(M0,2,sum);
    WPs0[b0,] = apply(M0,2,max);
    MPs0[b0,] = M0[1,]; 
  }
  alpha = 0.05;
  qS = apply(SPs0,2,function(x){quantile(x=x,1-alpha)});
  qW = apply(WPs0,2,function(x){quantile(x=x,1-alpha)});
  qB = apply(MPs0,2,function(x){quantile(x=x,1-alpha/(k-1))});
  return(list(qS=qS, qW=qW, qB=qB, MPs0=MPs0))
}