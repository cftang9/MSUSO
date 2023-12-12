#source("EGJ_USO_Library.r")
rR_samples = function(nv=c(6000,6000,6000),jh=c(1,-1),delta1=0.4,delta2=0.8){
  k = length(nv); 
  rR_samples = list()
  for(j in 1:k){
    rR_samples[[j]] = runif(nv[j])
  }
  for(j in k:2){
    s = j
    while(s>1){
      if(jh[s-1]>=0){
        rR_samples[[j]] = NSR(rR_samples[[j]],delta1=delta1*abs(jh[s-1]),delta2=delta2)
      }
      if(jh[s-1]< 0){
        rR_samples[[j]] = InvNSR(rR_samples[[j]],delta1=delta1*abs(jh[s-1]),delta2=delta2)
      }
      s = s-1; 
    }
  }
  return(rR_samples)
}
