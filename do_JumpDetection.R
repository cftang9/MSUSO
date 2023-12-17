library(xtable)
do_JumpDetection = function(nv = c(60,60,60), 
                            jh=c(1,0), delta2=0.8, 
                            JDCV, 
                            B=1000, 
                            Demo=1, Demo_b=250){
  k = length(nv); 
  cs = array(,k-1); 
  qU = JDCV$qU;  
  Jump.p0 = array(,c(B,k-1,3)); 
  Jump.ps = array(,c(B,k-1,3)); 
   
  #set.seed(102220221);
  set.seed(102220221 + sum(jh*10^((length(jh)-1):0)));   
  start.time = Sys.time()
  for(b in 1:B){
    X = rR_samples(nv=nv,jh=jh,delta2=delta2) 
    R = list(); 
    D = array(,c((k-1),3)); 
    M = array(,c((k-1),3)); 
    for(j in 1:(k-1)){
      cs[j] = sqrt((nv[j]*nv[j+1])/(nv[j]+nv[j+1])); 
      R[[j]] = Rmn(X[[j]],X[[j+1]]); 
      D[j,] = DP(R[[j]],nv[j],nv[j+1]); 
      M[j,] = MP(R[[j]],nv[j],nv[j+1]); 
      Jump.p0[b,j,] = c(D[j,]>qU)
    }
    
    JP1 = array(,c(k,k-1)); QP1 = array(,k); 
    JP2 = array(,c(k,k-1)); QP2 = array(,k); 
    JPs = array(,c(k,k-1)); QPs = array(,k); 

    constp1 = log(log(cs))*2/3; 
    constp2 = log(log(cs))*3/4; 
    constps = log(log(cs))*1; 
    
    delta.p1 = c(0,sort(D[,1]));  
    delta.p2 = c(0,sort(D[,2])); 
    delta.ps = c(0,sort(D[,3])); 
    
    for(j in 1:k){
      JP1[j,] = c(D[,1]>delta.p1[j]); 
      JP2[j,] = c(D[,2]>delta.p2[j]);
      JPs[j,] = c(D[,3]>delta.ps[j]);
      QP1[j] = sum( (1-JP1[j,])*(D[,1]/cs) + (JP1[j,])*((M[,1]+log(cs)*constp1)/cs));
      QP2[j] = sum( (1-JP2[j,])*(D[,2]/cs) + (JP2[j,])*((M[,2]+log(cs)*constp2)/cs));
      QPs[j] = sum( (1-JPs[j,])*(D[,3]/cs) + (JPs[j,])*((M[,3]+log(cs)*constps)/cs));
    }
    
    delta.p1.star = delta.p1[which.min(QP1)]
    delta.p2.star = delta.p2[which.min(QP2)]
    delta.ps.star = delta.ps[which.min(QPs)]
    
    Jump.ps[b,,1] = c(D[,1]>delta.p1.star);
    Jump.ps[b,,2] = c(D[,2]>delta.p2.star);
    Jump.ps[b,,3] = c(D[,3]>delta.ps.star); 
  }
  
  
  exact.jump = c(jh>0);
  
  exact.rate.p0 = array(0,c(1,3)); 
  TP.p0 = array(0,c(1,3)); 
  FP.p0 = array(0,c(1,3)); 
  
  exact.rate.ps = array(0,c(1,3)); 
  TP.ps = array(0,c(1,3)); 
  FP.ps = array(0,c(1,3)); 
  
  for(b in 1:B){
    exact.rate.p0[1] = all(Jump.p0[b,,1] == exact.jump)/B + exact.rate.p0[1]; 
    exact.rate.p0[2] = all(Jump.p0[b,,2] == exact.jump)/B + exact.rate.p0[2]; 
    exact.rate.p0[3] = all(Jump.p0[b,,3] == exact.jump)/B + exact.rate.p0[3]; 
    TP.p0 = exact.jump%*%Jump.p0[b,,]/B + TP.p0; 
    FP.p0 = (1-exact.jump)%*%(Jump.p0[b,,])/B + FP.p0; 
    
    
    exact.rate.ps[1] = all(Jump.ps[b,,1] == exact.jump)/B + exact.rate.ps[1]; 
    exact.rate.ps[2] = all(Jump.ps[b,,2] == exact.jump)/B + exact.rate.ps[2]; 
    exact.rate.ps[3] = all(Jump.ps[b,,3] == exact.jump)/B + exact.rate.ps[3]; 
    TP.ps =  exact.jump%*%Jump.ps[b,,]/B + TP.ps; 
    FP.ps =  (1-exact.jump)%*%(Jump.ps[b,,])/B + FP.ps; 
  }
  
  
  
  
  Jp10 = c(exact.rate.p0[1], TP.p0[1], FP.p0[1])
  Jp20 = c(exact.rate.p0[2], TP.p0[2], FP.p0[2])
  Jps0 = c(exact.rate.p0[3], TP.p0[3], FP.p0[3])
  
  Jp1s = c(exact.rate.ps[1], TP.ps[1], FP.ps[1])
  Jp2s = c(exact.rate.ps[2], TP.ps[2], FP.ps[2])
  Jpss = c(exact.rate.ps[3], TP.ps[3], FP.ps[3])
  
  Jp = data.frame( rbind(array(c(Jp10,Jp20,Jps0),c(1,9)), 
                   array(c(Jp1s,Jp2s,Jpss),c(1,9)))
                   )
  
  colnames(Jp) = c("C","TP","FP","C","TP","FP","C","TP","FP")
  rownames(Jp) = c("Jp0","Jps")
  #print(Jp)
  
  return(list( Jp10 = Jp10, 
               Jp20 = Jp20, 
               Jps0 = Jps0, 
               Jp1s = Jp1s, 
               Jp2s = Jp2s, 
               Jpss = Jpss, 
               Jptable = xtable(Jp,digits=3)
               ))
}



