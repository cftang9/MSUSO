do_EqualityTests <- function(nv = c(60,60,60), 
                             jh=c(1,-1), delta2=0.8, 
                             dists="R", 
                             eqCV, 
                             B=1000, 
                             MEB=1, 
                             DR=1, ms=c(4,7,10), 
                             Demo=1, Demo_b=250){
  set.seed(102220221); 
  k = length(nv); 
  
  UPs = array(,c(B,3)); TPs = array(,c(B,3)); 
  Cauchy_com = array(,c(B,3)); 
  Holm_pvalues = array(,c(B,3)); 
  Hochberg_pvalues = array(,c(B,3)); 
  Hommel_pvalues = array(,c(B,3)); 
  BH_pvalues = array(,c(B,3)); 
  BY_pvalues = array(,c(B,3)); 
  fdr_pvalues = array(,c(B,3)); 
  Bon_pvalues = array(,c(B,3)); 
  dMEB = array(,B); 
  DR_1 = array(,B); 
  DR_2 = array(,B); 
  DR_3 = array(,B); 
  
  start.time = Sys.time()
  for(b in 1:B){
    X = rR_samples(nv=nv,jh=jh,delta2=delta2) 
    #for(i in 1:k){
      # if(dists=="R"){
      #   X[[i]] = rMixUSO(n=nv[i],pl=param[i],delta1=0.8,delta2=0.8)
      # }
      # if(dists=="A"){
      #   X[[i]] = rMW(n=nv[i],a=1+exp(param));
      # }
      # if(dists=="NR"){
      #   X[[i]] = rMixNonUSO(n=nv[i],pl=param[i],delta1=0.8,delta2=0.8)
      # }
    #}
    R = list(); D = array(,c((k-1),3)); 
    for(j in 1:(k-1)){
      R[[j]] = Rmn(X[[j]],X[[j+1]]); 
      D[j,] = DP(R[[j]],nv[j],nv[j+1]); 
    }
    TPs[b,] = apply(D,2,sum);
    UPs[b,] = apply(D,2,max);

    pD = array(,c((k-1),3))
    for(p in 1:3){
      for(j in 1:(k-1)){
        pD[j,p] = mean(eqCV$DPs0[,p]>D[j,p]); 
      }
      Holm_pvalues[b,p] = min(p.adjust(pD[,p], method = "holm")); 
      Hochberg_pvalues[b,p] = min(p.adjust(pD[,p], method = "hochberg")); 
      Hommel_pvalues[b,p] = min(p.adjust(pD[,p], method = "hommel")); 
      BH_pvalues[b,p] = min(p.adjust(pD[,p], method = "BH")); 
      BY_pvalues[b,p] = min(p.adjust(pD[,p], method = "BY")); 
      fdr_pvalues[b,p] = min(p.adjust(pD[,p], method = "fdr")); 
      Bon_pvalues[b,p] = min(p.adjust(pD[,p], method = "bon")); 
      Cauchy_com[b,p] = mean( tan((0.5-pD[,p])*pi) ); 
    }
    
    if(MEB==1){
      temp = MEB_Eqd_USO(X_data=X)
      dMEB[b] = as.numeric(temp$Reject);
    }
    
    DR_1[b] = DR_Eqd_USO(X,m=ms[1])$LogEmpLikRatio; 
    DR_2[b] = DR_Eqd_USO(X,m=ms[2])$LogEmpLikRatio; 
    DR_3[b] = DR_Eqd_USO(X,m=ms[3])$LogEmpLikRatio; 
    
    
    ind = 1:b
    
    qT = eqCV$qT; qU = eqCV$qU; 

    if(b%%Demo_b==0 & Demo==1){
      print(start.time)
      print(b)
      int.time = (Sys.time() - start.time)/b
      print(Sys.time()+(B-b)*int.time)
      
      print(paste(mean(TPs[ind,1] > qT[1]), "&", 
                  mean(UPs[ind,1] > qU[1]), "&",
                  mean(Cauchy_com[ind,1] > qcauchy(0.95)), "&", 
                  mean(TPs[ind,2] > qT[2]), "&", 
                  mean(UPs[ind,2] > qU[2]), "&",
                  mean(Cauchy_com[ind,2] > qcauchy(0.95)), "&", 
                  mean(TPs[ind,3] > qT[3]), "&", 
                  mean(UPs[ind,3] > qU[3]), "&", 
                  mean(Cauchy_com[ind,3] > qcauchy(0.95)), "&", 
                  mean(DR_1[ind]>eqCV$qDR_1), "&",
                  mean(DR_2[ind]>eqCV$qDR_2), "&",
                  mean(DR_3[ind]>eqCV$qDR_3), "&",
                  mean(dMEB[ind])))
      print(dists)
      print(jh)
    }
  }
  
  Tab_pow = paste(mean(TPs[,1] > qT[1]), "&", 
                  mean(UPs[,1] > qU[1]), "&",
                  mean(TPs[,2] > qT[2]), "&", 
                  mean(UPs[,2] > qU[2]), "&",
                  mean(TPs[,3] > qT[3]), "&", 
                  mean(UPs[,3] > qU[3]), "&", 
                  mean(DR_1>eqCV$qDR_1), "&",
                  mean(DR_2>eqCV$qDR_2), "&",
                  mean(DR_3>eqCV$qDR_3), "&",
                  mean(dMEB)
  )
  
  Tab_pow_p = paste(mean(Cauchy_com[,1] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,1] < 0.05), "&", 
                     mean(Bon_pvalues[,1] < 0.05), "&", 
                     mean(Cauchy_com[,2] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,2] < 0.05), "&", 
                     mean(Bon_pvalues[,2] < 0.05), "&", 
                     mean(Cauchy_com[,3] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,3] < 0.05), "&", 
                     mean(Bon_pvalues[,3] < 0.05)
  )
  
  Tab_pow_1p = paste(mean(TPs[,1] > qT[1]), "&",
                     mean(UPs[,1] > qU[1]), "&",
                     mean(Cauchy_com[,1] > qcauchy(0.95)), "&",
                     mean(Holm_pvalues[,1] < 0.05), "&",
                     mean(Hochberg_pvalues[,1] < 0.05), "&",
                     mean(Hommel_pvalues[,1] < 0.05), "&",
                     mean(BH_pvalues[,1] < 0.05), "&",
                     mean(BY_pvalues[,1] < 0.05), "&",
                     mean(fdr_pvalues[,1] < 0.05), "&",
                     mean(Bon_pvalues[,1] < 0.05)
  )

  Tab_pow_2p = paste(mean(TPs[,2] > qT[2]), "&",
                     mean(UPs[,2] > qU[2]), "&",
                     mean(Cauchy_com[,2] > qcauchy(0.95)), "&",
                     mean(Holm_pvalues[,2] < 0.05), "&",
                     mean(Hochberg_pvalues[,2] < 0.05), "&",
                     mean(Hommel_pvalues[,2] < 0.05), "&",
                     mean(BH_pvalues[,2] < 0.05), "&",
                     mean(BY_pvalues[,2] < 0.05), "&",
                     mean(fdr_pvalues[,2] < 0.05), "&",
                     mean(Bon_pvalues[,2] < 0.05)
  )

  Tab_pow_sp = paste(mean(TPs[,3] > qT[3]), "&",
                     mean(UPs[,3] > qU[3]), "&",
                     mean(Cauchy_com[,3] > qcauchy(0.95)), "&",
                     mean(Holm_pvalues[,3] < 0.05), "&",
                     mean(Hochberg_pvalues[,3] < 0.05), "&",
                     mean(Hommel_pvalues[,3] < 0.05), "&",
                     mean(BH_pvalues[,3] < 0.05), "&",
                     mean(BY_pvalues[,3] < 0.05), "&",
                     mean(fdr_pvalues[,3] < 0.05), "&",
                     mean(Bon_pvalues[,3] < 0.05)
  )
  return(list(Tab_pow=Tab_pow, 
              Tab_pow_p=Tab_pow_p,
              Tab_pow_1p=Tab_pow_1p,
              Tab_pow_2p=Tab_pow_2p,
              Tab_pow_sp=Tab_pow_sp))
}
