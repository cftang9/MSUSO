do_GOFTests <- function(nv = c(60,60,60), 
                        jh=c(1,-1), delta2=0.8, 
                        dists="R", 
                        GOFCV, 
                        B=1000, 
                        MEB=1, 
                        DR=1, 
                        Demo=1, Demo_b=250){
  set.seed(102220221); 
  k = length(nv); 
  
  DSks = array(,c(B,3)); DWks = array(,c(B,3)); DBon = array(,c(B,3));
  Cauchy_com = array(,c(B,3)); 
  Holm_pvalues = array(,c(B,3)); 
  Hochberg_pvalues = array(,c(B,3)); 
  Hommel_pvalues = array(,c(B,3)); 
  BH_pvalues = array(,c(B,3)); 
  BY_pvalues = array(,c(B,3)); 
  fdr_pvalues = array(,c(B,3)); 
  Bon_pvalues = array(,c(B,3)); 
  dMEB = array(,B); 
  
  dDR_1 = array(,B); 
  dDR_2 = array(,B); 
  dDR_3 = array(,B); 
  
  start.time = Sys.time()
  for(b in 1:B){
    X = rR_samples(nv=nv,jh=jh,delta2=delta2) 
    # for(i in 1:k){
    #   if(dists=="R"){
    #     X[[i]] = rMixUSO(n=nv[i],pl=param[i],delta1=0.8,delta2=0.8)
    #   }
    #   if(dists=="A"){
    #     X[[i]] = rMW(n=nv[i],a=1+exp(param));
    #   }
    #   if(dists=="NR"){
    #     X[[i]] = rMixNonUSO(n=nv[i],pl=param[i],delta1=0.8,delta2=0.8)
    #   }
    # }
    R = list(); M = array(,c((k-1),3)); 
    for(j in 1:(k-1)){
      R[[j]] = Rmn(X[[j]],X[[j+1]]); 
      M[j,] = MP(R[[j]],nv[j],nv[j+1]); 
    }
    
    pM = array(,c((k-1),3))
    for(p in 1:3){
      for(j in 1:(k-1)){
        pM[j,p] = mean(GOFCV$MPs0[,p]>M[j,p]); 
      }
      Holm_pvalues[b,p] = min(p.adjust(pM[,p], method = "holm")); 
      Hochberg_pvalues[b,p] = min(p.adjust(pM[,p], method = "hochberg")); 
      Hommel_pvalues[b,p] = min(p.adjust(pM[,p], method = "hommel")); 
      BH_pvalues[b,p] = min(p.adjust(pM[,p], method = "BH")); 
      BY_pvalues[b,p] = min(p.adjust(pM[,p], method = "BY")); 
      fdr_pvalues[b,p] = min(p.adjust(pM[,p], method = "fdr")); 
      Bon_pvalues[b,p] = min(p.adjust(pM[,p], method = "bon")); 
      Cauchy_com[b,p] = mean( tan((0.5-pM[,p])*pi) ); 
    }
    
    temp = MUSOGOF(X)
    DSks[b,] = as.numeric(temp$Sks);
    DWks[b,] = as.numeric(temp$Wks);
    
    if(k<=6)DBon[b,] = as.numeric(temp$Bon);
    
    
    #qS = GOFCV$qS; qW = GOFCV$qW; 
    
    if(MEB==1){
      temp = MEB_GOF_USO(X_data=X)
      dMEB[b] = as.numeric(temp$Reject);
    }
    
    dDR_1[b] = as.numeric(DR_GOF_USO(X,m=4)$dDR); 
    dDR_2[b] = as.numeric(DR_GOF_USO(X,m=7)$dDR); 
    dDR_3[b] = as.numeric(DR_GOF_USO(X,m=10)$dDR); 
    
    
    if(b%%Demo_b==0 & Demo==1){
      print(start.time)
      print(b)
      int.time = (Sys.time() - start.time)/b
      print(Sys.time()+(B-b)*int.time)
      ind = 1:b
      print(paste(mean(DSks[ind,1]), "&", 
                  mean(DWks[ind,1]), "&",
                  mean(Cauchy_com[ind,1] > qcauchy(0.95)), "&", 
                  mean(DSks[ind,1]), "&", 
                  mean(DWks[ind,2]), "&",
                  mean(Cauchy_com[ind,2] > qcauchy(0.95)), "&", 
                  mean(DSks[ind,3]), "&", 
                  mean(DWks[ind,3]), "&", 
                  mean(Cauchy_com[ind,3] > qcauchy(0.95)), "&", 
                  mean(dDR_1[ind]), "&", 
                  mean(dDR_2[ind]), "&", 
                  mean(dDR_3[ind]), "&", 
                  mean(dMEB[ind])))
      print(jh)
    }
  }
  
  
  
  Tab_pow = paste(mean(DSks[,1]), "&", 
                  mean(DWks[,1]), "&",
                  mean(DSks[,2]), "&", 
                  mean(DWks[,2]), "&",
                  mean(DSks[,3]), "&", 
                  mean(DWks[,3]), "&", 
                  mean(dDR_1), "&", 
                  mean(dDR_2), "&", 
                  mean(dDR_3), "&", 
                  mean(dMEB))
  
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
  
  Tab_pow_1p = paste(mean(DSks[,1]), "&", 
                     mean(DWks[,1]), "&",
                     mean(Cauchy_com[,1] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,1] < 0.05), "&", 
                     mean(Hochberg_pvalues[,1] < 0.05), "&", 
                     mean(Hommel_pvalues[,1] < 0.05), "&", 
                     mean(BH_pvalues[,1] < 0.05), "&", 
                     mean(BY_pvalues[,1] < 0.05), "&", 
                     mean(fdr_pvalues[,1] < 0.05), "&", 
                     mean(Bon_pvalues[,1] < 0.05)
  )
  
  Tab_pow_2p = paste(mean(DSks[,2]), "&", 
                     mean(DWks[,2]), "&",
                     mean(Cauchy_com[,2] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,2] < 0.05), "&", 
                     mean(Hochberg_pvalues[,2] < 0.05), "&", 
                     mean(Hommel_pvalues[,2] < 0.05), "&", 
                     mean(BH_pvalues[,2] < 0.05), "&", 
                     mean(BY_pvalues[,2] < 0.05), "&", 
                     mean(fdr_pvalues[,2] < 0.05), "&", 
                     mean(Bon_pvalues[,2] < 0.05)
  )
  
  Tab_pow_sp = paste(mean(DSks[,3]), "&", 
                     mean(DWks[,3]), "&", 
                     mean(Cauchy_com[,3] > qcauchy(0.95)), "&", 
                     mean(Holm_pvalues[,3] < 0.05), "&", 
                     mean(Hochberg_pvalues[,3] < 0.05), "&",
                     mean(Hommel_pvalues[,3] < 0.05), "&",
                     mean(BH_pvalues[,3] < 0.05), "&",
                     mean(BY_pvalues[,3] < 0.05), "&",
                     mean(fdr_pvalues[,3] < 0.05), "&", 
                     mean(Bon_pvalues[,3] < 0.05)
  )
  return(list(Tab_pow=Tab_pow, Tab_pow_p=Tab_pow_p, Tab_pow_1p=Tab_pow_1p, Tab_pow_2p=Tab_pow_2p, Tab_pow_sp=Tab_pow_sp))
}
