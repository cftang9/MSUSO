source("ME-B.R")

k=3; n=60; 

start.time=Sys.time()
if(1==1){
  B = 3000; test = array(,B); 
  set.seed(102220220)
  for(b in 1:B){
    X_data = list()
    for(j in 1:k){
      X_data[[j]] <- runif(n);
    }
    
    test[b] = MEB_GOF_USO(X_data)$EL;
    if(b%%5==0){
      int.time = (Sys.time() - start.time)/b
      print(b)
      print(Sys.time()+int.time*(B-b))
    }
  }
  quantile(test,0.95)
}

if(1==0){
  save(test, 
       file=paste0("MEBGOFTS_k",k,"n",n,".Rdata"))
}