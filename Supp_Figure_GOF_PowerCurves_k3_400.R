rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

######## C1 n=400 #################3
nv = c(400,400,400); k = length(nv); 
alpha = 0.05; B = 1000; BB = 1000; 
Cases = 1:10; nC = length(Cases)
pSks_C1_400 = array(0,c(nC,3)); pWks_C1_400 = array(0,c(nC,3)); pBon_C1_400 = array(0,c(nC,3))
set.seed(0527202103)
for(l in Cases){
  start = Sys.time(); 
  for(b in 1:B){
    X1 = RIS(sort(runif(nv[1])), delta=1/2-l/20);
    X2 = sort(runif(nv[2])); 
    X3 = NSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4);
    
    X = list(X1, X2, X3); 
    temp = MGOFUSO(X)
    
    pSks_C1_400[l,] = pSks_C1_400[l,] + temp$decision.Skps/B
    pWks_C1_400[l,] = pWks_C1_400[l,] + temp$decision.Wkps/B
    pBon_C1_400[l,] = pBon_C1_400[l,] + temp$decision.Bon/B
    
    end = Sys.time(); 
    ave.interval = (end-start)/b
    if(b%%250==0){
      print(l)
      print(start); print("C1 n400")
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + (ave.interval*(B-b)) + ave.interval*B*(nC-l) ); 
    } 
  }
}

######## C3 n=400 #################3
nv = c(400,400,400); k = length(nv); 
alpha = 0.05; B = 1000; BB = 1000; 
Cases = 1:10; nC = length(Cases)
pSks_C3_400 = array(0,c(nC,3)); pWks_C3_400 = array(0,c(nC,3)); pBon_C3_400 = array(0,c(nC,3))
set.seed(0527202104)
for(l in Cases){
  start = Sys.time(); 
  for(b in 1:B){
    X1 = RIS(sort(runif(nv[1])), delta=1/2-l/20);
    X2 = sort(runif(nv[2])); 
    X3 = RLS(sort(runif(nv[3])), delta=1/2-l/20);
    
    X = list(X1, X2, X3); 
    temp = MGOFUSO(X)
    
    pSks_C3_400[l,] = pSks_C3_400[l,] + temp$decision.Skps/B
    pWks_C3_400[l,] = pWks_C3_400[l,] + temp$decision.Wkps/B
    pBon_C3_400[l,] = pBon_C3_400[l,] + temp$decision.Bon/B
    
    end = Sys.time(); 
    ave.interval = (end-start)/b
    if(b%%250==0){
      print(l)
      print(start); print("C3 n400")
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + (ave.interval*(B-b)) + ave.interval*B*(nC-l) ); 
    } 
  }
}

save(file="PCsk3n400.Rdata",
     pSks_C1_400, pWks_C1_400, pBon_C1_400,
     pSks_C3_400, pWks_C3_400, pBon_C3_400)

#pdf("Supp_Figure_GOF_PowerCurves_k3_400.pdf", width=19.5, heigh=13)

png("Supp_Figure_GOF_PowerCurves_k3_400.pdf", width=1872, heigh=1248)

par(mfrow=c(2,3))
par(mar=c(2,2,2,2)) 
plot((Cases-1)/10, pSks_C1_400[,1], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_400[,1],col="blue"); 
lines((Cases-1)/10, pBon_C1_400[,1],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C1_400[,2], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_400[,2],col="blue"); 
lines((Cases-1)/10, pBon_C1_400[,2],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C1_400[,3], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_400[,3],col="blue"); 
lines((Cases-1)/10, pBon_C1_400[,3],lty=2); 
abline(h=0.05,lty=3); 

plot((Cases-1)/10, pSks_C3_400[,1], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_400[,1],col="blue"); 
lines((Cases-1)/10, pBon_C3_400[,1],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C3_400[,2], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_400[,2],col="blue"); 
lines((Cases-1)/10, pBon_C3_400[,2],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C3_400[,3], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_400[,3],col="blue"); 
lines((Cases-1)/10, pBon_C3_400[,3],lty=2); 
abline(h=0.05,lty=3); 
par(mfcol=c(1,1))

dev.off()