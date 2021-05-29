rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

######## C1 n=200 #################3
nv = c(200,200,200); k = length(nv); 
alpha = 0.05; B = 1000; BB = 1000; 
Cases = 1:10; nC = length(Cases)
pSks_C1_200 = array(0,c(nC,3)); pWks_C1_200 = array(0,c(nC,3)); pBon_C1_200 = array(0,c(nC,3))
set.seed(0527202101)
for(l in Cases){
  start = Sys.time(); 
  for(b in 1:B){
    X1 = RIS(sort(runif(nv[1])), delta=1/2-l/20);
    X2 = sort(runif(nv[2])); 
    X3 = NSR(sort(runif(nv[3])),delta1=0.4,delta2=0.4);
    
    X = list(X1, X2, X3); 
    temp = MGOFUSO(X)
    
    pSks_C1_200[l,] = pSks_C1_200[l,] + temp$decision.Skps/B
    pWks_C1_200[l,] = pWks_C1_200[l,] + temp$decision.Wkps/B
    pBon_C1_200[l,] = pBon_C1_200[l,] + temp$decision.Bon/B
    
    end = Sys.time(); 
    ave.interval = (end-start)/b
    if(b%%250==0){
      print(l)
      print(start); print("C1 n200")
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + (ave.interval*(B-b)) + ave.interval*B*(nC-l) ); 
    } 
  }
}

######## C3 n=200 #################3
nv = c(200,200,200); k = length(nv); 
alpha = 0.05; B = 1000; BB = 1000; 
Cases = 1:10; nC = length(Cases)
pSks_C3_200 = array(0,c(nC,3)); pWks_C3_200 = array(0,c(nC,3)); pBon_C3_200 = array(0,c(nC,3))
set.seed(0527202102)
for(l in Cases){
  start = Sys.time(); 
  for(b in 1:B){
    X1 = RIS(sort(runif(nv[1])), delta=1/2-l/20);
    X2 = sort(runif(nv[2])); 
    X3 = RLS(sort(runif(nv[3])), delta=1/2-l/20);
    
    X = list(X1, X2, X3); 
    temp = MGOFUSO(X)
    
    pSks_C3_200[l,] = pSks_C3_200[l,] + temp$decision.Skps/B
    pWks_C3_200[l,] = pWks_C3_200[l,] + temp$decision.Wkps/B
    pBon_C3_200[l,] = pBon_C3_200[l,] + temp$decision.Bon/B
    
    end = Sys.time(); 
    ave.interval = (end-start)/b
    if(b%%250==0){
      print(l)
      print(start); print("C3 n200")
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + (ave.interval*(B-b)) + ave.interval*B*(nC-l) ); 
    } 
  }
}

save(file="PCsk3n200.Rdata",
     pSks_C1_200, pWks_C1_200, pBon_C1_200,
     pSks_C3_200, pWks_C3_200, pBon_C3_200)

#pdf("Figure_2_GOF_PowerCurves_k3.pdf", width=19.5, heigh=13)

png("Figure_2_GOF_PowerCurves_k3_200.png", width=1872, heigh=1248)

par(mfrow=c(2,3))
par(mar=c(2,2,2,2)) 
plot((Cases-1)/10, pSks_C1_200[,1], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_200[,1],col="blue"); 
lines((Cases-1)/10, pBon_C1_200[,1],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C1_200[,2], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_200[,2],col="blue"); 
lines((Cases-1)/10, pBon_C1_200[,2],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C1_200[,3], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C1_200[,3],col="blue"); 
lines((Cases-1)/10, pBon_C1_200[,3],lty=2); 
abline(h=0.05,lty=3); 

plot((Cases-1)/10, pSks_C3_200[,1], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_200[,1],col="blue"); 
lines((Cases-1)/10, pBon_C3_200[,1],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C3_200[,2], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_200[,2],col="blue"); 
lines((Cases-1)/10, pBon_C3_200[,2],lty=2); 
abline(h=0.05,lty=3); 
plot((Cases-1)/10, pSks_C3_200[,3], type="l",ylim=c(0,1),col="green"); 
lines((Cases-1)/10, pWks_C3_200[,3],col="blue"); 
lines((Cases-1)/10, pBon_C3_200[,3],lty=2); 
abline(h=0.05,lty=3); 
par(mfcol=c(1,1))

dev.off()