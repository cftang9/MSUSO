rm(list=ls(all=TRUE))
source("MUSOLibrary.R")
load("MFAP4.Rdata")

Rmn01 = Rmn(Data_MFAP4$F0, Data_MFAP4$F1);
Rmn12 = Rmn(Data_MFAP4$F1, Data_MFAP4$F2);
Rmn23 = Rmn(Data_MFAP4$F2, Data_MFAP4$F3);
Rmn34 = Rmn(Data_MFAP4$F3, Data_MFAP4$F4);

LSMRmn010 = LSMRmn0(Rmn01); LSMRmn011 = LSMRmn1(Rmn01);
LSMRmn120 = LSMRmn0(Rmn12); LSMRmn121 = LSMRmn1(Rmn12);
LSMRmn230 = LSMRmn0(Rmn23); LSMRmn231 = LSMRmn1(Rmn23);
LSMRmn340 = LSMRmn0(Rmn34); LSMRmn341 = LSMRmn1(Rmn34);

n0 = length(Data_MFAP4$F0); 
n1 = length(Data_MFAP4$F1); 
n2 = length(Data_MFAP4$F2); 
n3 = length(Data_MFAP4$F3); 
n4 = length(Data_MFAP4$F4); 

u01 = seq(0,1,by=1/n1)
u12 = seq(0,1,by=1/n2)
u23 = seq(0,1,by=1/n3)
u34 = seq(0,1,by=1/n4)

#pdf("Figure_3_MFAP4.pdf", width=13, heigh=13)

png("Figure_3_MFAP4.png", width=1248, heigh=1248)

par(mfrow=c(2,2))
par(mar=c(2,2,0.1,0.1))
### F0 vs F1
plot(u01,Rmn01,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.25)
lines(c(0,1),c(0,1),lty=3)
for(i in 1:(n1)){
  lines(c(u01[i],u01[i+1]), c(Rmn01[i+1],Rmn01[i+1]),col="blue")
  lines(c(u01[i],u01[i+1]), c(LSMRmn010[i],LSMRmn011[i]))
  polygon(c(u01[i],u01[i],u01[i+1],u01[i+1]), c(Rmn01[i+1],LSMRmn010[i],LSMRmn011[i],Rmn01[i+1]),col="yellow",border=FALSE)
}
legend("bottomright",legend=c(expression(hat(R)[1]),
                              expression(M~hat(R)[1]),
                              expression(R[0])),
       col=c("black","blue","black"),
       lty=c(1,1,3), title = "F0 vs F1")
### F1 vs F2
plot(u12,Rmn12,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.25)
lines(c(0,1),c(0,1),lty=3)
for(i in 1:(n2)){
  lines(c(u12[i],u12[i+1]), c(Rmn12[i+1],Rmn12[i+1]),col="blue")
  lines(c(u12[i],u12[i+1]), c(LSMRmn120[i],LSMRmn121[i]))
  polygon(c(u12[i],u12[i],u12[i+1],u12[i+1]), c(Rmn12[i+1],LSMRmn120[i],LSMRmn121[i],Rmn12[i+1]),col="yellow",border=FALSE)
}
legend("bottomright",legend=c(expression(hat(R)[2]),
                              expression(M~hat(R)[2]),
                              expression(R[0])),
       col=c("black","blue","black"),
       lty=c(1,1,3), title = "F1 vs F2")
### F2 vs F3
plot(u23,Rmn23,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.25)
lines(c(0,1),c(0,1),lty=3)
for(i in 1:(n3)){
  lines(c(u23[i],u23[i+1]), c(Rmn23[i+1],Rmn23[i+1]),col="blue")
  lines(c(u23[i],u23[i+1]), c(LSMRmn230[i],LSMRmn231[i]))
  polygon(c(u23[i],u23[i],u23[i+1],u23[i+1]), c(Rmn23[i+1],LSMRmn230[i],LSMRmn231[i],Rmn23[i+1]),col="yellow",border=FALSE)
}
legend("bottomright",legend=c(expression(hat(R)[3]),
                              expression(M~hat(R)[3]),
                              expression(R[0])),
       col=c("black","blue","black"),
       lty=c(1,1,3), title = "F2 vs F3")
### F3 vs F4
plot(u34,Rmn34,xlim=c(0,1),ylim=c(0,1),pch=20,cex=0.25)
lines(c(0,1),c(0,1),lty=3)
for(i in 1:(n4)){
  lines(c(u34[i],u34[i+1]), c(Rmn34[i+1],Rmn34[i+1]),col="blue")
  lines(c(u34[i],u34[i+1]), c(LSMRmn340[i],LSMRmn341[i]))
  polygon(c(u34[i],u34[i],u34[i+1],u34[i+1]), c(Rmn34[i+1],LSMRmn340[i],LSMRmn341[i],Rmn34[i+1]),col="yellow",border=FALSE)
}
legend("bottomright",legend=c(expression(hat(R)[4]),
                              expression(M~hat(R)[4]),
                              expression(R[0])),
       col=c("black","blue","black"),
       lty=c(1,1,3), title = "F3 vs F4")
par(mfrow=c(1,1))


dev.off()