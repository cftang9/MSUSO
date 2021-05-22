rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

u = seq(0,1,by=1/1000)
R1 = InvQBezier(u,delta=0.4,delta2=0.4); 
R2 = InvQBezier(u,delta=0.4,delta2=0.8); 
R3 = InvQBezier(u,delta=0.8,delta2=0.4); 
R4 = InvQBezier(u,delta=0.8,delta2=0.8); 

#pdf("Figure_1_ODCs.pdf", width=6, heigh=6)
png("Supp_Figure_1_ODCs.png", width=576, heigh=576)

par(mar=c(2,2,0.1,0.1))
plot(R1,u,type="l")
lines(R2,u)
lines(R3,u)
lines(R4,u)
lines(u,R1)
lines(c(0,1),c(0,1),lty=3)

text(0.8,0.84,expression(R[0]))
text(0.67,0.83,expression(R[1]))
text(0.45,0.7,expression(R[2]))
text(0.61,0.9,expression(R[3]))
text(0.25,0.85,expression(R[4]))
text(0.87,0.66,expression(R[1]^-1))

dev.off()
