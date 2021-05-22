rm(list=ls(all=TRUE))
source("MUSOLibrary.R")


rm(list=ls(all=TRUE))
library(xtable)
working.dir = getwd(); 
setwd(".."); parent = getwd(); 
setwd(working.dir)
source(paste0(parent,"/MUSOLibrary.R"))
source(paste0(parent,"/CorreCurves.R"))
#source(paste0(parent,"/AreaBasedLibrary.R"))

u = seq(0,1,by=1/1000)

R1 = InvQBezier(u,delta=0.4,delta2=0.4); 
R2 = InvQBezier(u,delta=0.4,delta2=0.8); 
R3 = InvQBezier(u,delta=0.8,delta2=0.4); 
R4 = InvQBezier(u,delta=0.8,delta2=0.8); 

#pdf("Figure_1_ODCs.pdf", width=13, heigh=6)

png("Figure_1_ODCs.png", width=1248, heigh=576)

par(mfrow=c(1,2))

par(mar=c(2,2,0.1,0.1))
plot(R1,u,type="l")
#lines(R2,u)
#lines(R3,u)
#lines(R4,u)
lines(u,R1)
lines(c(0,1),c(0,1),lty=3)

text(0.8,0.84,expression(R[0]))
text(0.65,0.84,expression(R[1]))
#text(0.45,0.7,expression(R[2]))
#text(0.61,0.9,expression(R[3]))
#text(0.25,0.85,expression(R[4]))
text(0.87,0.66,expression(R[1]^-1))



u = seq(0,1,1/100); 
plot(RIS(u,delta=1/2),u,type="l"); 
for(l in 1:9){
  lines(RIS(u,delta=1/2-l/20),u)
}

text(0.54,0.96, expression(paste(delta , " = 0" )))
arrows(0.48,0.93,0.12,0.57,length=0.1)
text(0.10,0.52, expression(paste(delta , " = 9" )))




par(mfrow=c(1,1))

dev.off()
