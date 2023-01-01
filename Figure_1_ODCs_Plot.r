rm(list = ls(all=TRUE))
source("EGJ_USO_Library.R")

#pdf(file="./Figure_1_ODC_Plot.pdf", width = 12, height = 6)
par(mfrow=c(1,2))
u = seq(0,1,by=0.01)
par(mar=c(2,2,.1,.1))
plot(RIS(u,delta=1/2),u,type="l"); 
text(0.54,0.96, expression(G[0]))
lines(u,RIS(u,delta=1/2),type="l"); 
text(0.91,0.51, expression(G[0]^{-1}))
lines(c(0,1),c(0,1))
text(0.64,0.70, expression(R[0]))

u = seq(0,1,1/100); 
plot(RIS(u,delta=1/2),u,type="l"); 
for(l in 1:9){
  lines(RIS(u,delta=1/2-l/20),u)
}

text(0.54,0.96, expression(paste(delta , " = 0" )))
arrows(0.48,0.93,0.12,0.57,length=0.1)
text(0.10,0.52, expression(paste(delta , " = 9" )))

#graphics.off()
