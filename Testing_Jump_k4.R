rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_JumpDetection.r")
source("JDCV.r")

n=60; k=4; nv = rep(n,k); 
JDCV0 = JD_CV(nv=nv); 

Test000 = do_JumpDetection(nv = nv, jh=c(0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test000

Test200 = do_JumpDetection(nv = nv, jh=c(0.2,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test200

Test400 = do_JumpDetection(nv = nv, jh=c(0.4,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test400

Test600 = do_JumpDetection(nv = nv, jh=c(0.6,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test600

Test800 = do_JumpDetection(nv = nv, jh=c(0.8,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test800

Testo00 = do_JumpDetection(nv = nv, jh=c(1.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Testo00


Test220 = do_JumpDetection(nv = nv, jh=c(0.2,0.2,0.0), JDCV=JDCV0, Demo_b=250); Test220

Test420 = do_JumpDetection(nv = nv, jh=c(0.4,0.2,0.0), JDCV=JDCV0, Demo_b=250); Test420

Test640 = do_JumpDetection(nv = nv, jh=c(0.6,0.4,0.0), JDCV=JDCV0, Demo_b=250); Test640

Test860 = do_JumpDetection(nv = nv, jh=c(0.8,0.6,0.0), JDCV=JDCV0, Demo_b=250); Test860

Testo80 = do_JumpDetection(nv = nv, jh=c(1.0,0.8,0.0), JDCV=JDCV0, Demo_b=250); Testo80

Testoo0 = do_JumpDetection(nv = nv, jh=c(1.0,1.0,0.0), JDCV=JDCV0, Demo_b=250); Testoo0

Test642 = do_JumpDetection(nv = nv, jh=c(0.6,0.4,0.2), JDCV=JDCV0, Demo_b=250); Test642

Test864 = do_JumpDetection(nv = nv, jh=c(0.8,0.6,0.4), JDCV=JDCV0, Demo_b=250); Test864

Testo86 = do_JumpDetection(nv = nv, jh=c(1.0,0.8,0.6), JDCV=JDCV0, Demo_b=250); Testo86

Testo88 = do_JumpDetection(nv = nv, jh=c(1.0,0.8,0.8), JDCV=JDCV0, Demo_b=250); Testo88

Testooo = do_JumpDetection(nv = nv, jh=c(1.0,1.0,1.0), JDCV=JDCV0, Demo_b=250); Testooo

if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestJump_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}

