rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_JumpDetection.r")
source("JDCV.r")

n=60; k=3; nv = rep(n,k); 
JDCV0 = JD_CV(nv=nv); 

Test00 = do_JumpDetection(nv = nv, jh=c(0.0,0.0), JDCV=JDCV0, Demo_b=250); Test00

Test20 = do_JumpDetection(nv = nv, jh=c(0.2,0.0), JDCV=JDCV0, Demo_b=250); Test20

Test40 = do_JumpDetection(nv = nv, jh=c(0.4,0.0), JDCV=JDCV0, Demo_b=250); Test40

Test60 = do_JumpDetection(nv = nv, jh=c(0.6,0.0), JDCV=JDCV0, Demo_b=250); Test60

Test80 = do_JumpDetection(nv = nv, jh=c(0.8,0.0), JDCV=JDCV0, Demo_b=250); Test80

Testo0 = do_JumpDetection(nv = nv, jh=c(1.0,0.0), JDCV=JDCV0, Demo_b=250); Testo0


Test22 = do_JumpDetection(nv = nv, jh=c(0.2,0.2), JDCV=JDCV0, Demo_b=250); Test22

Test42 = do_JumpDetection(nv = nv, jh=c(0.4,0.2), JDCV=JDCV0, Demo_b=250); Test42

Test64 = do_JumpDetection(nv = nv, jh=c(0.6,0.4), JDCV=JDCV0, Demo_b=250); Test64

Test86 = do_JumpDetection(nv = nv, jh=c(0.8,0.6), JDCV=JDCV0, Demo_b=250); Test86

Testo8 = do_JumpDetection(nv = nv, jh=c(1.0,0.8), JDCV=JDCV0, Demo_b=250); Testo8

Testoo = do_JumpDetection(nv = nv, jh=c(1.0,1.0), JDCV=JDCV0, Demo_b=250); Testoo


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestJump_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}

