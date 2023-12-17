rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_JumpDetection.r")
source("JDCV.r")

n=60; k=5; nv = rep(n,k); 
JDCV0 = JD_CV(nv=nv); 


Test0000 = do_JumpDetection(nv = nv, jh=c(0.0,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test0000

Test2000 = do_JumpDetection(nv = nv, jh=c(0.2,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test2000

Test4000 = do_JumpDetection(nv = nv, jh=c(0.4,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test4000

Test6000 = do_JumpDetection(nv = nv, jh=c(0.6,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test6000

Test8000 = do_JumpDetection(nv = nv, jh=c(0.8,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test8000

Testo000 = do_JumpDetection(nv = nv, jh=c(1.0,0.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Testo000


Test2200 = do_JumpDetection(nv = nv, jh=c(0.2,0.2,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test2200

Test4200 = do_JumpDetection(nv = nv, jh=c(0.4,0.2,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test4200

Test6400 = do_JumpDetection(nv = nv, jh=c(0.6,0.4,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test6400

Test8600 = do_JumpDetection(nv = nv, jh=c(0.8,0.6,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test8600

Test8800 = do_JumpDetection(nv = nv, jh=c(0.8,0.8,0.0,0.0), JDCV=JDCV0, Demo_b=250); Test8800

Testoo00 = do_JumpDetection(nv = nv, jh=c(1.0,1.0,0.0,0.0), JDCV=JDCV0, Demo_b=250); Testoo00

Test8640 = do_JumpDetection(nv = nv, jh=c(0.8,0.6,0.4,0.0), JDCV=JDCV0, Demo_b=250); Test8640

Testoo80 = do_JumpDetection(nv = nv, jh=c(1.0,1.0,0.8,0.0), JDCV=JDCV0, Demo_b=250); Testoo80

Test8642 = do_JumpDetection(nv = nv, jh=c(0.8,0.6,0.4,0.2), JDCV=JDCV0, Demo_b=250); Test8642

Testo864 = do_JumpDetection(nv = nv, jh=c(1.0,0.8,0.6,0.4), JDCV=JDCV0, Demo_b=250); Testo864

Testoo88 = do_JumpDetection(nv = nv, jh=c(1.0,1.0,0.8,0.8), JDCV=JDCV0, Demo_b=250); Testoo88

Testoooo = do_JumpDetection(nv = nv, jh=c(1.0,1.0,1.0,1.0), JDCV=JDCV0, Demo_b=250); Testoooo



if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestJump_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestJump_k",k,"n",n,".Rdata"))
}

