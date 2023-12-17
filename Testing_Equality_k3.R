rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_EqualityTests.r")
source("EqCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=3; nv = rep(n,k); ms=c(4,7,10); 
eqCV0 = EqTest_CV(nv=nv,ms=ms); 
MEB=0

Test00 = do_EqualityTests(nv = nv, jh=c(0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test00

Test20 = do_EqualityTests(nv = nv, jh=c(0.2,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test20

Test40 = do_EqualityTests(nv = nv, jh=c(0.4,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test40

Test60 = do_EqualityTests(nv = nv, jh=c(0.6,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test60

Test80 = do_EqualityTests(nv = nv, jh=c(0.8,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test80

Testo0 = do_EqualityTests(nv = nv, jh=c(1.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo0


Test22 = do_EqualityTests(nv = nv, jh=c(0.2,0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test22

Test42 = do_EqualityTests(nv = nv, jh=c(0.4,0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test42

Test64 = do_EqualityTests(nv = nv, jh=c(0.6,0.4), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test64

Test86 = do_EqualityTests(nv = nv, jh=c(0.8,0.6), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test86

Testo8 = do_EqualityTests(nv = nv, jh=c(1.0,0.8), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo8


Test0N2 = do_EqualityTests(nv = nv, jh=c(0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test0N2 

Test2N2 = do_EqualityTests(nv = nv, jh=c(0.2,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test2N2 

Test4N2 = do_EqualityTests(nv = nv, jh=c(0.4,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test4N2 

Test6N2 = do_EqualityTests(nv = nv, jh=c(0.6,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test6N2 


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestEq_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}

