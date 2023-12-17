rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_EqualityTests.r")
source("EqCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=4; nv = rep(n,k); ms=c(4,7,10); 
eqCV0 = EqTest_CV(nv=nv,ms=ms); 
MEB=1

Test000 = do_EqualityTests(nv = nv, jh=c(0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test000

Test200 = do_EqualityTests(nv = nv, jh=c(0.2,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test200

Test400 = do_EqualityTests(nv = nv, jh=c(0.4,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test400

Test600 = do_EqualityTests(nv = nv, jh=c(0.6,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test600

Test800 = do_EqualityTests(nv = nv, jh=c(0.8,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test800

Testo00 = do_EqualityTests(nv = nv, jh=c(1.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo00


Test220 = do_EqualityTests(nv = nv, jh=c(0.2,0.2,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test220

Test420 = do_EqualityTests(nv = nv, jh=c(0.4,0.2,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test420

Test640 = do_EqualityTests(nv = nv, jh=c(0.6,0.4,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test640

Test860 = do_EqualityTests(nv = nv, jh=c(0.8,0.6,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test860

Test642 = do_EqualityTests(nv = nv, jh=c(0.6,0.4,0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test642

Test864 = do_EqualityTests(nv = nv, jh=c(0.8,0.6,0.4), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test864


Test00N2 = do_EqualityTests(nv = nv, jh=c(0.0,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test00N2 

Test20N2 = do_EqualityTests(nv = nv, jh=c(0.2,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test20N2 

Test22N2 = do_EqualityTests(nv = nv, jh=c(0.2,0.2,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test22N2 

Test40N2 = do_EqualityTests(nv = nv, jh=c(0.4,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test40N2 

Test64N2 = do_EqualityTests(nv = nv, jh=c(0.6,0.4,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test64N2 

Testo8N2 = do_EqualityTests(nv = nv, jh=c(1.0,0.8,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo8N2 


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}

if(1==0){
  load(file=paste0("TestEq_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}
