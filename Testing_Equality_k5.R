rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_EqualityTests.r")
source("EqCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=5; nv = rep(n,k); ms=c(4,7,10); 
eqCV0 = EqTest_CV(nv=nv,ms=ms); 
MEB=1

Test0000 = do_EqualityTests(nv = nv, jh=c(0.0,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test0000

Test2000 = do_EqualityTests(nv = nv, jh=c(0.2,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test2000

Test4000 = do_EqualityTests(nv = nv, jh=c(0.4,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test4000

Test6000 = do_EqualityTests(nv = nv, jh=c(0.6,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test6000

Test8000 = do_EqualityTests(nv = nv, jh=c(0.8,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test8000

Testo000 = do_EqualityTests(nv = nv, jh=c(1.0,0.0,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo000


Test2200 = do_EqualityTests(nv = nv, jh=c(0.2,0.2,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test2200

Test4200 = do_EqualityTests(nv = nv, jh=c(0.4,0.2,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test4200

Test6400 = do_EqualityTests(nv = nv, jh=c(0.6,0.4,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test6400

Test8600 = do_EqualityTests(nv = nv, jh=c(0.8,0.6,0.0,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test8600

Test8640 = do_EqualityTests(nv = nv, jh=c(0.8,0.6,0.4,0.0), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test8640

Test8642 = do_EqualityTests(nv = nv, jh=c(0.8,0.6,0.4,0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test8642

Testo864 = do_EqualityTests(nv = nv, jh=c(1.0,0.8,0.6,0.4), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo864


Test000N2 = do_EqualityTests(nv = nv, jh=c(0.0,0.0,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test000N2 

Test200N2 = do_EqualityTests(nv = nv, jh=c(0.2,0.0,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test200N2 

Test220N2 = do_EqualityTests(nv = nv, jh=c(0.2,0.2,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test220N2 

Test420N2 = do_EqualityTests(nv = nv, jh=c(0.4,0.2,0.0,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test420N2 

Test642N2 = do_EqualityTests(nv = nv, jh=c(0.6,0.4,0.2,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Test642N2 

Testo86N2 = do_EqualityTests(nv = nv, jh=c(1.0,0.8,0.6,-0.2), eqCV=eqCV0, ms=ms, MEB=MEB, Demo_b=250); Testo86N2 


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}

if(1==0){
  load(file=paste0("TestEq_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestEq_k",k,"n",n,".Rdata"))
}
