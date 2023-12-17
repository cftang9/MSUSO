rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_GOFTests.r")
source("GOFCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=4; nv = rep(n,k); B=1000; Demo_b = 250; 
GOFCV0 = GOFTest_CV(nv=nv); 
MEB=1

Test000 = do_GOFTests(nv = nv, jh=c(0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test000

Test200 = do_GOFTests(nv = nv, jh=c(0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test200

Test400 = do_GOFTests(nv = nv, jh=c(0.4,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test400

Test220 = do_GOFTests(nv = nv, jh=c(0.2,0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test220

Test420 = do_GOFTests(nv = nv, jh=c(0.4,0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test420

Test222 = do_GOFTests(nv = nv, jh=c(0.2,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test222



TestN200 = do_GOFTests(nv = nv, jh=c(-0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN200

TestN400 = do_GOFTests(nv = nv, jh=c(-0.4,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN400

TestN600 = do_GOFTests(nv = nv, jh=c(-0.6,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN600

TestN800 = do_GOFTests(nv = nv, jh=c(-0.8,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN800

TestNo00 = do_GOFTests(nv = nv, jh=c(-1.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNo00




TestN202 = do_GOFTests(nv = nv, jh=c(-0.2,0.0,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN202

TestN402 = do_GOFTests(nv = nv, jh=c(-0.4,0.0,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN402

TestN2N20 = do_GOFTests(nv = nv, jh=c(-0.2,-0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2N20

TestN4N20 = do_GOFTests(nv = nv, jh=c(-0.4,-0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4N20

TestN6N4N2 = do_GOFTests(nv = nv, jh=c(-0.6,-0.4,-0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN6N4N2

TestN8N6N4 = do_GOFTests(nv = nv, jh=c(-0.8,-0.6,-0.4), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN8N6N4

TestNoN8N6 = do_GOFTests(nv = nv, jh=c(-1.0,-0.8,-0.6), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNoN8N6


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}

if(1==0){
  load(file=paste0("TestGOF_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}
