rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_GOFTests.r")
source("GOFCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=3; nv = rep(n,k); B=1000; Demo_b = 250; 
GOFCV0 = GOFTest_CV(nv=nv); 
MEB=1

Test00 = do_GOFTests(nv = nv, jh=c(0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test00

Test20 = do_GOFTests(nv = nv, jh=c(0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test20

Test40 = do_GOFTests(nv = nv, jh=c(0.4,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test40

Test42 = do_GOFTests(nv = nv, jh=c(0.4,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test42

Test64 = do_GOFTests(nv = nv, jh=c(0.6,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test64



TestN20 = do_GOFTests(nv = nv, jh=c(-0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN20

TestN40 = do_GOFTests(nv = nv, jh=c(-0.4,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN40

TestN60 = do_GOFTests(nv = nv, jh=c(-0.6,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN60

TestN80 = do_GOFTests(nv = nv, jh=c(-0.8,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN80

TestNo0 = do_GOFTests(nv = nv, jh=c(-1.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNo0



TestN22 = do_GOFTests(nv = nv, jh=c(-0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN22

TestN42 = do_GOFTests(nv = nv, jh=c(-0.4,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN42

TestN62 = do_GOFTests(nv = nv, jh=c(-0.6,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN62

TestN82 = do_GOFTests(nv = nv, jh=c(-0.8,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN82


TestN2N2 = do_GOFTests(nv = nv, jh=c(-0.2,-0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2N2

TestN4N2 = do_GOFTests(nv = nv, jh=c(-0.4,-0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4N2

TestN6N4 = do_GOFTests(nv = nv, jh=c(-0.6,-0.4), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN6N4

TestN8N6 = do_GOFTests(nv = nv, jh=c(-0.8,-0.6), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN8N6

TestNoN8 = do_GOFTests(nv = nv, jh=c(-1.0,-0.8), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNoN8


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestGOF_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}
