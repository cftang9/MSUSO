rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
source("rUSO_samples.r")
source("do_GOFTests.r")
source("GOFCV.r")
source("DR.r")
source("ME-B.r")

n=60; k=5; nv = rep(n,k); B=1000; Demo_b = 250; 
GOFCV0 = GOFTest_CV(nv=nv); 
MEB=1

Test0000 = do_GOFTests(nv = nv, jh=c(0.0,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test0000

Test2000 = do_GOFTests(nv = nv, jh=c(0.2,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test2000

Test4000 = do_GOFTests(nv = nv, jh=c(0.4,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test4000

Test2200 = do_GOFTests(nv = nv, jh=c(0.2,0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test2200

Test4200 = do_GOFTests(nv = nv, jh=c(0.4,0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test4200

Test2222 = do_GOFTests(nv = nv, jh=c(0.2,0.2,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); Test2222



TestN2000 = do_GOFTests(nv = nv, jh=c(-0.2,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2000

TestN4000 = do_GOFTests(nv = nv, jh=c(-0.4,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4000

TestN6000 = do_GOFTests(nv = nv, jh=c(-0.6,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN6000

TestN8000 = do_GOFTests(nv = nv, jh=c(-0.8,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN8000

TestNo000 = do_GOFTests(nv = nv, jh=c(-1.0,0.0,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNo000




TestN2002 = do_GOFTests(nv = nv, jh=c(-0.2,0.0,0.0,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2002

TestN4002 = do_GOFTests(nv = nv, jh=c(-0.4,0.0,0.0,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4002

TestN2N200 = do_GOFTests(nv = nv, jh=c(-0.2,-0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2N200

TestN4N200 = do_GOFTests(nv = nv, jh=c(-0.4,-0.2,0.0,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4N200

TestN6N4N20 = do_GOFTests(nv = nv, jh=c(-0.6,-0.4,-0.2,0.0), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN6N4N20

TestN8N6N4N2 = do_GOFTests(nv = nv, jh=c(-0.8,-0.6,-0.4,-0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN8N6N4N2

TestNoN8N6N4 = do_GOFTests(nv = nv, jh=c(-1.0,-0.8,-0.6,-0.4), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNoN8N6N4


TestN2022 = do_GOFTests(nv = nv, jh=c(-0.2,0.0,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2022

TestN4022 = do_GOFTests(nv = nv, jh=c(-0.4,0.0,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4022

TestN2N222 = do_GOFTests(nv = nv, jh=c(-0.2,-0.2,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN2N222

TestN4N222 = do_GOFTests(nv = nv, jh=c(-0.4,-0.2,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN4N222

TestN6N422 = do_GOFTests(nv = nv, jh=c(-0.6,-0.4,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN6N422

TestN8N622 = do_GOFTests(nv = nv, jh=c(-0.8,-0.6,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN8N622

TestNoN822 = do_GOFTests(nv = nv, jh=c(-1.0,-0.8,0.2,0.2), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestNoN822


TestN22rep = do_GOFTests(nv = nv, jh=c(rep(c(-0.2,0.2),2)), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN22rep

TestN42rep = do_GOFTests(nv = nv, jh=c(rep(c(-0.4,0.2),2)), GOFCV=GOFCV0, MEB=MEB, B=B, Demo_b=Demo_b); TestN42rep


if(1==0){
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}



if(1==0){
  load(file=paste0("TestGOF_k",k,"n",n,".Rdata"))
  varlist=ls()
  save(list=varlist[substring(varlist, 1, 4)=="Test"],file=paste0("TestGOF_k",k,"n",n,".Rdata"))
}
