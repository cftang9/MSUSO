rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
#source("MUSOLibrary.r")

n = 200
nv = c(n,n,n,n); k = length(nv); 
Table.row.200 = array(,c(7,9)); 
# Bon.cv.p1 = 0.6585602; # (n0=3000)
# Bon.cv.p2 = 0.7628116; # (n0=3000)
# Bon.cv.ps = 1.4588240; # (n0=3000)

B = 1000; 

jumps = c(0,0,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020221)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[1,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(1,0,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020222)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[2,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(0,1,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020223)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[3,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(1,1,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020224)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[4,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(-1,1,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020225)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[5,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(-1,0,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020226)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[6,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )


jumps = c(-1,-1,0); 
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(101020227)
start = Sys.time();
for(b in 1:B){
  X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0),
                X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1),  
                X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2),
                X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3))
  temp = MUSOGOF(X.data)
  power.Sks = power.Sks + temp$Sks/B;
  power.Wks = power.Wks + temp$Wks/B;
  power.Bon = power.Bon + temp$Bon/B;
  end = Sys.time(); 
  ave.interval = (end-start)/b
  if(b%%50==0){
    print(start); 
    print(power.Sks*B/b); 
    print(power.Wks*B/b); 
    print(power.Bon*B/b); 
    print(b); print(ave.interval); 
    print(Sys.time() + ave.interval*(B-b)); 
  }
}
Table.row.200[7,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                      power.Sks[2],power.Wks[2],power.Bon[2],   
                      power.Sks[3],power.Wks[3],power.Bon[3])
jumps
paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )

Table.row.200 = round(Table.row.200,3)
rn = 7; 
paste(Table.row.200[rn,1], "&", Table.row.200[rn,2], "&", Table.row.200[rn,3], "&", 
      Table.row.200[rn,4], "&", Table.row.200[rn,5], "&", Table.row.200[rn,6], "&", 
      Table.row.200[rn,7], "&", Table.row.200[rn,8], "&", Table.row.200[rn,9] )


