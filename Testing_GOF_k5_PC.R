rm(list = ls(all=TRUE))
source("EGJ_USO_Library.r")
#source("MUSOLibrary.r")

n = 60; 
nv = c(n,n,n,n,n); k = length(nv); 
Table.row.200 = array(,c(10,10)); 
# Bon.cv.p1 = 0.6585602; # (n0=3000)
# Bon.cv.p2 = 0.7628116; # (n0=3000)
# Bon.cv.ps = 1.4588240; # (n0=3000)

B = 1000; jumps = c(1,0,0,0); 
for(j in 0:9){
  power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
  set.seed(101020220+j)
  start = Sys.time();
  for(b in 1:B){
    X.data = list(X1 = rRILS_Jumps(n=nv[1],jumps=jumps,at=0,delta=(1-j/10)/2),
                  X2 = rRILS_Jumps(n=nv[2],jumps=jumps,at=1,delta=(1-j/10)/2),  
                  X3 = rRILS_Jumps(n=nv[3],jumps=jumps,at=2,delta=(1-j/10)/2), 
                  X4 = rRILS_Jumps(n=nv[4],jumps=jumps,at=3,delta=(1-j/10)/2),
                  X5 = rRILS_Jumps(n=nv[5],jumps=jumps,at=4,delta=(1-j/10)/2))
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
      print(Sys.time() + ave.interval*(B-b) + (9-j)*ave.interval*B); 
    }
  }
  Table.row.200[j+1,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   
                          power.Sks[2],power.Wks[2],power.Bon[2],   
                          power.Sks[3],power.Wks[3],power.Bon[3])
  jumps
}
Table.row.200 = round(Table.row.200,3)

paste(power.Sks[1], "&", power.Wks[1], "&", power.Bon[1], "&", 
      power.Sks[2], "&", power.Wks[2], "&", power.Bon[2], "&", 
      power.Sks[3], "&", power.Wks[3], "&", power.Bon[3] )

#save(file="k5r1000n60.R",Table.row.200)

rn = 7; 
paste(Table.row.200[rn,1], "&", Table.row.200[rn,2], "&", Table.row.200[rn,3], "&", 
      Table.row.200[rn,4], "&", Table.row.200[rn,5], "&", Table.row.200[rn,6], "&", 
      Table.row.200[rn,7], "&", Table.row.200[rn,8], "&", Table.row.200[rn,9] )






