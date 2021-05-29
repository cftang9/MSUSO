rm(list=ls(all=TRUE))
source("MUSOLibrary.R")


####################################################################### 
# n = 200

nv = c(200,200,200,200); 
B = 1000; BB = 1000; 
Table.GOF.K4.RR.200 = array(,c(7,9)); 
start = Sys.time();
for(j in 1:7){
  power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
  set.seed(0212202000+(j-1))
  for(b in 1:B){
    X.data = Gen.RR(n=nv,case=j)
    temp = MGOFUSO(X.data) 
    power.Sks = power.Sks + temp$decision.Skps/B;
    power.Wks = power.Wks + temp$decision.Wkps/B;
    power.Bon = power.Bon + temp$decision.Bon/B; 
    
    end = Sys.time(); 
    ave.interval = (end-start)/((j-1)*B+b)
    if(b%%100==0){
      print(start); print(c("RR",nv,j))
      print(round(c(power.Sks*B/b,power.Wks*B/b,power.Bon*B/b),3)); 
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + ave.interval*(B-b) + (7-j)*B*ave.interval); 
    } 
  }
  Table.GOF.K4.RR.200[j,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])
}


nv = c(200,200,200,200); 
B = 1000; BB = 1000; 
Table.GOF.K4.MN.200 = array(,c(7,9)); 
start = Sys.time();
for(j in 1:7){
  power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
  set.seed(0212202000+(j-1)*10)
  for(b in 1:B){
    X.data = Gen.MN(n=nv,case=j)
    temp = MGOFUSO(X.data) 
    power.Sks = power.Sks + temp$decision.Skps/B;
    power.Wks = power.Wks + temp$decision.Wkps/B;
    power.Bon = power.Bon + temp$decision.Bon/B; 
    
    end = Sys.time(); 
    ave.interval = (end-start)/((j-1)*B+b)
    if(b%%100==0){
      print(start); print(c("MN",nv,j))
      print(round(c(power.Sks*B/b,power.Wks*B/b,power.Bon*B/b),3)); 
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + ave.interval*(B-b) + (7-j)*B*ave.interval); 
    } 
  }
  Table.GOF.K4.MN.200[j,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])
}


#########################################################################3
# n = 400

nv = c(400,400,400,400);
B = 1000; BB = 1000; 
Table.GOF.K4.RR.400 = array(,c(7,9)); 
start = Sys.time();
for(j in 1:7){
  power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
  set.seed(0212202000+(j-1)*20)
  for(b in 1:B){
    X.data = Gen.RR(n=nv,case=j)
    temp = MGOFUSO(X.data) 
    power.Sks = power.Sks + temp$decision.Skps/B;
    power.Wks = power.Wks + temp$decision.Wkps/B;
    power.Bon = power.Bon + temp$decision.Bon/B; 
    
    end = Sys.time(); 
    ave.interval = (end-start)/((j-1)*B+b)
    if(b%%100==0){
      print(start); print(c("RR",nv,j))
      print(round(c(power.Sks*B/b,power.Wks*B/b,power.Bon*B/b),3)); 
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + ave.interval*(B-b) + (7-j)*B*ave.interval); 
    } 
  }
  Table.GOF.K4.RR.400[j,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])
}

nv = c(400,400,400,400);
B = 1000; BB = 1000; 
Table.GOF.K4.MN.400 = array(,c(7,9)); 
start = Sys.time();
for(j in 1:7){
  power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
  set.seed(0212202000+(j-1)*30)
  for(b in 1:B){
    X.data = Gen.MN(n=nv,case=j)
    temp = MGOFUSO(X.data) 
    power.Sks = power.Sks + temp$decision.Skps/B;
    power.Wks = power.Wks + temp$decision.Wkps/B;
    power.Bon = power.Bon + temp$decision.Bon/B; 
    
    end = Sys.time(); 
    ave.interval = (end-start)/((j-1)*B+b)
    if(b%%100==0){
      print(start); print(c("RR",nv,j))
      print(round(c(power.Sks*B/b,power.Wks*B/b,power.Bon*B/b),3)); 
      print(b); print(ave.interval); 
      print(Sys.time() + ave.interval*(B-b)); 
      print(Sys.time() + ave.interval*(B-b) + (7-j)*B*ave.interval); 
    } 
  }
  Table.GOF.K4.MN.400[j,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])
}

# save(file="SPsk4.Rdata", Table.GOF.K4.RR.200, Table.GOF.K4.RR.400, Table.GOF.K4.MN.200, Table.GOF.K4.MN.400)