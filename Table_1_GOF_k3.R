rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

### n = 200

nv = c(200,200,200); q = 0.2; 
B = 1000; BB = 1000; unif.n = 1000; L = 2000; alpha = 0.05; 
Table.row.200 = array(,c(7,9)); 

### G2 G2 G2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220200)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2,2)) 
  temp = MUSODecision(X.data)
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
Table.row.200[1,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G2 G2.6 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220201)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2.6,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.200[2,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G2 G2 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220202)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.200[3,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G2 G2.6 G3.2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220203)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2.6,3.2)) 
  temp = MUSODecision(X.data)
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
Table.row.200[4,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G3.2 G2.6 G3.2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220204)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,3.2)) 
  temp = MUSODecision(X.data)
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
Table.row.200[5,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G3.2 G2.6 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220205)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.200[6,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G3.2 G 2.6 G 2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220206)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,2)) 
  temp = MUSODecision(X.data)
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
Table.row.200[7,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

row.names(Table.row.200) <- c("G2=G2=G2", "G2<G2.6=G2.6", "G2=G2=G2.6", "G2<G2.6=G3.2", "G3.2!<G2.6<G3.2", "G3.2!<G2.6=G2.6", "G3.2!<G2.6!<G2")
colnames(Table.row.200) <- c("sp1","wp1","bp1","sp2","wp2","bp2","sps","wps","bps")
xtable(Table.row.200,digits=3)
kable(Table.row.200, caption = "size and power with mixture normal, k=3, n=200")



### n = 400
nv = c(400,400,400); q = 0.2; 
B = 1000; BB = 1000; unif.n = 1000; L = 2000; alpha = 0.05; 
Table.row.400 = array(,c(7,9)); 

### G2 G2 G2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220200)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2,2)) 
  temp = MUSODecision(X.data)
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
Table.row.400[1,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G2 G2.6 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220201)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2.6,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.400[2,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])


### G2 G2 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220202)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.400[3,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G2 G2.6 G3.2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220203)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(2,2.6,3.2)) 
  temp = MUSODecision(X.data)
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
Table.row.400[4,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])

### G3.2 G2.6 G3.2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220204)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,3.2)) 
  temp = MUSODecision(X.data)
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
Table.row.400[5,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])


### G3.2 G2.6 G2.6
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220205)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,2.6)) 
  temp = MUSODecision(X.data)
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
Table.row.400[6,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])


### G3.2 G2.6 G2
power.Sks = array(0,3); power.Wks = array(0,3); power.Bon = array(0,3); 
set.seed(021220206)
start = Sys.time();
for(b in 1:B){
  X.data = MixNormal(nv,rep(q,3), m=c(3.2,2.6,2)) 
  temp = MUSODecision(X.data)
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
Table.row.400[7,] = c(power.Sks[1],power.Wks[1],power.Bon[1],   power.Sks[2],power.Wks[2],power.Bon[2],   power.Sks[3],power.Wks[3],power.Bon[3])



row.names(Table.row.400) <- c("G2=G2=G2", "G2<G2.6=G2.6", "G2=G2=G2.6", "G2<G2.6=G3.2", "G3.2!<G2.6<G3.2", "G3.2!<G2.6=G2.6", "G3.2!<G2.6!<G2")
colnames(Table.row.400) <- c("sp1","wp1","bp1","sp2","wp2","bp2","sps","wps","bps")
xtable(Table.row.400,digits=3)
kable(Table.row.400, caption = "size and power with mixture normal, k=3, n=400")

