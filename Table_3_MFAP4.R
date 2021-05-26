rm(list=ls(all=TRUE))
source("MUSOLibrary.R")

load("MFAP4.Rdata")

set.seed(05222021)
MDDUSO(Data_MFAP4)

set.seed(05222021)
MGOFUSO(Data_MFAP4)


