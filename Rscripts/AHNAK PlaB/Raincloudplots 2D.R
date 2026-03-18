library(ggplot2)
require(plyr)
require(dplyr)
library(raincloudplots)
library(cowplot)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())

################################################################################
################################################################################
###############################functions########################################
################################################################################
################################################################################


importfile <- function (file,filefold,Data,inputs) {
  for(i in 1:length(file))
  {
    temp <- read.csv(file[i],header=FALSE)
    temp2 <- cbind(temp,t(replicate(dim(temp)[1], filefold[i,1:2])),inputs)
    Data <- rbind(Data,temp2)
  }
  return(Data)
}

importfile2 <- function (file,filefold,Data,inputs) {
  for(i in 1:length(file))
  {
    temp <- read.csv(file[i],header=FALSE)
    dist <- as.numeric(sqrt(temp[,2]^2+temp[,3]^2)*39.6825)
    temp2 <- cbind(temp[,1],rep(paste("1"),dim(temp)[1]),
                   rep(paste("1"),dim(temp)[1]), dist,
                   t(replicate(dim(temp)[1], filefold[i,1:2])),inputs)
    Data <- rbind(Data,temp2)
  }
  return(Data)
}

Distances <- function (Data,position) {
  if(position=="5mid")
  {
    Dist <- sqrt((Data$X1)^2+(Data$Y1)^2)
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],rep(paste("5mid"),dim(Data)[1]))
  }
  else if(position=="mid3")
  {
    Dist <- sqrt((Data$X2)^2+(Data$Y2)^2)
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],rep(paste("mid3"),dim(Data)[1]))
  }
  else if(position=="53")
  {
    Dist <- sqrt((Data$X1-Data$X2)^2+(Data$Y1-Data$Y2)^2)
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],rep(paste("53"),dim(Data)[1]))
  }
  
  names(Disttable) <- c( "No","X1","Y1","Value","Gene", "Field","Compartment")
  
  return(Disttable)
}

# %||%" <- function(a, b) {
#   if (!is.null(a)) a else b
# }

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}


################################################################################
################################################################################
################################################################################
################################################################################


cyt <- list.files(pattern = "Cytoplasmic Relative Coordinates.csv", recursive = TRUE)
nuc <- list.files(pattern = "Nuclear Relative Coordinates.csv", recursive = TRUE)
loc <- list.files(pattern = "Pixel shift.csv", recursive = TRUE)
mdn1 <- list.files(pattern = "Nuclear Distances.csv", recursive = TRUE)

split <- function(inputss){
  return(t(matrix(unlist(strsplit(dirname(inputss),"/")),ncol=length(strsplit(dirname(inputss),"/")))))
}

cytfold <- split(cyt)
nucfold <- split(nuc)
locfold <- split(loc)
mdn1fold <- split(loc)

Data1 <- data.frame()
Data2 <- data.frame()
Data3 <- data.frame()
Data5 <- data.frame()

Data1 <- importfile(cyt,cytfold,Data1,"Cytoplasm")
Data2 <- importfile(nuc,nucfold,Data2,"Nucleus")
Data3 <- importfile2(loc,locfold,Data3,"Coloc")
Data5 <- importfile2(mdn1,mdn1fold,Data5,"MDN1")

Data3$dist <-as.numeric(Data3$dist)
Data4 <- Data3

Data4[Data4$dist<100,] -> Data4

names(Data1) <- c( "No","X1","Y1","Z1","X2","Y2","Z2","Relative","Gene", "Field","Compartment")
names(Data2) <- c( "No","X1","Y1","Z1","X2","Y2","Z2","Relative","Gene", "Field","Compartment")
names(Data4) <- c( "No","X1","Y1","Value","Gene", "Field","Compartment")

Data5mid <- Distances(Data2,"5mid")
Datamid3 <- Distances(Data2,"mid3")
Data53 <- Distances(Data2,"53")

# Datanuc <- rbind(Data4,Data5mid,Datamid3,Data53)
Datanuc <- rbind(Data5mid,Datamid3,Data53)

Data5mid <- Distances(Data1,"5mid")
Datamid3 <- Distances(Data1,"mid3")
Data53 <- Distances(Data1,"53")

Datacyt <- rbind(Data5mid,Datamid3,Data53)


#######################################################################
#######################################################################
################## To plot different conditions #######################
#######################################################################
#######################################################################

# L <- c("Oddeven")
#L <- c("Oddeven2")
L <- c("AHNAK Untreated 5'-mid-3'")
L <- rbind(L,c("AHNAK Untreated 5'-mid-3'"))
#L <- rbind(L,c("Puro","Cytoplasm","DYNC1I2"))
L <- data.frame(L)
#L <- t(L)
colnames(L) <- c("Gene")

# Data <- Datanuc
Data <- Datanuc
Data %>% match_df(L) -> plotval
##Datanuc -> plotval


lens <- length(unique(Data$Compartment))+0.45
#orders <- c("Coloc")
orders <- unique(Data$Compartment)
dia_me <- ddply(plotval, .(Compartment),numcolwise(median))
dia_me$Value <- round(dia_me$Value,2)

dia_me <- dia_me[match(orders, dia_me$Compartment),]
orders <- rev(orders)

names <- names(plotval)

names[4] <- "y_axis"
names[5] <- "x_axis"

plotvals <- plotval
names(plotvals) <- names

dev.new()

ggplot(plotval,aes(x=Compartment,y=Value, fill = Compartment,
                   colour = Compartment))+
  geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2,
                   trim = FALSE,colour = "BLACK")+
  geom_point(position = position_jitter(width = .15), size = .35)+
  geom_boxplot(aes(x = Compartment, y = Value),outlier.shape = NA, 
               #position = position_nudge(x = 0.15, y = 0), 
               alpha = 0.3, width = .1, colour = "BLACK")+
  ylab('Distance (nm)')+
  xlab('Probe')+
  coord_flip()+
  theme_cowplot()+
  scale_x_discrete(limits=orders)+
  geom_text(data = dia_me,size=6, aes(label = Value, y=270, x=lens:1.25))+
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")
  #ggtitle("Figure 6: Change in Colour Palette")



