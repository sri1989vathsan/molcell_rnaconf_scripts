require(ggplot2)
require(plyr)
require(dplyr)

this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
rm(list = ls())

cyt <- list.files(pattern = "Mask colocalization.csv", recursive = TRUE)
# nuc <- list.files(pattern = "Nuclear Distances.csv", recursive = TRUE)
# loc <- list.files(pattern = "^Distances.csv", recursive = TRUE)
# cytcor <- list.files(pattern = "Cytoplasmic Relative Coordinates.csv", recursive = TRUE)
# nuccor <- list.files(pattern = "Nuclear Relative Coordinates.csv", recursive = TRUE)
# loccor <- list.files(pattern = "^Pixel shift.csv", recursive = TRUE)


split <- function(inputss){
  return(t(matrix(unlist(strsplit(dirname(inputss),"/")),ncol=length(strsplit(dirname(inputss),"/")))))
}

cytfold <- split(cyt)
# nucfold <- split(nuc)
# cytcorfold <- split(cytcor)
# nuccorfold <- split(nuccor)
# locfold <- split(loc)
# loccorfold <- split(loccor)

Data <- data.frame()
# Dataloc <- data.frame()
# Datacor <-data.frame()
# Dataloccor <- data.frame()

importfile <- function (file,filefold,Data,inputs) {
  for(i in 1:length(file))
  {
    temp <- read.csv(file[i],header=TRUE)
    temp2 <- cbind(temp,t(replicate(dim(temp)[1], filefold[i,1:2])),inputs)
    Data <- rbind(Data,temp2)
  }
  return(Data)
}

Data <- importfile(cyt,cytfold,Data,"Cytoplasm")
# Data <- importfile(nuc,nucfold,Data,"Nucleus")
# Dataloc <- importfile(loc,locfold,Dataloc,"Coloc")
# 
# Datacor <- importfile(cytcor,cytcorfold,Datacor,"Cytoplasm")
# Datacor <- importfile(nuccor,nuccorfold,Datacor,"Nucleus")
# Dataloccor <- importfile(loccor,loccorfold,Dataloccor,"Coloc")


names(Data) <- c("#","5p","mid","exon","5p-mid","5p-exon","mid-exon",
                 "5p-mid-exon","Dataset","Field","Compartment")

Data2nd <- Data[Data$Dataset==unique(Data$Dataset)[1],]
Datalast <-Data[Data$Dataset==unique(Data$Dataset)[2],]

Datas <- Data2nd

Data2 <- Datas[Datas$`5p`==1,]
allintron <- dim(Data2)[1]
intronexon <- dim(Data2[Data2$`5p-exon`>1,])[1]

perc <- intronexon*100/allintron


Dataslast <- Datalast
Data2last <- Dataslast[Dataslast$`5p`==1,]
allintronlast <- dim(Data2last)[1]
intronexonlast <- dim(Data2last[Data2last$`5p-exon`>1,])[1]

perclast <- intronexonlast*100/allintronlast

# names(Datacor) <- c("#","X1","Y1","Z1","X2","Y2","Z2","Reference","Gene","Field","Compartment")
# names(Dataloccor) <- c("#","X1","Y1","Z1","Gene","Field","Compartment")
# 
# Dataoldcor <- Datacor
# 
# ##Datanewcor <- Datacor
# ##Datanewcor$X1 <- Datacor$X1 - median(Datacor$X1)
# ##Datanewcor$Y1 <- Datacor$Y1 - median(Datacor$Y1)
# 
# ##Datanewcor$X2 <- Datacor$X2 - median(Datacor$X2)
# ##Datanewcor$Y2 <- Datacor$Y2 - median(Datacor$Y2)
# 
# Datacor %>% group_by(Gene) %>% mutate(X1 = X1-median(X1),X2 = X2-median(X2),
#                                       Y1 = Y1-median(Y1),Y2 = Y2-median(Y2))-> Datanewcor
# 
# 
# 
# Datacor <- Datanewcor
# # Data2 <- Data
# # distancexy <- 39.6825*sqrt(Data2[,2]^2+Data2[,3]^2)
# # distancexyz <- sqrt((39.6825*Data2[,2])^2+(39.6825*Data2[,3])^2+(180*Data2[,4])^2)
# # Data <- cbind(Data2[,1:4],distancexy,distancexyz,Data2[,5:7])
# 
# dist2d <- Datacor$'#'
# dist2d <- as.data.frame(dist2d)
# mid5p <- sqrt((Datacor$X1)^2+(Datacor$Y1)^2)
# dist2d$mid5p <- mid5p
# mid3p <- sqrt((Datacor$X2)^2+(Datacor$Y2)^2)
# dist2d$mid3p <- mid3p
# p53p <- sqrt((Datacor$X1-Datacor$X2)^2+(Datacor$Y1-Datacor$Y2)^2)
# dist2d$p53p <- p53p
# dist2d$Compartment <- Datacor$Compartment
# dist2d$Reference <- Datacor$Reference
# dist2d$Gene <- Datacor$Gene
# distlocx <- (Dataloccor$X1)^2
# distlocy <- (Dataloccor$Y1)^2
# distloc <- data.frame(Dataloccor[,1],39.6925*sqrt(distlocx+distlocy),Dataloccor[,5:7],Dataloccor[,7])
# Dist1 <- cbind(dist2d[,c(1:2,5:7)],data.frame(rep('mid5p',nrow(dist2d))))
# Dist2 <- cbind(dist2d[,c(1,3,5:7)],data.frame(rep('mid3p',nrow(dist2d))))
# Dist3 <- cbind(dist2d[,c(1,4,5:7)],data.frame(rep('p53p',nrow(dist2d))))
# names(Dist1) <- c("#","Distance","Compartment","Reference","Gene","Relative")
# names(Dist2) <- c("#","Distance","Compartment","Reference","Gene","Relative")
# names(Dist3) <- c("#","Distance","Compartment","Reference","Gene","Relative")
# names(distloc) <- c("#","Distance","Compartment","Reference","Gene","Relative")
# Dist <- rbind(Dist1,Dist2,Dist3,distloc)
# 
# names(Data) <- c( "No","X1","Y1","Z1","Value","Gene", "Field","Compartment")

#######################################################################
#######################################################################
################## To plot different conditions #######################
#######################################################################
#######################################################################

L <- c("Coloc","Coloc","Oddeven")
L <- rbind(L,c("CDK6 5'-mid-3'","mid5p","Nucleus"))
L <- rbind(L,c("CDK6 5'-mid-3'","mid3p","Nucleus"))
L <- rbind(L,c("CDK6 5'-mid-3'","p53p","Nucleus"))
L <- data.frame(L)
##L <- t(L)
colnames(L) <- c("Gene","Relative","Compartment")

L <- as.data.frame(L)
Dist %>% match_df(L) -> plotval

##Data -> plotval

##ist -> plotval
##lens <- length(unique(plotval$Relative))+0.25
lens <- length(L$Relative)+0.25
#orders <- c("Coloc")
orders <- unique(plotval$Relative)
dia_me <- ddply(plotval, .(Relative),numcolwise(median))
dia_me$Value <- round(dia_me$Distance,2)

dia_me <- dia_me[match(orders, dia_me$Relative),]
orders <- rev(orders)



dev.new()
ggplot(plotval,aes(x=Relative,y=Distance))+
  geom_violin(trim=FALSE,aes(fill = factor(Relative)))+#,scale="width")+
  geom_boxplot(width=0.05, fill="white",outlier.colour=NA)+
  labs(title="",x="", y = "Distances (nm)")+
  # facet_grid(. ~ Compartment, scales="free") +
  scale_y_continuous(limits =c(0,300))+#, breaks =seq(from = 0, to = 300, by = 50), minor_breaks =seq(from = 25, to = 325, by = 50)) +#,limits =c(-100,100))  +
  coord_flip()+
  scale_x_discrete(limits=orders)+#,
  geom_text(data = dia_me,size=6, aes(label = Value, y=270, x=lens:1.25))+
  
  #labels=c("Cyclo100"="Cyclo","Untreated"="Untreated","Coloc"="Coloc")) +
  # stat_summary(fun.data = n_fun, geom = "text",hjust=Inf)+
  # geom_text(data = data.frame(), aes(x = names(meds) , y = meds, label = paste("n =", meds))) +
  #scale_x_continuous(breaks =c(-100,0,100), minor_breaks =c(-50,50),limits = c(-100,100)) +
  theme(axis.line = element_line(colour = "black"),
        title = element_text(size=18),
        #plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #panel.grid.major = element_line(color ="black",size=2),
        #panel.grid.minor = element_line(color ="black"),
        # panel.border = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        strip.text.y = element_blank(),
        legend.position="none",
        axis.title.x = element_text(face="bold", size=18),
        # axis.text.x  = element_text(face="bold", size=15),
        axis.text.x  = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=18,margin=margin(0,+10,0,5)),
        axis.text.y  = element_text(face="bold", size=15),
        plot.margin=unit(c(-7,1,1,-10),"mm"))


#######################################################################
#######################################################################
#######################################################################
#######################################################################
ggplot(plotval,aes(x=Relative,y=Distance))+
  geom_violin(trim=FALSE,aes(fill = factor(Relative)))+
  geom_boxplot(width= 0.5, trim=FALSE,aes(fill = factor(Compartment)))+ 
  labs(title="Middle-3'",x="", y = "Distances (nm)")+
  # facet_grid(. ~ Compartment, scales="free") +
  scale_y_continuous(limits =c(0,300))+#, breaks =seq(from = 0, to = 300, by = 50), minor_breaks =seq(from = 25, to = 325, by = 50)) +#,limits =c(-100,100))  +
  coord_flip()+
  scale_x_discrete(limits=orders)+#,
  geom_text(data = dia_me,size=6, aes(label = Value, y=270, x=lens:1.25))+
  
  #labels=c("Cyclo100"="Cyclo","Untreated"="Untreated","Coloc"="Coloc")) +
  # stat_summary(fun.data = n_fun, geom = "text",hjust=Inf)+
  # geom_text(data = data.frame(), aes(x = names(meds) , y = meds, label = paste("n =", meds))) +
  #scale_x_continuous(breaks =c(-100,0,100), minor_breaks =c(-50,50),limits = c(-100,100)) +
  theme(axis.line = element_line(colour = "black"),
        title = element_text(size=18),
        #plot.title = element_text(hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #panel.grid.major = element_line(color ="black",size=2),
        #panel.grid.minor = element_line(color ="black"),
        # panel.border = element_blank(),
        panel.border = element_rect(colour = "black",fill=NA),
        panel.background = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(),
        strip.text.y = element_blank(),
        legend.position="none",
        axis.title.x = element_text(face="bold", size=18),
        # axis.text.x  = element_text(face="bold", size=15),
        axis.text.x  = element_text(face="bold", size=15),
        axis.title.y = element_text(face="bold", size=18,margin=margin(0,+10,0,5)),
        axis.text.y  = element_text(face="bold", size=15))
#plot.margin=unit(c(-7,1,1,-10),"mm"))
