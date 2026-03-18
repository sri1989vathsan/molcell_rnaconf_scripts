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

importfile3 <- function (file,filefold,Data,inputs) {
  for(i in 1:length(file))
  {
    temp <- read.csv(file[i],header=FALSE)
    temp2 <- cbind(temp,t(replicate(dim(temp)[1], filefold[i,1:2])),inputs)
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
                            Dist,Data[,9:10],Data$Compartment,rep(paste("5mid"),dim(Data)[1]))
  }
  else if(position=="mid3")
  {
    Dist <- sqrt((Data$X2)^2+(Data$Y2)^2)
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],Data$Compartment,rep(paste("mid3"),dim(Data)[1]))
  }
  else if(position=="53")
  {
    Dist <- sqrt((Data$X1-Data$X2)^2+(Data$Y1-Data$Y2)^2)
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],Data$Compartment,rep(paste("53"),dim(Data)[1]))
  }
  
  names(Disttable) <- c( "No","X1","Y1","Value","Gene", "Field","Compartment","Reference")
  
  return(Disttable)
}

Distances3D <- function (Data,position) {
  if(position=="5mid")
  {
    Dist <- Data$X2
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],rep(paste("5mid"),dim(Data)[1]))
  }
  else if(position=="mid3")
  {
    Dist <- Data$Y2
    Disttable <- data.frame(Data$No, rep(paste("100"),dim(Data)[1]),
                            rep(paste("100"),dim(Data)[1]), 
                            Dist,Data[,9:10],rep(paste("mid3"),dim(Data)[1]))
  }
  else if(position=="53")
  {
    Dist <- Data$Z2
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


cyt <- list.files(pattern = "Cytoplasmic Absolute Coordinates.csv", recursive = TRUE)
nuc <- list.files(pattern = "Nuclear Absolute Coordinates.csv", recursive = TRUE)

split <- function(inputss){
  return(t(matrix(unlist(strsplit(dirname(inputss),"/")),ncol=length(strsplit(dirname(inputss),"/")))))
}

cytfold <- split(cyt)
nucfold <- split(nuc)


Data1 <- data.frame()
Data2 <- data.frame()

Data1 <- importfile(cyt,cytfold,Data1,"Cytoplasm")
Data2 <- importfile(nuc,nucfold,Data2,"Nucleus")


names(Data1) <- c( "No","X1","Y1","Z1","X2","Y2","Z2","X3","Y3","Z3",
                   "Relative","Gene", "Field","Compartment")
names(Data2) <- c( "No","X1","Y1","Z1","X2","Y2","Z2","X3","Y3","Z3",
                   "Relative","Gene", "Field","Compartment")

Data <- rbind(Data1,Data2)


#######################################################################
#######################################################################
################## To plot different conditions #######################
#######################################################################
#######################################################################
# 
# L <- c("Oddeven","Coloc")
#L <- c("Oddeven2")
L <- c("AHNAK Untreated 5'-mid-3'","Nucleus")
L <- rbind(L,c("AHNAK PlaB 5'-mid-3'","Nucleus"))
#L <- rbind(L,c("Puro","Cytoplasm","GART"))
#L <- rbind(L,c("Puro","Cytoplasm","DYNC1I2"))
L <- data.frame(L)
#L <- t(L)
colnames(L) <- c("Gene","Compartment")
library(plyr)
# Data <- Datanuc
Data %>% match_df(L) -> plotval
##Datanuc -> plotval
detach("package:plyr", unload = TRUE)

################################################################################
################################################################################
################################################################################

temp1 <- data.frame(cbind(plotval$No,as.numeric(plotval$X1),as.numeric(plotval$Y1),
                          as.numeric(plotval$X2),as.numeric(plotval$Y2),
                          as.numeric(plotval$X3),as.numeric(plotval$Y3),
                          plotval$Gene,plotval$Compartment))

names(temp1) <-  c( "No","X1","Y1","X2","Y2","X3","Y3","Gene","Compartment")
temp1mean <- data.frame()
X <- (as.numeric(temp1$X1)+as.numeric(temp1$X2)+as.numeric(temp1$X3))/3
Y <- (as.numeric(temp1$Y1)+as.numeric(temp1$Y2)+as.numeric(temp1$Y3))/3

fivep <- data.frame(cbind(as.numeric(temp1$X2) - X,as.numeric(temp1$Y2) - Y,
                          1:length(X),
                          temp1$Gene,temp1$Compartment))
# fivep <- data.frame(cbind(as.numeric(temp1$Y1) - Y,1:length(Y),temp1$Gene,temp1$Compartment))

mid <- data.frame(cbind(as.numeric(temp1$X1) - X,as.numeric(temp1$Y1) - Y,
                        1:length(X),
                        temp1$Gene,temp1$Compartment))

threep <- data.frame(cbind(as.numeric(temp1$X3) - X,as.numeric(temp1$Y3) - Y,
                           1:length(X),
                           temp1$Gene,temp1$Compartment))

names(fivep) <- c("X","Y","Num","Condition","Compartment") -> names(mid)
names(threep) <- names(fivep)

Datasss <- data.frame(c("X"=c(),"Y"=c(),"Num"=c(),"Condition"=c(),"Compartment"=c()))

Datasss[fivep$Num,] <- fivep

plotvals <- gdata::interleave(fivep, mid, threep)

plotvals$X <- as.numeric(plotvals$X)*39.6825
plotvals$Y <- as.numeric(plotvals$Y)*39.6825

cols <- c(rep(paste("Red"),(sum(plotval$Gene==(unique(plotval$Gene)[1])))),
          rep(paste("Black"),(sum(plotval$Gene==(unique(plotval$Gene)[2])))))
cols1 <- c(rep(paste("Red"),(sum(plotval$Gene==(unique(plotval$Gene)[1])))))
cols2 <- c(rep(paste("Black"),(sum(plotval$Gene==(unique(plotval$Gene)[2])))))

################################################################################
################################################################################
################################################################################

# cols <- c(rep(paste("Black"),length()))

dev.new()
scatterPlot <- ggplot(plotvals,aes(X, Y,color=Num)) +
  geom_path(size = 0.5,lineend = "round")+
  # geom_line()+
  scale_color_manual(values=cols)+
  scale_y_continuous(limits =c(-200,200))+
  scale_x_continuous(limits =c(-200,200))+
  theme(axis.line = element_blank(),
        title = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #panel.grid.major = element_line(color ="black",size=2),
        #panel.grid.minor = element_line(color ="black"),
        #panel.border = element_blank(),
        # panel.border = element_rect(size=1),
        panel.background = element_blank(),
        axis.line.x = element_line(size=1),
        axis.line.y = element_line(size=1),
        strip.text.y = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        # axis.text.x  = element_text(face="bold", size=15),
        axis.text.x  = element_text(face="bold", size=15),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(face="bold", size=15))
scatterPlot

plotvals1 <- plotvals[plotvals$Condition==unique(plotvals$Condition)[1],]
plotvals2 <- plotvals[plotvals$Condition==unique(plotvals$Condition)[2],]




scatterPlot1 <- ggplot(plotvals1,aes(X, Y,color=Num)) +
  geom_path(size = 0.5,lineend = "round")+
  # geom_line()+
  scale_color_manual(values=cols1)+
  scale_y_continuous(limits =c(-200,200))+
  scale_x_continuous(limits =c(-200,200))+
  theme(axis.line = element_blank(),
        title = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #panel.grid.major = element_line(color ="black",size=2),
        #panel.grid.minor = element_line(color ="black"),
        #panel.border = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.line.x = element_line(size=1),
        axis.line.y = element_line(size=1),
        strip.text.y = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        # axis.text.x  = element_text(face="bold", size=15),
        axis.text.x  = element_text(face="bold", size=15),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(face="bold", size=15))
scatterPlot1

dev.new()
scatterPlot2 <- ggplot(plotvals2,aes(X, Y,color=Num)) +
  geom_path(size = 0.5,lineend = "round")+
  # geom_line()+
  scale_color_manual(values=cols2)+
  scale_y_continuous(limits =c(-200,200))+
  scale_x_continuous(limits =c(-200,200))+
  theme(axis.line = element_blank(),
        title = element_blank(),
        plot.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        #panel.grid.major = element_line(color ="black",size=2),
        #panel.grid.minor = element_line(color ="black"),
        #panel.border = element_blank(),
        # panel.border = element_rect(size=1),
        panel.background = element_blank(),
        axis.line.x = element_line(size=1),
        axis.line.y = element_line(size=1),
        strip.text.y = element_blank(),
        legend.position="none",
        axis.title.x = element_blank(),
        # axis.text.x  = element_text(face="bold", size=15),
        axis.text.x  = element_text(face="bold", size=15),
        axis.title.y = element_blank(),
        axis.text.y  = element_text(face="bold", size=15))
scatterPlot2
################################################################################
################################################################################
################################################################################


rog1 <- sqrt((as.numeric(fivep[fivep$Condition==unique(fivep$Condition)[1],1])^2+
             as.numeric(fivep[fivep$Condition==unique(fivep$Condition)[1],2])^2+
             as.numeric(mid[mid$Condition==unique(mid$Condition)[1],1])^2+
             as.numeric(mid[mid$Condition==unique(mid$Condition)[1],2])^2+
             as.numeric(threep[threep$Condition==unique(threep$Condition)[1],1])^2+
             as.numeric(threep[threep$Condition==unique(threep$Condition)[1],2])^2)/3)*39.6825

rog2 <- sqrt((as.numeric(fivep[fivep$Condition==unique(fivep$Condition)[2],1])^2+
                as.numeric(fivep[fivep$Condition==unique(fivep$Condition)[2],2])^2+
                as.numeric(mid[mid$Condition==unique(mid$Condition)[2],1])^2+
                as.numeric(mid[mid$Condition==unique(mid$Condition)[2],2])^2+
                as.numeric(threep[threep$Condition==unique(threep$Condition)[2],1])^2+
                as.numeric(threep[threep$Condition==unique(threep$Condition)[2],2])^2)/3)*39.6825

# 
# lens <- length(unique(Data$Reference))
# 
# ranges <- c(lens+0.45)
# 
# for (p in 1:(lens-1))
# {
#   q <- lens-p
#   temp <- c(q+0.6,q+0.4)
#   ranges <- c(ranges,temp)
#   print(ranges)
# }
# 
# orders <- unique(Data$Reference)
# dia_me <- plotval %>% group_by(Gene,Reference) %>% 
#   summarise(Value=median(Value))
# dia_me$Value <- round(dia_me$Value,2)
# 
# ord <- sapply(orders, function(value) which(dia_me$Reference == value))
# 
# dia_me <- dia_me[unlist(ord),]
# orders <- rev(orders)
# 
# 
# dev.new()
# ggplot(plotval,aes(x=Reference,y=Value, fill = Gene,
#                    colour = Gene))+
#   geom_flat_violin(position = position_nudge(x = .25, y = 0), 
#                    alpha = 0.4, trim = FALSE,colour = "BLACK")+
#   geom_point(position = position_jitterdodge(
#     jitter.width = 0.15,dodge.width = 0.35), size = .6,
#     alpha=0.7)+
#   geom_boxplot(aes(x = Reference, y = Value, fill = Gene),outlier.shape = NA, 
#                position = position_dodge(width= 0.35), 
#                alpha = 0.3, width = .1, colour = "BLACK")+
#   #geom_text(aes(colour=Gene, label = Value), show.legend = FALSE)+
#   ylab('Distance (nm)')+
#   xlab('Probe')+
#   coord_flip()+
#   theme_cowplot()+
#   ylim(0,500)+
#   scale_x_discrete(limits=orders)+
#   geom_text(data = dia_me, size=5, aes(label = Value, colour=Gene, 
#                                        group = Reference,
#                                        y=270, x=ranges), show.legend = FALSE)+
#   #guides(fill = FALSE, colour = FALSE) +
#   scale_colour_brewer(palette = "Dark2")+
#   scale_fill_brewer(palette = "Dark2")+
#   theme(legend.position = "bottom")
# 
# 
# 
# ################################################################################
# ################################################################################
# 
# 
# lens <- length(unique(Data$Compartment))+0.45
# #orders <- c("Coloc")
# orders <- unique(Data$Compartment)
# dia_me <- ddply(plotval, .(Compartment),numcolwise(median))
# dia_me$Value <- round(dia_me$Value,2)
# 
# dia_me <- dia_me[match(orders, dia_me$Compartment),]
# orders <- rev(orders)
# 
# names <- names(plotval)
# 
# names[4] <- "y_axis"
# names[5] <- "x_axis"
# 
# plotvals <- plotval
# names(plotvals) <- names
# 
# dev.new()
# 
# ggplot(plotval,aes(x=Compartment,y=Value, fill = Compartment,
#                    colour = Compartment))+
#   geom_flat_violin(position = position_nudge(x = .25, y = 0),adjust =2,
#                    trim = FALSE,colour = "BLACK")+
#   geom_point(position = position_jitter(width = .15), size = .35)+
#   geom_boxplot(aes(x = Compartment, y = Value),outlier.shape = NA, 
#                #position = position_nudge(x = 0.15, y = 0), 
#                alpha = 0.3, width = .1, colour = "BLACK")+
#   ylab('Distance (nm)')+
#   xlab('Probe')+
#   coord_flip()+
#   theme_cowplot()+
#   scale_x_discrete(limits=orders)+
#   geom_text(data = dia_me,size=6, aes(label = Value, y=270, x=lens:1.25))+
#   guides(fill = FALSE, colour = FALSE) +
#   scale_colour_brewer(palette = "Dark2")+
#   scale_fill_brewer(palette = "Dark2")
# #ggtitle("Figure 6: Change in Colour Palette")
# 
# ################################################################################
# ################################################################################
# ################################################################################
# 
# dev.new()
# ggplot(plotval,aes(x=Compartment,y=Value))+
#   geom_violin(trim=FALSE,aes(fill = factor(Compartment)))+#,scale="width")+
#   geom_boxplot(width=0.05, fill="white",outlier.colour=NA)+
#   labs(title="",x="", y = "Distances (nm)")+
#   # facet_grid(. ~ Compartment, scales="free") +
#   scale_y_continuous(limits =c(0,300))+#, breaks =seq(from = 0, to = 300, by = 50), minor_breaks =seq(from = 25, to = 325, by = 50)) +#,limits =c(-100,100))  +
#   coord_flip()+
#   scale_x_discrete(limits=orders)+#,
#   geom_text(data = dia_me,size=6, aes(label = Value, y=270, x=lens:1.25))+
#   #labels=c("Cyclo100"="Cyclo","Untreated"="Untreated","Coloc"="Coloc")) +
#   # stat_summary(fun.data = n_fun, geom = "text",hjust=Inf)+
#   # geom_text(data = data.frame(), aes(x = names(meds) , y = meds, label = paste("n =", meds))) +
#   #scale_x_continuous(breaks =c(-100,0,100), minor_breaks =c(-50,50),limits = c(-100,100)) +
#   theme(axis.line = element_line(colour = "black"),
#         title = element_text(size=18),
#         plot.title = element_text(hjust = 0.5),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         #panel.grid.major = element_line(color ="black",size=2),
#         #panel.grid.minor = element_line(color ="black"),
#         # panel.border = element_blank(),
#         panel.border = element_rect(colour = "black",fill=NA),
#         panel.background = element_blank(),
#         axis.line.x = element_line(),
#         axis.line.y = element_line(),
#         strip.text.y = element_blank(),
#         legend.position="none",
#         axis.title.x = element_text(face="bold", size=18),
#         # axis.text.x  = element_text(face="bold", size=15),
#         axis.text.x  = element_text(face="bold", size=15),
#         axis.title.y = element_text(face="bold", size=18,margin=margin(0,+10,0,5)),
#         axis.text.y  = element_text(face="bold", size=15),
#         plot.margin=unit(c(-7,1,1,-10),"mm"))
