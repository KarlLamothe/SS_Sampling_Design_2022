################################################################################
# Selecting sample sites for 2022 Silver Shiner surveys in 16 Mile Creek
# R Code written by Karl A. Lamothe - karl.lamothe@dfo-mpo.gc.ca
# Great Lakes Laboratory for Fisheries and Aquatic Sciences
# Historical data were collected under Fisheries and Oceans Canada permitting 
# Edited 2022-09-16  R Version: 4.2.1
################################################################################
# load libraries
library(pacman)      # download/load packages 'p_load'
p_load(unmarked)     # occupancy models
p_load(ggplot2)      # plotting
p_load(patchwork)    # multiple ggplots per figure
p_load(ape)          # Moran.I

#~~~~~~~~~~~~~~~~~~~~~~~#
# personal ggplot theme #
#~~~~~~~~~~~~~~~~~~~~~~~#
theme_me <- theme_bw() +
  theme(axis.title   = element_text(size= 11, family= "sans", colour= "black"),
        axis.text.x  = element_text(size= 10, family= "sans", colour= "black"),
        axis.text.y  = element_text(size= 10, family= "sans", colour= "black"),
        legend.title = element_text(size= 10, family= "sans", colour= "black"),
        legend.text  = element_text(size= 8,  family= "sans", colour= "black"),
        strip.text   = element_text(size= 11, family= "sans", colour= "black"),
        plot.title   = element_text(size= 11, family= "sans", colour= "black"),
        panel.border = element_rect(colour = "black", fill=NA),
        axis.ticks   = element_line(colour="black"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# load previously collected data #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
Habitat      <- read.csv("Data/Habitat.csv", header = T)
Depth_2022   <- read.csv("Data/Depth_Site_Selection.csv", header=T)
SSdata.Adult <- read.csv("Data/Adult_Occurrence.csv", header=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Testing for spatial autocorrelation in habitat variables using Moran I
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
colnames(Depth_2022)

# plotting depth histogram
sampleddepths<-ggplot(Depth_2022, aes(x=Mean.pool.depth))+
  geom_histogram(aes(y=..density..), color='white', bins=25)+
  geom_density(col = "red", lwd=3)+
  ggtitle("All surveyed sites")+
  labs(y='Density', x='Mean pool depth (m)') + ylim(0,4)+
  theme_me 
sampleddepths

# prep data for spatial autocorrelation calculation
colnames(Depth_2022)
Depth_2022<-Depth_2022[1:8] # select relevant columns
Depth_2022_2<-Depth_2022[complete.cases(Depth_2022), ] # remove data without coords
head(Depth_2022_2)

# Calculate geographic distances and inverse distance matrix
SS_dists <- as.matrix(dist(cbind(Depth_2022_2$Long, Depth_2022_2$Lat)))
SS.dists.inv <- 1/SS_dists
diag(SS.dists.inv) <- 0
SS.dists.inv[1:5, 1:5]

# Test for spatial autocorrelation
Moran.I(Depth_2022_2$Mean.pool.depth, SS.dists.inv, alternative="greater")

# crude map of sampling locations as a function of depth
ggplot(Depth_2022_2, aes(x=-Long, y=Lat, color=Mean.pool.depth))+
  geom_point()+
  labs(x='Longitude',y='Latitude',color='Depth')+
  theme_me+
  coord_fixed()+
  theme(legend.position = c(0.85,0.65),
        legend.background = element_blank())

################################################################################
########################### Occupancy models ###################################
###############################################################################
# Making the response variable data frame for the adult occupancy model
SS.frame2<-cbind.data.frame(Haul.1 = SSdata.Adult$Haul.1, 
                            Haul.2 = SSdata.Adult$Haul.2, 
                            Haul.3 = SSdata.Adult$Haul.3)

# Occupancy model data frame using 'unmarked'
SS.umf2<-unmarkedFrameOccu(y=SS.frame2, siteCovs = Habitat)
SS.umf2<-SS.umf2[-c(15),]
head(SS.umf2)

# run model
fm3 <- occu(~ 1 ~Depth, SS.umf2) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Collect data and make prediction #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Using measured mean depth from 2022 initial sampling
newdata   <- cbind.data.frame(Depth = Depth_2022$Mean.pool.depth)

# prediction using the occupancy model
occu <- cbind.data.frame(predict(fm3, type='state', newdata=newdata),
                         Depth=newdata$Depth)

# plot predictions as a function of Depth
occ.plt2<-ggplot(data=occu)+
  geom_ribbon(aes(x=Depth,ymin=lower, ymax=upper), fill="grey75")+
  geom_line(aes(x=Depth, y=Predicted), size=0.5)+
  labs(x='Depth (m)', y='Occupancy Probability')+
  ylim(0,1)+ 
  scale_x_continuous(limits=c(0.25,1),
                     breaks=c(0.25,0.50,0.75,1.00),
                     labels=function(x){sprintf("%.2f", x)})+
  theme_me
occ.plt2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Extract predictions
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
OccupancyEstimate <- occu$Predicted # predicted occupancy for measured depths
Depth_2022$Psi    <- OccupancyEstimate # include it in the original df

# select survey sites
set.seed(517438)
sites <- cbind.data.frame(Pool.ID=sample(Depth_2022$Pool.ID, size=100, 
                                         replace=F, prob=OccupancyEstimate))
Selected <- merge(x=sites, y=Depth_2022)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# making plot of depth for selected sites and comparing with full distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# histogram of selected depths
selectedsites<-ggplot(Selected, aes(x=Mean.pool.depth))+
  geom_histogram(aes(y=..density..), color='white', bins=20)+
  geom_density(col = "red", lwd=3)+
  ggtitle('Selected sites')+
  labs(y='Density', x='Mean pool depth (m)') + ylim(0,4)+ 
  theme_me 

sampleddepths+selectedsites
# in this plot you see that deeper sites were preferentially selected, indicated
# by the flattening of the density curve. However, its not an outstanding 
# difference, which is a result of the relationship between 
# occupancy and depth in combination of the depths sampled (see: occ.plt2)

# histogram of selected hydraulic head
ggplot(Selected, aes(x=HH))+
  geom_histogram(color='white', bins=20)+
  ggtitle('Selected sites')+
  labs(y='Frequency', x='Hydraulic head (mm)') + #ylim(0,4)+ 
  theme_me 

nrow(Depth_2022[Depth_2022$Mean.pool.depth>0.8,])
nrow(Depth_2022[Depth_2022$P1..max.>0.8,])
max(Depth_2022$P1..max.)
