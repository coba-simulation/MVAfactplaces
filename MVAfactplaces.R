# ------------------------------------------------------------------------------
# Book: MVA - Applied Multivariate Statistical Analysis
# ------------------------------------------------------------------------------
# Quantlet: MVAfactplaces
# ------------------------------------------------------------------------------
# Description: Factor analysis of the Places data set from the package tourr.
# Shows how to do applied factor analysis and visualization. Including selection
# of the number of factors, estimation, rotation and visualization. Factor
# Scores are also evaluated and plotted on a map.
# ------------------------------------------------------------------------------
# Usage: -
# ------------------------------------------------------------------------------
# Inputs: - 
# ------------------------------------------------------------------------------
# Output: Factor Analysis Results, Scree Plots and different other visualization
# material of scores, communalities as well as a spatial representation of the
# Factor Scores.
# ------------------------------------------------------------------------------
# Keywords: Factor Analysis, Factor Viszalization, Places Rated, Rotation,
# Factor Scores
# ------------------------------------------------------------------------------
# See also: SMSfactanal, SMSfactuscrime, SMSfactushealth, MVAfactcarm,
# MVAfacthous, SMSfactbank,SMSfactfood, SMSfacthletic	
# ------------------------------------------------------------------------------
# Author: Johannes Stoiber, Michael Lebacher 01/07/2016
# ------------------------------------------------------------------------------

#################################################################################
############### Loading and opening the needed Packages #########################
#################################################################################

#install.packages("tourr")         ;library("tourr") 
#install.packages("stats")         ;library("stats")
#install.packages("foreign")       ;library("foreign")
#install.packages("psych")         ;library("psych")
#install.packages("ggplot2")       ;library("ggplot2")
#install.packages("reshape2"       ;library("reshape2")
#install.packages("scatterplot3d") ;library("scatterplot3d")
#install.packages("stargazer")     ;library("stargazer")
#install.packages("GPArotation")   ;library("GPArotation")
#install.packages("rworldmap")     ;library("rworldmap")

#################################################################################
########################### Data and Description ################################
#################################################################################

# The dataset is taken from the Places Rated Almanac, by Richard Boyer and David 
# Savageau, copyrighted and published by Rand McNally. The nine rating criteria 
# used by Places Rated Almanac are:
# Climate and Terrain
# Housing
# Health Care and Environment
# Crime
# Transportation
# Education
# Arts
# Recreation
# Economics
# For all but two of the above criteria, the higher the score, the better.
# For Housing and Crime, the lower the score the better.
# Additionally spatial data and information on population is included in the 
# data set.

#################################################################################
########################### Preparation #########################################
#################################################################################

data<-places[,-c(10,11,12,14)] # drop non-metric data
names.vec <- names(data) # save the names in a vector
names.vec[2] <- "housing" # renames housingcost in housing
names(data) <- names.vec # gives the names to the dataset

#################################################################################
########################### Suitabilit of the Dataset ###########################
#################################################################################

stargazer(KMO(data)$MSAi,type="text") # Drop varialbes with KMO<0.5
data <- data[,-9] # Drops econ
stargazer(KMO(data)$MSAi,type="text") # Now all variables satisfy: KMO>0.5

#################################################################################
########################### On the Number of Factors  ###########################
#################################################################################

scree(data,fa=F) # Elbow Criteria: 1 Factor, Eigenvalues: 2-3 Factors

#################################################################################
############################ Factor Analysis ####################################
#################################################################################

# Defines a function that compares the result for a specifict variable and
# a specifict factor across different estimation methods and rotations
fun<-function(x,y){                                # x gives the variable and y gives the factor
  Rot<-c("quartimax","varimax","oblimin","promax") # Possible rotations
  Est<-c("minres", "ml", "gls", "wls", "pa")       # Possible estimation methods
  Mat<-matrix(0,ncol=length(Rot),nrow=length(Est))
  rownames(Mat)<-Est; colnames(Mat)<-Rot
  for (i in 1:length(Est)) {
        for (j in 1:length(Rot))
      {
      Mat[i,j]<-fa(data,3,rotate=Rot[j],fm=Est[i])$loadings[x,y]       
      }
  }
    a<-rownames(fa(data,3)$loadings)[x]
    g<-function(x){max(x)-min(x)}
    spaest<-apply(Mat,2,g)
    sparot<-apply(Mat,1,g)
    Mat<-cbind(Mat,sparot)
    Mat<-rbind(Mat,c(spaest,0))
    rownames(Mat)[length(Est)+1]<-"spaest"
    return(stargazer(Mat,type="text",title=paste("Variable",a,"Factor",y)))
}
# function can be applied to finde suitable estimation and rotation

# Analysis with three Factors and different rotations
fact_none<-fa(data,3,rotate="none")           ;fact_none
fact_varimax<-fa(data,3,rotate="biquartimin") ;fact_varimax
fact_oblim<-fa(data,3,rotate="oblimin")       ;fact_oblim
fact_promax<-fa(data,3,rotate="promax")       ;fact_promax

# Pick one rotation for example variamax and check the loadings
fvs_varimax<-fact_varimax$Structure 

# Plot the Factor Loadings with Barplots
par(mfrow=c(1 ,3))
loadings1<-fact_varimax$loadings[,1]
loadings2<-fact_varimax$loadings[,2]
loadings3<-fact_varimax$loadings[,3]
col1<-1-as.numeric(is.na(fvs_varimax[,1]))
col2<-1-as.numeric(is.na(fvs_varimax[,2]))
col3<-1-as.numeric(is.na(fvs_varimax[,3]))
barplot(loadings1, horiz = T,main="Factor 1",col=col1*2)
barplot(loadings2, horiz = T,main="Factor 2",col=col2*2)
barplot(loadings3, horiz = T,main="Factor 3",col=col3*2)

# Plot the Factor Loadings in a directed graph scheme
par(mfrow=c(1 ,1))
factor.plot(fact_varimax, cut=0.5)
fa.diagram(fact_varimax,cut=0.5)

# Plot the Factor Loadings in a 3D-Plot
scores_varimax <- factor.scores(x=data, f=fact_varimax) # Extract the Scores 

f_loadings<-fact_varimax$loadings
order.loadings <- cbind(f_loadings[,1], f_loadings[,2], f_loadings[,3])
colnames(order.loadings) <- c("MR1", "MR2", "MR3")
vec.col <- c(1,2,4,1,1,3,4,1,4) # 1 loads on no factor, 4 on factor 1, 2 on factor 2, 3 on factor 3

s3d <- scatterplot3d(order.loadings,
              main="3D factor loadings",
              color=vec.col,
              lwd=2,
              pch=16,
              cex.lab = 1,
              cex.axis = 1,
              type="h", lty.hplot=2,
              xlab="MR 1", # corresponds to factor 1 (order.loadings[,1])
              ylab="MR 2", # corresponds to factor 2 (order.loadings[,2])
              zlab="MR 3"  # corresponds to factor 2 
              )
s3d.coords <- s3d$xyz.convert(order.loadings) # convert 3D coords to 2D projection
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
     labels=row.names(order.loadings),               # text to plot
     cex=0.8, pos=4)           # shrink text 50% and place to right of points)


# Communalities
fact_biquartimin<-fa(data,3,rotate="biquartimin");fact_biquartimin
factor_biq <- fact_biquartimin$loadings
h2 <- fact_biquartimin$communality
load_com <- cbind(factor_biq, h2)

# communality by hand:
g <- function(x){sum(x^2)}
communalities <- apply(factor_biq, 1, g)

# output for LaTeX
stargazer(load_com)

# sort communaliteis and show them in a plot, including total communality
com.sort <- sort(h2, decreasing = TRUE)
total <- sum(com.sort)/9
com.sort <- c(com.sort, total)
names(com.sort)[10] <- "total"
plot(com.sort,
     main ="Communalities",
     xlim = c(0, 11),
     ylim = c(0.1, 1),
     ylab = "explained variance",
     xaxt = "n",
     xlab = ''
)
abline(v=9.5, lty=2)
text(cbind(seq(1:10),com.sort),
     labels =  names(com.sort),
     cex=0.8, pos=1) 
stargazer(com.sort)


########################################################################################################
###########################        Doing everything again         ######################################
###########################            Use 2 Factors              ######################################
########################################################################################################

# Factor analysis with 2 factors and different rotations
fact_none<-fa(data,2,rotate="none");fact_none
fact_varimax<-fa(data,2,rotate="varimax");fact_varimax
fact_oblim<-fa(data,2,rotate="oblim");fact_oblim
fact_promax<-fa(data,2,rotate="promax");fact_promax


# We have a closer look at fact_varimax
fvs_varimax<-fact_varimax$Structure 

# Plot the Factor Loadings with Barplots
par(mfrow=c(1 ,2)) 
loadings1<-fact_varimax$loadings[,1]
loadings2<-fact_varimax$loadings[,2]
col1<-1-as.numeric(is.na(fvs_varimax[,1]))
col2<-1-as.numeric(is.na(fvs_varimax[,2]))
barplot(loadings1, horiz = T,main="Factor 1",col=col1*2)
barplot(loadings2, horiz = T,main="Factor 2",col=col2*2)

# Plot the Factor Loadings as a directed graph
par(mfrow=c(1 ,1))
fa.diagram(fact_varimax,cut=0.5)

# 2D plot visualizing factors and loadings
f_loadings<-fact_varimax$loadings
order.loadings <- cbind(f_loadings[,1], f_loadings[,2])
colnames(order.loadings) <- c("MR1", "MR2")
vec.col <- c(1,2,4,1,1,4,4,2) # 1 loads on no factor, 4 on factor 1, 2 on factor 2
plot(order.loadings,
     main = "2D factor loadings",
     col=vec.col,
     lwd = 2,
     pch = 16,
     cex.lab = 1,
     cex.axis = 1,
     xlab = "MR 1", # corresponds to order.loadings[,1]
     ylab = "MR2"   # corresponds to order.loadings[,2]
     )
abline(h=0.5, col=2) # corresponds to factor 2
abline(v=0.5, col=4) # corresponds to factor 1
text(f_loadings,
     labels=row.names(order.loadings),
     cex=0.8, pos=4) 
     )

# Communalities
fact_varimax<-fa(data,2,rotate="varimax");fact_varimax
h2.2 <- fact_varimax$communality
load_com.2 <- cbind(fact_varimax$loadings, h2.2)

# communality by hand:
com.2 <- apply(fact_varimax$loadings, 1, g)

# output for LaTeX
colnames(load_com.2)[3]<- "h2"
stargazer(load_com.2)
# sort communaliteis and show them in a plot, including total communality
com.sort <- sort(h2.2, decreasing = TRUE)
total <- sum(com.sort)/9
com.sort <- c(com.sort, total)
names(com.sort)[10] <- "total"
plot(com.sort,
     main ="Communalities",
     xlim = c(0, 11),
     ylim = c(0.1, 1),
     ylab = "explained variance",
     xaxt = "n",
     xlab = ''
)
abline(v=9.5, lty=2)
text(cbind(seq(1:10),com.sort),
     labels =  names(com.sort),
     cex=0.8, pos=1) 

########################################################################################################
########################### Spatial Visualiztion  #####################################################
########################################################################################################

scores_varimax<- factor.scores(x=data, f=fact_varimax)$scores # 2 factors
par(mfrow=c(1,1))
score_1<-scores_varimax[,1] # Select the first Factor Score
rbPal <- colorRampPalette(c('yellow','red')) # Coloring between yellow and red by intensity
cols <- rbPal(100)[as.numeric(cut(score_1,breaks = 100))]

newmap <- getMap(resolution = "low") # load in the world as a map
plot(newmap, xlim = c(-124, -65), ylim = c(27, 48.88), asp = 1,main="Geographic distribution of Factor Scores (I)") # Plot north america
coord<-cbind(places$long, places$lat) # load in longiture and latitude of the data
coord<-subset(coord, coord[,1]>-125)
points(coord, asp=T,bg=cols,pch=21,lwd=2) # gives the points

score_2<-scores_varimax[,2] # select the second Factor Score
rbPal <- colorRampPalette(c('yellow','red')) # Coloring between yellow and red by intensity
cols <- rbPal(6)[as.numeric(cut(score_2,breaks = 6))] # Coloring between yellow and red by intensity
plot(newmap, xlim = c(-124, -65), ylim = c(27, 48.88), asp = 1,main="Geographic distribution of Factor Scores (II)")# Plot north america
coord<-cbind(places$long, places$lat)# load in longiture and latitude of the data
coord<-subset(coord, coord[,1]>-125)
points(coord, asp=T,bg=cols,pch=21,lwd=2) # gives the points



