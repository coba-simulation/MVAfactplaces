#################################################################################
##                           Structure                                         ##
#################################################################################

#-------------------------------------------------------------------------------#
##  * 1   Preperation
##  * 1.1 Install and load packages
##  * 1.2 Load data and data description
##  * 1.3 Suitability of the dataset
##  * 1.4 Decide on number of factors
##  * 2   Factor Analysis (with 3 factors)
##  * 2.1 Methods and rotations
##  * 2.2 Visualization
##  * 3   Validation
##  * 3.1 Communalites
##  * 4   Repeat steps 2 and 3 using 1 factors
##  * 5   Spatial Visualization for 1 factors
#-------------------------------------------------------------------------------#

#################################################################################
##                     1   Preperation                                         ##  
##                     1.1 Installl and load packages                          ##
#################################################################################


# clear console
rm(list=ls())

# install and load packages
libraries = c("tourr", "stats","foreign", "psych", "ggplot2", "reshape2", 
              "scatterplot3d", "stargazer", "GPArotation", "rworldmap")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
}
)
lapply(libraries, library, quietly = TRUE, character.only = TRUE)



#################################################################################
##                  1.2 Load data and data description                         ##
#################################################################################

## load data 
head(places)

#-------------------------------------------------------------------------------#

# The dataset is taken from the Places Rated Almanac, by Richard Boyer and David 
# Savageau, copyrighted and published by Rand McNally. The nine rating criteria 
# used by Places Rated Almanac are:
#
# * Climate and Terrain
# * Housing
# * Health Care and Environment
# * Crime
# * Transportation
# * Education
# * Arts
# * Recreation
# * Economics
#
# For all but two of the above criteria, the higher the score, the better.
# For Housing and Crime, the lower the score the better.
# Additionally spatial data and information on population is included in the 
# data set.
#-------------------------------------------------------------------------------#



#################################################################################
##                         1.3 Suitabilit of the Dataset                       ##
#################################################################################

data        = places[,-c(10,11,12,14)]  # drop non-metric data
n.vec       = names(data)               # save the names in a vector
n.vec[2]    = "housing"                 # renames housingcost in housing (for later use)
names(data) = n.vec                     # gives the names to the dataset

stargazer(KMO(data)$MSAi,type="text")   # Drop varialbes with KMO<0.5
data = data[,-9]                        # Drops econ
stargazer(KMO(data)$MSAi,type="text")   # Now all variables satisfy: KMO>0.5

#################################################################################
##                   1.4 Decide on number of factors                           ##
#################################################################################

scree(data,fa=F) # Elbow Criteria: 1 Factor, Eigenvalues: 2-3 Factors

#################################################################################
##                           2   Factor Analysis                               ##
##                           2.1 Method and rotations                          ##
#################################################################################

# Analysis with three Factors and different rotations
fact_none = fa(data,3,rotate="none")        ;fact_none
fact_biqu = fa(data,3,rotate="biquartimin") ;fact_biqu
fact_obli = fa(data,3,rotate="oblimin")     ;fact_obli
fact_prom = fa(data,3,rotate="promax")      ;fact_prom

################################################################################
##                         2.2. Visualization                                 ##
################################################################################


# Pick one rotation (e.g. biquartimin) and check the loadings
fvs = fact_biqu$Structure 
fvs[abs(fvs)<0.5] = NA


#------------------------Barplot with laodings---------------------------------#

#  Plot the Factor Loadings with Barplots
par(mfrow=c(1 ,3))
loadings1 = fact_biqu$loadings[,1]
loadings2 = fact_biqu$loadings[,2]
loadings3 = fact_biqu$loadings[,3]
col1      = 1-as.numeric(is.na(fvs[,1]))
col2      = 1-as.numeric(is.na(fvs[,2]))
col3      = 1-as.numeric(is.na(fvs[,3]))
barplot(loadings1, horiz = T,main="Factor 1",col=col1*2)
barplot(loadings2, horiz = T,main="Factor 2",col=col2*2)
barplot(loadings3, horiz = T,main="Factor 3",col=col3*2)

#-----------------------Directed graph scheme-----------------------------------#

#  Plot the Factor Loadings in a directed graph scheme
par(mfrow=c(1 ,1))
factor.plot(fact_biqu, cut=0.5)
fa.diagram(fact_biqu,cut=0.5)

#------------------------3-D Plot of the 3 factors------------------------------#

# Plot the Factor Loadings in a 3D-Plot
f_loadings     = fact_biqu$loadings
order.loadings = cbind(f_loadings[,1], f_loadings[,2], f_loadings[,3])
colnames(order.loadings) = c("MR1", "MR2", "MR3")
vec.col = c(1,2,4,1,1,3,4,1,4) # 1 loads on no factor, 4 on factor 1, 2 on factor 2, 3 on factor 3

s3d  = scatterplot3d(order.loadings,
                     main     = "3D factor loadings",
                     color    = vec.col,
                     lwd      = 2,
                     pch      = 16,
                     cex.lab  = 1,
                     cex.axis = 1,
                     type     = "h",
                     lty.hplot= 2,
                     xlab     = "MR 1", # corresponds to factor 1 
                     ylab     = "MR 2", # corresponds to factor 2 
                     zlab     = "MR 3"  # corresponds to factor 3
)
s3d.coords = s3d$xyz.convert(order.loadings) # convert 3D coords to 2D projection
text(s3d.coords$x, s3d.coords$y,             # x and y coordinates
     labels=row.names(order.loadings),       # text to plot
     cex=0.8, pos=4)                         # shrink text 80% and place to right of points)

################################################################################
##                           3.  Validation                                   ##
##                           3.1 Communalities                                ##
################################################################################

# Communalities
h2       = fact_biqu$communality

# sort communalities in decreasing order and plot them, including total communality
com.sort = sort(h2, decreasing = TRUE)
total    = sum(com.sort)/9
com.sort = c(com.sort, total)
names(com.sort)[10] = "total"
plot(com.sort,
     main = "Communalities",
     xlim = c(0, 11),
     ylim = c(0.1, 1),
     ylab = "explained variance",
     xaxt = "n",
     xlab = ''
)
# vertical dashed line to seperate total communality value from the individual communalities
abline(v=9.5, lty=2)
text(cbind(seq(1:10),com.sort),
     labels =  names(com.sort),
     cex=0.8, pos=1) 


#################################################################################
##                    4   Repeat step 2 and 3 for 1 factors                    ##
#################################################################################

## Rotation and estimation
# Factor analysis with 1 factors using varimax roation
fact = fa(data,1);fact
# Information about factor loadings and communalities
tab = cbind(fact$Structure, fact$communality)
colnames(tab) = c("MR1", "h2")

fvs  = fact$Structure 
fvs[abs(fvs)<0.5] = NA

## Visualization
#------------------------Barplot with laodings---------------------------------#
# Plot the Factor Loadings with Barplots
loadings = fact$loadings[,1]
col      = 1-as.numeric(is.na(fvs[,1]))
barplot(loadings, horiz = T,main="Factor 1",col=col*2)

#-----------------------Directed graph scheme-----------------------------------#
# Plot the Factor Loadings as a directed graph
fa.diagram(fact,cut=0.5)


## Validation
# Communalities
h2.2       = fact$communality

# sort communaliteis and show them in a plot, including total communality
com.sort = sort(h2.2, decreasing = TRUE)
total    = sum(com.sort)/9
com.sort = c(com.sort, total)
names(com.sort)[10] = "total"
plot(com.sort,
     main = "Communalities",
     xlim = c(0, 11),
     ylim = c(-0.1, 1),
     ylab = "explained variance",
     xaxt = "n",
     xlab = ''
)
abline(v=9.5, lty=2)
text(cbind(seq(1:10),com.sort),
     labels =  names(com.sort),
     cex=0.8, pos=1) 

#################################################################################
##                5   Spatial visualiztion for 2 factors                       ##
#################################################################################

## load in the world as a map
newmap  = getMap(resolution = "low")                        

## extract the scores
score = factor.scores(x=data, f=fact)$scores 

#--------------------------Scores factor 1--------------------------------------#
plot(score)                                                  # most betweeen -1 and 2
score[score>2] = 2                                           # for illustration reasons
rbPal   = colorRampPalette(c('yellow','red'))                # Coloring between yellow and red by intensity
cols    = rbPal(10)[as.numeric(cut(score_1,breaks = 10))]
plot(newmap,
     xlim = c(-124, -65),
     ylim = c(27, 48.88),
     asp  = 1,
     main = "Geographic distribution of Factor Scores (I)"   # Plot north america
)
coord = cbind(places$long, places$lat)                       # load in longiture and latitude of the data
coord = subset(coord, coord[,1]>-125)
points(coord, asp=T,bg=cols,pch=21,lwd=2)                    # gives the points


