#https://www.statology.org/k-means-clustering-in-r/

#Things to do
#determine optimum number of clusters based on wsse and silhouette width
#Cluster analysis of sh parameter using kmean, k mediod, and hierarchical method
#make inferences based on clusters(correlate the cluster result with indivual parameter value,
#total cumulative loss, elevation, change in housing price, other social demographics)

#load libraries
library(tidyverse)
library(Hmisc)
library(corrplot)
library(MESS) #to use diag.panel = panel.hist
#set working directory
setwd("C:/NRE5605/project")

#import sh_parameters
input <- read.csv("SH_Parameters.csv")
input <- input %>% 
  select(-census_tract)

#Creating the scatter plot matrix
pairs(input,
      upper.panel = NULL,         # Disabling the upper panel
      diag.panel = panel.hist,
      cex.labels = 1.2, font.labels = 1.2, main = "Pair plot of SH parameters")  # Adding the histograms

#correlation matrix that shows ranked correlation
cormat <- cor(input, method = "spearman") #make a correlation matrix
corrplot(cormat, method = "number", type = "lower", outline = T, order = "AOE") #make a corrplot
#function to calcualte wsse
wsse <- function(k) { #k is the number of clusters
  kmeans(input, k, nstart = 1000)$tot.withinss  # k-mean function, data is your input data
}

k.values <- 2:10 #gives a range for number of clusters
wsse.values <- map_dbl(k.values, wsse) #maps all the clusters numbers and gives corresponding wsse
plot(wsse.values, type = "o")
#From the wsse plot, it looks the elbow point is 3, so going for 3 no. of clusters

#function to calculate silhouette width

