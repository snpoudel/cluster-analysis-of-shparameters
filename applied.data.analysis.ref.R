###the first part to get your started on processing the NetCDF data
rm(list=ls(all=T))
#install.packages('ncdf4') # most likely you will need to install this.
library('ncdf4')
setwd('/Users/rfeng/Documents/UConn/2021Spring/homework/HW3')
ncFile<-nc_open('ts_Amon_GISS-E2-1-G_historical_r1i1p1f1_gn_195101-200012.nc')
LonIdx <- which( ncFile$dim$lon$vals > 230 & ncFile$dim$lon$vals < 300)  # find indexes a range of longitudes
LatIdx <- which( ncFile$dim$lat$vals > 20. & ncFile$dim$lat$vals < 55.) # find indexes a range of latitudes
ts <- ncvar_get( ncFile, 'ts')[LonIdx, LatIdx,] # using the indexes to select ts
ts_ann=apply(ts,c(1,2), mean) # calculate temporal averages of ts - the third dimension
lat = ncvar_get( ncFile, 'lat')[LatIdx] #find lat corresponding to latIdx
lon = ncvar_get( ncFile, 'lon')[LonIdx] #find lon corresponding to lonIdx
ts_ann = (ts_ann - mean(ts_ann))/sd(ts_ann) # standardize the variable, ts is on different scale than the difference of ts
#######now read the prediction of future ts
ncFile<-nc_open('ts_Amon_GISS-E2-1-G-CC_esm-ssp585_r1i1p1f1_gn_200501-205012.nc')
ts_early <- ncvar_get( ncFile, 'ts')[LonIdx, LatIdx,1:240]
ncFile<-nc_open('ts_Amon_GISS-E2-1-G-CC_esm-ssp585_r1i1p1f1_gn_205101-210012.nc')
ts_end <- ncvar_get( ncFile, 'ts')[LonIdx, LatIdx,361:600]
diff_ts = apply(ts_end,c(1,2), mean) - apply(ts_early,c(1,2), mean)
diff_ts = (diff_ts - mean(diff_ts))/sd(diff_ts)
nlength = dim(diff_ts)[1]*dim(diff_ts)[2]  # calculate the total number of lat/lon raster grid
#install.packages('pracma')
library('pracma')
ts_ann = Reshape(ts_ann, nlength,1) # convert the 2D array into 1D 
diff_ts = Reshape(diff_ts, nlength, 1)
ts_data = data.frame(ts_ann, diff_ts) # combine two variables to form a data.frame
##########let's make a WSSE, average silhouette width function to help select the numbers of clusters
set.seed(123)
#install.packages('cluster')
library(cluster)
wsse <- function(k) {
  kmeans(ts_data, k, nstart = 1000)$tot.withinss  # k-mean function
}
avg_sil <- function(kk) {
  hm <- agnes(ts_data, diss=F, method="ward") # construct the tree with ward's method
  labels=cutree(as.hclust(hm), k=kk)  # kk is the input number of clusters
  ss <- silhouette(labels, dist(ts_data)) 
  return(mean(ss[, 3]))
}
########now we will use these function to calculate a WSSE and avg silhouette width for a range of cluster numbers
k.values <- 2:15
#install.packages('purrr')
library('purrr')  # purrr package contains the map function to execute functions with a range of input values (k,values)
avg_sil_values <- map_dbl(k.values, avg_sil) # map_dbl execute the function avg_sil with inputs of k.values from 2 to 15
 <- map_dbl(k.values, wsse) # same as the above
#install.packages('mclust')
library(mclust) #mclust package contains the Gaussian Mixture model functions
###Use Mclust function to produce a whole range of solutions for cluster number 1 to 15
ts.mclust=Mclust(ts_data, G = 2:15)
###check out which cluster option has the biggest BIC value
maxBIC = apply(ts.mclust$BIC,1,max)
######from here down, you are in charge of coding, good luck! 
######now make the three panel plot with avg_sil_values, wsse_values, and maxBIC 
#make a decision on the optimal numbers of clusters
######complete question 2
## generate cluster labels based on three different methods
par(mfrow=c(1,3))
title( "Different method of determining the number of cluster")
plot(k.values,wsse_values,type='b', xlab="Number of Clusters (K)", ylab= "Total within sum of Square",cex.lab=1.5,col="blue", cex=1.5, cex.axis=1.5)
plot(k.values,avg_sil_values,type='b',xlab="Number of Clusters (K)", ylab= "Average Silhouette Width",cex.lab=1.5, col="blue",cex=1.5, cex.axis=1.5)
plot(k.values,maxBIC, type='b',xlab="Number of Clusters (K)", ylab= "Gap Statistics",cex.lab=1.5, col="blue",cex=1.5, cex.axis=1.5)



par(mfrow=c(1,3))

pm = pam(ts_data, 8, diss = F)
pm.labels = pm$cluster
plot(ts_ann, diff_ts, col=rainbow(10)[pm.labels],main="Partition Around Medoid (PAM) Method", xlab="Temperature (K)", ylab= "Change in Temperature",cex.lab=1.5, cex=1.5, cex.axis=1.5)
hm=agnes(ts_data, diss=F, method="ward")
hm.labels=cutree(as.hclust(hm), k=7)
plot(ts_ann, diff_ts, col=rainbow(10)[hm.labels],main= "Hierarchical method ", xlab="Temperature (K)", ylab= "Change in Temperature",cex.lab=1.5, cex=1.5, cex.axis=1.5)
ts.mclust=Mclust(ts_data, G = 9)
mm.labels = ts.mclust$classification
plot(ts_ann, diff_ts, col=rainbow(10)[mm.labels],main="Model-based Method ", xlab="Temperature (K)", ylab= "Change in Temperature",cex.lab=1.5, cex=1.5, cex.axis=1.5)
#title("Scatterplots from different methods", line =-2, outer = TRUE, cex=2, font=10)
mtext("Scatterplots from different methods", side =3, line = -1.4, outer = TRUE, font=2)

##restore the raw values - the following is to get you started to answer question 2
ts_ann=Reshape(apply(ts,c(1,2), mean), nlength, 1) 
diff_ts = apply(ts_end,c(1,2), mean) - apply(ts_early,c(1,2), mean)
diff_ts = Reshape(diff_ts, nlength, 1)
####calculate group average

par(mfrow=c(1,3))

#tapply() applies a function or operation on subset of the vector broken down by a given factor variable.
### calculate the averages for each cluster center
tann_avg=tapply(ts_ann, pm.labels, mean)-273.15
tdiff_avg=tapply(diff_ts, pm.labels, mean)
## make a data frame by combining coordinates for annual mean T and change of T
pm.cluster_center = data.frame(tann_avg, tdiff_avg)
## plot barplot: barplot function only allows for matrix as input; it also thinks each row is a variable instead of each column
## so we have to transpose the matrix
barplot(t(as.matrix(pm.cluster_center)),xlab="Number of Clusters",ylab = "Susceptibility",legend = c("Median of temperature", "Median of Temperature change"), beside=T, col=c("aquamarine3","coral"))

tann_avg=tapply(ts_ann, hm.labels, mean)-273.15
tdiff_avg=tapply(diff_ts, hm.labels, mean)
## make a data frame by combining coordinates for annual mean T and change of T
pm.cluster_center = data.frame(tann_avg, tdiff_avg)
## plot barplot: barplot function only allows for matrix as input; it also thinks each row is a variable instead of each column
## so we have to transpose the matrix
barplot(t(as.matrix(pm.cluster_center)),xlab="Number of Clusters",ylab = "Susceptibility",legend = c("Median of temperature", "Median of Temperature change"), beside=T, col=c("aquamarine3","coral"))

tann_avg=tapply(ts_ann, mm.labels, mean)-273.15
tdiff_avg=tapply(diff_ts, mm.labels, mean)
## make a data frame by combining coordinates for annual mean T and change of T
pm.cluster_center = data.frame(tann_avg, tdiff_avg)
## plot barplot: barplot function only allows for matrix as input; it also thinks each row is a variable instead of each column
## so we have to transpose the matrix
barplot(t(as.matrix(pm.cluster_center)),xlab="Number of Clusters",ylab = "Susceptibility",legend = c("Median of temperature", "Median of Temperature change"), beside=T, col=c("aquamarine3","coral"))
#head(pm.cluster_center)

