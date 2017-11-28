#Self-exciting point process (SEPP) for crime hotspot prediction
#----------------------------------
#By: monsuru
#----------------------------------
#In the installation: R must be version 3.1.1; CRAN Window must be London 1 (not http one)
#Note that the package used here "etasFLP"

#install.packages("etasFLP")
#install.packages("maptools")
#install.packages("ggplot2")
#install.packages("sp")

library("etasFLP")
library(maptools)
library(ggplot2)
library(sp)

#------------------------------------------
#Inputs data sets
#------------------------------------------

#import spatial grid units e.g. see script "Creating a spatial grid system.....R"
#-----------------------------------
grids <- readOGR(dsn="//ds.leeds.ac.uk/staff/staff7/geomad/WORK/etasFLP_SEPP", layer="spatial_grid_system")
#-----------------------------------

#to check the CRS of the grid system is in WGS84, otherwise transform
#proj4string(grids)

#import historical crime data sets (as in the format below)
#Crime record description:
#Burglary crime record of South-Chicago area
#Data range: 01 Mar. 2011 - 6 Jan. 2012

#-----------------------------------
crime_R <- read.table(file="burglary_record.csv", sep=",", head=T)  
#-----------------------------------

#Data Preview:
#     time       lat      long        z   magn  
#1   	2011-03-01 41.83410 -87.65090   2  	4.5    
#2   	2011-03-01 41.76172 -87.62762   2  	4.5    
#3   	2011-03-01 41.84009 -87.61767   2  	4.5    
#4   	2011-03-01 41.78694 -87.61594   2  	4.5    
#5   	2011-03-01 41.77375 -87.60719   2  	4.5    
#6   	2011-03-01 41.77295 -87.60590   2  	4.5  
#.	    .		  .	      .    	  .	 .
#.	    .		  .	      .    	  .	 .	
#1253 2012-04-02 41.74371 -87.620819  2  	4.5  

#comments: The 'z' and 'magn1' 
#Ensure that the time format is as shown above.

#RUN THE PROGRAM:
####--------------------------------------------------------------####
	Hotspot_surface(s_grids=grids , 	 				
					train_endDate="2011-09-28", 		
						train_startDate="2011-03-01", 
							coverage_perc=10, 	
								plot=TRUE)		
####--------------------------------------------------------------####											####


#-----------------------------------------
# SEPP function
#-----------------------------------------

Hotspot_surface <- function(s_grids = grids , 
					train_endDate="2011-09-28", 
						train_startDate="2011-03-01", 
							coverage_perc=20, 
								plot=NULL){ ##

endDate <- as.Date(train_endDate) #the current date
startDate <- as.Date(train_startDate) #the start date of training dataset

#-----------------------------------------
#Note: the longer the training dataset, the more accurate the prediction. 
#This however, comes with a price of enormous computational time, especially when the training data contains too many records.
#Based on experience, the algorithm will take approximately 2 minutes to run a training dataset of 2,000 records 
#Thus, I restrict the training data length to 2,000 records in this script (see line )
#Plot=NULL (TRUE/FALSE) #to plot the hotspot surface
#-----------------------------------------

crime_R <- crime_R[which(as.Date(crime_R[,1])<= as.Date(endDate)),]


copy_s_grids <- NULL 

#-------------------------------------------------------

final_D <- NULL

#extract the training dataset
train_subset_1 <- which(as.Date(crime_R[,1])<(startDate))
train_subset_2 <- which(as.Date(crime_R[,1])<=(endDate))

train_subset_diff <- crime_R[setdiff(train_subset_2, train_subset_1),]
train_subset_diff <- train_subset_diff[order(train_subset_diff[,1]),]

#to limit the number of records to 2,000 in order to reduce the computational time to approx. 2 minutes
how_roww <- nrow(train_subset_diff)
if(how_roww > 2000){
		train_subset_diff <- train_subset_diff[which(as.Date(train_subset_diff[,1])>=as.Date(train_subset_diff[(how_roww-2000),1])),]
	}
rownames(train_subset_diff)<- 1:nrow(train_subset_diff)

#creating a point location ID
SN <- row.names(train_subset_diff)


#-----------------------------------------------
#convert dates of training dataset to "seconds": this is required by the etasFLP algorithm
crime_R2<-NULL
for(i in 1:nrow(train_subset_diff)){#111 i=1
	conv_R <- as.numeric(as.POSIXct(train_subset_diff[i,1]))
	crime_R2 <- c(crime_R2,conv_R)
	}#111
#---------------------------------------
#now generate dataset to be used in the training etasFLP
final_D <- train_subset_diff[,2:5]
final_D <- cbind(crime_R2,final_D)
colnames(final_D)<- c("time","lat","long","z","magn1")

t1<-Sys.time() #capture training start time
#--------------------------------------
#etasFLP algorithm
result_N <- etasclass(cat.orig = final_D, magn.threshold = 3, magn.threshold.back = 4.5,
		k0 = 0.005, c = 0.005, p = 1.01, a = 1.05, gamma = 0.6, d = 1.1, q = 1.52, 
		params.ind = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,TRUE, TRUE), 
		declustering = TRUE, thinning = TRUE, flp = FALSE,
		ndeclust = 15, onlytime = FALSE, is.backconstant = FALSE,
		description = "etas flp", sectoday = TRUE, usenlm = TRUE,
		epsmax = 0.1)
#------------------------------------------------
#The risk values are computed at each unique point location; stored as "result_N$l"
#------------------------------------------------
t2<-Sys.time() #capture training end time
tt <- t2-t1
tt #time spent on training

#Now, append each risk values with their corresponding point location (record).
final_risk_Info <- cbind(train_subset_diff, result_N$l)

#sorting the locations in descending order of risk values
#final_risk_Info <- final_risk_Info[(order(-as.numeric(final_risk_Info[,6]))),]

#To generate the hotspot surface using the spatial grid system
#Convert the point records to spatial (.shp)
#Overlay the spatial grid system on the study area
#Sum up the risk values of all point location that fall within each grid unit. 

# Setting existing coordinate as lat-long system
cord.dec = SpatialPoints(cbind(final_risk_Info$long, final_risk_Info$lat), proj4string = CRS("+proj=longlat"))

# Transforming coordinate to UTM (which is the coordinate system of the grid system based 
# on the "Creating spatial unit system ....R" script": EPSG=32748 for WGS=84, UTM Zone=48M,# Southern Hemisphere)

point.UTM <- spTransform(cord.dec, (s_grids@proj4string))

#plot(s_grids)
#plot(point.UTM, add=TRUE)
#plot(s_grids[17,], col="red", add=TRUE)
#plot(point.UTM_risk[1187,],col="blue", add=TRUE)

#appending the attributes to points
point.UTM_risk = SpatialPointsDataFrame(point.UTM, final_risk_Info)

copy_s_grids <- s_grids
copy_s_grids$riskV <- NA

for(i in 1:length(copy_s_grids)){#1 i<-1
this_g_unit <- sum(point.UTM_risk@data[!is.na(point.UTM_risk %over% s_grids[i,]),][,6])
copy_s_grids$riskV[i] <- this_g_unit
}#1


#-----------------------------------
#plot the hotspot surface based on the risk values?
#-----------------------------------
if(plot==TRUE){#
stop_Plotting <- round(length(copy_s_grids@data$riskV)*(coverage_perc/100), digits=0)
plot(copy_s_grids)
j=1
	order_P <- order(-copy_s_grids@data$riskV)
		for(j in 1:stop_Plotting){#1
			plot(copy_s_grids[order_P[j],], col="red", add=TRUE)
			
		}#2
}

} ##









