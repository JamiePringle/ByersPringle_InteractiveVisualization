#this code defines function to calculate the average location of a series of
#points on a sphere. It also calculates the square root of the sum of the square
#distances, which is equivalent to the standard deviation of their distances.

library('geosphere')
library('testit')

sumSquareDist <-function(lonLatPnt,lonVec,latVec) {
  #this function returns the sum of distances squared between lonLatPnt 
  #(two element vector, c(longitude,latitude) in degrees)
  #degrees) and all the points in the vectors lonVec and latVec. Units are km^2
  
  lonPnt=lonLatPnt[1]
  latPnt=lonLatPnt[2]
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  lonLatMat=matrix(0,nrow=length(lonVec),ncol=2)
  lonLatMat[,1]=lonVec
  lonLatMat[,2]=latVec
  
  dist=distHaversine(lonLatMat,c(lonPnt,latPnt),r=6378.137)
  
  sumSqDist=sum(dist^2)
  
  #for debugging
  #print(c('Point is',c(lonPnt,latPnt),sumSqDist))
  
  return(sumSqDist)
}

sumSquareDist_weighted <-function(lonLatPnt,lonVec,latVec,numVec) {
  #same as sumSquareDist() but with each (lon,lat) pair weighted by
  #numVec, which must be the same length. 
  
  lonPnt=lonLatPnt[1]
  latPnt=lonLatPnt[2]
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  lonLatMat=matrix(0,nrow=length(lonVec),ncol=2)
  lonLatMat[,1]=lonVec
  lonLatMat[,2]=latVec
  
  dist=distHaversine(lonLatMat,c(lonPnt,latPnt),r=6378.137)
  
  sumSqDist=sum((dist^2)*numVec) #does this work?
  
  #for debugging
  #print(c('Point is',c(lonPnt,latPnt),sumSqDist))
  
  return(sumSqDist)
}

findMeanAndSTD<-function(lonVec,latVec){
  #this function takes two vectors of equal size, the first of longitudes,
  #the second of latitudes, and returns a vector of the point that minimizes
  #the sum of square distances to all the points (the "mean" location) and 
  #the sqrt(sum of square distances) to all the points from this minimum. The
  #later is equivalent to a standard deviation.
  #
  #It does this by minimizing the function sumSquareDist as a function of location of 
  #the mean point to find the mean point
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  
  #find the lon and lat to start the optimization at by taking mean 
  #of lat and lon vectors
  
  #latStart=mean(latVec['lat'], na.rm = TRUE)  CJMP lets talk why this was a good start, but not the right start
  #lonStart=mean(lonVec['lon'], na.rm = TRUE)  CJMP lets talk why this was a good start, but not the right start
  latStart=mean(latVec, na.rm = TRUE)   
  lonStart=mean(lonVec, na.rm = TRUE)   
  #print(latStart)
  #print(lonStart)
  #for syntax of optim, refer to 
  #https://stackoverflow.com/questions/24623488/how-do-i-use-a-function-with-parameters-in-optim-in-r
  result=optim(par=c(lonStart,latStart),fn=sumSquareDist,lonVec=lonVec,latVec=latVec, lower = c(-180,-90), upper = c(180,90), method = "L-BFGS-B")
  
  #CJMP altered to return standard deviation, not sum of squares, by normalizing sum of squares by N-1
  N=length(latVec)
  #print(paste('RETURNING result$value,N',result$value,N))
  theResult=data.frame(lon=result$par[1],lat=result$par[2],std=sqrt(result$value/(N-1))) #CJMP changed how answer is returned
  return(theResult)
  
}

findMeanAndSTD_weighted<-function(lonVec,latVec,numVec){
  #this is the same as findMeanAndSTD() except the points are weighted by numVec
  #this accounts for the case where there are different numbers of particles at
  #each (lon,lat) pair, as is the case in our connectivity data
  
  assert('length of lonVec and latVec should be equal',length(lonVec)==length(latVec))
  
  #find the lon and lat to start the optimization at by taking mean 
  #of lat and lon vectors
  
  #latStart=mean(latVec['lat'], na.rm = TRUE)  CJMP lets talk why this was a good start, but not the right start
  #lonStart=mean(lonVec['lon'], na.rm = TRUE)  CJMP lets talk why this was a good start, but not the right start
  latStart=mean(latVec, na.rm = TRUE)   
  lonStart=mean(lonVec, na.rm = TRUE)   
  #print(latStart)
  #print(lonStart)
  #for syntax of optim, refer to 
  #https://stackoverflow.com/questions/24623488/how-do-i-use-a-function-with-parameters-in-optim-in-r
  result=optim(par=c(lonStart,latStart),fn=sumSquareDist_weighted,lonVec=lonVec,latVec=latVec,numVec=numVec,
               lower = c(-180,-90), upper = c(180,90), method = "L-BFGS-B")
  
  #CJMP altered to return standard deviation, not sum of squares, by normalizing sum of squares by N-1
  N=sum(numVec)
  #print(paste('RETURNING result$value,N',result$value,N))
  theResult=data.frame(lon=result$par[1],lat=result$par[2],std=sqrt(result$value/(N-1))) #CJMP changed how answer is returned
  return(theResult)
  
}

if (FALSE) {
  #write some testing code for sumSquareDist
  latVec<-c(-1.0,1.0,0.0)
  lonVec<-c(10.0,10.0,11.0)
  latPnt<-0.0
  lonPnt<-10.0
  sumSqDist<-sumSquareDist(c(lonPnt,latPnt),lonVec,latVec)
  print(paste('The sumSquareDist in meters is',sumSqDist,'km^2. Aswers should be 3*(111^2)=36963ish'))
  print(' ')
  
  #now test findMeanAndStd for three points in a line
  latVec<-c(-1.0,1.0,0.0)
  lonVec<-c(10.0,10.0,10.0)
  result=findMeanAndSTD(lonVec,latVec)
  print(paste(result$lon,result$lat,result$std,'answer should be 10.0E,0.0N, sqrt((111^2+111^2))/2=111ish'))
  print(' ')
  
  print('now test weighted code')
  #write some testing code for sumSquareDist
  latVec<-c(-1.0,1.0,0.0)
  lonVec<-c(10.0,10.0,11.0)
  numVec<-c(1.0,1.0,1.0)
  latPnt<-0.0
  lonPnt<-10.0
  sumSqDist<-sumSquareDist_weighted(c(lonPnt,latPnt),lonVec,latVec,numVec)
  print(paste('The sumSquareDist in meters is',sumSqDist,'km^2. Aswers should be 3*(111^2)=36963ish'))
  print(' ')
  
  #now test findMeanAndStd for three points in a line
  latVec<-c(-1.0,1.0,0.0)
  lonVec<-c(10.0,10.0,10.0)
  numVec<-c(1.0,1.0,1.0)
  result=findMeanAndSTD_weighted(lonVec,latVec,numVec)
  print(paste(result$lon,result$lat,result$std,'answer should be 10.0E,0.0N, sqrt((111^2+0.0+111^2)/2)=111ish'))
  print(' ')
  
  latVec<-c(-1.0,1.0,0.0)
  lonVec<-c(10.0,10.0,10.0)
  numVec<-c(100.0,100.0,100.0)
  result=findMeanAndSTD_weighted(lonVec,latVec,numVec)
  print(paste(result$lon,result$lat,result$std,'answer should be 10.0E,0.0N, sqrt((100*111^2+100*0.0+100*111^2)/299)=91ish'))
  print(' ')
  
}


