#set up libraries
source('connectivityUtilities.R')
library(tictoc)
source('averagePlaceOnSphere.R') #code to average on a sphere and calculate Ldiff
future::plan(multisession) #set up parallel calculation


#And get data for several months and several regions, trim points close to land, and combine

 for (minPLD in seq(2,30,2)) {
   for (season in c(2)){
     for (depth in c(1)) {
            
      #data frame to accumulate data in
      Eall<-data.frame()
      
      for (regionName in c('theAmericas','AsiaPacific','EuropeAfricaMiddleEast')) { 
      #for (regionName in c('theAmericas')) {  #FOR DEBUGGING ONLY
        print(paste('working on',regionName))
        
        #define depth and year/climatology and vertical behavior
        depthName=paste('_',as.character(depth),'m',sep='')
        year<-'climatology'
        
        if (depth==1) {
          verticalBehavior<-'starts'
          verticalBehavior<-'fixed' #choose which one to do!
        }else {
          verticalBehavior<-'fixed'
        }
        #minPLD<-14; 
        PLDname=paste('PLD',as.character(minPLD),sep=''); maxPLD<-minPLD
        
        #get data
        tic('getting data')
        if (season==2) {
          timeName='AprMayJun_'
          month<-4; E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-5; E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-6; E3<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
        } else if (season==4) {
          timeName='OctNovDec_'
          month<-10; E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-11; E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-12; E3<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
        } else if (season==1) {
          timeName='JanFebMar_'
          month<-1; E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-2; E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-3; E3<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
        }else if (season==3) {
          timeName='JulAugSep_'
          month<-7; E1<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-8; E2<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
          month<-9; E3<-getConnectivityData(regionName,depth,year,verticalBehavior,month,minPLD,maxPLD)
        }
        toc()
        
        #the distance from the nearest solid ground depends on your depth
        #so lets use the appropriate one
        if (depth==1) {
          thisGridDist=gridDist
        } else if (depth==10) {
          thisGridDist=gridDist_10m
        } else if (depth==20) {
          thisGridDist=gridDist_20m
        } else if (depth==40) {
          thisGridDist=gridDist_40m
        }
        
        #trim to land
        tic(paste('have trimmed',regionName))
        trimTo=TRUE
        E1<-subsetConnectivity_byGridValue(E1,thisGridDist,0.0,2.1,trimTo=trimTo)
        E2<-subsetConnectivity_byGridValue(E2,thisGridDist,0.0,2.1,trimTo=trimTo)
        E3<-subsetConnectivity_byGridValue(E3,thisGridDist,0.0,2.1,trimTo=trimTo)
        toc() 
        
        tic(paste('have combined',regionName))
        E<-combineConnectivityData(E1,E2)
        E<-combineConnectivityData(E,E3)
        toc()
        
        if (nrow(Eall)==0) {
          Eall<-E
        }
        else {
          tic(paste('have combined',regionName,'with other regions'))
          Eall<-combineConnectivityData(Eall,E)
          toc()
        }
      }
      
      #And now lets compute Ladv and Ldiff, and save a subset with Ladv, Ldiff, and mean ending points of data.
      #rename as Eall as E
      E<-Eall
      
      #add lat lon
      E<-addLatLon(E)
      
      #calculate fraction returned and add to data
      whatLength<-function(a){
        return(sum(a$numTo))
        print(a$numTo)
      }
      numReturn<-apply(E,1,whatLength)
      E$fracReturn<-numReturn/E$numLaunched
      
      #now calculate Ladv and Ldiff
      Esub<-E#[1:10000,] #subset for testing
      tic('calculating Ladv, Ldiff, calc took')
      
      #the following code iterates over each row of Esub (usually Esub<-E), and calculates
      #Ladv and Ldiff and returns them as a data.frame(). If you use future_pmap_dfr(), the 
      #calculation is done in parallel on all cores. The details of how this works can be 
      #found at https://blog.az.sg/posts/map-and-walk/
      Edists<-future_pmap_dfr(Esub,function(...){
        Erow<-tibble(...)
        
        if (length(Erow$lonFrom)>1) {
          lonToVec<-Erow$lonTo ; lonFrom<-Erow$lonFrom[[1]]
          latToVec<-Erow$latTo ; latFrom<-Erow$latFrom[[1]]
          numToVec<-Erow$numTo
          
          #print(paste('STARTING lonToVec',lonFrom,length(lonToVec)))
          
          #calculate mean location of *To locations
          #meanPlaceAndStd<-findMeanAndSTD(lonToVec,latToVec) #old, incorrect form, that did not do weighting.
          meanPlaceAndStd<-findMeanAndSTD_weighted(lonToVec,latToVec,numToVec)
          latMean<-meanPlaceAndStd$lat
          lonMean<-meanPlaceAndStd$lon
          Ldiff<-meanPlaceAndStd$std
          
          #debugging code
          #print(paste('BARF',lonFrom,latFrom,lonMean,latMean,Ldiff))
          
          #now calculate mean distace from starting *From point to mean *To Point
          Ldist=distHaversine(c(lonFrom,latFrom),c(lonMean,latMean),r=6378.137)
          
          theResult<-data.frame(Ladv=Ldist,Ldiff=Ldiff,lonMeanTo=lonMean,latMeanTo=latMean)
        } else {
          theResult<-data.frame(Ladv=NaN,Ldiff=NaN)
        }
      })
      toc()
      
      #now, lets add this to the existing dataFrame
      #NOTE WELL: this will be very wrong if Esub was not E
      E<-cbind(Esub,Edists)
      
      #write full data to file, and save subset as CSV for plotting with python
      #saveRDS(E,paste(timeName,PLDname,'_',verticalBehavior,depthName,'_withLadvLdiff.RDS',sep=''))
      
      Esubset<-E[c('lonFrom','latFrom','Ladv','Ldiff','fracReturn','lonMeanTo','latMeanTo')]
      write.csv(Esubset,file=gzfile(paste('dataFiles/',timeName,PLDname,'_',verticalBehavior,depthName,'_withLadvLdiff.csv.gz',sep='')),
                row.names=FALSE)
      
    }
  }
}