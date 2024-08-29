#this script plots the results calculated by 01_Get_and_Trim_data.R
#on a regional map. This can be used to recreate the maps in the paper,
#and extend them to your applications. 

#first, set up the graphics environment and load map data
#a good guide to making these maps (on which this is based)
#can be found at https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
library("plotly")
library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")

#lets get a nice colormap
library(viridis)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

#now load one of the files for a given PLD, behavior, and season
fileInName<-'dataFiles/AprMayJun_PLD4_fixed_1m_withLadvLdiff.csv.gz'
dispersalData<-read.csv(fileInName)

#lets choose what region to plot
lonMin=-85.0
lonMax=-70.0
latMin=20.0
latMax=34.9

#now, by default, it will plot all the points in the world. We only
#want to plot the data near the Gulf of Mexico, so lets trim the data
#file
dispersalData<-dispersalData[dispersalData$latFrom>=latMin,]
dispersalData<-dispersalData[dispersalData$latFrom<=latMax,]
dispersalData<-dispersalData[dispersalData$lonFrom>=lonMin,]
dispersalData<-dispersalData[dispersalData$lonFrom<=lonMax,]

#lets limit the plot to the region around the Gulf of Mexico
#remember, we need to print to show the results
#plot fracReturn, Ladv and Ldiff in different plots
myMap<-ggplot(data = world) +
  geom_sf()+ggtitle(fileInName)+
  geom_point(data=dispersalData,mapping=aes(x=lonFrom,y=latFrom,color=fracReturn))+
  scale_color_viridis_c(option = "plasma")+
  coord_sf(xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), expand = FALSE)
jnk<-ggplotly(myMap)
show(jnk)


myMap<-ggplot(data = world) +
  geom_sf()+ggtitle(fileInName)+
  geom_point(data=dispersalData,mapping=aes(x=lonFrom,y=latFrom,color=Ladv))+
  scale_color_viridis_c(option = "plasma")+
  coord_sf(xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), expand = FALSE)
jnk<-ggplotly(myMap)
show(jnk)


myMap<-ggplot(data = world) +
  geom_sf()+ggtitle(fileInName)+
  geom_point(data=dispersalData,mapping=aes(x=lonFrom,y=latFrom,color=Ldiff))+
  scale_color_viridis_c(option = "plasma")+
  coord_sf(xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), expand = FALSE)
jnk<-ggplotly(myMap)
show(jnk)

#now, lets compute the physical adversity
dispersalData$lnRc<-dispersalData$Ladv^2/(2*dispersalData$Ldiff^2)-log(dispersalData$fracReturn)
myMap<-ggplot(data = world) +
  geom_sf()+ggtitle(fileInName)+
  geom_point(data=dispersalData,mapping=aes(x=lonFrom,y=latFrom,color=lnRc))+
  scale_color_viridis_c(option = "plasma")+
  coord_sf(xlim = c(lonMin, lonMax), ylim = c(latMin, latMax), expand = FALSE)
jnk<-ggplotly(myMap)
show(jnk)
