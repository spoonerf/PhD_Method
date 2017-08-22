install.packages("dismo")
library(dismo)
library(demoniche)
library(biomod2)

capra <-gbif("capra", "ibex", geo=FALSE)
capgeo <- subset(capra, !is.na(lon) & !is.na(lat) & (lon!="NA" & lat !="NA") | year!="NA") 
dups <- duplicated(capgeo[, c("lon", "lat")])
capg <-capgeo[!dups, ]
capg2 <- capg[capg$lon > 0 & capg$lon<25 & capg$lat > 43 & capg$year>=1985, ] 
capg2$presence<-rep(1)
capg2$ID<-1:nrow(capg2)
capc<-as.matrix(capg2[ , c( "lon","lat", "presence")])
capc<-matrix(capc[complete.cases(capc)], ncol=3)
xy<-as.matrix(capc[,c(1,2)])
df<-data.frame(capc[,3])
sp<-SpatialPointsDataFrame(coords=xy, data=df)
x <- circles(sp, d=50000, lonlat=TRUE)
pol <- polygons(x)
plot(pol)
points(sp)

samp1 <- spsample(pol, 250, type='random', iter=25)
points(samp1, col="red")

e<-extent(sp)+4 #adding a buffer to the extent, should be dispersal distance of species ideally
year<-1950
bio<-stack(paste(wd, "/Bioclim/",year,"_bioclim_variable_stack.tif", sep=""))  
mask<-bio[[1]]
set.seed(1963)

bg<-randomPoints(mask, 500)

bg2 <- randomPoints(mask, 50, ext=e)


par(mfrow=c(1,2))
plot(!is.na(mask), legend=FALSE)
points(bg, cex=0.5)

cells <- cellFromXY(mask, samp1)
length(cells)
cells <- unique(cells)
length(cells)

xy <- xyFromCell(mask, cells)
plot(pol, axes=TRUE)
points(xy, cex=0.75, pch=20, col='blue')
spxy <- SpatialPoints(xy, proj4string=CRS('+proj=longlat +datum=WGS84'))
o <- over(spxy, geometry(x))
xyInside <- xy[!is.na(o), ]

v <- extract(mask, x@polygons, cellnumbers=T)
v <- do.call(rbind, v)
v <- unique(v[,1])
head(v)

m <- mask
m[]<-NA
m[v]<-1
plot(m, ext=extent(x@polygons)+1)
plot(x@polygons, add=T)
points(sp)

lf<-list.files(paste(wd, "/Bioclim/", sep=""))

first<-which(grepl("1985", lf))
last<-which(grepl("2016", lf))
all_years<-stack(paste(wd, "/Bioclim/",lf[first:last], sep=""))

bios<-seq(1,nlayers(all_years), by=19)

for (i in 1:19){
  layers<-bios
  bio_layer<-mean(all_years[[layers]])
  writeRaster(bio_layer, paste(wd, "/Bioclim/Bio_",i,"_1985_2016_average.tif",sep=""))
  bios<-bios+1
  print(layers)
}

bio_layer_pred<-c(5,6,7,13,19) # 5 = max temp warmest month , 6 = min temp coldest month , 7 = temp annual range ,
                                    # 10 = mean temp of warmest quarter, 11 = mean temp coldest quarter
                                    #13 = precipitation of wettest month, 19 = precipitatin of coldest quarter

predictors<-stack(paste(wd, "/Bioclim/Bio_", bio_layer_pred,"_1985_2016_average.tif",sep="" ))

plot(predictors)
predictors_alps<-crop(predictors, e)

presvals <- extract(predictors_alps, sp)

set.seed(0)
backgr <- randomPoints(predictors, 500)
absvals <- extract(predictors, backgr)

pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
head(sdmdata)
pairs(sdmdata[,2:ncol(sdmdata)], cex=0.1, fig=TRUE)



