install.packages("dismo")
library(dismo)
library(demoniche)
library(biomod2)


wd<- getwd()
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

set.seed(0)
group<-kfold(sp, 4)


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

bio_layer_pred<-c(5,6,13,15,19) # 5 = max temp warmest month , 6 = min temp coldest month , 7 = temp annual range ,
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

colnames(sdmdata)<-c("pb", "bio5", "bio6", "bio13","bio15" ,"bio19")
pairs(sdmdata[,2:ncol(sdmdata)], cex=0.1, fig=TRUE)

m1 <- glm(pb ~ bio5+bio6+bio13 +bio15 + bio19, data=sdmdata)
summary(m1)

years<-1950:2016
predict_stack<-stack(paste(wd, "/Bioclim/",years,"_bioclim_variable_stack.tif", sep=""))
#predict_alps<-crop(predict_stack, e)
bio_layer_pred<-c(5,6,13,15,19)

#glm predictions for alpine ibex
for (i in 1:length(years)){
  predict_stack_year<-predict_stack[[bio_layer_pred]]
  predict_alps_year<-crop(predict_stack_year, e)
  names(predict_alps_year)<-c("bio5", "bio6", "bio13", "bio15", "bio19")
  pm1<-predict(predict_alps_year, m1)
  writeRaster(pm1, paste("D:/Fiona/Git_Method/Git_Method/Alp_SDMs/GLM/glm_",years[i],"capra_ibex.tif", sep=""), overwrite=T)
  bio_layer_pred<-bio_layer_pred+19
  print(i)
  }

predict_alps_year<-predict_alps[[bio_layer_pred]]
names(predict_alps_1950)<-c("bio5", "bio6", "bio7", "bio13", "bio19")

pm1<-predict(predict_alps_1950, m1)
plot(pm1)


colnames(presvals)<-c('bio5', 'bio6', 'bio7', 'bio13', 'bio19')
bc <- bioclim(presvals[,c('bio5', 'bio6', 'bio7', 'bio13', 'bio19')])
class(bc)

pbc<-predict(predict_alps_1950, bc)
plot(pbc)

years<-1950:2016
predict_stack<-stack(paste(wd, "/Bioclim/",years,"_bioclim_variable_stack.tif", sep=""))
#predict_alps<-crop(predict_stack, e)
bio_layer_pred<-c(5,6,13,15,19)


for (i in 1:length(years)){
  predict_stack_year<-predict_stack[[bio_layer_pred]]
  predict_alps_year<-crop(predict_stack_year, e)
  names(predict_alps_year)<-c("bio5", "bio6", "bio7", "bio13", "bio19")
  pbc<-predict(predict_alps_year, bc)
  writeRaster(pbc, paste("D:/Fiona/Git_Method/Git_Method/Alp_SDMs/Bioclim/bioclim_",years[i],"capra_ibex.tif", sep=""), overwrite=T)
  bio_layer_pred<-bio_layer_pred+19
  print(i)
}









