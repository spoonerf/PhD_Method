library(raster)
setwd("C:/Users/Fiona/Desktop/alpine_gif")
alpine<-brick("Alpine_Ibex_Projection.tif")

year<-seq(1950:2016)

png(file="%04d_alpine.png", width=1440, height=720)
for (i in 1:nlayers(alpine)) {   
   plot(alpine[[i]],axes=FALSE,legend=FALSE, main=i+1949)
  #writeRaster(alpine[[i]], paste(i+1949, "_alpine.png", sep=""))
  }
dev.off()

system('"C:\\Program Files\\ImageMagick-7.0.7-Q16\\convert.exe\" -delay 20 -loop 1 *.png alpine_ibex.gif')

file.remove(list.files(pattern=".png"))
