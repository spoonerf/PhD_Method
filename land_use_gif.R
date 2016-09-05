library(raster)

primf<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/states.nc", varname="primf")
primn<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/states.nc", varname="primn")
#decadal
pastr<-brick("C:/Users/Fiona/Documents/PhD/PhD_Method/states.nc", varname="pastr")

year<-seq(850,2015,10)
lyr<-seq(1,1166,10)

png(file="%04d_primf.png", width=1440, height=720)
for (i in lyr) {    #1:length(year)

  prim<-primf[[i]]+primn[[i]]  
  plot(prim,axes=FALSE,legend=FALSE, main=i+849)
}
dev.off()


system('"C:\\Program Files\\ImageMagick-7.0.2-Q16\\convert.exe\" -delay 20 -loop 1 *.png prim_both_dec_850c.gif')

file.remove(list.files(pattern=".png"))


#annual
primf2<-primf[[1091:1166]]
year<-1940:2015

png(file="%04d_primf.png", width=1440, height=720)
for (i in 1:nlayers(primf)) {   
  
  plot(primf2[[i]],axes=FALSE,legend=FALSE, main=year[i])
}
dev.off()

system('"C:\\Program Files\\ImageMagick-7.0.2-Q16\\convert.exe\" -delay 20 -loop 0 *.png primf_annual_1940.gif')

file.remove(list.files(pattern=".png"))


##decadal to 1850 then annual

lyr1<-seq(1,991,30)
lyr2<-seq(1001,1166,2)

#year<-c(year1,year2)
lyr<-c(lyr1,lyr2)

png(file="%04d_primf.png", width=1440, height=720)
for (i in lyr) {    #1:length(year)
  prim<-(primf[[i]]+primn[[i]])*100 
  plot(prim,axes=F,legend=T, main=paste("Percentage Cover of Primary Habitat - Year:",i+849,sep=" "), ylim=c(-55,90), xlim=c(-180,180))
}
dev.off()

system('"C:\\Program Files\\ImageMagick-7.0.2-Q16\\convert.exe\" -delay 20 -loop 0 *.png prim_both_dec_and_ann_short.gif')

file.remove(list.files(pattern=".png"))

####pasture decadel to 1850 then annual

lyr1<-seq(1,991,10)
lyr2<-seq(1001,1166,1)

#year<-c(year1,year2)
lyr<-c(lyr1,lyr2)

png(file="%04d_pastr.png", width=1440, height=720)
for (i in lyr) {    #1:length(year)
  past<-pastr[[i]]*100 
  plot(past,axes=F,legend=T, main=paste("Percentage Cover of Pasture - Year:",i+849,sep=" "), ylim=c(-55,90), xlim=c(-180,180))
}
dev.off()

system('"C:\\Program Files\\ImageMagick-7.0.2-Q16\\convert.exe\" -delay 20 -loop 0 *.png pastr_both_dec_and_ann_leg_twit.gif')

file.remove(list.files(pattern=".png"))







