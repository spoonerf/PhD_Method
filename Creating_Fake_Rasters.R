example<-raster( paste(wd,"/Alp_SDMs/Ensembles_2006_2016/pres_abs_true_skill_stat_weighted_ensemble_sdm_1950.tif", sep=""))

values(example)[values(example) == 0] <- 1
plot(example)
sub<-round(cellStats(example, "sum")/(length(1950:2016)+1))
ncell(example)
     
for (i in 1950:2016){                  
    r<-sampleRandom(example, size=sub, asRaster=T, na.rm=T)
    r[is.na(r)] <- 0
    plot(r)
    example<-example - r
    writeRaster(example, paste(wd, "/Alp_SDMs/Fake_Landscape/Fake_Landscape_", i, ".tif", sep=""), overwrite=T)
    example[example == 0] <- NA
    plot(example)
    print(cellStats(example, "sum"))
}

