LPI<-read.csv("LPI_pops_20160523_edited.csv")
LPI<-subset(LPI, Specific_location==1 )

xy<-cbind(LPI$Longitude, LPI$Latitude)

xy<-cbind(df2$Longitude, df2$Latitude)

xy<-cbind(df_bird$Longitude, df_bird$Latitude)

xy<-cbind(df_mammal$Longitude, df_mammal$Latitude)

googEl <- function(locs)  {
  require(RJSONIO)
  locstring <- paste(do.call(paste, list(locs[, 2], locs[, 1], sep=',')),
                     collapse='|')
  u <- sprintf('http://maps.googleapis.com/maps/api/elevation/json?locations=%s&sensor=false',
               locstring)
  res <- fromJSON(u)
  out <- t(sapply(res[[1]], function(x) {
    c(x[['location']]['lat'], x[['location']]['lng'], 
      x['elevation'], x['resolution']) 
  }))    
  rownames(out) <- rownames(locs)
  return(out)
}


alt<-googEl(xy)

