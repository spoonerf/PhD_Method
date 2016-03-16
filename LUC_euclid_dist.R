
LPILU<-read.csv("Global_Land_Use_All_LPI.csv") 
#should probablyy get rid of those with high "other" category as may be reflective of crappy data

LPI<-read.csv("LPI_populations_IP_fishedit_20140310_nonconf.csv")

LPIsp<-subset(LPI, Specific_location==1 & System !="Marine" & Class != "Actinopterygii"& Class != "Cephalaspidomorphi" )

doDist = function(sp_name) {
  spid2 = subset(LPIsp, ID == sp_name)   #subsetting the population data by each population
  spid = spid2[,63:118]                     #subsetting only the dates
  colnames(spid)<-1950:2005              #renaming the date column names as R doesn't like numbered column names
  lu_id=subset(LPILU, ID == sp_name)  #subsetting the climate data by each population
  
  id<-spid2$ID
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  
  if (sum(!is.na(spidt)) > 0) {
  Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  minyr<-matrix(subset(lu_id, Year==min(Year))[,c(6:11)],nrow=1)
  maxyr<-matrix(subset(lu_id, Year==max(Year))[,c(6:11)],nrow=1)
  
  mat<-rbind(minyr, maxyr)
  euc<-dist(mat, method = "euclidean")
  euc<-as.numeric(euc)  
  } else{
    euc<-NA
  }
  euc_df<-data.frame(id,euc)
  print(euc_df)  
  return(euc_df)
}

all_df_list <- lapply(unique(LPIsp$ID), doDist)

all_matrix <- matrix(unlist(all_df_list), ncol=2, byrow=TRUE)
mean_df <- data.frame(all_matrix)
colnames(mean_df) <- c("ID", "LUC_dist")

write.csv(mean_df, "LUC_distance_all.csv")




