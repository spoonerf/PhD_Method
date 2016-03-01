setwd("D:/Fiona/Git_Method/Git_Method")
getwd()
library(plyr)
library(taRifx)
library(broom)

# LPI_EDScores<-read.csv("2016_01_05_LMEModel_data.csv")
# 
# LPI_EDScores<-LPI_EDScores[,c(2,16,20,22,25,27:34,47:75,223,245:247)]
# # 
# sens<-read.csv("Sensitivity_Land_Use.csv")
# clim_sens<-read.csv("Climate_Sensitivity.csv")
# hot<-read.csv("hot_days_pop.csv")
# cold<-read.csv("cold_days_pop.csv")
# 
# LPI_EDScores2<-merge(LPI_EDScores, sens, by="ID")
# LPI_EDScores2<-merge(LPI_EDScores2, clim_sens, by="ID")
# LPI_EDScores2<-merge(LPI_EDScores2, hot, by="ID")
# LPI_EDScores2<-merge(LPI_EDScores2, cold, by="ID")
# 
# LPI<- subset(LPI_EDScores3, !is.na(change_rate_25)&!is.na(mean_slope) &System != "Marine"&Include10 == "Yes")
# 
# 
# nrow(LPI)
# 
# max<-read.csv("all_EU_daily_max_temp.csv")
# min<-read.csv("all_EU_daily_min_temp.csv")
# pcp<-read.csv("all_EU_daily_pcp_temp.csv")
# 
# pcp$binomial<-as.character(pcp$binomial)
# 
# pcp$id<-data.frame(matrix(unlist(strsplit(pcp$binomial, split='_', fixed=TRUE)), ncol=3, byrow=T))[,1]
# 
# colnames(max)[4:8]<-c("Lon", "Lat", "Year", "Month", "Day")
# colnames(min)[4:8]<-c("Lon", "Lat", "Year", "Month", "Day")
# colnames(pcp)[4:8]<-c("Lon", "Lat", "Year", "Month", "Day")
# 
# max$date <- as.Date( paste(max$Year,max$Month,max$Day , sep = "." )  , format = "%Y.%m.%d" )
# min$date <- as.Date( paste(min$Year,min$Month,min$Day , sep = "." )  , format = "%Y.%m.%d" )
# pcp$date <- as.Date( paste(pcp$Year,pcp$Month,pcp$Day , sep = "." )  , format = "%Y.%m.%d" )
# 




LPI_EDScores<-read.csv("D:/Fiona/Git_Method/Git_Method/LPI_populations_IP_fishedit_20140310_nonconf.csv")

max<-read.csv("Daily_Prec_All_EU_form.csv")

max$Date<-as.Date(max$ind, "%Y-%m-%d")
max$Year<-format(max$Date, "%Y")
colnames(max)[2]<-"id"

LPI<-LPI_EDScores[LPI_EDScores$ID %in% max$id,] 


doMax = function(sp_name) {
  spid2 = subset(LPI, ID == sp_name)   #subsetting the population data by each population
  spid = spid2[,63:125]                     #subsetting only the dates
  colnames(spid)<-1950:2012              #renaming the date column names as R doesn't like numbered column names
  climid=subset(max, id == sp_name)  #subsetting the climate data by each population

  name<-spid2$Binomial
  id<-spid2$ID
  name_id<-paste(name, id, sep="_") #creating id for naming files of plots
  Date<-as.numeric(colnames(spid))
  spidt<-destring(t(spid))
  
  Year<-Date[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  Population<-spidt[min(which(!is.na(spidt))):max(which(!is.na(spidt)))]
  
  Max_ch<-na.omit(climid[climid$Year %in% Year, ])
  Max_val<-Max_ch[,5]  
  Max_date<-Max_ch[,7]  
  
   if (sum(is.nan(Max_val))!=length(Max_val)){     #checking that there are values extracted for this population - i.e. that it is on land
    
    lm_max<-lm(Max_val~Max_date)
    lm_max_df<-tidy(lm_max)[2,]  
    res_max<-residuals(lm_max)
    Year_df<-format(Max_date,'%Y')
    res_df<-data.frame(res_max, Year_df)

    year_res<-ddply(res_df, "Year_df", summarise, var_res=var(res_max), mean_res= mean(res_max), range_res=range(res_max)[2] - range(res_max)[1])
    var_res<-mean(year_res$var_res)
    range_var<-range(year_res$var_res)[2] - range(year_res$var_res)[1]    #range of the annual variance - calculated from daily measurements
    mean_an_mean<-mean(year_res$mean_res)
    mean_an_range<-mean(year_res$range_res)
    
    
    max_df<-cbind(id,lm_max_df, var_res, range_var, mean_an_mean, mean_an_range)
    
  } else{
    
    max_df<-matrix(c(id,NA,NA,NA,NA,NA,NA,NA), nrow=1, ncol=10)
    colnames(max_df)<-c("id", "term", "estimate", "std.error", "statistic", "p.value", "range_var","var_res", "mean_an_mean", "mean_an_range")
    max_df<-data.frame(max_df)
  }

  print(max_df)
  return(max_df)
}

all_df_list <- lapply(unique(LPI$ID), doMax)

all_matrix <- matrix(unlist(all_df_list), ncol=10, byrow=TRUE)
min_df2 <- data.frame(all_matrix)
colnames(min_df2) <- c("ID", "Term","Estimate","SE","Statistic","p.val", "var_res","range_var", "mean_an_mean", "mean_an_range")

min_df2$mean_an_range<-as.numeric(as.character(min_df2$mean_an_range))
min_df2$var_res<-as.numeric(as.character(min_df2$var_res))
min_df2$range_var<-as.numeric(as.character(min_df2$range_var))
min_df2$mean_an_mean<-as.numeric(as.character(min_df2$mean_an_mean))
min_df2$mean_an_range<-as.numeric(as.character(min_df2$mean_an_range))


write.csv(min_df2, "prec_residuals_all_EU.csv")


crex<-subset(min_df2, ID==327)

