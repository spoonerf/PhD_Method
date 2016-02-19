LPI_EDScores<-read.csv("2016_01_05_LMEModel_data.csv")
Inc5MammalsTemp <- subset(LPI_EDScores,!is.na(mean_slope_se)& !is.na(change_rate)&!is.na(mean_slope) &System != "Marine"&Include10 == "Yes"&r_sq2 >= 0.5 )
nrow(Inc5MammalsTemp)
years<-Inc5MammalsTemp[,c(2,97:161)]

first=numeric(0)

for (i in 1:length(years$ID)){
  
  first_year<-which(!is.na(years[i,c(2:66)]))[1] + 1949
  print(first_year)
  first<-rbind(first, first_year)
  }

Inc5MammalsTemp$first_year<-first

Inc5MammalsTemp<-  Inc5MammalsTemp[order(Inc5MammalsTemp$first_year),]   #min year only starts at 1970 whereas data starts at 1950
#,-Inc5MammalsTemp$length_time

year_plot<-Inc5MammalsTemp[,c(97:161)]

yr<-data.matrix(year_plot)

yr[yr != "NA"]<-1
yr[is.na(yr)]<-0

yrr<-raster(yr)
plot(yrr)

