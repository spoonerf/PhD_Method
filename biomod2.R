install.packages("biomod2",repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("slam",repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("spocc", dependencies = T,repos=c("http://rstudio.org/_packages", "http://cran.rstudio.com"))
install.packages("stringi")
library(biomod2)
library(spocc)

sloth<-occ("acridotheres tristis", from="gbif")
df<-occ2df(sloth)

myRespName = "Bradypusvariegatus"


myRespXY<-data.frame(df$longitude, df$latitude)
myRespXY<-unique(na.omit(myRespXY))

myResp<-data.frame(rep(1, length(myRespXY$df.longitude)))

myResp<-data.frame(myRespXY, myResp)

smp_size<-round(nrow(myResp)*0.7)
train_ind <- sample(seq_len(nrow(myResp)), size = smp_size)

train <- myResp[train_ind, ]
test <- myResp[-train_ind, ]

xy<-data.frame(train[,c(1:2)])
dt<-data.frame(train[,3])

myResp<-SpatialPointsDataFrame(xy,dt)

xye<-data.frame(test[,c(1:2)])
dte<-data.frame(test[,3])

evalResp<-SpatialPointsDataFrame(xye,dte)


env <- getData("worldclim", var="bio", res=5)

myExpl<-env[[c(3,4,7,11,12)]]

myBiomodData<-BIOMOD_FormatingData(resp.var = myResp, expl.var = myExpl, 
                     resp.xy = myRespXY, resp.name= myRespName,
                     PA.nb.rep = 1,
                     PA.nb.absences = 20, PA.strategy = "sre")
                     
#                     eval.resp.var = dte, eval.resp.xy = xye)

plot(myBiomodData)

