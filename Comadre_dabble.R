dir<-getwd()
load(paste(dir, "COMADRE_v.1.0.0.RData", sep="/"))

tempMetadata<-subset(comadre$metadata, SpeciesAccepted=="Ovis_aries")

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]   #MatA is matrix pop model, can be split into U, F and C

#first matrix in list is the mean of others?

#tempMatmean<-(tempMat[[1]][[1]]+ tempMat[[2]][[1]]+tempMat[[3]][[1]]+tempMat[[4]][[1]]+tempMat[[5]][[1]]+tempMat[[6]][[1]]+tempMat[[7]][[1]])/length(keep)

MatList<-list(tempMat[[1]][[1]], tempMat[[2]][[1]], tempMat[[3]][[1]],tempMat[[4]][[1]],tempMat[[5]][[1]],tempMat[[6]][[1]],tempMat[[7]][[1]])
AllMat<-unlist(MatList)
Matcol<-matrix(AllMat, ncol=7)


output<-data.frame(lambdas = rep(NA, length(tempMat)), damps= rep (NA, length(tempMat)))

require(popbio)

for (i in 1:length(tempMat)){
  tryCatch({
    output$lambdas[i]<-max(Re(eigen(tempMat[[i]]$matA)$value))
    output$damps[i]<-damping.ratio(tempMat[[i]]$matA)
    print(paste("Species number ", i, ": ", tempMetadata$SpeciesAuthor[i], sep= ""))
  }, error = function(e){})
}

par(mfrow =c(1,2))
hist(log(output$lambdas), xlab="Log population growth rate", col="gold", main=)
  abline(v=0, col="red", lwd=4, lty=3)

  

  
  
####################################################################################################
####################################################################################################
#plot a world map wild carnivores in the northern hemisphere  
  
tempMetadata<-subset(comadre$metadata, MatrixComposite == "Mean" & Order =="Carnivora" & MatrixCaptivity == "W" &LatNS =="N" &SurvivalIssue < 1 & MatrixSplit == "Divided" & MatrixFec == "Yes")

keep<-as.numeric(rownames(tempMetadata))

tempMat<-comadre$mat[keep]

output<-data.frame(lambdas =rep(NA, length(tempMat)))

#create dummy variables to convert geographic information to be plotted in map

tempMetadata$LAT=NA
tempMetadata$LON=NA

for (i in 1:dim(tempMetadata)[1]){
  if(tempMetadata$LonWE[i]=="W"){ tempMetadata$LonDeg[i]<- -tempMetadata$LonDeg[i]}
  tempMetadata$LAT[i]<-tempMetadata$LatDeg[i] + tempMetadata$LatMin[i]/60 + tempMetadata$LatSec[i]/3600
  tempMetadata$LON[i]<-tempMetadata$LonDeg[i] + tempMetadata$LonMin[i]/60 + tempMetadata$LonSec[i]/3600
}

tempMetadata$lambdas<-NA

require(popbio)

for (i in 1:length(tempMat)){
  tryCatch({
    tempMetadata$lambdas[i]<-max(Re(eigen(tempMat[[i]]$matA)$value))
    print(paste("Species number ", i, ": ", tempMetadata$SpeciesAuthor[i], sep=""))
  }, error = function(e){})
}

#creating colour codes for lambdas, blue >1, white =1, red <1

tempMetadata$LambdaCol<-NA
tempMetadata$LambdaCol[which(tempMetadata$lambdas > 1)] ="Blue"
tempMetadata$LambdaCol[which(tempMetadata$lambdas == 1)] ="White"
tempMetadata$LambdaCol[which(tempMetadata$lambdas < 1)] ="Red"

require(maps)
require(scales)

par(mfrow =c(1,1))
map("world", col="white", fill=T, xlim=c(-175, -20), ylim=c(20,60), border="black", mar=rep(0,4))
points(jitter(tempMetadata$LON, amount=.6), jitter(tempMetadata$LAT, amount=.6), col="black", bg=alpha(tempMetadata$LambdaCol,.5), ylim=c(-90,90), xlim=c(-180,180),type="p",pch=21,cex=1)
















