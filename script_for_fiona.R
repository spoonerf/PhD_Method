ibex<-read.csv("LPI_pops_20160523_edited.csv")
ibex<-ibex[ibex$ID == 539 | ibex$ID == 10692|ibex$ID == 10694|ibex$ID == 10695|ibex$ID == 10696|
             ibex$ID == 10710 |ibex$ID == 10713|ibex$ID == 10714|ibex$ID == 10717 |ibex$ID == 10718,]

pop_datab <- ibex[,c(65:130)]

#539
pop_539<-c(2767, 2873, 3367, 3541, 3370, 3674, 3822, 2803, 3222, 3537, 3266, 3513, 3822, 
           3203,3431,3382, 2746, 3158, 3116, 3590, 3754, 2410, 3084, 3042, 3234,
           3412, 3255, 3187, 3362, 3740, 3914, 4283, 4303, 4694, 4631, 4790, 4754, 4991, 4136,
           4360, 3998, 3581, 3701, 3674, 3632)
years_539<-1956:2000

#10692
pop_10692<-c(158,128,135,139,134,116,128,130,159,189,211,242,249,246,254,256,251,235,234,
             222,204,211,232,256,270,277,304,275,287,268,256,253,247,241,NA,250)
years_10692<-1950:1985
df_10692<-data.frame(pop_10692, years_10692)


pop_10692_gam<-mgcv:::gam(pop_10692 ~ s(years_10692, k=round(length(years_10692)/2)))
pred_10692<-predict(pop_10692_gam,df_10692,type="response",se=FALSE) 
pop_10692[which(is.na(pop_10692))]<-pred_10692[which(is.na(pop_10692))]

#10694
pop_10694<-c(45,36,34,34,36,38,37,44,44,48,51,58,63,61,58,54,51,51,51,54,51,47,46,44,44,42,33,30,26,26,25,21,21,21,23)
years_10694<-1950:1984

#10695 

pop_10695<-c(211,237,269,299,323,355,377,398,415,437,448,439,432,424,415,408,399,392,384,375,377,
             385,397,409,421,432,444,454,468,471,454)
years_10695<-1950:1980

#10696

pop_10696<-c(16,13,13,13,13,11,13,18,25,32,39,51,54,56,65,65,66,70,78,82,92,111,127,146,153,156,
             166,172,185,206,215,218,241,246,265)
years_10696<-1950:1984


#10710

pop_10710<-c(14,28,NA,NA,NA,NA,63,88,92,NA,NA,NA,NA,185,210,284,341)
years_10710<-1989:2005
df_10710<-data.frame(pop_10710, years_10710)

pop_10710_gam<-mgcv:::gam(pop_10710 ~ s(years_10710, k=round(length(years_10710)/2)))
pred_10710<-predict(pop_10710_gam,df_10710,type="response",se=FALSE) 
pop_10710[which(is.na(pop_10710))]<-pred_10710[which(is.na(pop_10710))]


#10713

pop_10713<-c(193,173,183,220,230,248,258,267,277,248,248,272,291,286,251,221,216,223,231,228,145,163,
             179,194,176,157,178,148,178,169,169,182,193,205,225,233,214,240,279,229)
years_10713<-1950:1989


#10714

pop_10714<-c(8,15,26,30,NA,31,32,47,50,58,64,66,72,75,90)

years_10714<-1950:1964
df_10714<-data.frame(pop_10714, years_10714)

pop_10714_gam<-mgcv:::gam(pop_10714 ~ s(years_10714, k=round(length(years_10714)/2)))
pred_10714<-predict(pop_10714_gam,df_10714,type="response",se=FALSE) 
pop_10714[which(is.na(pop_10714))]<-pred_10714[which(is.na(pop_10714))]

#10717
pop_10717<-c(23,32,34,36,42,44,46,55,60,73)
years_10717<-1955:1964

#10718
pop_10718<-c(652,353,387,464,481,515,626,523,531,694,762,531,821,855,787,889,
             847,761,864,1034,NA,949,957,1017,1094,NA,1076,1042,923,735,649,658,
             700,880,1067,999,1007,1212,1434,1519,1553,1596,1468,1280,1254,1152,1254,1263,
             1237,1220,1237,1075,801)
years_10718<-1950:2002
df_10718<-data.frame(pop_10718, years_10718)

pop_10718_gam<-mgcv:::gam(pop_10718 ~ s(years_10718, k=round(length(years_10718)/2)))
pred_10718<-predict(pop_10718_gam,df_10718,type="response",se=FALSE) 
pop_10718[which(is.na(pop_10718))]<-pred_10718[which(is.na(pop_10718))]


library(filzbach)
library(plotrix)
# transitions says which

pop_id<-pop_10718

transitions=rep(TRUE,length(pop_id))# this says which data point to use for model fitting
#transitions[c(3:5)]=FALSE # set to false years where elephants left KNP because of a drought
loglike_var=function(r,k,mu,sigma){
  res=matrix(nrow=nobs-1)
  prob=matrix(nrow=nobs-1)
  for (te in 2:nobs){
    if (transitions[te-1]){
      res[te-1]=(pop[te]/pop[te-1])-1-(r*(1-pop[te-1]/k))
      prob[te-1]=dnorm(res[te-1],mean=mu,sd=sigma)
      loglike=sum(log(prob),na.rm=TRUE)
    }
  }
  return(loglike)
} 

filz.par_var=list(r=c(-0.05,0.25,0.1,0,-1,1), # 1 st argument lower bound, 2nd upper bound, 3rd best guess. these are the only thing you need to change
                  k=c(300,2400,900,0,-1,1),
                  mu=c(-0.1,0.1,0,0,-1,1),
                  sigma=c(0,0.6,0.2,0,-1,1)) 


pop=pop_id # set training dataset to kruger elephant data in the no-culling period
t=length(pop_id)
nobs=t
No=pop_id[1] # starting population to the first data point
parameters_var=matrix(nrow=0,ncol=4)


# run MonteCarlo Markov Chain parameter optimization with 50,000 it burn-in and 40,000 sampling using the likelihood function
# above and the input parameters in filx.par_var. sample 1 in 10 parameters from posterior distribution
# run the sampler 5 times and rbind it all.
for (i in 1:5){
  parameters_var=rbind(parameters_var,runMCMC(50000,40000,loglike_var,nobs,filz.par_var,thinning=10)) 
}
summary(parameters_var) # growth rate, carrying capacity, mean and std of noise



colMeans(parameters_var)

hist(parameters_var[,1])
hist(parameters_var[,2])
hist(parameters_var[,3])
hist(parameters_var[,4])


quantile(parameters_var[,1],probs=c(0.025,0.5,0.975))




k_539<-parameters_var[,2]
k_10692<-parameters_var[,2]
k_10694<-parameters_var[,2]
k_10695<-parameters_var[,2]
k_10696<-parameters_var[,2]
k_10710<-parameters_var[,2]
k_10713<-parameters_var[,2]
k_10714<-parameters_var[,2]
k_10717<-parameters_var[,2]
k_10718<-parameters_var[,2]


k_539m<-mean(k_539)
k_10692m<-mean(k_10692)
k_10694m<-mean(k_10694)
k_10695m<-mean(k_10695)
k_10696m<-mean(k_10696)
k_10710m<-mean(k_10710)
k_10713m<-mean(k_10713)
k_10714m<-mean(k_10714)
k_10717m<-mean(k_10717)
k_10718m<-mean(k_10718)

ID<-c(539,10692,10694,10695,10696,10710,10713,10714,10717,10718)
ks<-c(k_539m,k_10692m, k_10694m, k_10695m, k_10696m, k_10710m, k_10713m, k_10714m, k_10717m, k_10718m)


dfk<-data.frame(ID, ks)

write.csv(dfk, "carrying_capacity_log_likelihood_capra_ibex.csv")
