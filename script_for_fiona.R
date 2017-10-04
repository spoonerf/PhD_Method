KN_nocull <- c(7900, 8064, 8320, 8371, 8869, 9152, 8356, 9276, 10459, 11500, 11672,11483,12427,13050,12930,13573,13750,14454,16571)
KN_cull <- c(6000, 6400, 7600, 8200, 8700, 7800, 7500, 7800, 7600, 7400, 7000, 7500, 7300, 7250, 7200, 8000, 8500, 8300, 6800, 7400, 6900, 7200, 7300, 7200, 7250, 7400, 7600)
KN_full <- c(KN_cull, KN_nocull)
years=1967:2012
years_cull=1967:1993
years_nocull=1994:2012


west_alps<- c(16,25,34,42,50,59,69,81,92,102,114,130,152,183,225,278,336) #results from western alps gam,
years<-1989:2005


library(filzbach)
library(plotrix)
# transitions says which

#transitions=rep(TRUE,length(KN_nocull))# this says which data point to use for model fitting
transitions=rep(TRUE,length(west_alps))# this says which data point to use for model fitting
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

#original
filz.par_var=list(r=c(0,0.5,0.25,0,-1,1), # 1 st argument lower bound, 2nd upper bound, 3rd best guess. these are the only thing you need to change
                  k=c(40000,70000,55000,0,-1,1),
                  mu=c(-0.1,0.1,0,0,-1,1),
                  sigma=c(0,0.6,0.2,0,-1,1)) 
#original^^

filz.par_var=list(r=c(-0.2,0.56,0.21,0,-1,1), # 1 st argument lower bound, 2nd upper bound, 3rd best guess. these are the only thing you need to change
                  k=c(200,800,500,0,-1,1),
                  mu=c(-0.1,0.1,0,0,-1,1),
                  sigma=c(0,0.6,0.2,0,-1,1)) 

# pop=KN_nocull # set training dataset to kruger elephant data in the no-culling period
# t=length(KN_nocull)
# nobs=t
# No=KN_nocull[1] # starting population to the first data point
# parameters_var=matrix(nrow=0,ncol=4)

pop=west_alps # set training dataset to kruger elephant data in the no-culling period
t=length(west_alps)
nobs=t
No=west_alps[1] # starting population to the first data point
parameters_var=matrix(nrow=0,ncol=4)



# run MonteCarlo Markov Chain parameter optimization with 50,000 it burn-in and 40,000 sampling using the likelihood function
# above and the input parameters in filx.par_var. sample 1 in 10 parameters from posterior distribution
# run the sampler 5 times and rbind it all.
for (i in 1:5){
  parameters_var=rbind(parameters_var,runMCMC(50000,40000,loglike_var,nobs,filz.par_var,thinning=10)) 
}

summary(parameters_var) # growth rate, carrying capacity, mean and std of noise
