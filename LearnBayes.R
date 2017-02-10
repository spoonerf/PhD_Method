install.packages("LearnBayes")
library(LearnBayes)

tstatistic=function(x,y){
  m=length(x)
  n=length(y)
  sp=sqrt(((m-1)*sd(x)^2)+(n-1)*sd(y)^2/m+n-2)
  t.stat=(mean(x)-mean(y))/(sp*sqrt(1/m+1/n))
  return(t.stat)
}


alpha=0.1   #the stated level of significance
m=10
n=10
N=10000
n.reject=0
for (i in 1:N){
  x=rnorm(m, mean=0, sd=1)
  y=rnorm(n, mean=0, sd=1)
  t.stat=tstatistic(x,y)
  if (abs(t.stat)>qt(1-alpha/2,n+m-2))
    n.reject=n.reject+1
  }

true.sig.level=n.reject/N   # the true level of significance - the proportion of times the t statistic passes a certain threshold
true.sig.level

###################plotting the t statistic

m=10
n=10
my.tsimulations=function()
tstatistic(rnorm(m,mean=10,sd=2), rexp(n,rate=1/10))        #storing the t statistic, comparing two distributions

tstat.vector=replicate(10000,my.tsimulations())             #replicating the sampling and comparison 10000 times

plot(density(tstat.vector),xlim=c(-5,8), ylim=c(0,0.7), lwd=3)
curve(dt(x,df=18), add=TRUE)
legend(4, .3, c("exact", "t(18)"), lwd=c(3,1))


######Exercises Chapter 1

data("studentdata")

hist(studentdata$Dvds)
summary(studentdata$Dvds)
table(studentdata$Dvds)
barplot(table(studentdata$Dvds))

#2

output<-boxplot(studentdata$Height ~ studentdata$Gender)
output

#3
plot(WakeUp ~ ToSleep, data=studentdata)
lm(WakeUp~ToSleep, data=studentdata)
abline(lm(WakeUp~ToSleep, data=studentdata))

#4 Performance of the traditional confidence interval for a proportion
binomial.conf.interval=function(y,n){
  z=qnorm(.95)
  phat=y/n
  se=sqrt(phat*(1-phat)/n)
  return(c(phat-z*se, phat+z*se))
}

?rbinom
n<-20
p<-0.05
y<-rbinom(n, 20 ,p)

conf<-binomial.conf.interval(y, 20)

confmat<-matrix(conf, ncol=2)

counter<-0

for (i in 1:length(confmat[,1])){
  
  if(confmat[i,2]> p &confmat[i,1] < p){
    counter<-counter+1
  }
  
}

prob<-counter/n
prob

#5

alpha=0.1   #the stated level of significance
m=10
n=10
N=10000
n.reject=0
for (i in 1:N){
  x=rnorm(m, mean=0, sd=1)
  y=rnorm(n, mean=0, sd=1)
  t.stat=tstatistic(x,y)
  if (abs(t.stat)>qt(1-alpha/2,n+m-2))
    n.reject=n.reject+1
}

true.sig.level=n.reject/N   # the true level of significance - the proportion of times the t statistic passes a certain threshold
true.sig.level


mcfun<-function(n,p,m){
  counter<-0
  true_sigs<-numeric(m)
  
  for (i in 1:m){
    
    y<-rbinom(n, 20 ,p)
    conf<-binomial.conf.interval(y, 20)
    
    confmat<-matrix(conf, ncol=2)
    
    counter<-0
    
    for (j in 1:length(confmat[,1])){
      
      if(confmat[j,2]> p &confmat[j,1] < p){
        counter<-counter+1
      }
    
      true_sig<-counter/n
      true_sigs<-true_sig
  
      }
  return(true_sigs)
  }
}

mcfun(10,0.05,1000)

mcfun(25,0.05,1000)

mcfun(100,0.05,1000)
  
mcfun(10,0.25,1000)

mcfun(25,0.25,1000)

mcfun(100,0.25,1000)

mcfun(10,0.5,1000)

mcfun(25,0.5,1000)

mcfun(100,0.5,1000)

##########Chapter 2 - Introduction to Bayesian Thinking

p=seq(0.05, 0.95, by=0.1)

prior=c(1,5.2,8,7.2,4.6,2.1,0.7,0.1,0,0)
prior=prior/sum(prior)
plot(p,prior,type="h",ylab="Prior Probability")

data<-c(11,16)   #11 successes and 16 failures - 11 with >8hrs sleep
post<-pdisc(p,prior,data)
round(cbind(p,prior,post),2)

library(lattice)
PRIOR<-data.frame("prior",p,prior)
POST<-data.frame("posterior",p,post)
names(PRIOR)<-c("Type", "P", "Probability")
names(POST)<-c("Type", "P", "Probability")
data<-rbind(PRIOR, POST)
xyplot(Probability~P|Type, data=data, layout=c(1,2), type="h",lwd=3,col="black")


###Using a Beta Prior - for continuous probabilities

quantile2=list(p=0.9, x=0.5)  #90% sure the probability is less than 0.5 
quantile1=list(p=0.5, x=0.3)  #believe that the median value is 0.3, 50% above and 50% below 
beta.select(quantile1, quantile2)  #to find the beta density values 

##plotting the prior, likelihood and posterior as curves

a<-3.26
b<-7.19
s<-11
f<-16

#density function of beta distribution
curve(dbeta(x, a+s, b+f), from=0, to=1,
      xlab="p", ylab="Density", lty=1, lwd=4)
curve(dbeta(x, s+1, f+1), add=TRUE, lty=2,lwd=4)
curve(dbeta(x,a,b), add=TRUE, lty=3, lwd=4)
legend(0.7,4,c("Prior", "Likelihood", "Posterior"),
       lty=c(3,2,1), lwd=c(3,3,3))

#probability of heavy sleepers is greater than 0.5?

1-pbeta(0.5, a+s, b+f)

qbeta(c(0.05, 0.95),a+s, b+f)







