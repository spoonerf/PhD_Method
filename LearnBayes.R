##Chapter 1
install.packages("LearnBayes")
library(LearnBayes)


data(studentdata)
studentdata[1,]
table(studentdata$Drink)
barplot(table(studentdata$Drink), xlab="Drink", ylab="Count")

hours.of.sleep<-studentdata$WakeUp - studentdata$ToSleep

summary(hours.of.sleep)

hist(hours.of.sleep, main="")

boxplot(hours.of.sleep ~ studentdata$Gender, ylab="Hours of Sleep")

female.Haircut<- studentdata$Haircut[studentdata$Gender =="female"]
male.Haircut<- studentdata$Haircut[studentdata$Gender =="male"]

summary(female.Haircut)
summary(male.Haircut)

plot(jitter(studentdata$ToSleep), jitter(hours.of.sleep))

fit<-lm(hours.of.sleep ~ studentdata$ToSleep)
abline(fit)



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

###################plotting the t statistic

m=10
n=10
my.tsimulations=function()
tstatistic(rnorm(m,mean=10,sd=2), rexp(n,rate=1/10))        #storing the t statistic, comparing two distributions

tstat.vector=replicate(10000,my.tsimulations())             #replicating the sampling and comparison 10000 times

plot(density(tstat.vector),xlim=c(-5,8), ylim=c(0,0.7), lwd=3)
curve(dt(x,df=18), add=TRUE)
legend(4, .3, c("exact", "t(18)"), lwd=c(3,1))


