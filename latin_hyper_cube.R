install.packages("lhs")
library(lhs)



med_dens<-log(70.6)
max_dens<-log(443)

lsd<-(max_dens - med_dens)/3

X <- randomLHS(93, 4)
X[,1] <- qnorm(X[,1], 0.25, 0.08333)
X[,2] <- qnorm(X[,2], 0.25, 0.08333)
X[,3] <- qnorm(X[,3], 0.25, 0.08333)
X[,4] <- qlnorm(X[,4], med_dens ,lsd)

X



