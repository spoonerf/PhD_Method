dt<-read.csv("inla_test.csv")

source("rsquaredglmm.R")

library(lme4) 

m0<-lmer(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+(1|Binomial)+(1|loc_id),data=dt, REML=F)

summary(m0)
library(INLA)

formula<-formula(lambda_mean ~ change_rate_scale+mean_slope_scale+change_rate_scale:mean_slope_scale+Bodymass_scale+ f(loc_id, model="iid") + f(Binomial, model="iid"))		

inla2<-inla(formula,data=dt,
            family="gaussian",
            Ntrials=rep(1,nrow(dt)),
            verbose = TRUE,
            control.predictor=list(compute=TRUE))

tmp = inla2$summary.fitted.values[1:nrow(dt), 'mean']-dt$lambda_mean
plot(density(tmp))
xvals = seq(-1, 1, length.out=1000)
lines(xvals, dnorm(xvals, mean=mean(tmp), sd=sd(tmp)), col='green')


summary(inla2)

str(inla2)


############Survival Analysis

data(Leuk)
sapply(Leuk, summary)


##building the mesh - required for any spatial models
loc<-cbind(Leuk$xcoord, Leuk$ycoord)
bnd1<-inla.nonconvex.hull(loc, convex = 0.05)
bnd2<-inla.nonconvex.hull(loc, convex = 0.25)
mesh<- inla.mesh.2d(loc, boundary = list(bnd1, bnd2), max.edge = c(0.05, 0.2), cutoff=0.005)

##projector matrix is obtained with

A<-inla.spde.make.A(mesh, loc)

## The spde model is created with 

spde<- inla.spde2.matern(mesh, alpha=2)

##The model formula including the intercept, covariates and the SPDE model is:

formula <- inla.surv(time, cens) ~ 0 + a0 +
  sex + age + wbc + tpi + 
  f(spatial, model=spde)
 #response variables are time and cens

stk <- inla.stack(data=list(time=Leuk$time, cens=Leuk$cens),
                  A=list(A, 1), 
                  effect=list(
                    list(spatial=1:spde$n.spde), 
                    data.frame(a0=1, Leuk[,-c(1:4)])))


r <- inla(formula, family="weibull", data=inla.stack.data(stk),control.predictor=list(A=inla.stack.A(stk)))

#The intercept and the covariate effects can be extracted with
round(r$summary.fix, 4)

#and the hyperparameters with

round(r$summary.hy, 3)

##We visualize the spatial effect into the map. The map of the districts is also available into the INLA package. 
##First, we define a projection from the mesh into a grid

r0 <- diff(range(bbox(nwEngland)[1,]))/diff(range(bbox(nwEngland)[2,])) 
prj <- inla.mesh.projector(mesh, xlim=bbox(nwEngland)[1,], 
                           ylim=bbox(nwEngland)[2,], 
                           dims=c(200*r0, 200))

#then we interpolate it and assign NA for grid points not inside the map

m.spat <- inla.mesh.project(prj, r$summary.ran$spatial$mean) 
sd.spat <- inla.mesh.project(prj, r$summary.ran$spatial$sd) 
ov <- over(SpatialPoints(prj$lattice$loc), nwEngland)
sd.spat[is.na(ov)] <- m.spat[is.na(ov)] <- NA

# The posterior mean and standard deviation are in Figure 5.1. 
#As a result, the spatial effect has continuous variation along the 
#region, rather than constant inside each district.

library(fields)

par(mfrow=c(1,2), mar=c(0,0,0,0))
image.plot(x=prj$x, y=prj$y, z=m.spat, asp=1, 
           xlab='', ylab='', axes=FALSE, horizontal=TRUE) 
plot(nwEngland, add=TRUE)
image.plot(x=prj$x, y=prj$y, z=sd.spat, asp=1, 
           xlab='', ylab='', axes=FALSE, horizontal=TRUE)
plot(nwEngland, add=TRUE)






