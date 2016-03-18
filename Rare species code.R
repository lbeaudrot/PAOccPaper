 source("helper functions WPI and rare species analysis.R")

#Parameters for a species that is not too rare
sites=60
days=15
psi=0.1
p=0.1
phi=0.2
gamma=0.1
nyears=5

 #function to generate some data
 f.data.generator<-function(sites,days,psi,p,phi,gamma,nyears) {
   #first year of data
   y1<-matrix(NA,nr=sites,nc=days)
   #generate the expected occupancies
   z<-rbinom(sites,1,psi)
   #generate the observations
   for(i in 1:sites)
     y1[i,]<-rbinom(days,1,z[i]*p)
   #subsequent years
   #three dimensional matrix to store the results
   yk<-array(NA,dim=c(sites,days,nyears))
   yk[,,1]<-y1
   for(k in 2:nyears){
     #generate the deterministic part of the model
     occ<-apply(yk[,,k-1],1,max,na.rm=T)
     z<-rbinom(sites,1,occ*phi+(1-occ)*gamma)
     #generate the observations
     for(i in 1:sites)
       yk[i,,k]<-rbinom(days,1,z[i]*p)
     
   }  
   yk
 } 
 
notRare<-f.data.generator(sites,days,psi,p,phi,gamma,nyears)

#estimate detection probability from the data
#ndetections<-apply(notRare,c(1,3),sum,na.rm=T)
#nsamples<-apply(!is.na(notRare),c(1,3),sum)
#detProb<-ndetections/nsamples
#indx<-which(detProb>0)
#muPrior<-mean(detProb[indx])
#thauPrior<-1/(sd(detProb[indx]))^2

#use these estimates as priors in the JAGS model
jags.data <- list(y = notRare, nsite = sites, nrep = days, nyear = nyears)

# Initial values
initial <- apply(notRare, c(1,3), max, na.rm = TRUE)
inits <- function(){ list(z = initial)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "lambda") 


# MCMC settings
ni <- 5000
nt <- 4
nb <- 500
nc <- 3
require(R2jags)
out <- jags(jags.data, inits, params, "Dynocc-jags-priorOnDetProb.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out, dig = 2)

#compare estimates
summ<-out$BUGSoutput$summary
indx<-which(rownames(summ)=="psi[1]")

plot(1:nyears,summ[indx:(indx+nyears-1),1],t='b',ylim=c(0,1));lines(summ[indx:(indx+nyears-1),3],lty=2);lines(summ[indx:(indx+nyears-1),7],lty=2)
lines(f.cal.psis(psi,phi,gamma,nyears),col='red',lty=2)

#do it with a species with half the detection probability

sites=60
days=15
psi=0.1
p=0.05
phi=0.2
gamma=0.1
nyears=3
notRare<-f.data.generator(sites,days,psi,p,phi,gamma,nyears)

#estimate detection probability from the data
ndetections<-apply(notRare,c(1,3),sum,na.rm=T)
nsamples<-apply(!is.na(notRare),c(1,3),sum)
detProb<-ndetections/nsamples
indx<-which(detProb>0)
muPrior<-mean(detProb[indx])
thauPrior<-1/(sd(detProb[indx]))^2

#use these estimates as priors in the JAGS model
jags.data <- list(y = notRare, nsite = sites, nrep = days, nyear = nyears,mu=muPrior,thau=thauPrior)

# Initial values
initial <- apply(notRare, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}

# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "lambda") 


# MCMC settings
ni <- 5000
nt <- 4
nb <- 500
nc <- 3
require(R2jags)
out <- jags(jags.data, inits, params, "Dynocc-jags-priorOnDetProb.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out, dig = 2)

#compare estimates
summ<-out$BUGSoutput$summary
indx<-which(rownames(summ)=="psi[1]")

plot(1:nyears,summ[indx:(indx+nyears-1),1],t='b',ylim=c(0,1));lines(summ[indx:(indx+nyears-1),3],lty=2);lines(summ[indx:(indx+nyears-1),7],lty=2)
lines(f.cal.psis(psi,phi,gamma,nyears),col='red',lty=2)

#do it with a species with very low  detection probability

sites=60
days=15
psi=0.1
p=0.01
phi=0.2
gamma=0.1
nyears=3
notRare<-f.data.generator(sites,days,psi,p,phi,gamma,nyears)

#estimate detection probability from the data
ndetections<-apply(notRare,c(1,3),sum,na.rm=T)
nsamples<-apply(!is.na(notRare),c(1,3),sum)
detProb<-ndetections/nsamples
indx<-which(detProb>0)
(muPrior<-mean(detProb[indx]))
(thauPrior<-1/(sd(detProb[indx]))^2)

#use these estimates as priors in the JAGS model
jags.data <- list(y = notRare, nsite = sites, nrep = days, nyear = nyears)

# Initial values
initial <- apply(notRare, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}
inits1 <- list(z = initial,p=0.5,phi=0.1)
 inits2 <- list(z = initial,p=.2,phi=.5)
 inits3 <- list(z = initial, .RNG.seed=8)
 
# Parameters monitored
params <- c("psi", "phi", "gamma", "p", "lambda") 


# MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 100
require(R2jags)
out <- jags(jags.data, inits, params, "Dynocc-jags-priorOnDetProb.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out, dig = 2)

#compare estimates
summ<-out$BUGSoutput$summary
indx<-which(rownames(summ)=="psi[1]")

plot(1:nyears,summ[indx:(indx+nyears-1),1],t='b',ylim=c(0,1));lines(summ[indx:(indx+nyears-1),3],lty=2);lines(summ[indx:(indx+nyears-1),7],lty=2)
lines(f.cal.psis(psi,phi,gamma,nyears),col='red',lty=2)

#calculate the expected modes for psi(1), psi(2) and psi(3)
#extraction the simulation results matrix
qwe<-out$BUGSoutput$sims.array
#psi's are the 7th, 8th and 9th variables

#for psi1
psi1<-qwe[,,7]
modesDist<-(apply(psi1,2,f.mode))
plot(density(modesDist))
abline(v=quantile(modesDist,probs=c(0.025,0.975)))
abline(v=f.mode(modesDist),lwd=3)
abline(v=0.1)
abline(v=median(modesDist),lty=2)
(median(modesDist))

#Try another approach - use the prior for both the detection probability and the first year occupancy
obs.psi1<-sum(apply(notRare[,,1],1,max,na.rm=T))/sites
#Use this value as a prior for first year occupancy in the jags model
jags.data <- list(y = notRare, nsite = sites, nrep = days, nyear = nyears,mu=muPrior,thau=1000,mupsi=obs.psi1,thaupsi=100)

# MCMC settings
ni <- 10000
nt <- 2
nb <- 500
nc <- 5

require(R2jags)
out <- jags(jags.data, inits, params, "Dynocc-jags-priorOnDetProb-psiYear1.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out, dig = 2)

summ<-out$BUGSoutput$summary
indx<-which(rownames(summ)=="psi[1]")

plot(1:nyears,summ[indx:(indx+nyears-1),1],t='b',ylim=c(0,1));lines(summ[indx:(indx+nyears-1),3],lty=2);lines(summ[indx:(indx+nyears-1),7],lty=2)
lines(f.cal.psis(psi,phi,gamma,nyears),col='red',lty=2)

qwe<-out$BUGSoutput$sims.array
#psi's are the 7th, 8th and 9th variables

#for psi1
psi1<-qwe[,,7]
plot(density(psi1))
abline(v=0.1)
abline(v=mean(psi1),lty=2)
abline(v=quantile(psi1,c(0.025,0.975)),lty=3)

#for p
p<-qwe[,,5]
plot(density(p))
abline(v=0.01)
abline(v=mean(p),lty=2)
abline(v=quantile(p,c(0.025,0.975)),lty=3)
abline(v=f.mode(p),col='red')
f.mode(p)
median(p)
abline(v=median(p),col='green')

#go for a very VERY RARE species
sites=60
days=15
psi=0.01
p=0.01
phi=0.2
gamma=0.1
nyears=3
notRare<-f.data.generator(sites,days,psi,p,phi,gamma,nyears)
#estimate occupancy in year1 from the data
(obs.psi1<-sum(apply(notRare[,,1],1,max,na.rm=T))/sites)
#estimate detection probability from the data
ndetections<-apply(notRare,c(1,3),sum,na.rm=T)
nsamples<-apply(!is.na(notRare),c(1,3),sum)
detProb<-ndetections/nsamples
indx<-which(detProb>0)
(muPrior<-mean(detProb[indx]))
(thauPrior<-1/(sd(detProb[indx]))^2)

jags.data <- list(y = notRare, nsite = sites, nrep = days, nyear = nyears,mu=muPrior,thau=thauPrior,mupsi=obs.psi1,thaupsi=100)

# Initial values
initial <- apply(notRare, c(1,3), max, na.rm = TRUE)
#tmp[tmp=="-Inf"]<-NA # remove the -Inf's that result when the camera trap was stolen
inits <- function(){ list(z = initial)}

# MCMC settings
ni <- 1000
nt <- 2
nb <- 500
nc <- 5

require(R2jags)
out <- jags(jags.data, inits, params, "Dynocc-jags-priorOnDetProb-psiYear1.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb, working.directory = getwd())
print(out, dig = 2)

summ<-out$BUGSoutput$summary
indx<-which(rownames(summ)=="psi[1]")

plot(1:nyears,summ[indx:(indx+nyears-1),1],t='b',ylim=c(0,1));lines(summ[indx:(indx+nyears-1),3],lty=2);lines(summ[indx:(indx+nyears-1),7],lty=2)
lines(f.cal.psis(psi,phi,gamma,nyears),col='red',lty=2)

qwe<-out$BUGSoutput$sims.array
#psi's are the 7th, 8th and 9th variables

#for psi1
psi1<-qwe[,,7]
plot(density(psi1))
abline(v=0.01)
abline(v=mean(psi1),lty=2)
abline(v=f.mode(psi1),col="red")
abline(v=median(psi1),col="green")
abline(v=quantile(psi1,c(0.025,0.975)),lty=3)

#for p
p<-qwe[,,5]
plot(density(p))
abline(v=0.01)
abline(v=mean(p),lty=2)
abline(v=quantile(p,c(0.025,0.975)),lty=3)
abline(v=f.mode(p),col='red')
f.mode(p)
median(p)
abline(v=median(p),col='green')

#Major conclusion: When data is sparse,  the mode of the posterior distribution of psi and p seems to be a better
#estimator of the real values compared to the median (next best) or the mean (worse)
#also using the observed values of psi and p as priors seems to help
#NEXT STEPS: Setup a simulation for various combinations of p and psi where data sets are generated and then
#compared the recovered mean, median and modes of the posterior distributions with the actual true values.
#NEED to do a workflow of this process

#test out ggmcmc package (to plot diagnostics from JAGS output)
# see http://xavier-fim.net/packages/ggmcmc/ for tutorial
m.jags <- jags.model("Dynocc-jags-priorOnDetProb-psiYear1.txt", data = jags.data, inits=inits, n.adapt = 5000, n.chains = nc)
S <- coda.samples(m.jags, params, n.iter = 10000)
D<-ggs(S,parallel=FALSE)
ggmcmc(D)
ggs_running(D)