#Test code

data.generator<-function(points,days,psi,p,phi,gamma,years) {
        require(unmarked)
        #first year of data
        y1<-matrix(NA,nr=points,nc=days)
        #generate the expected occupancies
        z<-rbinom(points,1,psi)
        #generate the observations
        for(i in 1:points)
                y1[i,]<-rbinom(days,1,z[i]*p)
        #subsequent years
        #three dimensional matrix to store the results
        yk<-array(NA,dim=c(points,days,years))
        yk[,,1]<-y1
        for(k in 2:years){
                #generate the deterministic part of the model
                occ<-apply(yk[,,k-1],1,max,na.rm=T)
                z<-rbinom(points,1,occ*phi+(1-occ)*gamma)
                #generate the observations
                for(i in 1:points)
                        yk[i,,k]<-rbinom(days,1,z[i]*p)
                
        }  
        # convert results to a two dimensional matrix for colext
        yk <- matrix(yk, points, days*years)
        #ny <- matrix(as.character(1:years),nrow(yk), years, byrow=TRUE)
        # store data in an unmarkedMultFrame object so it is ready for analysis
        yUMF <- unmarkedMultFrame(y = yk, numPrimary = years)
        yUMF
} 

set.seed(400)
data <- data.generator(points = 100,days = 30, psi = 0.9, p = 0.2, phi = 0.9, 
                       gamma = 0, years = 5)

#Look at the first year of data
summary(data)
model <- colext(psiformula = ~1,gammaformula = ~1,epsilonformula = ~1, 
                pformula = ~1, data = data, method = "BFGS",se = FALSE)
summary(model)
backTransform(model,type=c("psi"))
backTransform(model,type=c("col"))
backTransform(model,type=c("ext"))
backTransform(model,type=c("det"))

year <- 1:5
y <- round(t(smoothed(model))*30)
y <- y[,c(2,1)]
lrmodel <- glm(y ~ year, family = "binomial")
lrmodel$fitted.values
summary(lrmodel)

# do a linear regression
yl <- smoothed(model)[2,]
lmodel <- lm(yl ~ year)
plot(smoothed(model)[2,])
lines(lrmodel$fitted.values)
lines(lmodel$fitted.values,lty=2)

data3d <- array(data@y,dim = c(30,30,5))
obs <- apply(apply(data3d,c(1,3),max),2,sum)/dim(data3d)[1]
mod <- as.numeric(projected(model)[2,])

plot(1:5, obs, type="p"); lines(1:5, mod) ; lines(lmodel,lty=2)

simulatetrend <- function(points = 10,days = 30, psi = 0.5, p = 0.2, phi = 0.2, gamma = 0, years = 5, nsim = 5) {
        #store results
        results <- matrix(NA, nr = nsim, nc = years)
        logisticRes <- matrix(NA, nr = nsim, nc = 2)
        linearRes <- matrix(NA, nr = nsim, nc = 2)
        gm_lambda <- numeric()
        time <- 1:years
        for (i in 1:nsim) {
                #generate data and store projected results
                data <- data.generator(points,days,psi,p,phi,gamma,years)
                model <- colext(~1, ~1, ~1, ~1, data = data, method = "BFGS",se = FALSE)
                results[i,] <- as.numeric(smoothed(model)[2,])
                
                #Calculate geometric mean of lambda
                lambda <- results[i,2:years]/results[i,1:(years-1)]
                gm_lambda[i] <- -(gm_mean(lambda) - 1)
                # Calculate logistic regression
                y <- round(t(smoothed(model))*points)
                #lrmodel <- glm(y ~ time, family = "binomial")
                #logisticRes[i,] <- lrmodel$coeff
                
                # Calculate a linear regression
                #y <- smoothed(model)[2,]
                #lmodel <- lm (y ~ time)
                #linearRes[i,] <- lmodel$coefficients
        }
        list(gmlambda = gm_lambda,logistic = logisticRes, linear = linearRes,results = results)
}

system.time(test <- simulatetrend(points = 50, psi = 0.9, gamma = 0, phi = 0.8, nsim = 100))
require(ggplot2)
res <- data.frame(test$results, id = 1:100)
names(res) <- c(1:(ncol(res)-1),"id")
res <- melt(res,"id")
names(res)[2] <- "year"
ggplot(res, aes(x=year,y=value, group=id))+geom_line(alpha=0.2)

slopes <- as.data.frame(test$gmlambda)
names(slopes) <- c("slope")#ggplot(slopes,aes(x=exp(X2))) + geom_histogram()
ggplot(slopes,aes(x=slope)) + geom_histogram() + geom_density()

gm_mean = function(x, na.rm=TRUE){
        exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
        if(any(x < 0, na.rm = TRUE)){
                return(NaN)
        }
        if(zero.propagate){
                if(any(x == 0, na.rm = TRUE)){
                        return(0)
                }
                exp(mean(log(x), na.rm = na.rm))
        } else {
                exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
        }
}
library(parallel)
library(foreach)

library(doParallel)



simulatetrend <- function(points = 10,days = 30, psi = 0.5, p = 0.2, phi = 0.2, gamma = 0, years = 5, nsim = 5) {
        #store results
        #gm_lambda <- numeric()
        time <- 1:years
        results <- foreach(i = 1:nsim, .combine = rbind, .packages = "unmarked", .export = c("data.generator","gm_mean")) %dopar% {
                #generate data and store projected results
                data <- data.generator(points,days,psi,p,phi,gamma,years)
                model <- colext(~1, ~1, ~1, ~1, data = data, method = "BFGS",se = FALSE)
                timeseries <- as.numeric(smoothed(model)[2,])
                
                #Calculate geometric mean of lambda
                lambda <- timeseries[2:years]/timeseries[1:(years-1)]
                lambda <- -(gm_mean(lambda) - 1)
                c(timeseries,lambda)
        }
        results
        #list(gmlambda = gm_lambda, results = results)
}
# Calculate the number of cores
no_cores <- detectCores() - 1

cl<-makeCluster(no_cores)
registerDoParallel(cl)

test <- simulatetrend(points = 10, psi = 0.9, gamma = 0, phi = 0.8, nsim = 100)

res <- data.frame(test)
slopes <- data.frame(slopes = test[,ncol(res)])
res <- res[,-ncol(res)]
res <- data.frame(res, id = 1:nrow(res))
names(res) <- c(1:(ncol(res)-1), "id")

res <- melt(res,c("id"))
names(res)[2] <- "year"
ggplot(res, aes(x=year,y=value, group=id))+geom_line(alpha=0.2)
ggplot(res, aes(x=year, y= value)) + geom_violin() + geom_line(aes(group=id),alpha=0.2)
ggplot(slopes,aes(x=slopes)) + geom_histogram() + geom_density()
