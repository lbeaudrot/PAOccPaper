############ TEAM CT DYNAMIC OCCUPANCY SIMULATION STUDY ##############

# Load functions for simulating data
source(file="dynoFunctions.R")

# Generate list of input values for simulation study where each item in the list has a vector of input values for
# function(points, days, psi,p ,phi ,gamma, years, nsim)
pts <- c(10, 20, 30, 60, 90, 120)
days <- c(15, 30, 45, 60)
psi <- c(0.5)
p <- c(0.08, 0.2, 0.5)
phi <- c(0.99, 0.95, 0.9, 0.85)
g <- 0
yrs <- c(12)
nsim <- c(500)

# Create object with all combinations of input parameters
mat <- expand.grid(pts, days, psi, p, phi, g, yrs, nsim)
colnames(mat) <- c("pts", "days", "psi", "p", "phi", "g", "yrs", "nsim")

# Create a string of seeds so that exact results can be obtained with parallel processing
# Seed needs to be set for each core to ensure repeatability; write file for records
#set.seed(400)
#seeds <- round(runif(nsim, 0, 10000), 0)
#write.csv(seeds, file="seeds_2016-06-22.csv")

# Setting up the machine for parallel processing
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoParallel(cl)
#clusterSetRNGStream(cl=cl, iseed=400)

# Run loop to generate simulations for each combination of parameters in object "mat"
result <- list()
for(m in 1:dim(mat)[1]){
#for(m in 240:dim(mat)[1]){
#for(m in 20:288){
    print(m)
    test <- simulatetrend(points=mat[m,1], days=mat[m,2], psi=mat[m,3], p=mat[m,4], phi=mat[m,5], gamma=mat[m,6], years=mat[m,7], nsim=mat[m,8])
    result[[m]] <- test
    names(result)[m] <- paste("input values", paste(mat[m,], collapse=", "))
}

hold <- vector()
  for(m in 1:dim(mat)[1]){
    hold[m] <- paste("input values", paste(mat[m,], collapse=", "))
  }


stopCluster(cl)

write.csv(melt.list(result), file="result.csv", row.names=TRUE)

# Below is code to combine different output files when melted and written to a file
#test <- read.csv(file="result_1to238.csv") 
#test2 <- read.csv(file="result_239.csv") 
#test3 <- read.csv(file="result_240to288.csv") 
#test4 <- rbind(test, test2, test3)

#val <- sprintf("%03d", 1:500)
#ID <- rep(val, dim(test4)[1]/length(val))
#test4 <- data.frame(test4, ID)
#test4$L1 <- factor(test4$L1)
#abcd <- cast(test4, ID ~ X2 ~ L1)
#test5 <- alply(.data=abcd, .margins=3, .fun="[")
#result <- test5
#hold2 <- levels(test4$L1)

# Change "mat" order of phi to correct order of phi found in levels of L1.
PHIextract <- strsplit(as.character(levels(test4$L1)), split=",")
reordered <- ldply(PHIextract, .fun="[")
names(reordered) <- c("pts", "days", "psi", "p", "phi", "g", "yrs", "nsim")
reordered$days <- as.numeric(reordered$days)
reordered$psi <- as.numeric(reordered$psi)
reordered$p <- as.numeric(reordered$p)
reordered$phi <- as.numeric(reordered$phi)

label <- vector()
  for(m in 1:dim(reordered)[1]){
    label[m] <- paste(reordered[m,], collapse=", ")
  }


# Use phi^t to calculate expected value to be used below to establish confidence intervals
phi.all <- data.frame(mat$phi^1, mat$phi^2, mat$phi^3, mat$phi^4, 
                      mat$phi^5, mat$phi^6, mat$phi^7, mat$phi^8, 
                      mat$phi^9, mat$phi^10, mat$phi^11, mat$phi^12, 
                      mat$phi^13, mat$phi^14)

# Use phi^t to calculate expected value to be used below to establish confidence intervals
phi.all <- data.frame(reordered$phi^1, reordered$phi^2, reordered$phi^3, reordered$phi^4, 
                      reordered$phi^5, reordered$phi^6, reordered$phi^7, reordered$phi^8, 
                      reordered$phi^9, reordered$phi^10, reordered$phi^11, reordered$phi^12, 
                      reordered$phi^13, reordered$phi^14)

# Run second loop to calculate changes for each set of simulation values

result2 <- list()
z.first <- vector()
y.first <- vector()

for(i in 1:length(result)){
#for(i in 240:length(result)){
    print(i)
    res <- data.frame(result[[i]])
    slopes <- data.frame(slopes=res[,ncol(res)])
    changes <- data.frame(res[,(ncol(res)/2+1):(ncol(res)-1)])
    res <- res[,-c((ncol(res)/2+1):(ncol(res)))]
    res <- data.frame(res, id = 1:nrow(res))
    names(res) <- c(1:(ncol(res)-1), "id")
    res <- melt(res,c("id"), variable_name = "year")
    result2[[i]] <- res
    
    # Need to correct the names given to result2 to match the correct order 
    names(result2)[i] <- names(result)[i] # only if running loop 2 directly after loop 1, otherwise names order differs
    #names(result2)[i] <- label[i] # only if running loop 2 after casting data to recreate result
    
    names(changes) <- c(1:ncol(changes))
    confint <- apply(changes,2,quantile,c(0.1,0.9))
    z <- which(round(confint[1,],2) > 0) 
    z.first[i] <- ifelse(length(z)==0,0,min(z))
    #y <- which(round(confint[1,],2) > phi.all[i,]) # y.first[i] <- the first year that confint contains phi, with phi^n
    y <- which(confint[1,] < phi.all[i,] & confint[2,] > phi.all[i,])
    y.first[i] <- ifelse(length(y)==0,0,min(y))
}


#det.year <- data.frame(y.first=y.first, z.first=z.first, mat, input=names(result)) # only if running loop 2 directly after loop 1, otherwise names order differs
det.year <- data.frame(y.first=y.first, z.first=z.first, reordered) # only if running loop 2 after casting data to recreate result
#colnames(det.year) <- c("y.first", "z.first", "pts", "days", "psi", "p", "phi", "g", "yrs", "nsim")
write.csv(det.year, file="det.year.csv", row.names=TRUE)
write.csv(melt.list(result2), file="result2.csv", row.names=TRUE)


