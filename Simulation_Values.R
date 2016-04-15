# Generate list of input values for simulation study where each item in the list has a vector of input values for
# function(points,days,psi,p,phi,gamma,years)

# Define parameter values of interest
pts <- c(10, 20, 30, 60, 90, 120)
days <- c(15, 30, 45, 60)
psi <- c(0.5)
p <- c(0.08, 0.2, 0.5)
phi <- c(0.99, 0.95, 0.9, 0.85)
g <- 0
yrs <- c(10)

mat <- expand.grid(pts, days, psi, p, phi, g, yrs)
colnames(mat) <- c("pts", "days", "psi", "p", "phi", "g", "yrs")
mat.list <- alply(.data=mat, .margins=1, fun="[")


# If we are running the same version of R, set.seed should produce similar output, although if R installation differed, they could be different
# see https://stat.ethz.ch/pipermail/r-help/2005-September/079391.html

# If parallelizing, we should NOT use set.seed, but should instead use the function clusterSetRNGStreatm() from the parallel package
# see http://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/

# It doesn't matter what number we use for the seed
# see http://www.r-bloggers.com/what-are-the-most-common-rng-seeds-used-in-r-scripts-on-github/
# see also http://www.r-bloggers.com/a-look-at-random-seeds-in-r-or-85-why-cant-you-be-more-like-548/