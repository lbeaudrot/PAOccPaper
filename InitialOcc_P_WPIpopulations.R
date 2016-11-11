# Determine the % population change for the 511 populations in the WPI paper using PowerSensor! input, Initial Occ and P values

# Load WPI results and initial occupancy probabilites calculated in "SimulationStudy_InitialOcc_P.R" in "WPI Paper" project
data <- read.csv(file="WPI_results_InitialOcc.year1.csv")

# Load detection probabilities from Eric (file=mean_detection_prob_v1_5_Oct_2016.csv)
p.data <- read.csv(file="mean_detection_prob_v1_5_Oct_2016.csv")
det.p <- data.frame(species_id = p.data$species_id, site_id = p.data$site_id, mdp=p.data$mdp)

# Create bins for analysis of detection probabilities, similar to what was done for initial occupancy
bins <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
p_bins <- .bincode(det.p$mdp, bins, right=TRUE)
det.p <- data.frame(det.p , p_bins)

# Merge initial occ and p probabilites with population result table BASED ON BOTH SPECIES ID AND SITE ID
merged <- merge(data, det.p, by=c("site_id", "species_id"))

# Create values for p_bins and year1_bins that will match with initial occ and detection probabilities in sims. 
# Categories are lower bounds (i.e. if p=0.1, then the detection probability was greater than or equal to 0.1, but less than 0.2)
p_cat <- ifelse(merged$p_bins==2, 0.1, 
       ifelse(merged$p_bins==3, 0.2,
       ifelse(merged$p_bins==4, 0.3,
       ifelse(merged$p_bins==5, 0.4,
       ifelse(merged$p_bins==6, 0.5,
       ifelse(merged$p_bins==7, 0.6,
       ifelse(merged$p_bins==8, 0.7,
       ifelse(merged$p_bins==9, 0.8,
       ifelse(merged$p_bins==10, 0.9, 0)))))))))

year1_cat <- ifelse(merged$p_bins==2, 0.1, 
       ifelse(merged$year1_bins==3, 0.2,
       ifelse(merged$year1_bins==4, 0.3,
       ifelse(merged$year1_bins==5, 0.4,
       ifelse(merged$year1_bins==6, 0.5,
       ifelse(merged$year1_bins==7, 0.6,
       ifelse(merged$year1_bins==8, 0.7,
       ifelse(merged$year1_bins==9, 0.8,
       ifelse(merged$year1_bins==10, 0.9, 0)))))))))

# Create a single column with p and initial occ values to match on.
tomatch <- paste(p_cat, year1_cat, sep=".")

merged <- data.frame(merged, p_cat, year1_cat, tomatch)

# Load output of simulation loop with values to extract
sims <- read.csv(file="det.year_2016-10-08.csv")
# Reduce to only include appropriate # of days and pts for each site
sims <- sims[,2:13]
sims <- sims[sims$days==30,]
sims <- sims[sims$pts==60,] # double check # of camera traps reported per site in WPI paper - treat PSH separately (30 CTs)






######## GENERALIZE CODE STARTING HERE TO EXTRACT PHI VALUE OF CHOICE (i.e. 0.99, 0.95, 0.90, 0.85)
#Alternatively, melt and cast data so that z.first for multipe psi values are in a table


# Limit to 15% declines to begin with
sims <- sims[sims$phi==0.85,]

Tomatch <- paste(sims$p,sims$psi1, sep=".")
sims <- data.frame(sims, Tomatch)

# 81 of the 511 populations in the WPI paper have initial occ and p values that we modeled in the similation study
merged[merged$tomatch %in% sims$Tomatch,]

# Assign percent decline in "sims" using merged$p_bins, and merged$year1bins (i.e. Initial occupancy probability) in "merged"
# Assign percent decline in "sims" using Tomatch, and tomatch (i.e. Detection . Initial occupancy probability in "merged") 
test <- merged[merged$tomatch %in% sims$Tomatch,]
test2 <- sims[match(test$tomatch,sims$Tomatch),]
test3 <- data.frame(test, test2)
test4 <- test3[test3$z.first80>0,]

data.frame(Site=test4$sitecode, Species=test4$bin, Status=test4$NewRare80, phi=test4$phi, Y1=test4$z.first80, Guild=test4$guild, RLS=test4$rls)

# Now need to connect the # of years of CT data per site in the WPI paper to identify minimum level of decline detected
