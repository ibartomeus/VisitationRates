#This script uses data extracte DiTrani et al 2024 to recreate a 
#plausible distribution of raw data.

#Original data
wm <- read.csv(file = "data/Watermelon_Visitation.csv")
head(wm)

#Reconstrubt plausible distrinution
#ensuring non negative values...
l <- -1
while(any(unlist(l)<0)){
  l <- list()
  for(i in 1:13){
    l[[i]] <- rnorm(n = wm$n[i], mean = wm$seeds[i], sd = wm$SE[i]*sqrt(wm$n[i]))
  }
}

#Recreate the dataset
dat <- data.frame(visits = rep(wm$visits, wm$n), 
                  species = rep(wm$species, wm$n),
                  seeds = round(unlist(l)))
dat

#Save it. Data is added to the main VisitationRates_Data.csv
write.csv(dat, "data/watermelon_expnded.csv")
