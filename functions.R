#functions needed

V_transect_to_flower <- function(V_transect, flw_x_m2, lifespan){
  V_flowers <- V_transect * (100/flw_x_m2) * lifespan
  V_flowers
}

#V_transect_to_flower(V_transect = 50, flw_x_m2 = 100, lifespan = 8)

plot_visits <- function(a, b, c, from_ = 0, to_ = 1, add_ = FALSE, col_ = 2){
  a2 <- (b*a)/100
  b2 <- b-a2
  curve(a2+b2*(1-exp(-c*x)), from = from_, to = to_, add = add_, col = col_, ylim = c(0, b), las = 1)
}

#plot_visits(a = 30, b = 200, c = 10, add = FALSE)

fit_data <- function(seedset, visitation, a_start = 0, b_start = 10, c_start = 0.5,
                     simplify = "params"){
  nlmod3 <- nls(seedset ~  ((b*a)/100) + (b-((b*a)/100)) * (1-exp(-c*visitation)), 
                start = list(a = a_start, b = b_start,c = c_start),
                control= nls.control(maxiter = 1000))
  #summary(nlmod3)
  if(simplify == "params"){
    out <- coef(nlmod3)
    }
  if(simplify == "stats"){
    out <- summary(nlmod3)$coefficients
  }
  if(simplify == "model"){
    out <- nlmod3
  }
  out
}

calculate_visits0 <- function(a, b, SVD, loss = 0.1, to_ = 100, 
                             plot_ = FALSE, col_ = 2){
  #set a
  if(!is.numeric(a)){"a must be numeric"}
  #set b
  if(!is.numeric(b)){"b must be numeric"}
  #force from_ = 0 because the way this is written.
  from_ = 0
  #create a data.frame of accumulated pollen dep per visits
  visits <- c(from_:to_)
  poldep <- rep(NA, length(visits))
  poldep[1] <- (b*a)/100 #0 visits = %of seed set without visits.
  poldep[2] <- SVD #1 visit (assuming selfing rates are not discounted form the empirical calculation)
  for(i in 3:to_){
    lost <- (1-loss)^(i-2)
    poldep[i] <- SVD*lost
  }
  poldep[to_+1] <- SVD*((1-loss)^(to_-2)) #+1 visit
  #plot(poldep, las = 1)
  #cumulative:
  poldep2 <- rep(NA, length(visits))
  for(j in from_:to_){
    poldep2[j+1] <- sum(poldep[1:(j+1)])
  }
  #plot(poldep2, las = 1)
  #fit c
  nlmod3 <- nls(poldep2 ~  ((b*a)/100) + (b-((b*a)/100)) * (1-exp(-c*visits)), start = list(c = 0.2),
                control= nls.control(maxiter = 1000))
  #summary(nlmod3)
  c <- coef(nlmod3)
  #plot it
  if(plot_ == TRUE){
    a2 <- (b*a)/100
    b2 <- b-a2
    plot(poldep2 ~ visits, las = 1,  xlab = "visits",
         ylab = "seed set", xlim = c(0,to_), ylim = c(0, max(c(poldep2, b))))
    curve(a2+b2*(1-exp(-c*x)), from = from_, to = to_, add = TRUE, 
          col = col_, ylim = c(0, max(c(poldep2, b))), las = 1)
  }
  c
}

#A way to ensure loss (the diminishing returns function) aproximates b is to calculate

loss <- function(a, b, SVD){
  SVD/(b-(b*a/100))
}

#calculate_visits0(a = 20, b = 100, SVD = 50, loss = 0.5, to_ = 100, 
 #                plot_ = TRUE, col_ = 2)
#loss(a = 20, b = 100, SVD = 50)  
#calculate_visits0(a = 20, b = 100, SVD = 50, loss = 0.625, to_ = 100, 
 #                plot_ = TRUE, col_ = 2)

#to automatically calculate loss:

calculate_visits <- function(a, b, SVD, to_ = 100, 
                             plot_ = FALSE, col_ = 2){
  baseline <- (b * a) / 100   # Calculate the baseline (autonomous selfing)
  target_gap <- b - baseline  # Calculate the "Target Gap" pollinators need to fill
  SVD <- SVD - baseline #discounting selfing from the empiricall estimation.
  calculated_loss <- SVD / target_gap   # Calculate the perfect loss value
  # Safety check: loss cannot be > 1 or the math breaks
  if(calculated_loss > 1) {
    stop("SVD is higher than the remaining gap; total seed set will exceed 'b'.")
  }
  #create a dataframe to store accumulated pollen
  visits <- c(0:to_)
  poldep <- rep(NA, length(visits))
  poldep[1] <- baseline
  for(i in 2:length(poldep)){
    # This simulates the decay towards the limit b
    poldep[i] <- SVD * (1 - calculated_loss)^(i - 2)
  }
  poldep2 <- cumsum(poldep)
  #fit c
  #calculate start values
  c_start <- -log(1 - (SVD / (b - baseline)))
  if(is.na(c_start) | is.infinite(c_start)) c_start <- 0.1
  nlmod3 <- nls(poldep2 ~  ((b*a)/100) + (b-((b*a)/100)) * (1-exp(-c*visits)), 
                start = list(c = c_start),
                control= nls.control(maxiter = 2000, minFactor = 1/2048, warnOnly = TRUE))
  #summary(nlmod3)
  c <- coef(nlmod3)
  c_se <- summary(nlmod3)$coefficients[2]
  #plot it
  if(plot_ == TRUE){
    a2 <- (b*a)/100
    b2 <- b-a2
    plot(poldep2 ~ visits, las = 1,  xlab = "visits",
         ylab = "seed set", xlim = c(0,to_), ylim = c(0, max(c(poldep2, b))))
    curve(a2+b2*(1-exp(-c*x)), from = from_, to = to_, add = TRUE, 
          col = col_, ylim = c(0, max(c(poldep2, b))), las = 1)
  }
  return(list(optimal_loss = calculated_loss, c_parameter = c, c_se = c_se))
}

#calculate_visits(a = 20, b = 100, SVD = 50, to_ = 100, 
 #                            plot_ = TRUE, col_ = 2)
  
#Next function can give you visitation rate to achieve 95% of seed set.

calculate_required_visits <- function(c_val, b, a, target_percent = 0.95) {
  #Calculate the baseline (autonomous selfing)
  baseline <- (b * a) / 100
  #Define the target seed set (e.g., 95% of b)
  target_seeds <- b * target_percent
  #Check if baseline already exceeds target
  if (baseline >= target_seeds) {
    return(0) # Zero visits needed if selfing covers it
  }
  #solving for x (visitation)
  #Rearranging: (target - baseline) / (b - baseline) = 1 - exp(-c*x)
  numerator <- target_seeds - baseline
  denominator <- b - baseline
  visits_needed <- as.numeric(-log(1 - (numerator / denominator)) / c_val)
  return(visits_needed)
}

#my_c <- calculate_visits(a=20, b=100, SVD=35)
#calculate_required_visits(c_val = my_c$c_parameter, b = 100, a = 20)
#abline(v=4.8)

pseudeR2 <- function(model, seedset){
  RSS <- sum(residuals(model)^2)
  TSS <- sum((seedset - mean(seedset))^2)
  pseudo_R2 <- 1 - (RSS / TSS) 
  pseudo_R2
}
