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

fit_data <- function(seedset, visitation, a_start = 0, b_start = 10, c_start = 0.5){
  nlmod3 <- nls(seedset ~  ((b*a)/100) + (b-((b*a)/100)) * (1-exp(-c*visitation)), 
                start = list(a = a_start, b = b_start,c = c_start),
                control= nls.control(maxiter = 1000))
  #summary(nlmod3)
  c <- coef(nlmod3)
  c
}

calculate_visits <- function(a, b, SVD, loss = 0.1, from_ = 0, to_ = 100, 
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
  poldep[1] <- a #0 visits
  poldep[2] <- SVD #1 vist
  for(i in 3:to_){
    lost <- (1-loss)^(i-2)
    poldep[i] <- SVD*lost
  }
  poldep[to_+1] <- SVD*((1-loss)^(to_-2)) #+1 vist
  #plot(poldep, las = 1)
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

#calculate_visits(a = 0, b = 100, SVD = 10, loss = 0.1, from_ = 0, to_ = 100, 
#                             plot_ = TRUE, col_ = 2)
  
#calculate_visits(a = 20, b = 100, SVD = 10, loss = 0.1, from_ = 0, to_ = 100, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 20, b = 100, SVD = 30, loss = 0.1, from_ = 0, to_ = 100, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 20, b = 100, SVD = 30, loss = 0.3, from_ = 0, to_ = 100, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 20, b = 90, SVD = 6, loss = 0.1, from_ = 0, to_ = 100, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 20, b = 90, SVD = 5, loss = 0.1, from_ = 0, to_ = 100, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 0, b = 10, SVD = 6, loss = 0.5, from_ = 0, to_ = 10, 
#                 plot_ = TRUE, col_ = 2)

#calculate_visits(a = 0, b = 10, SVD = 6, loss = 0.6, from_ = 0, to_ = 10, 
#                 plot_ = TRUE, col_ = 2)

#Next function can give you visitation rate to achibe 95% of seed set.
#Next function can make elasticity analysis
