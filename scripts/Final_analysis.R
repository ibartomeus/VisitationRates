#Final analysis

#Source functions. Should be made into a package.
# install.packages("remotes")
#remotes::install_github("ibartomeus/VisitationR")
library(VisitationR)

#Fig1-----
c <- calculate_visits(a = 20, b = 100, SVD = 40, to_ = 15, 
                      plot_ = TRUE, col_ = 2)
ovr <- calculate_required_visits(a = 20, b = 100, c = c$c_parameter)
abline(h = 100, lty = 2, col = "brown")
lines(y = c(-2,95), x = c(ovr,ovr), lty = 2, col = "blue")
lines(y = c(-2,40), x = c(1,1), lty = 2)
points(x = 1, y = 40, pch = 19)
text(x = 2, y = 20, "a = selfing rate", col = "darkgrey", cex = 0.7)
text(x = 4.3, y = 40, "SVE = single visit efficiency", cex = 0.7)
text(x = 6, y = 58, "c = diminishing returns param.", col = "red", cex = 0.7)
text(x = 2, y = 95, "b = maximum seed set", col = "brown", cex = 0.7)
text(x = 12, y = 80, "OVR = optimal visits", col = "blue", cex = 0.7)

#Read Data----
#read observational data
dat <- read.csv("data/VisitationRates_Data.csv")

#First turn BB weight into seeds follwing Bernardine and Vance equation

WtoS <- function(w){
  ifelse(w == 0, 0, (w-0.09)/0.06)
}

summary(SHighbushBlueberry$PollinationValue)

#make a dataset per crop
#unique(dat$Crop)
#LimaBean <- subset(dat, Crop == "Lima.bean")
Strawberry <- subset(dat, Crop == "Strawberry") 
SHighbushBlueberry <- subset(dat, Crop == "Southern.Highbush.Blueberry")  
SHighbushBlueberry$PollinationValue <- WtoS(SHighbushBlueberry$PollinationValue)
NHighbushBlueberry <- subset(dat, Crop == "Northern.Highbush.Blueberry")
NHighbushBlueberry$PollinationValue <- WtoS(NHighbushBlueberry$PollinationValue)
Rabbiteye.Blueberry <- subset(dat, Crop == "Rabbiteye.Blueberry")
Rabbiteye.Blueberry$PollinationValue <- WtoS(Rabbiteye.Blueberry$PollinationValue)
Apple <- subset(dat, Crop == "Apple")
Watermelon <- subset(dat, Crop == "Watermelon")

#mean(subset(Watermelon, Watermelon$Number.of.visits == 30)[,"PollinationValue"])

#Fig 3-----

#AIC and R2 in Real_crops.qmd

#apple
#fit models
model_nls <- fit_data(seedset = Apple$PollinationValue, visitation = Apple$Number.of.visits, simplify = "model")
model_lm <- lm(Apple$PollinationValue ~ Apple$Number.of.visits)
model_cuadratic <- lm(PollinationValue ~ Number.of.visits + I(Number.of.visits^2), data = Apple)

#Extract a, b and c values
mf <- fit_data(seedset = Apple$PollinationValue, visitation = Apple$Number.of.visits, simplify = "params")

#Plot
par(mfrow=c(2,2))
plot(Apple$PollinationValue ~ jitter(Apple$Number.of.visits), las = 1, xlab = "Visits", ylab = "Seed set", col = "grey", main = "Apple")
plot_visits(mf["a"], mf["b"], mf["c"], to_ = max(Apple$Number.of.visits), add_ = TRUE, col = "#E69F00")
abline(model_lm, col = "#56B4E9")
curve(predict(model_cuadratic, newdata = data.frame(Number.of.visits = x)), add = TRUE, col = "#009E73")

#Strawberry
#Strawberry <- subset(Strawberry, Pollinator %in% c("Apis mellifera", #"Andrena.spp", "Ceratina.spp"))

#fit models
model_nls <- fit_data(seedset = Strawberry$PollinationValue, visitation = Strawberry$Number.of.visits, simplify = "model",
                      a_start = 10, b_start = 100, c_start = 0.5)
model_lm <- lm(Strawberry$PollinationValue ~ Strawberry$Number.of.visits)
model_cuadratic <- lm(PollinationValue ~ Number.of.visits + I(Number.of.visits^2), data = Strawberry)

#Extract a, b and c values
mf <- fit_data(seedset = Strawberry$PollinationValue, visitation = Strawberry$Number.of.visits, simplify = "params", a_start = 10, b_start = 100, c_start = 0.5)

#Plot
plot(Strawberry$PollinationValue ~ jitter(Strawberry$Number.of.visits), las = 1, xlab = "Visits", ylab = "Seed set", col = "grey", xlim = c(0,20), main = "Strawberry")
plot_visits(mf["a"], mf["b"], mf["c"], to_ = max(Strawberry$Number.of.visits), add_ = TRUE, col = "#E69F00")
plot_visits(mf["a"], mf["b"], mf["c"], to_ = 30, add_ = TRUE, col = "#E69F00")
abline(model_lm, col = "#56B4E9")
curve(predict(model_cuadratic, newdata = data.frame(Number.of.visits = x)), 
      add = TRUE, col = "#009E73")

#Blueberry
SHighbushBlueberry <- subset(SHighbushBlueberry, PollinationMetric == "Fruit.weight")

#fit models
model_nls <- fit_data(seedset = SHighbushBlueberry$PollinationValue, 
                      visitation = SHighbushBlueberry$Number.of.visits, simplify = "model")
model_lm <- lm(SHighbushBlueberry$PollinationValue ~ SHighbushBlueberry$Number.of.visits)
model_cuadratic <- lm(PollinationValue ~ Number.of.visits + I(Number.of.visits^2), data = SHighbushBlueberry)

#Extract a, b and c values
mf <- fit_data(seedset = SHighbushBlueberry$PollinationValue, 
               visitation = SHighbushBlueberry$Number.of.visits, simplify = "params")

#Plot
plot(SHighbushBlueberry$PollinationValue ~ jitter(SHighbushBlueberry$Number.of.visits), 
     las = 1, xlab = "Visits", ylab = "Seed set", col = "grey", main = "Blueberry")
plot_visits(mf["a"], mf["b"], mf["c"], to_ = max(SHighbushBlueberry$Number.of.visits), 
            add_ = TRUE, col = "#E69F00")
abline(model_lm, col = "#56B4E9")
curve(predict(model_cuadratic, newdata = data.frame(Number.of.visits = x)), add = TRUE, col = "#009E73")

#Watermelon
#fit models
Watermelon <- subset(Watermelon, Watermelon$Number.of.visits != 30)

model_nls <- fit_data(seedset = Watermelon$PollinationValue, 
                      visitation = Watermelon$Number.of.visits, simplify = "model", a_start = 75, b_start = 300)
model_lm <- lm(Watermelon$PollinationValue ~ Watermelon$Number.of.visits)
model_cuadratic <- lm(PollinationValue ~ Number.of.visits + I(Number.of.visits^2), data = Watermelon)

#Extract a, b and c values
mf <- fit_data(seedset = Watermelon$PollinationValue, 
               visitation = Watermelon$Number.of.visits, simplify = "params", 
               a_start = 75, b_start = 300)

#Plot
plot(Watermelon$PollinationValue ~ jitter(Watermelon$Number.of.visits), 
     las = 1, xlab = "Visits", ylab = "Seed set", xlim = c(0,30), col = "grey", main = "Watermelon")
plot_visits(mf["a"], mf["b"], mf["c"], to_ = 30, 
            add_ = TRUE, col = "#E69F00")
abline(model_lm, col = "#56B4E9")
curve(predict(model_cuadratic, newdata = data.frame(Number.of.visits = x)), add = TRUE, col = "#009E73")
#Add legend with colors.
par(mfrow=c(1,1))

#Fig 3----

#read data
dat <- read.csv("data/VisitationRates_ParameterData.csv")
#head(dat)
dat <- dat[-24,] #removing cheak peas, which is wierd as a = b = SV
#rename params
dat$a <- dat$a_seed_set
dat$b <- dat$b_ovules

#param a (literature) need to be in %
dat$a_emp <- (dat$a / dat$b)*100

#param a as in Sipoa (already in %)
dat$a_siopa_seed_set <- dat$b * (1 - dat$a_siopa)
#remove NA's
#dat <- dat[-which(is.na(dat$a)),]


#SV we can take the average
dat$SV <- rowMeans(dat[,c("Honeybee", "Bumblebee", 
                          "Other.wild.bee", "Hoverfly", "Other")], na.rm = TRUE)
#But some are wierd
dat$SV[which(dat$Crop == "cherry")] <- 0.25 #Using only wild bees.
dat$SV[which(dat$Crop == "strawberry")] <- 150 #Using max efficiency...

#plot siopa vs SV
#par(mfrow=c(1,2))
plot(x = dat$a_siopa_seed_set, y = dat$SV, xlim = c(0,950), ylim = c(0,950), xlab = "a (Siopa)", ylab = "SVE", las = 1)
abline(0, 1)

plot(x = dat$a_siopa_seed_set, y = dat$SV, xlim = c(0,100), ylim = c(0,100), xlab = "a (Siopa)", ylab = "SVE", las = 1)
abline(0, 1)
#par(mfrow=c(1,1))

#length(which(dat$a_siopa_seed_set > dat$SV))
#11 crops show larger selfing rates than SVE...

#create final a
dat$a <- dat$a_emp

#output data frame
out <- data.frame(crop = rep(NA,nrow(dat)), 
                  family = rep(NA,nrow(dat)),
                  a = rep(NA,nrow(dat)), 
                  b = rep(NA,nrow(dat)), 
                  c = rep(NA,nrow(dat)), 
                  SVE = rep(NA,nrow(dat)), 
                  OVR = rep(NA,nrow(dat)))
#loop it
dat$Crop[23] <- "lupin"
for(i in 1:nrow(dat)){
  #extract c
  params <- calculate_visits(a = dat$a[i], b = dat$b[i], SVD = dat$SV[i], 
                             plot_ = FALSE)
  #warning can be dissmissed.
  #And visits needed.
  OVR_model <- calculate_required_visits(a = dat$a[i], b = dat$b[i], c_val = params$c_parameter, target_percent = 0.8) 
  #store the data
  out$crop[i] <- dat$Crop[i]
  out$family[i] <- dat$family[i]
  out$a[i] <- dat$a[i]
  out$b[i] <- dat$b[i]
  out$c[i] <- params$c_parameter
  out$SVE[i] <- dat$SV[i]
  out$OVR[i] <- OVR_model
}

par(mfrow=c(5,5), mar=c(3, 4, 1, 1))
for(i in 1:25){
  plot_visits(a = out$a[i], b = out$b[i], c = out$c[i], 
              xlab = "", ylab = "",
              to_ = 20, main = out$crop[i])
  abline(v = out$OVR[i], lty = 2, col = "grey")
}
par(mfrow=c(1,1), mar=c(4, 4, 1.5, 0.5))
mtext('Visits', side=1, line=3, at=6, cex = 0.9)
mtext('Seeds', side=2, line=3, at=20, cex = 0.9)


par(mfrow=c(1,1), mar=c(7, 4.2, 1.5, 0.5))
#bymedian <- with(out, reorder(family, -OVR, median))
#boxplot(out$OVR ~ bymedian, las = 2, xlab = "", ylab = "OVR")
bymedian <- with(out, reorder(crop, -OVR, median))
boxplot(out$OVR ~ bymedian, las = 2, xlab = "", ylab = "Optimal visitation rate",
        cex.axis = 0.75)
par(mfrow=c(1,1), mar=c(4, 4, 1.5, 0.5))
abline(h = 5, lty = 2)
abline(h = 10, lty = 2)
#parameter influence? To SUp Mat?
#scatter.smooth(out$OVR ~ out$a, las = 1)
#scatter.smooth(out$OVR ~ out$b, las = 1)
#scatter.smooth(out$OVR ~ out$c, las = 1)
#scatter.smooth(out$OVR ~ out$SVE, las = 1)
#It's a combination, but c is equivalent to OVR 


#Fig4----

#load params
all_data <- read.csv("data/parameters_real_crops.csv")
#head(all_data)
all_data <- all_data[-c(5,6),] #removing BB varieties

#Detour: Obtain params from Fig 1
#all_data$crop #"apple", "strawberry", "lima bean","SHighbushBlueberry", "Watermelon" 
#Extracted from raw data in Real_crops
spec_dat <- data.frame(a_specific = c(0.08, 64.6, 0.45, 0.35, 0),
                       b_specific = c(10, 199, 3, 1.9, 338),
                       SVE_specific = c(5.6, 74.1, 0.54, 0.84, 72.8))

spec_dat$a <- (spec_dat$a_specific / spec_dat$b_specific)*100
out2 <- data.frame(a = rep(NA,nrow(spec_dat)), 
                   b = rep(NA,nrow(spec_dat)), 
                   c = rep(NA,nrow(spec_dat)), 
                   SVE = rep(NA,nrow(spec_dat)), 
                   OVR = rep(NA,nrow(spec_dat)))
for(i in 1:nrow(spec_dat)){
  #extract c
  params <- calculate_visits(a = spec_dat$a[i], b = spec_dat$b_specific[i], SVD = spec_dat$SVE_specific[i], 
                             plot_ = FALSE)
  #And visits needed.
  OVR_model <- calculate_required_visits(a = spec_dat$a[i], b = spec_dat$b[i], c_val = params$c_parameter, target_percent = 0.80) 
  #store the data
  out2$a[i] <- spec_dat$a[i]
  out2$b[i] <- spec_dat$b_specific[i]
  out2$c[i] <- params$c_parameter
  out2$SVE[i] <- spec_dat$SVE_specific[i]
  out2$OVR[i] <- OVR_model
}
par(mfrow=c(2,2), mar=c(4, 4, 1.5, 0.5))
plot(all_data$a[-3], out2$a[-3], 
     ylim = c(0,max(all_data$a,out2$a)+5), 
     xlim = c(0,max(all_data$a,out2$a)+5),  
     xlab = "Obs_a", ylab = "Model_a",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)
plot(all_data$b[-3], out2$b[-3], 
     ylim = c(0,max(all_data$b,out2$b)+5), 
     xlim = c(0,max(all_data$b,out2$b)+5), 
     xlab = "Obs_b", ylab = "Model_b",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)

plot(all_data$c[-3], all_data$c_model[-3], 
     ylim = c(0,max(all_data$c,out2$c)+1), 
     xlim = c(0,max(all_data$c,out2$c)+1),  
     xlab = "Obs_c", ylab = "Model_c",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)
plot(all_data$OVR[-3], out2$OVR[-3], 
     #ylim = c(0,max(all_data$OVR,out2$OVR)+5), 
     #xlim = c(0,max(all_data$OVR,out2$OVR)+5),   
     ylim = c(0,20), 
     xlim = c(0,20),  
     xlab = "Obs_OVR", ylab = "Model_OVR", 
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)
par(mfrow=c(1,1), mar=c(4, 4, 1.5, 0.5))
#######END detour

#Sup mat S1
par(mfrow=c(2,2), mar=c(4, 4, 1, 0.5))
plot(all_data$a[-3], all_data$a_model[-3], 
     ylim = c(0,40), 
     xlim = c(0,40),  
     xlab = "Obs_a", ylab = "Model_a",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)

plot(all_data$b[-3], all_data$b_model[-3], 
     ylim = c(0,max(all_data$b,all_data$b_model)+5), 
     xlim = c(0,max(all_data$b,all_data$b_model)+5), 
     xlab = "Obs_b", ylab = "Model_b",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)

plot(all_data$c[-3], all_data$c_model[-3], 
     ylim = c(0,max(all_data$c,all_data$c_model)), 
     xlim = c(0,max(all_data$c,all_data$c_model)),  
     xlab = "Obs_c", ylab = "Model_c",
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)
plot(all_data$OVR[-3], all_data$OVR_model[-3], 
     ylim = c(0,max(all_data$OVR,all_data$OVR_model)), 
     xlim = c(0,max(all_data$OVR,all_data$OVR_model)),   
     xlab = "Obs_OVR", ylab = "Model_OVR", 
     col = c("green", "red", "blue","darkgreen"), las = 1)
abline(0, 1)
#Not too satisfactory... specially strawberry, but also watermelon and lima bean
par(mfrow=c(1,1), mar=c(4, 4, 1.5, 0.5))


#Fig 5-------

#read data
dat <- read.csv("data/VisitationRates_Data.csv")
#make a dataset per crop
SHighbushBlueberry <- subset(dat, Crop == "Southern.Highbush.Blueberry")  
SHighbushBlueberry <- subset(SHighbushBlueberry, PollinationMetric == "Fruit.weight")
SHighbushBlueberry$PollinationValue <- WtoS(SHighbushBlueberry$PollinationValue)
NHighbushBlueberry <- subset(dat, Crop == "Northern.Highbush.Blueberry")
NHighbushBlueberry <- subset(NHighbushBlueberry, PollinationMetric == "Fruit.weight")
NHighbushBlueberry$PollinationValue <- WtoS(NHighbushBlueberry$PollinationValue)
Rabbiteye.Blueberry <- subset(dat, Crop == "Rabbiteye.Blueberry")
Rabbiteye.Blueberry <- subset(Rabbiteye.Blueberry, PollinationMetric == "Fruit.weight")
Rabbiteye.Blueberry$PollinationValue <- WtoS(Rabbiteye.Blueberry$PollinationValue)

#We can use quantiles, but mean of 0 and and max visits, or simply min and max, might also work.

#For SHB
SHB_a_real <- quantile(subset(SHighbushBlueberry, 
                              Number.of.visits == 0)$PollinationValue)[2]
#max(SHighbushBlueberry$Number.of.visits)
SHB_b <- quantile(subset(SHighbushBlueberry, 
                         Number.of.visits == 15)$PollinationValue)[4]
SHB_SVE <- mean(subset(SHighbushBlueberry, 
                       Number.of.visits == 1)$PollinationValue)
#We need a in %
SHB_a <- (SHB_a_real/ SHB_b)*100

#For NHB
NHB_a_real <- quantile(subset(NHighbushBlueberry, 
                              Number.of.visits == 0)$PollinationValue)[2]
#max(NHighbushBlueberry$Number.of.visits)
NHB_b <- quantile(subset(NHighbushBlueberry, 
                         Number.of.visits == 15)$PollinationValue)[4]
NHB_SVE <- mean(subset(NHighbushBlueberry, 
                       Number.of.visits == 1)$PollinationValue)
#We need a in %
NHB_a <- (NHB_a_real/ NHB_b)*100

#For REB
REB_a_real <- quantile(subset(Rabbiteye.Blueberry, 
                              Number.of.visits == 0)$PollinationValue)[2]
#max(Rabbiteye.Blueberry$Number.of.visits)
REB_b <- quantile(subset(Rabbiteye.Blueberry, 
                         Number.of.visits == 10)$PollinationValue)[4]
REB_SVE <- mean(subset(Rabbiteye.Blueberry, 
                       Number.of.visits == 1)$PollinationValue)
#We need a in %
REB_a <- (REB_a_real/ REB_b)*100

#calculate param c
SHB_c <- calculate_visits(a = SHB_a, b = SHB_b, SVD = SHB_SVE, to_ = 15)$c_parameter
NHB_c <- calculate_visits(a = NHB_a, b = NHB_b, SVD = NHB_SVE, to_ = 15)$c_parameter
REB_c <- calculate_visits(a = REB_a, b = REB_b, SVD = REB_SVE, to_ = 15)$c_parameter

#now we can plot all three
par(mfrow=c(1,2))
plot(SHighbushBlueberry$PollinationValue ~ jitter(SHighbushBlueberry$Number.of.visits), 
     las = 1, xlab = "Visits", ylab = "Seed set", type = "n", main = "Different blueberry crops",
     ylim = c(0,50), xlim = c(0,15))
plot_visits(SHB_a, SHB_b, SHB_c, to_ = 15, 
            add_ = TRUE, col = "#223366")
plot_visits(NHB_a, NHB_b, NHB_c, to_ = 15, 
            add_ = TRUE, col = "#A3B18A")
plot_visits(REB_a, REB_b, REB_c, to_ = 15, 
            add_ = TRUE, col = "#B9D1EA")

#Now different pollinator compositions.
#For SHB
SHB_a_real <- quantile(subset(SHighbushBlueberry, 
                              Number.of.visits == 0)$PollinationValue)[2]
#We need a in %
SHB_b <- quantile(subset(SHighbushBlueberry, 
                         Number.of.visits == 10)$PollinationValue)[4]
SHB_a <- (SHB_a_real/ SHB_b)*100

#pollinators available
#unique(SHighbushBlueberry$Pollinator)
SHBmix <- subset(SHighbushBlueberry, Pollinator == "pollinator.mix")
SHBhb <- subset(SHighbushBlueberry, Pollinator == "Apis.mellifera")
SHBbb <- subset(SHighbushBlueberry, Pollinator == "Tetragonula.carbonaria")

#calculate different SVD values
#SHB_SVE_mix <- mean(subset(SHBmix, 
#                           Number.of.visits == 2)$PollinationValue)/2 #diviing by 2 penalizes a lot... see below
SHB_SVE_hb <- mean(subset(SHBhb, 
                          Number.of.visits == 1)$PollinationValue)
SHB_SVE_bb <- mean(subset(SHBbb, 
                          Number.of.visits == 1)$PollinationValue)
SHB_SVE_mix <- (SHB_SVE_hb + SHB_SVE_bb)/2


#calculate param c
SHB_c_mix <- calculate_visits(a = SHB_a, b = SHB_b, SVD = SHB_SVE_mix, to_ = 15)$c_parameter
SHB_c_hb <- calculate_visits(a = SHB_a, b = SHB_b, SVD = SHB_SVE_hb, to_ = 15)$c_parameter
SHB_c_bb <- calculate_visits(a = SHB_a, b = SHB_b, SVD = SHB_SVE_bb, to_ = 15)$c_parameter
#now we can plot all three
plot(SHighbushBlueberry$PollinationValue ~ jitter(SHighbushBlueberry$Number.of.visits), 
     las = 1, xlab = "Visits", ylab = "Seed set", type = "n", main = "Pollinator composition",
     ylim = c(0,30), xlim = c(0,10))
plot_visits(SHB_a, SHB_b, SHB_c_mix, to_ = 10, 
            add_ = TRUE, col = "#D55E00")
plot_visits(SHB_a, SHB_b, SHB_c_hb, to_ = 10, 
            add_ = TRUE, col = "#E69F00")
plot_visits(SHB_a, SHB_b, SHB_c_bb, to_ = 10, 
            add_ = TRUE, col = "darkred")

par(mfrow=c(1,1))

