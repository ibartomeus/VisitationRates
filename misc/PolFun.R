#load Lucas data

d <- read.table("Other/Sapes/Database_S1.txt", header = TRUE, dec = ",")
str(d)

#check ingolf data
levels(d$system_ID)

dd <- subset(d, system_ID == "Strawberry")
plot(dd$fruit_set ~ dd$visits_wild_insects, main = "strawberry")
abline(lm(dd$fruit_set ~ dd$visits_wild_insects))


#Sapes

plot(d$fruit_set ~ d$visits_wild_insects, col = d$system_ID)
scatter.smooth(d$fruit_set ~ d$visits_wild_insects+d2$visits_honey_bees)
scatter.smooth(d$fruit_set ~ d$visits_wild_insects)

unique(d$system_ID)

dd <- subset(d, is.na(d$pollen_dep) == TRUE)
str(dd)
unique(dd$system_ID)

d2 <- subset(d, system_ID %in% c("Buckwheat_A", "Buckwheat_B", 
                                "Red_clover", "Tumip_rape", "Springrape", "Strawberry", "Sunflower"))
scatter.smooth(d2$fruit_set ~ d2$visits_wild_insects+d2$visits_honey_bees, col = d$system_ID)
scatter.smooth(d2$fruit_set ~ d2$visits_wild_insects, col = d2$system_ID)

d2 <- subset(d, system_ID %in% c("Buckwheat_A", "Buckwheat_B"))
d2 <- subset(d, system_ID %in% c("Red_clover"))
#d2 <- subset(d, system_ID %in% c("Tumip_rape"))
d2 <- subset(d, system_ID %in% c("Springrape"))
d2 <- subset(d, system_ID %in% c("Strawberry"))
d2 <- subset(d, system_ID %in% c("Sunflower"))
scatter.smooth(d2$fruit_set ~ d2$visits_wild_insects+d2$visits_honey_bees, col = d2$system_ID)
scatter.smooth(d2$fruit_set ~ d2$visits_wild_insects, col = d2$system_ID)

#define a function where 

#from morris, using the same idea, but different notation:
#y is benefit, x is the number of visits by all taxa combined, 
#p1 (a) is the y-intercept (i.e., reproductive success for unvisited [V 1⁄4 0] flowers),
#p2 (b) is the asymptote (i.e., B when V 1⁄4 ‘), 
#and p3 (c) governs the rate of approach to the asymptote. 

# a = % of fruit set witout pollinators
# b = Max yield attained on the region
# c = Shape of the curve?

pol <- function(a, b, c, from_ = 0, to_ = 1, add_ = FALSE, col_ = 2){
  a2 <- (b*a)/100
  b2 <- b-a2
  curve(a2+b2*(1-exp(-c*x)), from = from_, to = to_, add = add_, col = col_, ylim = c(0, b))
}

pol(a = 30, b = 200, c = 10, add = FALSE)
pol(a = 30, b = 200, c = 0.4, add = FALSE)
for (i in 2:10){
  pol(a = 30, b = 200, c = i, add = TRUE)
}

pol(a = 0, b = 200, c = 5, add = FALSE)
for (i in seq(10,100,10)){
  pol(a = i, b = 200, c = 5, add = TRUE)
}

pol(a = 30, b = 1000, c = 5, add = FALSE)
for (i in seq(150,1000,50)){
  pol(a = 30, b = i, c = 5, add = TRUE)
}


pol2 <- function(a, b, c, from_ = 0, to_ = 1, add_ = FALSE, col_ = 2){
  a2 <- (b*a)/100
  curve(b /(1 + (b/a2 -1) * exp(-c*x)), from = from_, to = to_, add = add_, col = col_, ylim = c(0, b))
}

pol2(a = 30, b = 200, c = 10, add = FALSE)
for (i in 2:10){
  pol2(a = 30, b = 200, c = i, add = TRUE)
}

pol2(a = 0, b = 200, c = 5, add = FALSE)
for (i in seq(10,100,10)){
  pol2(a = i, b = 200, c = 5, add = TRUE)
}

pol2(a = 30, b = 1000, c = 5, add = FALSE)
for (i in seq(150,1000,50)){
  pol2(a = 30, b = i, c = 5, add = TRUE)
}


#fit lucas data via nls


y <- d$fruit_set + 3
x <- d$visits_wild_insects + d$visits_honey_bees + 3
x <- d$visits_wild_insects + 3 #used
y <- d2$fruit_set + 3
x <- d2$visits_wild_insects + d2$visits_honey_bees + 3 #used
x <- d2$visits_wild_insects + 3 
scatter.smooth(y ~ x, xlab = "visits", ylab = "fruit set")
nlmod3 <- nls(y ~  P1 + (P2-P1) * (1-exp(-P3*x)), start = list(P1 = 1, P2 = 3, P3 = 1),
              control= nls.control(maxiter = 1000))
summary(nlmod3)
AIC(nlmod3)
pol(a = (1.4*100)/4.4, b =4.4, c = 0.25, from = 0, to = 8, add = TRUE, col = "red") #all
pol(a = (1.9*100)/3.5, b = 3.5, c = 0.4, from = 0, to = 7, add = TRUE, col = "green") #Swedish crops
#pol(a = (-5.7*100)/3.2, b = 3.2, c = 1.4, from = 0, to = 7, add = TRUE, col = "green") #Swedish crops
pol(a = (0.59*100)/2.8, b = 2.8, c = 0.66, from = 1, to = 5, add = TRUE, col = "blue") #redclover
#pol(a = (1.74*100)/3.69, b = 3.69, c = 0.36, from = 1, to = 6, add = TRUE, col = "blue") #redclover
nlmod4 <- nls(y ~  P2 / (1 + (P2/P1 - 1) * exp(-P3*x)), start = list(P1 = 1, P2 = 3, P3 = 1),
              control= nls.control(maxiter = 1000))
summary(nlmod4)
AIC(nlmod4)
pol2(a = (1.7*100)/4.1, b = 4.1, c = 0.45, from = 0, to = 8, add = TRUE, col = "red") #all
pol2(a = (1.9*100)/3.4, b = 3.4, c = 0.65, from = 0, to = 7, add = TRUE, col = "green") #Swedish crops
pol2(a = (0.56*100)/2.68, b = 2.68, c = 1.5, from = 0, to = 5, add = TRUE, col = "blue") #redclover

#plot that parametrs for a fake crop:

pol(a = 30, b = 200, c = 0.4, from = 0, to = 10, col = 2)
#pol2(a = 30, b = 200, c = 0.65, add = TRUE, from = 0, to = 10, col = 4)

pol(a = 0, b = 200, c = 0.4, add = TRUE, from = 0, to = 10, col = 4)
#pol2(a = 0.1, b = 200, c = 0.65, add = TRUE, from = 0, to = 10, col = 3)
#pol2(a = 0.1, b = 200, c = 4, add = TRUE, from = 0, to = 10, col = 3)

pol(a = 80, b = 200, c = 0.4, add = TRUE, from = 0, to = 10, col = 4)

#NOTE POL2 is not working because A is not the intercept!

#pol(a = 0, b = 200, c = 0.25, add = TRUE, from = 0, to = 10, col = 4)
#pol(a = 0, b = 200, c = 0.66, add = TRUE, from = 0, to = 10, col = 4)

#add SE
pol(a = 30, b = 200, c = 0, add = TRUE, from = 0, to = 10, col = 1)
pol(a = 30, b = 200, c = 0.9, add = TRUE, from = 0, to = 10, col = 1)

pol(a = 0, b = 200, c = 0, add = TRUE, from = 0, to = 10, col = 1)
pol(a = 0, b = 200, c = 0.9, add = TRUE, from = 0, to = 10, col = 1)


#list % of selfing + max yield for main crops
#read alarm data

alarm <- read.csv("CompleteData.csv", header = TRUE)
str(alarm)

a2 <- subset(alarm, treat == "Net")
m <- tapply(a2$yield, a2$crop, mean, na.rm = TRUE)

a2 <- subset(alarm, treat == "Open")
m2 <- tapply(a2$yield, a2$crop, mean, na.rm = TRUE)

(m/m2)*100