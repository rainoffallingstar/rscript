age <- c(56,46,39,49,62,55,23,47,57,63,38,55)
weight <- c(70,56,89,79,62,45,53,67,47,93,78,65)
mean(weight)
sd(weight)
cor(age,weight)
plot(age,weight)
getwd() #Lists the current working directory.
q()
ls()#Lists the objects in the current workspace
rm(age)#Removes (deletes) one or more objects
setwd("/home/rstudio/R in action/")
options()
options(digits = 3)
x <- runif(20)
summary(x)
hist(x)
save.image("histofx") #save the current workspace not the pictures
install.packages("gclus")
installed.packages("gclus")
library("gclus")
# 卸载包
remove.packages("")