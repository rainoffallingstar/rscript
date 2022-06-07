manager <- c(1,2,3,4,5)
date <- c("10/24/08","10/28/08","10/1/08","10/12/08","5/1/09")
country <- c("US","US","UK","UK","UK")
gender <- c("M","F","F","M","F")
age <- c(32,45,25,39,99)
q1 <- c(5,3,3,3,2)
q2 <- c(4,5,5,3,2)
q3 <- c(5,2,5,4,1)
q4 <- c(5,5,5,NA,2)
q5 <- c(5,5,2,NA,1)
leadership <- data.frame(manager,date,country,
                         gender,age,q1,q2,q3,q4,q5)
leadership
mydata<-data.frame(x1 = c(2, 2, 6, 4), 
                   x2 = c(3, 4, 2, 8))
mydata <- transform(mydata,
                    sumx=x1+x2,
                    meanx=(x1+x2)/2)
mydata
leadership <- within(leadership,{
                       agecat <- NA
                       agecat[age>75] <- "elder"
                       agecat[age>=55 & age<=75] <- "middle aged"
                       agecat[age<55] <- "young"})
leadership
leadership <- within(leadership,{
                     agecat <- factor(agecat,order=TRUE,levels=c("young","middle aged","elder"))})
# in this situation,you can not change the fator form outside the data frame
agecat <- factor(agecat,order=TRUE,levels=c("young","middle aged","elder"))

library(plyr)
leadership <- rename(leadership,
                     c(manager="managerID",date="testDate"))
is.na(leadership[,6:10])
na.omit(leadership) # 删除带有na的个体
datefromat <- "%m/%d/%y"
leadership$testDate <- as.Date(leadership$testDate,datefromat)
Sys.Date()
newdata <- leadership[order(leadership$age),]
total <- merge(dataframeA, dataframeB, by="ID")
#merges dataframeA and dataframeB by ID. Similarly
myvars <- names(leadership) %in% c("q3", "q4") # 应注意逻辑运算符号 %in%
newdata <- leadership[!myvars] 

newdata <- leadership[1:3,]
newdata <- leadership[leadership$gender == "M" &
                        leadership$age >30]
startdate <- as.Date("2009-01-01")
enddate <- as.Date("2009-10-31")
newdate <- leadership[which(leadership$testDate >= startdate &
                             leadership$testDate <= enddate)]
