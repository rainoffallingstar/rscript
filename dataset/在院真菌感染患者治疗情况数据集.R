age <- c(56,46,39,49,62,55,23,47,57,63,38,55)
weight <- c(70,56,89,79,62,45,53,67,47,93,78,65)
sex <- c("male","male","male","male","female","female","male","male","male","female","female","female")
patientID <- c(1,2,3,4,5,6,7,8,9,10,11,12)
treament <- c("kabofenjing","kabofenjing","kabofenjing","kabofenjing","kabofenjing","kabofenjing","fuliconzo","fuliconzo","fuliconzo","fuliconzo","fuliconzo","fuliconzo")
Sevendaystatue <- c("excellent","improve","improve","excellent","excellent","improve","improve","poor","improve","improve","poor","improve")
# you can add a number(0:9)after the varity but not pre-it
fourteenstatue <- c("excellent","excellent","excellent","excellent","excellent","excellent","excellent","improve","excellent","excellent","improve","excellent")
treatment <- factor(treament,levels=c(1,2),labels=c("kabofenjing","fuliconzo")) 
sex <- factor(sex, levels=c(1, 2), labels=c("male", "female"))
Sevendaystatue <- c("excellent","improve","improve","excellent","excellent","improve","improve","poor","improve","improve","poor","improve")
Sevendaystatue <- factor(Sevendaystatue,order=TRUE,levels=c("poor","improve","excellent"))
fourteenstatue <- factor(fourteenstatue,order=TRUE,levels=c("poor","improve","excellent"))
patientdata <- data.frame(patientID,age,weight,sex,treament,Sevendaystatue,fourteenstatue,row.names=patientID)
str(patientdata)

#lists
title <- "在院真菌感染患者治疗情况"
mylist <- list(title=title,age=age,treament,Sevendaystatue,patientdata)
mylist
str(mylist)
