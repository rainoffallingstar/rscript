# vector
age <- c(56,46,39,49,62,55,23,47,57,63,38,55)
age[1]
age[c(2,6)] # you can't use age[2,6]
age[3:6]
sex <- c("male","male","male","male","female","female","male","male","male","female","female","female")

# matrices
y <- matrix(1:20,nrow = 5,ncol = 4)
cells <- c(1,26,24,68)
rnames <- c("R1","R2")
cnames <- c("c1","c2")
mymatrix <- matrix(cells,nrow = 2,ncol = 2,byrow = TRUE,dimnames=list(rnames,cnames))
# byrow = true,the vector will fill the matrix by row
mymatrix2 <- matrix(cells,nrow = 2,ncol = 2,byrow = FALSE,dimnames=list(rnames,cnames))
mymatrix2[1,] # mymatrix2[1,]not use the colnames or rownames
mymatrix2[1,1]

# array
dim1 <- c("A1","A2")
dim2 <- c("B1","B2","B3")
dim3 <- c("C1","C2","C3","C4")
z <- array(1:24,c(2,3,4),dimnames=list(dim1,dim2,dim3))
z
# you can't use byrow in array

patientID <- c(1,2,3,4,5,6,7,8,9,10,11,12)
treament <- c("kabofenjing","kabofenjing","kabofenjing","kabofenjing","kabofenjing","kabofenjing","fuliconzo","fuliconzo","fuliconzo","fuliconzo","fuliconzo","fuliconzo")
Sevendaystatue <- c("excellent","improve","improve","excellent","excellent","improve","improve","poor","improve","improve","poor","improve")
# you can add a number(0:9)after the varity but not pre-it
fourteenstatue <- c("excellent","excellent","excellent","excellent","excellent","excellent","excellent","improve","excellent","excellent","improve","excellent")
patientdata <- data.frame(patientID,age,treament,Sevendaystatue,fourteenstatue)
patientdata
patientdata[3:5]
patientdata$age
table(patientdata$treament,patientdata$Sevendaystatue)
summary(patientdata$age)
plot(patientdata$treament,patientdata$Sevendaystatue)

# the following codes won't work,because   The  limitations  with  this  approach  are  evident  when  more  than  one  object  can have the same name.
attach(patientdata)
  summary(treatment)
  plot(treament,Sevendaystatue)
detach(patientdata)

with(patientdata,{
  print(summary(patientdata$treatment))
  plot(patientdata$treament,patientdata$Sevendaystatue)
})

#factor

save.image("rinaction")

load("rinaction")

patientdata2 <- data.frame(patientID,age,treament,Sevendaystatue,fourteenstatue,row.names=patientID)

# factor define all the factors
treament <- factor(treatment,levels=c(1,2),labels=c("kabofenjing","fuliconzo")) 
sex <- factor(sex, levels=c(1, 2), labels=c("male", "female"))
Sevendaystatue <- c("excellent","improve","improve","excellent","excellent","improve","improve","poor","improve","improve","poor","improve")
Sevendaystatue <- factor(Sevendaystatue,order=TRUE,levels=c("poor","improve","excellent"))
fourteenstatue <- factor(fourteenstatue,order=TRUE,levels=c("poor","improve","excellent"))
str(patientdata)
summary(patientdata)

#lists
title <- "在院真菌感染患者治疗情况"
mylist <- list(title=title,age=age,treament,Sevendaystatue)
mylist
mylist[1]

install.packages("xlsx")
