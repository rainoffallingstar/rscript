# 本部分包括logistic回归、lda、qda及knn
library(ISLR)
library(MASS)
attach(Smarket)
names(Smarket)
str(Smarket)
summary(Smarket)
cor(Smarket[,-9])
plot(Smarket[,-9])

# logistic 回归

train=(Year<2005)
Smarket.2005=Smarket[!train,]
Direction.2005=Direction[!train]

glm.fit=glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume,data=Smarket,
            famliy=binomial,subset = train) # don't work?
summary(glm.fit)

glm.probs=predict(glm.fit,Smarket.2005,type="response")
glm.probs[1:10]
glm.pred=rep("down",1250)
glm.pred(glm.probs>.5)="up"
table(glm.pred,Direction.2005)

# LDA
lda.fit=lda(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume,data=Smarket,
            subset = train)
summary(lda.fit)
lda.fit
lda.pred=predict(lda.fit,Smarket.2005)
lda.pred 
names(lda.pred)
lda.class <- lda.pred$class
table(lda.class,Direction.2005)

#QDA
qda.fit=qda(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume,data=Smarket,
            subset = train)
summary(qda.fit)
qda.fit
qda.pred=predict(qda.fit,Smarket.2005)
qda.pred 
names(qda.pred)
qda.class <- qda.pred$class
table(qda.class,Direction.2005)

# KNN
library(class)
train.X <- cbind(Lag1,Lag2)[train,]
test.X <- cbind(Lag1,Lag2)[!train,]
train.Direction <- Direction[train]
set.seed(1)
knn.pred=knn(train.X,test.X,train.Direction,k=10)
table(knn.pred,Direction.2005)
  