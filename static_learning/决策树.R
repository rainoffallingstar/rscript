
# 决策树
library(tree)
library(ISLR)
attach(Carseats)
High=ifelse(Sales<=8,"No","Yes")
Carseats= data.frame(Carseats,High)
tree.catseats=tree(High~.-Sales,Carseats)
summary(tree.catseats)
plot(tree.catseats)
text(tree.catseats,pretty = 0)

# 测试集的引入

set.seed(2)
train=sample(1:nrow(Carseats),200)
carseats.test=Carseats[-train,]
High.test = High[-train]
tree.catseats=tree(High~.-Sales,Carseats,subset=train)
tree.pred=predict(tree.catseats,carseats.test,type = "class")
table(tree.pred,High.test)

# 剪枝
set.seed(3)
cv.carseats=cv.tree(tree.catseats,FUN=prune.misclass)
names(cv.carseats)

par(mfrow=c(1,2))
plot(cv.carseats$size,cv.carseats$dev,type="b")
plot(cv.carseats$k,cv.carseats$dev,type="b")

prune.carseats=prune.misclass(tree.catseats,best=9)
plot(prune.carseats)
text(prune.carseats,pretty=0)
tree.pred=predict(prune.carseats,carseats.test,type="b")

# 回归树

library(MASS)
set.seed(1)
train=sample(1:nrow(Boston),nrow(Boston)/2)
Boston.test=Boston[-train,"medv"]
tree.boston=tree(medv~.,Boston,subset=train)
tree.pred=predict(tree.boston,Boston.test,type = "class")
table(tree.pred,Boston.test)

# 剪枝
set.seed(3)
cv.boston=cv.tree(tree.boston)
names(cv.boston)

par(mfrow=c(1,2))
plot(cv.boston$size,cv.boston$dev,type="b")
plot(cv.boston$k,cv.boston$dev,type="b")

prune.boston=prune.tree(tree.boston,best=5)
plot(prune.boston)
text(prune.boston,pretty=0)
tree.pred=predict(prune.boston,carseats.test,type="b")
plot(tree.pred,Boston.test)

# 袋装法与随机森林

library(randomForest)
set.seed(1)
bag.boston=randomForest(medv~.,data=Boston,subset=train,
                        mtry=13,importance=TRUE)
bag.boston

pred.bag = predict(bag.boston,newdata=Boston[-train,])
plot(pred.bag,Boston.test)
abline(0,1)
mean((ped.bag-Boston.test)^2)
importance(bag.boston)

# 提升法
library(gbm)
set.seed(1)
boost.boston=gbm(medv~.,data=Boston[train,],
                 distribution="gaussian",n.trees=5000,interaction.depth=4)
summary(boost.boston)


