library(MASS)
library(ISLR)
# 定量分析时采用boston数据集
fix(Boston)
attach(Boston)
names(Boston)
lm.fit=lm(medv~lstat , data=Boston) 
lm.fit #模型具体信息
summary(lm.fit) #模型p值及R2
confint(lm.fit) # 置信区间 use confint not confit
predict(lm.fit , data.frame (lstat=(c (5 , 10 , 15))) , 
         interval="confidence") # 预测值及相应置信区间
# 以下作图

plot(lstat,medv)
abline(lm.fit) # use abline not abine
# 使用多个变量时
lm.fit2=lm(medv~lstat+age,data=Boston)
# 引入非线性项时
lm.fit1=lm(medv~lstat+I(age^2),data=Boston)
summary(lm.fit2)
# use all varities 
lm.fit3=lm(medv~.,data = Boston)
# 在lm.fit3的基础上进行向后选择
lm.fit4=update(lm.fit3,~.-age)
# 在模型间进行对比比较
anova(lm.fit3,lm.fit4)
detach(Boston)

# issues：3D散点图

# 利用定性变量进行回归时，系统会自动创造编码
attach(Carseats)
names(Carseats)
lm.fit=lm(Sales~.,data=Carseats)
summary(lm.fit)
contrasts(ShelveLoc)
