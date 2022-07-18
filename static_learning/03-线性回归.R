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

# 最优子集选择、向前选择、向后选择
library(leaps)
regfit.full=regsubsets(medv~. , data=Boston)
regfit.bwd=regsubsets(medv~. , data=Boston,
                      method = "backward")
regfit.fwd=regsubsets(medv~. , data=Boston,
                      method = "forward")
summary(regfit.bwd)
coef(regfit.bwd,8)
coef(regit.fwd,8)
coef(regfit.full,8)

# 岭回归和lasso

# 数据矩阵化

x=model.matrix(crim~.,data = Boston)[,-14] #将响应变量从x变量中剔除
y= medv

library(glmnet)
grid=10^seq(10,-2,length=100)
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)
plot(ridge.mod)
lasso.mod=glmnet(x,y,alpha=1,lambda=grid)
plot(lasso.mod)

# 主成分回归与偏最小二乘回归
library(pls)
set.seed(2)
pcr.fit=pcr(medv~. , data=Boston,scale=TRUE,
            validation="CV")
summary(pcr.fit)
validationplot(pcr.fit,val.type="MSEP")
pls.fit=plsr(medv~. , data=Boston,scale=TRUE,
            validation="CV")
summary(pls.fit)
validationplot(pls.fit,val.type="MSEP")

# 回归样条
library(aplines)
fit= lm(medv~bs(age,knots=c(25,)).,data = Boston)
# 自然样条
fit= lm(medv~ns(age,df=4).,data = Boston)
# 光滑样条
fit=smooth.spline(age,wage,df=4)
#局部回归
fit=loess(wage~age,span=.2,data=Wage)

#GAM
#略





