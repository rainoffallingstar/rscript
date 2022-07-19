# 主成分分析
states=row.names(USArrests)
states
names(USArrests)
apply(USArrests,2,mean)
pr.out=prcomp(USArrests,scale=TRUE)
names(pr.out)
biplot(pr.out,scale = 0)
# KNN

# 系统聚类

library(ISLR)
nci.labs=NCI60$labs
nci.data=NCI60$data
dim(nci.data)
table(nci.labs)
pr.out=prcomp(nci.data,scale=TRUE)
cols=function(vec){
  + cols=rainbow(length(unique(vec)))
  + return(cols[as.numeric(as.factor(vec))])
  }
par(mfrow=c(1,2))
plot(pr.out$x[,1:2],col=cols(nci.labs),pch=19,
              xlab = "z1",ylab = "z2")
plot(pr.out$x[,1:3],col=cols(nci.labs),pch=19,
     xlab = "z1",ylab = "z3")
