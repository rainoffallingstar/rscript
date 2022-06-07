jpeg("fig1.jpg") # don't forget the format .jpg
 attach(mtcars)
 plot(wt,mpg)
 abline(lm(mpg~wt)) # use ~ but not "-" 辅助线
 title("regression of mpg on weights")
 detach(mtcars)
dev.off()
 
dosage <- c(20,30,40,45,60)
drugA <- c(16,20,27,40,60)
drugB <- c(15,18,25,31,40)
plot(dosage,drugA,type = "e") # work when type=b/c

opar <- par(no.readonly = TRUE)
par(fig=c(0, 0.8, 0, 0.8))  # 设置画布大小
par(lty=2,pch=17) #lty设置线条类型，pch设置符号类型
plot(dosage,drugA,type = "b",lwd=3,cex=2,col="red",
     main="clinical trials for drugA",col.main="black",
     sub="this is a hypothetical data",col.sub="black",
     xlab="dosage",ylab="drug response",col.lab="black",cex.lab=.75,
     xlim=c(0,60),ylim=c(0,70))
lines(dosage,drugB,type = "b",pch=17,lty=2,col="blue")
install.packages("Hmisc")
install.packages(c("lattice","survival","Formula","ggplot2"))
library(Hmisc)
minor.tick(nx=3,ny=3,tick.ratio=0.5) #do not miss rato with ratio
legend("topleft",inset=.05,title = "drug type",c("A","B"),
       lty=c(1,2),pch=c(15,17),col=c("red","blue"))
#xlim ylim 设置坐标轴，xlab设置坐标轴名称
#type设置点类型,lwd设置线宽，cex设置点放大比例
par(opar)
install.packages("RColorBrewer")
library(RColorBrewer)
n <- 7
mycolors <- brewer.pal(n,"Set1") #use Set not set
barplot(rep(1,n),col=mycolors)
par(fig=c(0, 0.8, 0.2, 1), new=TRUE) 
boxplot(dosage,axes=FALSE,horizontal=TRUE)
