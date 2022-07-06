library(survival)
data(ovarian)
str(ovarian)
colevel <- c(1,2)
ovarian$resid.ds <- factor(ovarian$resid.ds,
                           levels = colevel,
                           labels = c("no","yes"))
ovarian$rx <- factor(ovarian$rx,
                     levels = colevel,
                     labels = c("A","B"))
ovarian$ecog.ps <- factor(ovarian$ecog.ps,
                          levels = colevel,
                          labels = c("good","bad"))
ovarian$agegr <- cut(ovarian$age,
                     breaks = c(0,50,75),
                     labels = c("<=50",">50"))
surv.obj <- Surv(time = ovarian$futime,event = ovarian$fustat)
surv.obj
surv.all <- survfit(surv.obj ~ 1)
summary(surv.all,censored = TRUE)
plot(surv.all,mark.time = TRUE)

# 对于rx进行生存分析并比较差异性
surv.treat <- survfit(surv.obj ~ rx,data = ovarian)
summary(surv.treat)
library(survminer)
ggsurvplot(surv.treat,data=ovarian,pval =TRUE)
survdiff(suv.obj ~ rx,data = ovarian)
