# R with pubmed
library(RISmed)
search_topic <- c('PAH-CHD')
search_query <- EUtilsSummary(search_topic,db="pubmed", retmax=10000,datetype='pdat', mindate=2010, maxdate=2022)
EUtilsSummary()
查查看下检索内容
summary(search_query)
看下这些文献的Id
QueryId(search_query)
获取检索结果
records<- EUtilsGet(search_query)
pubmed_data <- data.frame('Title'=ArticleTitle(records),
                          'Year'=YearAccepted(records),
                          'journal'=ISOAbbreviation(records))
head(pubmed_data,1)
pubmed_data[1:3,1]
write.csv(pubmed_data,file='PAH-CHD.csv')
分析文章情况
y <- YearPubmed(EUtilsGet(search_query))
可视化一下
library(ggplot2)
date()
count <- table(y)
count <- as.data.frame(count)
names(count)<-c("Year", "Counts")
library(RColorBrewer)
library(ggsci)
ggplot(data=count, aes(x=Year, y=Counts,fill=Year)) +
  geom_bar(stat="identity", width=0.5)+
  labs(y = "Number of articles",title="PubMed articles containing PAH-CHD"
  ) + theme_bw() + scale_fill_manual(values = colorRampPalette(brewer.pal(10, "Accent"))(10)) +
  theme(legend.position="bottom")