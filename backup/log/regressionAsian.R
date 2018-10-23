#!/usr/bin/env Rscript


library(gplots)
library(MASS)
library(MatchIt)
library(dummies)


data=read.csv("/home/acquireassets/Desktop/newdata.csv")



data1=data[data$stateicp%in%c("California","New York","Massachusetts","Maryland","Michigan","Ohio",
"Mississippi","Oklahoma","Texas","Florida","North Carolina","Virginia","Georgia","Tennessee",
"Missouri","Wisconsin","Vermont","Washington","Alabama","Arizona","Illinois",
"Idaho","Kansas","Mississippi","Nebraska","Oklahoma","South Carolina","South Dakota","Utah","Wyoming"),]



data1=data1[data1$year==2012|data1$year==2016,]



data1$time=as.numeric(data1$year!=2012)


data1$group=as.numeric(data1$stateicp%in%c("Alabama","Florida","Georgia","Idaho","Kansas","Mississippi","Missouri","Nebraska","North Carolina","Oklahoma","South Carolina","South Dakota","Tennessee","Texas","Utah","Virginia","Wisconsin","Wyoming"))

data1$group1=abs(1-data1$group)



data2=data1[data1$incearn>1,]
data2$outcome=log(data2$incearn)

data2=data2[which((data2$race=="Japanese") | (data2$race=="Chinese") | (data2$race=="Other Asian or Pacific Islander")),]

 model1=rlm(outcome~costelec+costgas+costwatr+costfuel+valueh+perwt+pernum+as.factor(sex)+
as.numeric(age)+as.factor(marst)+
as.factor(hcovany)+as.factor(empstat)
+as.factor(classwkrd)+as.factor(year)+as.factor(group1)+year:group1,data=data2)




good=summary(model1)
dd = data.frame(good$coefficients) 
dd$p.value =  2*pt(abs(dd$t.value),good$df[2], lower.tail=FALSE)      
 dd
