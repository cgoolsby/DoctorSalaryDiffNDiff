data1=data[data$stateicp%in%c("California","New York","Massachusetts","Maryland","Michigan","Ohio",
"Mississippi","Oklahoma","Texas","Florida","North Carolina","Virginia","Georgia","Tennessee",
"Missouri","Wisconsin","Vermont","Washington","Alabama","Arizona","Illinois",
"Idaho","Kansas","Mississippi","Nebraska","Oklahoma","South Carolina","South Dakota","Utah","Wyoming"),]


data1=data1[data1$year>=2008,]

data1=data1[data1$year<=2016,]

data1$group=as.numeric(data1$stateicp%in%c("Maine","Alabama","Florida","Georgia","Idaho","Kansas","Mississippi","Missouri","Nebraska","North Carolina","Oklahoma","South Carolina","South Dakota","Tennessee","Texas","Utah","Virginia","Wisconsin","Wyoming"))

data1$group1=abs(1-data1$group)


data2=data1[data1$incearn>1,]
data2$outcome=log(data2$incearn)

#model1=lm(outcome~costelec+costgas+costwatr+costfuel+valueh+perwt+pernum+as.factor(sex)+age+as.factor(marst)+
#as.factor(race)+as.factor(hcovany)+as.factor(educ)+as.factor(empstat)
#+as.factor(classwkrd)+incinvst+as.factor(time)+as.factor(group1)+time:group1,data=data2)

#Select Specific Dataset, i.e. race, gender etc
#data2=data2[data2$raced=="Black/African American/Negro",] 
data2=data2[data2$raced=="White",] 

data3=data2[data2$group1==0,]

data4=data2[data2$group1==1,]
par(mfrow=c(1,1))
#plotmeans(outcome~year,data=data3)
#plotmeans(outcome~year,data=data4)

comparison=aggregate(incearn ~ year, data3, mean)
comparison1=aggregate(incearn ~ year, data4, mean)

png('plot.png')
#plot(comparison1[,1],comparison1[,2],type="l",ylim=c(184000,250000),col="blue",xlab="year",ylab="Mean income")
plot(comparison1[,1],comparison1[,2],type="l",col="blue",xlab="year",ylab="Mean income")

lines(comparison[,1],comparison[,2], col=c("darkgreen"))


aggregate(incearn~ group1, data2, mean)

dev.off()



