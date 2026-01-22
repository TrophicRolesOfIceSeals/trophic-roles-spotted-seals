### take a look at wt-length relationships

Data <- read.csv('./data/Lth-Wt-Energy-Combined.csv')

### SandLance

SandLance = Data[which(Data$Species=="Pacific sand lance"),]

PacificLance = SandLance[which(SandLance$Latitude<65.818889),]
ArcticLance = SandLance[which(SandLance$Latitude>65.818889),]

PacificLM = lm(PacificLance$logwt ~ PacificLance$loglth)
ArcticLM = lm(ArcticLance$logwt ~ ArcticLance$loglth)

plot(PacificLance$loglth,PacificLance$logwt)
points(ArcticLance$loglth,ArcticLance$logwt,col="blue")

lines(PacificLance$loglth,fitted.values(PacificLM))
lines(ArcticLance$loglth,fitted.values(ArcticLM),col="blue")

loglth = c(15:23)/10
lines(loglth,-6.3871+3.18601*loglth,col="red")
lines(loglth,-6.0376+3.2793*loglth,col="green")

summary(PacificLM)
summary(ArcticLM)

#just base off of combined data

CombinedLM = lm(SandLance$logwt ~ SandLance$loglth)
summary(CombinedLM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(CombinedLM$residuals))))
bias10

SL_FL_to_wetwt <- function(FL){1.008*10^(-6.191 + 3.3545*log10(FL))
}


### Capelin

Capelin = Data[which(Data$Species=="Pacific capelin"),]
Capelin$logwt=as.numeric(Capelin$logwt)
Capelin = Capelin[-which(is.na(Capelin$logwt)),]

Pacific = Capelin[which(Capelin$Latitude<65.818889),]
Arctic = Capelin[which(Capelin$Latitude>65.818889),]

plot(Pacific$loglth,Pacific$logwt)
points(Arctic$loglth,Arctic$logwt,col="blue")

#these look the same, just run 1 analysis... restrict
# to mm > 50 since these fish don't fit linear model very
# well and smallest capelin in spotted seal stomach was 
# around 100 mm TL
Capelin50 = Capelin[which(Capelin$Length..mm.>50),]
Capelin50[,"logwt"]=as.numeric(Capelin50[,"logwt"])
CapelinLM = lm(Capelin50$logwt~Capelin50$loglth)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(CapelinLM$residuals))))
bias10


Capelin_FL_to_wetwt = function(FL){
  1.013*10^(-6.314 + 3.518*log10(FL))
}


### Saffron cod

SC = Data[which(Data$Species=="Saffron cod"),]

Bering = SC[which(SC$Latitude<65.818889),]
Arctic = SC[which(SC$Latitude>65.818889),]
plot(Bering$loglth,Bering$logwt)
points(Arctic$loglth,Arctic$logwt,col="blue")

#perhaps a tiny indication that Bering sea fish weight less 
# at shorter lengths and more and longer lengths than than 
#their Arctic brethren, but not so much that I think it's worth dealing
#with 

Summer = SC[which(SC$Season=="Summer"),]
Winter = SC[which(SC$Season=="Winter"),]
plot(Summer$loglth,Summer$logwt)
points(Winter$loglth,Winter$logwt,col="blue")

#seasons look combinable too

#combined linear model

SC = SC[-which(SC$logwt==''),]
SC_LM <- lm(SC$logwt ~ SC$loglth)
summary(SC_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(SC_LM$residuals))))
bias10


#note this is is total length, NOT fork length
SC_TL_to_wetwt <- function(TL){
  1.022*10^(-5.168+3.023*log10(TL))
}

#double check fit on abs scale
plot(SC$Length..mm.,SC$Whole.Body.Wet.Weight..g.)
TL = c(10:415)
lines(SC_TL_to_wetwt(TL))

### Arctic cod

#looks like sample sizes are really small, except for Arctic/Summer
#where n = 2507. Probably best just to combine into one relationship.

AC = Data[which(Data$Species=="Arctic cod"),]

AC_LM <- lm(AC$logwt ~ AC$loglth)
summary(AC_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(AC_LM$residuals))))
bias10

#note this is again total length, NOT fork length
AC_TL_to_wetwt <- function(TL){
  1.016*10^(-5.5436+3.1804*log10(TL))
}

#double check fit on abs scale
plot(AC$Length..mm.,AC$Whole.Body.Wet.Weight..g.)
TL = c(10:415)
lines(TL,AC_TL_to_wetwt(TL))


###  P Cod

#NELD has greater sample sizes (n=1030) than other potential sources 

#almost all samples come from the Bering in the summer... so probably
#just worth combining all into one equation

PC = Data[which(Data$Species=="Pacific cod"),]
PC_LM <- lm(PC$logwt ~ PC$loglth)
summary(PC_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(PC_LM$residuals))))
bias10


#note this is again total length, NOT fork length
PC_TL_to_wetwt <- function(TL){
  1.013*10^(-5.43+3.172*log10(TL))
}

#double check fit on abs scale
plot(PC$Length..mm.,PC$Whole.Body.Wet.Weight..g.)
TL = c(10:900)
lines(TL,PC_TL_to_wetwt(TL))


### Walleye pollock - looks like vast majority of samples are Bering/Summer


WP = Data[which(Data$Species=="Walleye pollock"),]

#look at arctic vs bering
Bering = WP[which(WP$Latitude<65.818889),]
Arctic = WP[which(WP$Latitude>65.818889),]
plot(Bering$loglth,Bering$logwt)
points(Arctic$loglth,Arctic$logwt,col="blue")

#there doesn't look like a reason to separate based on N/S of Bering Strait
# BUT... small fish lood really weird!!

# it looks like our smallest otolith length in a spotted seal stomach is 3.32 mm
# which translates to a 66.6 mm FL (using functions developed and saved in
# "otolith_regressions.RData)

# By contrast, the 'weird' relationships start with log10 length < 1.4 or so,
# which translates into < 25 mm FL.  So, we should probably just leave out these really small fish, since
# they'll mess up the regression and are much smaller than the range of otolith
# lengths we're interested in.  
#

WP = WP[which(WP$loglth>1.5),]

WP = WP[-which(WP$logwt==''),]  #take out missing values

WP_LM <- lm(WP$logwt ~ WP$loglth)
summary(WP_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(WP_LM$residuals))))
bias10

WP_FL_to_wetwt <- function(FL){
  1.012*10^(-5.27+3.047*log10(FL))
}

#double check fit on abs scale
plot(WP$Length..mm.,WP$Whole.Body.Wet.Weight..g.)
FL = c(32:700)
lines(FL,WP_FL_to_wetwt(FL))


### Rainbow smelt

# we have more data than the other potential source (Harvey 2000) - 

RS = Data[which(Data$Species=="Rainbow smelt"),]

#look at arctic vs bering
Bering = RS[which(RS$Latitude<65.818889),]
Arctic = RS[which(RS$Latitude>65.818889),]
plot(Bering$loglth,Bering$logwt)
points(Arctic$loglth,Arctic$logwt,col="blue")

#these look very combinable! 
#sample size = 119

RS_LM <- lm(RS$logwt ~ RS$loglth)
summary(RS_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(RS_LM$residuals))))
bias10


RS_FL_to_wetwt <- function(FL){
  1.018*10^(-5.887+3.339*log10(FL))
}

#double check fit on abs scale
plot(RS$Length..mm.,RS$Whole.Body.Wet.Weight..g.)
FL = c(32:300)
lines(FL,RS_FL_to_wetwt(FL))


### Pacific herring

# we have more data than other potential sources 
# there are only n=10 samples from north of Bering Strait... so will be best to combine


PH = Data[which(Data$Species=="Pacific herring"),]
PH = PH[-which(PH$logwt==''),]  #take out missing values
PH = PH[-which(PH$logwt=='#NUM!'),]  #take out missing values

#sample size = 1919

PH_LM <- lm(PH$logwt ~ PH$loglth)
summary(PH_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(1000000,0,sqrt(var(PH_LM$residuals))))
bias10

PH_FL_to_wetwt <- function(FL){
  1.009*10^(-5.596+3.26*log10(FL))
}

#double check fit on abs scale
plot(PH$Length..mm.,PH$Whole.Body.Wet.Weight..g.)
FL = c(32:400)
lines(FL,PH_FL_to_wetwt(FL),col="blue")

### Yellowfin sole

#looking at regression coefficients on Google Drive, 
#there's not much difference between N or S of Bering Strait - so suggest combining again

YS = Data[which(Data$Species=="Yellowfin sole"),]
YS = YS[-which(YS$logwt=='#NUM!'),]  #take out missing values

#sample size = 19842

YS_LM <- lm(YS$logwt ~ YS$loglth)
summary(YS_LM)

#for log10-normal bias correction use monte carlo
bias10 = mean(10^rnorm(10000000,0,sqrt(var(YS_LM$residuals))))
bias10

#note: TL not FL
YS_TL_to_wetwt <- function(TL){
  1.008^(-5.33+3.161*log10(TL))
}

#double check fit on abs scale
plot(YS$Length..mm.,YS$Whole.Body.Wet.Weight..g.)
TL = c(32:500)
lines(TL,YS_TL_to_wetwt(TL),col="blue")

#looks like we might be underpredicting weights of very large fish but this probably isn't an
# issue for the range of sizes of fish that spotted seals are eating


#### Export functions
save(YS_TL_to_wetwt,RS_FL_to_wetwt,PH_FL_to_wetwt,WP_FL_to_wetwt,
     PC_TL_to_wetwt,AC_TL_to_wetwt,SC_TL_to_wetwt,Capelin_FL_to_wetwt,
     SL_FL_to_wetwt,file="./data/weight_length_regressions.RData")
