#model with intercept, uncertainty intercept, slope, uncertainty slope, number of bottles, r2 adjusted, probability of significance
#model of linear regression applied, Being S1 the region under study, the control treatment K, and the light condition L
#install.packages("Hmisc")
#install.packages("multcomp")
#install.packages("seacarb")
#install.packages("itsadug")
library("Hmisc")
library("corrplot")
library("ggplot2")
library("matrixStats")
library("SuppDists")##Steel for multiple non-parametric comparisonsmultiples
library("kSamples")
library("MASS")
library("dunn.test")
library("multcomp")
library("sandwich")
library("seacarb")
library("oce")
library("gsw")
library("SolveSAPHE")
library("gridExtra")

source("/home/evaristo/Documentos/MOC2/Fun.R")

modelKL <- function(x,y,j=12)
    {
        model <- lm(y~x)
        sumario <- summary(model)
        S1KL <- data.frame(rep(c("Exp1"),j),rep(c("Control"),j),rep(c("Light"),j),x,y)
        colnames(S1KL) <- c("Experiment","Treatment","Illuminati","Time","Oxygen")
        interceptoS1KL <- sumario$coefficients[1,1]
        uinterceptoS1KL <- sumario$coefficients[1,2]
        slopeS1KL <- sumario$coefficients[2,1]
        uslopeS1KL <- sumario$coefficients[2,2]
        nS1KL <- length(na.omit(y))
        adjr2S1KL <- sumario$adj.r.squared
        pS1KL <-pf(sumario$fstatistic[1],sumario$fstatistic[2],sumario$fstatistic[3],lower.tail=F)
        pS1KL <- if(pS1KL<0.01) {pS1KL=0.01} else {as.numeric(pS1KL)}
        rectaS1KL <- c(interceptoS1KL,uinterceptoS1KL,slopeS1KL,uslopeS1KL,nS1KL,adjr2S1KL,pS1KL)
        print(rectaS1KL)
    }

PointsTimes2448 <- function(x,y)
{
    x24 <- x[1:8]
    y24 <- y[1:8]
    x48 <- x[5:12]
    y48 <- y[5:12]
    data.frame(x24,y24,x48,y48)
}

#Exp 3#### Este es el que se encuentra más al W
#AMZR Control Light
#el recta...48 es la recta del periodo 24 - 48
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <- c(197.94,197.27,197.75,197.47,195.06,194.57,194.56,194.26,191.60,192.22,190.73,194.53)
punts <- PointsTimes2448(x,y)
rectaAMZRKL <- modelKL(x,y,j=12)
rectaAMZRKL24 <- modelKL(punts[,1],punts[,2],j=8)
rectaAMZRKL48 <- modelKL(punts[,3],punts[,4],j=8)
#revision 2023
x0 <- mean(y[c(1,4)])
CI0 <- 2.776*sd(y[c(1,4)])/sqrt(4-1)
x24 <- mean(y[c(5,8)])
CI24 <- 2.776*sd(y[c(5,8)])/sqrt(4-1)
x48 <- mean(y[c(9,12)])
CI48 <- 2.776*sd(y[c(9,12)])/sqrt(4-1)
mean_CI_time <- c(x0,CI0,x24,CI24,x48,CI48)
dif <- c(x24-x0,x48-x0,x48-x24)/24
mean(c(dif))

#AMZR Control Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <- c(NA,197.18,197.36,197.25,194.19,196.50,195.80,194.53,191.40,191.90,191.35,188.69)
punts <- PointsTimes2448(x,y)
rectaAMZRKD <- modelKL(x,y,j=12)
rectaAMZRKD24 <- modelKL(punts[,1],punts[,2],j=8)
rectaAMZRKD48 <- modelKL(punts[,3],punts[,4],j=8)
#revision 2023
x0 <- mean(y[c(1,4)],na.rm=TRUE)
CI0 <- 2.776*sd(y[c(1,4)],na.rm=TRUE)/sqrt(4-1)
x24 <- mean(y[c(5,8)],na.rm=TRUE)
CI24 <- 2.776*sd(y[c(5,8)],na.rm=TRUE)/sqrt(4-1)
x48 <- mean(y[c(9,12)],na.rm=TRUE)
CI48 <- 2.776*sd(y[c(9,12)],na.rm=TRUE)/sqrt(4-1)
mean_CI_time <- c(x0,CI0,x24,CI24,x48,CI48)
dif <- c(x24-x0,x48-x0,x48-x24)/24
mean(c(dif))


#AMZR High CO2 Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(197.96,197.40,198.16,197.91,195.46,193.13,195.47,195.54,193.61,195.58,191.91,NA)
punts <- PointsTimes2448(x,y)
rectaAMZRCO2L <- modelKL(x,y,j=12)
rectaAMZRCO2L24 <- modelKL(punts[,1],punts[,2],j=8)
rectaAMZRCO2L48 <- modelKL(punts[,3],punts[,4],j=8)

#AMZR High CO2 Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(198.48,198.53,197.71,197.18,188.47,189.68,194.40,190.69,185.05,187.94,186.59,188.93)
punts <- PointsTimes2448(x,y)
rectaAMZRCO2D <- modelKL(x,y,j=12)
rectaAMZRCO2D24 <- modelKL(punts[,1],punts[,2],j=8)
rectaAMZRCO2D48 <- modelKL(punts[,3],punts[,4],j=8)

#Exp 2#### 
#SAM Control Light, SAMKL
x <- c(0,0,0,0,24,24,24,24,47,47,47,47)
y <- c(205.40,205.72,NA,206.21,203.83,203.92,202.90,202.82,196.34,198.02,201.65,202.40)
punts <- PointsTimes2448(x,y)
rectaSAMKL <- modelKL(x,y,j=12)
rectaSAMKL24 <- modelKL(punts[,1],punts[,2],j=8)
rectaSAMKL48 <- modelKL(punts[,3],punts[,4],j=8)

#SAM Control Dark (oxygen in micM O2)
x <- c(0,0,0,0,24,24,24,24,47,47,47,47)
y <- c(205.36,205.29,205.28,205.43,203.62,202.67,NA,200.31,200.12,201.58,NA,199.40)
punts <- PointsTimes2448(x,y)
rectaSAMKD <- modelKL(x,y,j=12)
rectaSAMKD24 <- modelKL(punts[,1],punts[,2],j=8)
rectaSAMKD48 <- modelKL(punts[,3],punts[,4],j=8)

#SAM High CO2 light
x <- c(0,0,0,0,24,24,24,24,47,47,47,47)
y <- c(204.26,NA,204.32,NA,199.74,200.09,201.75,201.54,199.34,197.02,199.33,197.99)
punts <- PointsTimes2448(x,y)
rectaSAMCO2L <- modelKL(x,y,j=12)
rectaSAMCO2L24 <- modelKL(punts[,1],punts[,2],j=8)
rectaSAMCO2L48 <- modelKL(punts[,3],punts[,4],j=8)

#SAM HighCO2 Dark
x <- c(0,0,0,0,24,24,24,24,47,47,47,47)
y <- c(206.43,203.34,198.79,203.51,194.11,195.41,196.71,194.65,198.24,193.16,191.68,NA)
punts <- PointsTimes2448(x,y)
rectaSAMCO2D <- modelKL(x,y,j=12)
rectaSAMCO2D24 <- modelKL(punts[,1],punts[,2],j=8)
rectaSAMCO2D48 <- modelKL(punts[,3],punts[,4],j=8)

#Exp 4#### 
#CAT Control Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(205.80,204.87,204.19,204.31,200.55,202.82,202.38,201.54,200.10,201.12,198.89,201.99)
punts <- PointsTimes2448(x,y)
rectaCATKL <- modelKL(x,y,j=12)
rectaCATKL24 <- modelKL(punts[,1],punts[,2],j=8)
rectaCATKL48 <- modelKL(punts[,3],punts[,4],j=8)
plot(y~x)

#CAT Control Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(205.18,203.23,203.89,203.32,200.80,201.51,200.69,200.61,192.43,197.36,195.86,196.36)
punts <- PointsTimes2448(x,y)
rectaCATKD <- modelKL(x,y,j=12)
rectaCATKD24 <- modelKL(punts[,1],punts[,2],j=8)
rectaCATKD48 <- modelKL(punts[,3],punts[,4],j=8)
points(y~x,pch=19)

#CAT High CO2 Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(202.89,203.40,203.45,203.06,200.67,199.48,200.19,200.95,196.58,196.98,197.84,195.91)
punts <- PointsTimes2448(x,y)
rectaCATCO2L <- modelKL(x,y,j=12)
rectaCATCO2L24 <- modelKL(punts[,1],punts[,2],j=8)
rectaCATCO2L48 <- modelKL(punts[,3],punts[,4],j=8)

#CAT High CO2 Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(203.85,203.15,203.40,203.13,200.53,200.33,199.83,199.45,195.80,194.97,196.09,196.41)
punts <- PointsTimes2448(x,y)
rectaCATCO2D <- modelKL(x,y,j=12)
rectaCATCO2D24 <- modelKL(punts[,1],punts[,2],j=8)
rectaCATCO2D48 <- modelKL(punts[,3],punts[,4],j=8)

#Exp 5#### 
#MAR Control Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(217.78,221.08,223.23,224.02,218.95,219.55,220.16,NA,216.24,217.73,216.44,217.48)
punts <- PointsTimes2448(x,y)
rectaMARKL <- modelKL(x,y,j=12)
rectaMARKL24 <- modelKL(punts[,1],punts[,2],j=8)
rectaMARKL48 <- modelKL(punts[,3],punts[,4],j=8)

#MAR Control Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(222.71,222.41,222.25,222.85,219.29,218.50,218.64,218.34,217.61,215.23,NA,220.12)
punts <- PointsTimes2448(x,y)
rectaMARKD <- modelKL(x,y,j=12)
rectaMARKD24 <- modelKL(punts[,1],punts[,2],j=8)
rectaMARKD48 <- modelKL(punts[,3],punts[,4],j=8)

#MAR High CO2 Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(222.89,222.70,223.33,223.60,223.44,217.16,220.09,217.48,216.99,215.41,217.45,217.09)
rectaMARCO2L <- modelKL(x,y,j=12)
rectaMARCO2L24 <- modelKL(punts[,1],punts[,2],j=8)
rectaMARCO2L48 <- modelKL(punts[,3],punts[,4],j=8)

#MAR High CO2 Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(NA,222.98,223.28,223.38,217.51,217.15,216.32,217.28,215.04,214.27,213.67,216.20)
rectaMARCO2D <- modelKL(x,y,j=12)
rectaMARCO2D24 <- modelKL(punts[,1],punts[,2],j=8)
rectaMARCO2D48 <- modelKL(punts[,3],punts[,4],j=8)

#Exp 6
#GD Control Light
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(200.87,192.37,205.42,191.06,193.91,193.49,NA,192.32,190.95,190.35,186.23,190.95)
punts <- PointsTimes2448(x,y)
rectaGDKL <- modelKL(x,y,j=12)
rectaGDKL24 <- modelKL(punts[,1],punts[,2],j=8)
rectaGDKL48 <- modelKL(punts[,3],punts[,4],j=8)

#GD Control Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(198.97,198.56,198.34,198.46,193.92,NA,193.26,194.36,188.51,190.49,191.14,191.86)
punts <- PointsTimes2448(x,y)
rectaGDKD <- modelKL(x,y,j=12)
rectaGDKD24 <- modelKL(punts[,1],punts[,2],j=8)
rectaGDKD48 <- modelKL(punts[,3],punts[,4],j=8)

#GD High CO2 ligth
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <-c(202.58,201.47,201.12,201.20,196.23,195.84,195.50,195.29,194.92,190.27,192.08,191.56)
rectaGDCO2L <- modelKL(x,y,j=12)
rectaGDCO2L24 <- modelKL(punts[,1],punts[,2],j=8)
rectaGDCO2L48 <- modelKL(punts[,3],punts[,4],j=8)

#GD High CO2 Dark
x <- c(0,0,0,0,24,24,24,24,48,48,48,48)
y <- c(200.49,NA,200.57,203.51,194.48,195.75,194.25,196.60,192.71,191.30,191.21,193.08)
rectaGDCO2D <- modelKL(x,y,j=12)
rectaGDCO2D24 <- modelKL(punts[,1],punts[,2],j=8)
rectaGDCO2D48 <- modelKL(punts[,3],punts[,4],j=8)

#### Tabla respiraciones

AMZR <- cbind(rep("AMZR",4),c("Control","Control","High CO2","High CO2"),c("Light","Dark"))
SAM <- cbind(rep("SAM",4),c("Control","Control","High CO2","High CO2"),c("Light","Dark"))
CAT <- cbind(rep("CAT",4),c("Control","Control","High CO2","High CO2"),c("Light","Dark"))
MAR <- cbind(rep("MAR",4),c("Control","Control","High CO2","High CO2"),c("Light","Dark"))
GD <- cbind(rep("GD",4),c("Control","Control","High CO2","High CO2"),c("Light","Dark"))
rbind(AMZR,SAM,CAT,MAR,GD)

Totselstemps <- rbind(rectaAMZRKL,rectaAMZRKD,rectaAMZRCO2L,rectaAMZRCO2D,rectaSAMKL,rectaSAMKD,rectaSAMCO2L,rectaSAMCO2D,rectaCATKL,rectaCATKD,rectaCATCO2L,rectaCATCO2D,rectaMARKL,rectaMARKD,rectaMARCO2L,rectaMARCO2D,rectaGDKL,rectaGDKD,rectaGDCO2L,rectaGDCO2D)

SuppTable3 <- data.frame(rbind(AMZR,SAM,CAT,MAR,GD),Totselstemps)
colnames(SuppTable3) <-c("Watermass","Treatment","NCP_R","Intercept","uIntercept","Slope","uSlope","n","Adj.R2","p")

SuppTable3[,4:10] <- round(SuppTable3[,4:10],2)


Zero24 <- rbind(rectaAMZRKL24,rectaAMZRKD24,rectaAMZRCO2L24,rectaAMZRCO2D24,rectaSAMKL24,rectaSAMKD24,rectaSAMCO2L24,rectaSAMCO2D24,rectaCATKL24,rectaCATKD24,rectaCATCO2L24,rectaCATCO2D24,rectaMARKL24,rectaMARKD24,rectaMARCO2L24,rectaMARCO2D24,rectaGDKL24,rectaGDKD24,rectaGDCO2L24,rectaGDCO2D24)

SuppTable3b <- data.frame(rbind(AMZR,SAM,CAT,MAR,GD),Zero24)
colnames(SuppTable3b) <-c("Watermass","Treatment","NCP_R","Intercept","uIntercept","Slope","uSlope","n","Adj.R2","p")
SuppTable3b[,4:10] <- round(SuppTable3b[,4:10],2)

Veinticuatro48 <- rbind(rectaAMZRKL48,rectaAMZRKD48,rectaAMZRCO2L48,rectaAMZRCO2D48,rectaSAMKL48,rectaSAMKD48,rectaSAMCO2L48,rectaSAMCO2D48,rectaCATKL48,rectaCATKD48,rectaCATCO2L48,rectaCATCO2D48,rectaMARKL48,rectaMARKD48,rectaMARCO2L48,rectaMARCO2D48,rectaGDKL48,rectaGDKD48,rectaGDCO2L48,rectaGDCO2D48)

SuppTable3c <- data.frame(rbind(AMZR,SAM,CAT,MAR,GD),Veinticuatro48)
colnames(SuppTable3c) <-c("Watermass","Treatment","NCP_R","Intercept","uIntercept","Slope","uSlope","n","Adj.R2","p")
SuppTable3c[,4:10] <- round(SuppTable3c[,4:10],2)

#Esto nos muestra que las diferencias entre las pendientes estimadas entre 0-24 y todos los puntos no son muy diferentes. 
mean(SuppTable3$Slope)
sd(SuppTable3$Slope)
mean(SuppTable3b$Slope)
sd(SuppTable3b$Slope)
mean(SuppTable3c$Slope)
sd(SuppTable3c$Slope)

Table3Tot <- data.frame(rbind(SuppTable3,SuppTable3b,SuppTable3c))
times_slope <- rep(c("0-24-48","0-24","24-48"),times=1,each=20)
Table3Tot2 <- data.frame(cbind(times_slope,Table3Tot))
AMZR<- Table3Tot2[Table3Tot2$Watermass=="AMZR",]
SAM<- Table3Tot2[Table3Tot2$Watermass=="SAM",]
CAT<- Table3Tot2[Table3Tot2$Watermass=="CAT",]
MAR<- Table3Tot2[Table3Tot2$Watermass=="MAR",]
GD<- Table3Tot2[Table3Tot2$Watermass=="GD",]

SuppTable3end <- Table3Tot2
##Tabla 3 buena!!
write.csv(SuppTable3end,"SuppTable3.csv",row.names=F)

#Si la distribución de los slopes fuese normal, los siguientes modelos muestran que la relación entre todos los puntos y lo 0-24 o 24-48 no son significativamente distintas, dado que si le quitamos un offset de 1 (pendiente distinta de 1), la relación pasa por el intercepto (no es significativamente distinto de 0 el intercepto) y la pendiente no es significativa. 

barup <- function(x,y,yse)
{epsilon=0.002
for(i in 1:length(yse))
    {segments(x[i],y[i]+yse[i],x[i],y[i]-0, col="gray80")
     segments(x[i]-epsilon, y[i]+yse[i] , x[i]+epsilon, y[i]+yse[i], col="gray80")        
    }}

bardown <- function(x,y,yse)
{epsilon=0.002
for(i in 1:length(xse))
    {segments(x[i],y[i]-yse[i],x[i],y[i], col="gray80")
     segments(x[i]-epsilon, y[i]-yse[i] , x[i]+epsilon, y[i]-yse[i], col="gray80")        
    }}

baright <-function(x,y,xse) 
{epsilon=0.002
for(i in 1:length(yse))
    {segments(x[i]+xse[i],y[i],x[i],y[i], col="gray80")
     segments(x[i]+xse[i], y[i]+epsilon,x[i]+xse[i], y[i]-epsilon, col="gray80")
    }}

barleft <-function(x,y,xse) 
{epsilon=0.002
for(i in 1:length(yse))
    {segments(x[i]-xse[i],y[i],x[i],y[i], col="gray80")
     segments(x[i]-xse[i], y[i]+epsilon,x[i]-xse[i], y[i]-epsilon, col="gray80")
    }}


tapply(SuppTable3end$Adj.R2,SuppTable3end$times_slope,mean)##entre 0-48 hay pocos valores significativos!


####Fig 3A
pdf(file="/home/evaristo/Documentos/MOC2/REV/ScientiaMarina/Figure3A.pdf",width=8,height=8)

par(mfrow=c(1,1))
par(mar=c(6,6.5,4,2)+0.1)
plot(SuppTable3b$Slope~SuppTable3$Slope, xlim=c(-0.35,0),ylim=c(-0.35,0),xlab="",ylab="",yaxt="n",cex=1.5,pch=19)
box(lwd=2.5)
mtext(bquote(bold(""*A*"")),cex=2, side=3, line=-2.5,las=1, adj=0.07, padj=0)
mtext(bquote(bold(""*Metabolism[0-24-48]*" ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")),side=c(1), line=c(4),cex=1.5)
mtext(bquote(bold(""*Metabolism[0-24|24-48]*" ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")),side=c(2), line=c(4),cex=1.5)
axis(side=2,las=2)
model <- lm(y~x)
x <- SuppTable3$Slope
y <- SuppTable3b$Slope
model <- lm(y~x)
ablinelimitMOC(model,"black")
confintcurv2secol(model,"black")### Formula de Fun para lineas de 2 se
xse <-SuppTable3$uSlope
yse <- SuppTable3b$uSlope
barup(x,y,yse)
baright(x,y,xse)
fit2 <- lm(SuppTable3b$Slope~SuppTable3$Slope +offset(1*SuppTable3$Slope))
x <- SuppTable3$Slope
y <- SuppTable3c$Slope
model2 <- lm(y~x)
fit3 <- lm(SuppTable3c$Slope~SuppTable3$Slope +offset(1*SuppTable3$Slope))
ablinelimitMOC(model2,"gray50")
confintcurv2secol(model2,"gray50")### Formula de Fun para lineas de 2 se
xse <-SuppTable3$uSlope
yse <- SuppTable3c$uSlope
bardown(x,y,yse)
barleft(x,y,xse)
points(SuppTable3b$Slope~SuppTable3$Slope, xlim=c(-0.35,0),ylim=c(-0.35,0),col="black",cex=1.5,pch=19)
points(SuppTable3c$Slope~SuppTable3$Slope, xlim=c(-0.35,0),ylim=c(-0.35,0),col="gray50",cex=1.5,pch=19)
summary(model)
summary(fit2)
summary(model2)
summary(fit3)

dev.off()

slopes <- c(SuppTable3$Slope,SuppTable3b$Slope,SuppTable3c$Slope)
se <- c(SuppTable3$uSlope,SuppTable3b$uSlope,SuppTable3c$uSlope)
tiempos <- as.factor(rep(c("0-24-48","0-24","24-48"),each=20))
tiempos_slopes<- data.frame(SuppTable3$Watermass,SuppTable3$Treatment,SuppTable3$NCP_R,tiempos,slopes,se)
colnames(tiempos_slopes)[1] <- "Watermass"
colnames(tiempos_slopes)[2] <- "Treatment"
colnames(tiempos_slopes)[3] <- "Light"
tiempos_slopes$Treatment <- as.factor(tiempos_slopes$Treatment)
#tiempos_slopes <- tiempos_slopes[!tiempos_slopes$tiempos=="24-48",]

slopes <- tiempos_slopes$slopes
tiempos <- tiempos_slopes$tiempos
Light <- tiempos_slopes$Light
Treatment <- tiempos_slopes$Treatment

model <- aov(slopes~tiempos)
summary(model)

model2 <- lm(slopes~tiempos)
summary(model2)
model3 <- lm(slopes~tiempos+Light)
summary(model3)
model4 <- lm(slopes~tiempos+Treatment)
summary(model4)
model5 <- lm(slopes~tiempos+Light+Treatment)
summary(model5)

anova(model2,model3,model4,model5)
AIC(model2,model3,model4,model5)

## si separamos los tiempos 0-24 y 24-48 tampoco salen diferencias significativas de pendientes
TukeyHSD(model, conf.level=.95)##existen diferencias significativas entre 0-24 y 24-48, pero esto asume normalidad de las variables y esas variables no son normales
shapiro.test(slopes)#sin embargo, la distribución no es normal
kruskal.test(slopes~tiempos)#Por tanto las diferencias entre pendientes detectadas no son significativas
pairwise.wilcox.test(slopes, tiempos, p.adjust="bonferroni")
#test max-t test para diferencias entre medias bajo heteroscedascitidad (Herberich et al. 2010)
model_glht <- glht(model, mcp(tiempos= "Tukey"), vcov = vcovHC)
summary(model_glht)#Hay ligeras diferencias entre 0-24 y 24-48, pero no són significativas!!!!!!!!!!
## si separamos los tiempos 0-24 y 24-48 tampoco salen diferencias significativas de pendientes
tiempos_slopes2 <- tiempos_slopes[21:60,]
wilcox.test(tiempos_slopes2$slopes ~ tiempos_slopes2$tiempos, correct=FALSE)
pairwise.wilcox.test(tiempos_slopes2$slopes, tiempos_slopes2$tiempos, p.adjust="bonferroni")
## Es decir que aunque las pendientes entre 24 y 48 tienen un slope promedio menos negativo que las pendientes entre 0-24 y 0-24-48, no podemos decir que esa diferencia sea significativa. Resumiendo, las tasas de consumo de oxígeno fueron las mísmas en los tres periodos de tiempo.

model6 <- aov(model5)
summary(model6)
model6_glht <- glht(model6, mcp(tiempos= "Tukey"), vcov = vcovHC)
summary(model6_glht)#Hay ligeras diferencias entre 0-24 y 24-48, pero no són significativas!!!!!!!!!!
model7_glht <- glht(model6, mcp(Treatment= "Tukey"), vcov = vcovHC)
summary(model7_glht)#Hay ligeras diferencias entre CO2 y Ctrl, pero no són significativas!!!!!!!!!!

tapply(tiempos_slopes$slopes,list(tiempos_slopes$tiempos,tiempos_slopes$Treatment,tiempos_slopes$Light),mean)
tapply(tiempos_slopes$slopes,list(tiempos_slopes$tiempos,tiempos_slopes$Light),mean)

###Quitando datos 24-48
tiempos_slopes <- tiempos_slopes[!tiempos_slopes$tiempos=="24-48",]
slopes <- tiempos_slopes$slopes
tiempos <- tiempos_slopes$tiempos
Light <- tiempos_slopes$Light
Treatment <- tiempos_slopes$Treatment
Watermass <- tiempos_slopes$Watermass

model2 <- lm(slopes~Light+Treatment)
summary(model2)## Este modelo explicaria un 33% de la variabilidad, existiendo un efecto de la luz y de la acidificación.

model3 <- lm(slopes~Watermass*Light+Treatment)
summary(model3)## Este modelo explicaria un 41% de la variabilidad, existiendo un efecto de la luz y de la acidificación. Sin embargo, no es significativvamente distinto del modelo con las dos interacciones, tiene un AIC menor y nos dice que hay diferencias entre luz y oscuridad que és lo que uno espera de un experimento de respiración.
model4 <- lm(slopes~Watermass*Light*Treatment)
summary(model4)#explicaria un 66% de la variabilidad. pero no es distinto del anterior y tiene una AIC mayor
AIC(model3,model4)



#######
### resumen de los slopes. Diferencias entre tasas de CR, NPP i GPP
resumen_slopes <- tapply(tiempos_slopes$slopes,list(tiempos_slopes$Watermass,tiempos_slopes$Treatment,tiempos_slopes$Light,tiempos_slopes$tiempos),mean)
resumen_se <- tapply(tiempos_slopes$se,list(tiempos_slopes$Watermass,tiempos_slopes$Treatment,tiempos_slopes$Light,tiempos_slopes$tiempos),mean)

R <- resumen_slopes[,,"Dark",]
Rse <- resumen_se[,,"Dark",]
NCP <- resumen_slopes[,,"Light",]
NCPse <- resumen_se[,,"Light",]
GPP <- NCP-R##Gonzalez et al. 2008
GPPse <- round(sqrt(NCPse^2+Rse^2),2)### Calculada la propagación de errores con esta fórmula!####

####Fig 3B
pdf(file="/home/evaristo/Documentos/MOC2/REV/ScientiaMarina/Figure3B.pdf",width=8,height=8)

RControl <- R[,"Control",]
RCO2 <- R[,"High CO2",]
par(mar=c(6,6.5,4,2)+0.1)
plot(RControl~RCO2,ylim=c(-0.35,.25),xlim=c(-0.35,.25),xlab="",ylab="",yaxt="n",cex=1.5,pch=19)
box(lwd=2.5)
mtext(bquote(bold(""*B*"")),cex=2, side=3, line=-2.5,las=1, adj=0.07, padj=0)
mtext(bquote(bold(""*Metabolism[CO[2]]*" ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")),side=c(1), line=c(4),cex=1.5)
mtext(bquote(bold(""*Metabolism*" ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")),side=c(2), line=c(4),cex=1.5)
axis(side=2,las=2)
modelR <- lm(as.vector(RControl)~as.vector(RCO2))
abline(0,1,lty=2)
#abline(modelR)
text(RCO2,RControl,rownames(RCO2),cex=0.50,pos=2,col="black")
xse <-Rse[,"High CO2",]
yse <- Rse[,"Control",]
bardown(RCO2,RControl,yse)
barleft(RCO2,RControl,xse)
points(RControl~RCO2,ylim=c(-0.35,.25),xlim=c(-0.35,.25),xlab="",ylab="",yaxt="n",cex=1.5,pch=19)

NCPControl <- NCP[,"Control",]
NCPCO2 <- NCP[,"High CO2",]
points(NCPControl~NCPCO2,ylim=c(-0.35,.3),xlim=c(-0.35,.25),col="green",cex=1.5,pch=19)
modelNCP <- lm(as.vector(NCPControl)~as.vector(NCPCO2))
#abline(modelNCP,col="green")
text(NCPCO2,NCPControl,rownames(NCPCO2),cex=0.50,pos=1,col="green")
xse <-NCPse[,"High CO2",]
yse <- NCPse[,"Control",]
bardown(NCPCO2,NCPControl,yse)
barleft(NCPCO2,NCPControl,xse)
points(NCPControl~NCPCO2,ylim=c(-0.35,.3),xlim=c(-0.35,.25),col="green",cex=1.5,pch=19)

GPPControl <- GPP[,"Control",]
GPPCO2 <- GPP[,"High CO2",]
points(GPPControl~GPPCO2,ylim=c(-0.35,.3),xlim=c(-0.35,.25),col="blue",cex=1.5,pch=19)
modelGPP <- lm(as.vector(GPPControl)~as.vector(GPPCO2))
#abline(modelGPP,col="blue")
text(GPPCO2,GPPControl,rownames(GPPCO2),cex=0.50,pos=1,col="blue")
xse <-GPPse[,"High CO2",]
yse <- GPPse[,"Control",]
barup(GPPCO2,GPPControl,yse)
baright(GPPCO2,GPPControl,xse)
points(GPPControl~GPPCO2,ylim=c(-0.35,.3),xlim=c(-0.35,.25),col="blue",cex=1.5,pch=19)

dev.off()

c(mean(as.vector(RControl)),mean(as.vector(RCO2)))
wilcox.test(RCO2, RControl, correct=FALSE)
c(mean(as.vector(NCPControl)),mean(as.vector(NCPCO2)))
wilcox.test(NCPCO2, NCPControl, correct=FALSE)
c(mean(as.vector(GPPControl)),mean(as.vector(GPPCO2)))
wilcox.test(GPPCO2, GPPControl, correct=FALSE)
## Es decir no hay diferencias entre NPC y GPP debido a la adición de CO2, pero si hay diferencias de respiración

tapply(AMZR$Slope,list(AMZR$Watermass,AMZR$Treatment,AMZR$NCP_R),mean)
tapply(SAM$Slope,list(SAM$Watermass,SAM$Treatment,SAM$NCP_R),mean)
tapply(CAT$Slope,list(CAT$Watermass,CAT$Treatment,CAT$NCP_R),mean)
tapply(MAR$Slope,list(MAR$Watermass,MAR$Treatment,MAR$NCP_R),mean)
tapply(GD$Slope,list(GD$Watermass,GD$Treatment,GD$NCP_R),mean)



###plot de la distribución de las tasas de consumo/producción neta de oxígeno
library("plyr")

slopes <- c(SuppTable3$Slope,SuppTable3b$Slope,SuppTable3c$Slope)
tiempos <- as.factor(rep(c("0-24-48","0-24","24-48"),each=20))
tiempos_slopes<- data.frame(SuppTable3$Treatment,tiempos,slopes)
colnames(tiempos_slopes)[1] <- "Treatment"
O2rate <- ddply(tiempos_slopes, "tiempos", summarise, grp.mean=mean(slopes))
head(O2rate)

p<-ggplot(tiempos_slopes, aes(x=slopes, fill = tiempos,color=tiempos))+
    geom_histogram(aes(y=..density..), position="identity", breaks=seq(-0.3,0,by=0.01))+
    geom_density(alpha=0.7)+
    geom_vline(data=O2rate, aes(xintercept=grp.mean, color=tiempos),linetype="dashed",lwd=1)+
    scale_color_manual(values=c("black", "gray30", "gray80"),name="Period (h)")+
    scale_fill_manual(values=c("black", "gray30", "gray80"),name="")+
    labs(title="",x=expression(bold("Oxygen metabolic rate ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")), y = "Density")+
    theme_classic()+
    theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.margin = margin(1,1,1,1, "cm"))
p <- p+guides(fill=guide_legend(title="Period (h)"))+geom_text(x=-0.30, y=20, label="All",colour="black",size=4)


slopes <- c(SuppTable3$Slope,SuppTable3b$Slope,SuppTable3c$Slope)
tiempos <- as.factor(rep(c("0-24-48","0-24","24-48"),each=20))
tiempos_slopes<- data.frame(SuppTable3$Treatment,tiempos,slopes)
colnames(tiempos_slopes)[1] <- "Treatment"
tiempos_slopes <- tiempos_slopes[tiempos_slopes$Treatment=="Control",]
O2rate <- ddply(tiempos_slopes, "tiempos", summarise, grp.mean=mean(slopes))
head(O2rate)

p2<-ggplot(tiempos_slopes, aes(x=slopes, fill = tiempos,color=tiempos))+
    geom_histogram(aes(y=..density..), position="identity", breaks=seq(-0.3,0,by=0.01))+
    geom_density(alpha=0.7)+
    geom_vline(data=O2rate, aes(xintercept=grp.mean, color=tiempos),linetype="dashed",lwd=1)+
    scale_color_manual(values=c("black", "gray30", "gray80"),name="Period (h)")+
    scale_fill_manual(values=c("black", "gray30", "gray80"),name="")+
    labs(title="",x=expression(bold("Oxygen metabolic rate ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")), y = "Density")+
    theme_classic()+
    theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.margin = margin(1,1,1,1, "cm"))
p2 <- p2+guides(fill=guide_legend(title="Period (h)"))+geom_text(x=-0.28, y=29, label="Control",colour="black",size=4)

slopes <- c(SuppTable3$Slope,SuppTable3b$Slope,SuppTable3c$Slope)
tiempos <- as.factor(rep(c("0-24-48","0-24","24-48"),each=20))
tiempos_slopes<- data.frame(SuppTable3$Treatment,tiempos,slopes)
colnames(tiempos_slopes)[1] <- "Treatment"
tiempos_slopes<- tiempos_slopes[tiempos_slopes$Treatment=="High CO2",]
O2rate <- ddply(tiempos_slopes, "tiempos", summarise, grp.mean=mean(slopes))
head(O2rate)

p3<-ggplot(tiempos_slopes, aes(x=slopes, fill = tiempos,color=tiempos))+
    geom_histogram(aes(y=..density..), position="identity", breaks=seq(-0.3,0,by=0.01))+
    geom_density(alpha=0.7)+
    geom_vline(data=O2rate, aes(xintercept=grp.mean, color=tiempos),linetype="dashed",lwd=1)+
    scale_color_manual(values=c("black", "gray30", "gray80"),name="Period (h)")+
    scale_fill_manual(values=c("black", "gray30", "gray80"),name="")+
    labs(title="",x=expression(bold("Oxygen metabolic rate ("*mu*mol*" "*O[2]*" "*L^-1*" "*h^-1*")")), y = "Density")+
    theme_classic()+
    theme(axis.text.y   = element_text(size=14),
        axis.text.x   = element_text(size=14),
        axis.title.y  = element_text(size=14),
        axis.title.x  = element_text(size=14),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        plot.margin = margin(1,1,1,1, "cm"))
p3 <- p3+guides(fill=guide_legend(title="Period (h)"))+geom_text(x=-0.30, y=29, label=expression("CO"[2]),colour="black",size=4)


####Fig TaxDist
pdf(file="/home/evaristo/Documentos/MOC2/REV/ScientiaMarina/FigureTaxDist.pdf",width=8,height=8)
grid.arrange(p2,p3,nrow=2)
dev.off()

#AMZR Control Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
AMZRKL <- c(197.94,197.27,197.75,197.47,195.06,194.57,194.56,194.26,191.60,192.22,190.73,194.53)
#AMZR Control Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
AMZRKD <- c(NA,197.18,197.36,197.25,194.19,196.50,195.80,194.53,191.40,191.90,191.35,188.69)
#AMZR High CO2 Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
AMZRCO2L <-c(197.96,197.40,198.16,197.91,195.46,193.13,195.47,195.54,193.61,195.58,191.91,NA)
#AMZR High CO2 Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
AMZRCO2D <-c(198.48,198.53,197.71,197.18,188.47,189.68,194.40,190.69,185.05,187.94,186.59,188.93)

#SAM Control Light, SAMKL
Time <- c(0,0,0,0,24,24,24,24,47,47,47,47)
SAMKL <- c(205.40,205.72,NA,206.21,203.83,203.92,202.90,202.82,196.34,198.02,201.65,202.40)
#SAM Control Dark (oxygen in micM O2)
Time <- c(0,0,0,0,24,24,24,24,47,47,47,47)
SAMKD <- c(205.36,205.29,205.28,205.43,203.62,202.67,NA,200.31,200.12,201.58,NA,199.40)
#SAM High CO2 light
Time <- c(0,0,0,0,24,24,24,24,47,47,47,47)
SAMCO2L <- c(204.26,NA,204.32,NA,199.74,200.09,201.75,201.54,199.34,197.02,199.33,197.99)
#SAM HighCO2 Dark
Time <- c(0,0,0,0,24,24,24,24,47,47,47,47)
SAMCO2D <- c(206.43,203.34,198.79,203.51,194.11,195.41,196.71,194.65,198.24,193.16,191.68,NA)

#CAT Control Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
CATKL <-c(205.80,204.87,204.19,204.31,200.55,202.82,202.38,201.54,200.10,201.12,198.89,201.99)
#CAT Control Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
CATKD <-c(205.18,203.23,203.89,203.32,200.80,201.51,200.69,200.61,192.43,197.36,195.86,196.36)
#CAT High CO2 Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
CATCO2L <-c(202.89,203.40,203.45,203.06,200.67,199.48,200.19,200.95,196.58,196.98,197.84,195.91)
#CAT High CO2 Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
CATCO2D <-c(203.85,203.15,203.40,203.13,200.53,200.33,199.83,199.45,195.80,194.97,196.09,196.41)

#MAR Control Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
MARKL <-c(217.78,221.08,223.23,224.02,218.95,219.55,220.16,NA,216.24,217.73,216.44,217.48)
#MAR Control Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
MARKD <-c(222.71,222.41,222.25,222.85,219.29,218.50,218.64,218.34,217.61,215.23,NA,220.12)
#MAR High CO2 Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
MARCO2L <-c(222.89,222.70,223.33,223.60,223.44,217.16,220.09,217.48,216.99,215.41,217.45,217.09)
#MAR High CO2 Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
MARCO2D <-c(NA,222.98,223.28,223.38,217.51,217.15,216.32,217.28,215.04,214.27,213.67,216.20)

#GD Control Light
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
GDKL <-c(200.87,192.37,205.42,191.06,193.91,193.49,NA,192.32,190.95,190.35,186.23,190.95)
#GD Control Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
GDKD <-c(198.97,198.56,198.34,198.46,193.92,NA,193.26,194.36,188.51,190.49,191.14,191.86)
#GD High CO2 ligth
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
GDCO2L <-c(202.58,201.47,201.12,201.20,196.23,195.84,195.50,195.29,194.92,190.27,192.08,191.56)
#GD High CO2 Dark
Time <- c(0,0,0,0,24,24,24,24,48,48,48,48)
GDCO2D <- c(200.49,NA,200.57,203.51,194.48,195.75,194.25,196.60,192.71,191.30,191.21,193.08)

Stations <- rep(c("AMZR","SAM","CAT","MAR","GD"),each=12,times=4)
Time <- c(rep(c(0,24,48),each=4),rep(c(0,24,47),each=4),rep(c(0,24,48),each=4),rep(c(0,24,48),each=4),rep(c(0,24,48),each=4))
Time2 <- rep(Time,times=4)
Treatment <- rep(c("K","CO2"),each=120)# K = Control, CO2 high CO2
Light <- rep(c("L","D"),each=60,times=2)# L = light, D = Dark
Oxygen <- c(AMZRKL,
            SAMKL,
            CATKL,
            MARKL,
            GDKL,
            AMZRKD,
            SAMKD,
            CATKD,
            MARKD,
            GDKD,
            AMZRCO2L,
            SAMCO2L,
            CATCO2L,
            MARCO2L,
            GDCO2L,
            AMZRCO2D,
            SAMCO2D,
            CATCO2D,
            MARCO2D,
            GDCO2D)

Metabolism <- data.frame(Stations,Time2,Treatment,Light,Oxygen)# L = light, D = Dark
Metabolism[order(Metabolism$Stations),]
Metabolism<-na.omit(Metabolism)
Metabolism$Stations <- as.factor(Metabolism$Stations)
Metabolism$Light <- as.factor(Metabolism$Light)
Metabolism$Treatment <- as.factor(Metabolism$Treatment)
  
#plot lattice de todo
#Gráfica ##1
library(lattice)

xyplot(Metabolism$Oxygen~Metabolism$Time2 | factor(Metabolism$Treatment,levels=c("K","CO2")):factor(Metabolism$Light,levels=c("L","D")):factor(Metabolism$Stations,levels=c("AMZR","SAM","CAT","MAR","GD")),as.table=TRUE,group=factor(Metabolism$Treatment):factor(Metabolism$Light),xlab=expression(bold("Time (h)")),ylab=expression(bold("Oxygen ("*mu*mol*" "*O[2]*" "*L^-1*")")),par.settings=list(layout.widths=list(left.padding=4, right.padding=1),layout.heights=list(top.padding=1, bottom.padding=5),superpose.symbol = list(pch = 19, cex = 1.5,col=c("red","pink","blue","cyan"))),panel = function(x,y, ...)
{
    panel.xyplot(x, y, ...)
    panel.lmline(x, y,col.line = "gray20",lty=2,...)
},
strip=strip.custom(factor.levels=c("Ctrl L AMZR","Ctrl L SAM","Ctrl L CAT","Ctrl L MAR","Ctrl L GD","Ctrl D AMZR","Ctrl D SAM","Ctrl D CAT","Ctrl D MAR","Ctrl D GD",expression("C"*O[2]*" L AMZR"),expression("C"*O[2]*" L SAM"),expression("C"*O[2]*" L CAT"),expression("C"*O[2]*" L MAR"),expression("C"*O[2]*" L GD"),expression("C"*O[2]*" D AMZR"),expression("C"*O[2]*" D SAM"),expression("C"*O[2]*" D CAT"),expression("C"*O[2]*" D MAR"),expression("C"*O[2]*" D GD"))))

### Modelo lineal de todo
model <- lm(Metabolism$Oxygen~Metabolism$Time2+Metabolism$Stations*Metabolism$Light*Metabolism$Treatment)
model1 <-lm(Metabolism$Oxygen~Metabolism$Time2*Metabolism$Stations*Metabolism$Light*Metabolism$Treatment)
model2 <-lm(Metabolism$Oxygen~Metabolism$Time2)
summary(model)
summary(model1)
summary(model2)
plot(model1,which=c(1),col=1,add.smooth=FALSE)
AIC(model,model1,model2)
anova(model,model1,model2)

op <- par(mfrow=c(2,2),mar=c(4,4,2,2))
plot(model1,which=c(1),col=1,add.smooth=FALSE)
plot(as.numeric(resid(model1))~as.factor(model1$model[,3]))
plot(as.numeric(resid(model1))~as.factor(model1$model[,5]))
plot(as.numeric(resid(model1))~as.factor(model1$model[,4]))

hist(as.numeric(resid(model1)))
bartlett.test(as.numeric(resid(model1)),as.factor(model1$model[,3]))## hay dif significativas en la varianza del modelo para las diferentes estaciones 
bartlett.test(as.numeric(resid(model1)),as.factor(model1$model[,4]))# no hay diferencias en el caso de la luz
bartlett.test(as.numeric(resid(model1)),as.factor(model1$model[,5]))# no hay diferencias en el caso del tratamiento

## Diferencias en las conc. de oxígeno inicial segun las diferentes estaciones para K
OxygenZero <- Metabolism[Metabolism$Time2==0,]
KOxygenZero <- OxygenZero[OxygenZero$Treatment=="K",c(1,5)]
KOxygenZero[,1] <-as.factor(KOxygenZero[,1]) 
kruskal.test(KOxygenZero$Oxygen~KOxygenZero$Stations)
Kmodel <- aov(Oxygen~Stations, data=KOxygenZero)
TukeyHSD(Kmodel)

#test max-t test para diferencias entre medias bajo heteroscedascitidad (Herberich et al. 2010)
tapply(OxygenZero$Oxygen,list(OxygenZero$Treatment,OxygenZero$Station),mean)
Kmodel_glht <- glht(Kmodel, mcp(Stations= "Tukey"), vcov = vcovHC)
summary(Kmodel_glht)#No hubo diferencias en el oxígeno inicial de los K  entre GD y AMZR

CO2OxygenZero <- OxygenZero[OxygenZero$Treatment=="CO2",c(1,5)]
CO2OxygenZero[,1] <-as.factor(CO2OxygenZero[,1]) 
CO2model <- aov(Oxygen~Stations, data=CO2OxygenZero)
CO2model_glht <- glht(CO2model, mcp(Stations= "Tukey"), vcov = vcovHC)
summary(CO2model_glht)#No hubo diferencias en los tratamientos CO2 entre el oxígeno inicial de SAM i CAT o SAM i DG 

## no hubo diferencias significativas de oxígeno inicial debido al burbujeo de CO2
Metabolism <- data.frame(Stations,Time2,Treatment,Light,Oxygen)# L = light, D = Dark
Metabolism[order(Metabolism$Stations),]
plot(Metabolism$Oxygen[Metabolism$Treatment=="CO2"]~Metabolism$Oxygen[Metabolism$Treatment=="K"])
model <- lm(Metabolism$Oxygen[Metabolism$Treatment=="CO2"]~Metabolism$Oxygen[Metabolism$Treatment=="K"])
abline(model,col="red")
summary(model)

####
##Resum valors oxygen inicial al treatment K i CO2
OxygenZero <- Metabolism[Metabolism$Time2==0,]
TreatmentLight <- paste(OxygenZero$Treatment,OxygenZero$Light)
OxygenZero <- data.frame(OxygenZero[,c(1:3)],TreatmentLight,OxygenZero[,c(4,5)])
OxygenZero <- na.omit(OxygenZero)
tapply(OxygenZero$Oxygen,list(OxygenZero$TreatmentLight,OxygenZero$Stations),mean)
hist(OxygenZero$Oxygen)
shapiro.test(OxygenZero$Oxygen)
model <- aov(OxygenZero$Oxygen~OxygenZero$Stations*OxygenZero$Treatment*OxygenZero$Light)
summary(model)
par(mfrow = c(2, 2))
plot(model,  add.smooth = FALSE)## hay violaciones de la homogeneidad, por eso salen outliers. Aplicaremos gls



library(nlme)## como se viola la homogeneidad hacemos un gls
y <-OxygenZero$Oxygen
x <- as.factor(OxygenZero$Stations)
w <-as.factor(OxygenZero$Light)
z <-as.factor(OxygenZero$Treatment)

####################### Permutación para ver si K y CO2 tienen la misma distribución, en principio la distribución del análisis deberia ser normal, pero como hay pocas muestras no tiene porqué dar significativo siempre
### Establecemos una permutación y miramos si las diferencias entre las muestras iniciales siguen normal, para ello separamnos tratamientos K o CO2, no así las de luz oscuridad por ser réplicas y el tratamiento no suponer ninguna manipulación que pueda añadir o extraer gas a las muestras iniciales.
perm = function(z, y) {
  # turn x to a character, easier to deal with
  x = as.character(z)
  # shuffle the x values:
  x_shuff = sample(z)
  # calculate the mean of each group:
  x_bar_add = mean(y[x_shuff == "CO2"])
  x_bar_ctl = mean(y[x_shuff == "K"])
  # calculate the difference:
  x_bar_add - x_bar_ctl
}
z <-OxygenZero[,3]
y <-OxygenZero[,6]     
Dnull = replicate(n = 5000, expr = perm(z, y))
Dobs = mean(y[z == "CO2"]) - mean(y[z == "K"])
hist(Dnull, col = "grey")
abline(v = Dobs, col = "blue", lwd = 3, lty = 2)
shapiro.test(Dnull)## la distribución de las diferencias dentre K y CO2 no es distinta de una distribución normal, por tanto las distribuciones iniciales de K y CO2 són las mismas== se puede aplicar test lineales por el teorema central del limite 
mean(abs(Dnull) <= Dobs)## H0 es cierta y, por tanto, no hay diferencias entre las dos poblaciones de diferencias, lo mismo de antes.

y <-OxygenZero$Oxygen
x <- as.factor(OxygenZero$Stations)
w <-as.factor(OxygenZero$Light)
z <-as.factor(OxygenZero$Treatment)

f1 <- formula(y~x*w*z)
M0 <-gls(f1,weights=varIdent(form=~1|x*z*w),method="REML",data=OxygenZero)
## se han corregido las violaciones de la homogeneidad haciendo varianzas distintas para cada interacción
summary(M0)
anova(M0)
M01 <- update(M0,.~.-w:z)
summary(M01)
anova(M01)
M02 <- update(M01,.~.-x:w:z)
summary(M02)
anova(M02)
M03 <- update(M02,.~.-x:w)
summary(M03)
M04 <- update(M03,.~.-w)
summary(M04)
anova(M04)
anova(M01,M02,M03,M04)
model <- aov(M04)
par(mfrow=c(2,2))
plot(model)


#esto nos dice que hay diferencias significativas entre estaciones, tratamientos y una interacción entre estaciones y tratamientos.
## La adición de CO2 extrajo -0.598 micromoles O2 de las muestras iniciales control p<0.05
## En CAT se adicionó a las muestras iniciales +1.78 micromoles O2 
## En GD se desgasificó de las muestras iniciales -2.6640 micromoles O2 
## En SAM se adicionó a las muestras iniciales +1.65 micromoles O2 
resultadoestacion <- data.frame(tapply(OxygenZero$Oxygen,list(OxygenZero$Station,OxygenZero$Treatment),mean))
difestacion <- resultadoestacion[,2]-resultadoestacion[,1]
difmedidopromedio <- mean(resultadoestacion[,2]-resultadoestacion[,1])


colorPalette <- c("hotpink","red","steelblue1","navy")
p <- ggplot(OxygenZero , aes(x=Stations,y=Oxygen,color=TreatmentLight))
p2 <- p+geom_boxplot(position=position_dodge(0.8))+geom_jitter(position=position_dodge(0.8))+scale_x_discrete(limits=c("AMZR","SAM","CAT","MAR","GD"))+
    theme(
            panel.background=element_rect(fill="white", colour="black", size=1),
            panel.grid=element_blank(),
            axis.text.x = element_text(size=rel(2)),
            axis.text.y = element_text(size=rel(2)),
            axis.title=element_text(size=rel(2), face="bold"),
            #legend.title = element_text("Treatment",size=rel(1.5)),
            legend.text = element_text(size=rel(1.5)),legend.position = "top",legend.key=element_rect(fill="white",colour="white"),plot.margin=margin(1,1,1,1,unit="cm"))+xlab("\n Experiment")+ylab(expression(bold("Oxygen ("*mu*mol*" "*O[2]*" "*L^-1*")")))+ ylim(190,225)+scale_color_manual(name="",values=c(colorPalette),labels=c("CO\u2082 L","CO\u2082 D","Ctrl L","Ctrl D"))#para poner letra dentro del plot + annotate("text", label='bold("A")',size=10, x=1, y=225,parse=TRUE) 
p2

####Fig 2
png(file="/home/evaristo/Documentos/MOC2/REV/ScientiaMarina/Figure2.png",width=1000,height=1000)
p2
dev.off()

#En resumen, las muestras vienen de dos distribuciones identicas pero no normales (permutación).  Al hacer las graficas con los datos originales vemos que hay heterogeneidad. Al aplicar un test gls y corregir para la heterogeineidad nos salen diferencias significativas, entre: estaciones, entre control y CO2, y entre estaciones y control o CO2 (M04).

## time courses
## time courses
## time courses
Metabolism <- data.frame(Stations,Time2,Treatment,Light,Oxygen)# L = light, D = Dark
TreatmentLight <- as.factor(rep(c("KL","KD","CO\u2082L","CO\u2082D"),each=60))
Metabolism <- data.frame(Metabolism[,-5],TreatmentLight,Metabolism[,5])
colnames(Metabolism)[c(5,6)] <- c("TreatmentLight","Oxygen")
Metabolism <- na.omit(Metabolism)


plotregtrozosK <- function(AMZRtotK,minimum,maximum,z='bold("A")')
{
df1 <- AMZRtotK
colorPalette <- c("red","pink","navy","steelblue1")

q1 <- ggplot(df1,aes(x=Time2,y=Oxygen,color=df1$TreatmentLight))+geom_point(size=3)+scale_color_manual(name="",values=c(colorPalette),breaks=unique(levels(TreatmentLight)),labels=c("CO\u2082_D","CO\u2082_L","Control_D","Control_L"))+
                    theme(
            panel.background=element_rect(fill="white", colour="black", size=0.75),
            panel.grid.major=element_line(size=0,"black",linetype=2),
            axis.text.x = element_text(size=rel(2)),
            axis.text.y = element_text(size=rel(2)),
            axis.title=element_text(size=rel(0), face="bold"),
            #legend.title = element_text("Treatment",size=rel(1.5)),
            legend.position="none",plot.margin=margin(1,1,1,1,unit="cm"))+
                                    annotate("text", label=z,size=4, x=5, y=minimum+0.5,parse=TRUE)+ ylim(minimum,maximum)

df2 <- df1[df1$Time2<=48&df1$Light=="D",]
df3 <- df1[df1$Time2<=48&df1$Light=="L",]
df4 <- df1[df1$Time2<=24&df1$Light=="D",]
df5 <- df1[df1$Time2>=24&df1$Light=="D",]
df6 <- df1[df1$Time2<=24&df1$Light=="L",]
df7 <- df1[df1$Time2>=24&df1$Light=="L",]

plot2 <- q1+geom_smooth(data = df2,aes(x=Time2,y=Oxygen),color="blue",fill="dodgerblue4",method="glm")+geom_smooth(data = df3,aes(x=Time2,y=Oxygen),color="cyan",fill="aquamarine",method="glm")+geom_smooth(data = df4,aes(x=Time2,y=Oxygen),se=FALSE,color="blue",linetype=c("dashed"),method="glm")+geom_smooth(data = df5,aes(x=Time2,y=Oxygen),se=FALSE,color="blue",linetype=c("dashed"),method="glm")+geom_smooth(data = df6,aes(x=Time2,y=Oxygen),se=FALSE,color="cyan",linetype=c("dashed"),method="glm")+geom_smooth(data = df7,aes(x=Time2,y=Oxygen),se=FALSE,color="cyan",linetype=c("dashed"),method="glm")

return(plot2)
}

plotregtrozosCO2 <- function(AMZRtotCO2,minimum,maximum,z='bold("B")')
{
df1 <- AMZRtotCO2
colorPalette <- c("red","pink","navy","steelblue1")

q1 <- ggplot(df1,aes(x=Time2,y=Oxygen,color=df1$TreatmentLight))+geom_point(size=3)+scale_color_manual(name="",values=c(colorPalette),breaks=unique(levels(TreatmentLight)),labels=c("CO\u2082_D","CO\u2082_L","Control_D","Control_L"))+
                    theme(
            panel.background=element_rect(fill="white", colour="black", size=0.75),
            panel.grid.major=element_line(size=0,"black",linetype=2),
            axis.text.x = element_text(size=rel(2)),
            axis.text.y = element_text(size=rel(2)),
            axis.title=element_text(size=rel(0), face="bold"),
            #legend.title = element_text("Treatment",size=rel(1.5)),
            legend.position="none",plot.margin=margin(1,1,1,1,unit="cm"))+
                                    annotate("text", label=z,size=4, x=5, y=minimum+0.5,parse=TRUE)+ ylim(minimum,maximum)

df2 <- df1[df1$Time2<=48&df1$Light=="D",]
df3 <- df1[df1$Time2<=48&df1$Light=="L",]
df4 <- df1[df1$Time2<=24&df1$Light=="D",]
df5 <- df1[df1$Time2>=24&df1$Light=="D",]
df6 <- df1[df1$Time2<=24&df1$Light=="L",]
df7 <- df1[df1$Time2>=24&df1$Light=="L",]


plot2 <- q1+geom_smooth(data = df2,aes(x=Time2,y=Oxygen),color="red",fill="tomato",method="glm")+geom_smooth(data = df3,aes(x=Time2,y=Oxygen),color="deeppink2",fill="pink",method="glm")+geom_smooth(data = df4,aes(x=Time2,y=Oxygen),se=FALSE,color="red",linetype=c("dashed"),method="glm")+geom_smooth(data = df5,aes(x=Time2,y=Oxygen),se=FALSE,color="red",linetype=c("dashed"),method="glm")+geom_smooth(data = df6,aes(x=Time2,y=Oxygen),se=FALSE,color="deeppink2",linetype=c("dashed"),method="glm")+geom_smooth(data = df7,aes(x=Time2,y=Oxygen),se=FALSE,color="deeppink2",linetype=c("dashed"),method="glm")

return(plot2)
}

#if(!require(devtools)) install.packages("devtools")
#devtools::install_github("kassambara/ggpubr")
library(ggpubr)

#Plot AMZR Control
OxygenAMZR <-Metabolism[Metabolism$Stations=="AMZR",]
AMZRtot <- OxygenAMZR[OxygenAMZR$Time2<=48,]# Este paso se puede evitar pero se ha definido así la función
maximum <- max(AMZRtot$Oxygen)
minimum <- min(AMZRtot$Oxygen)
AMZRtotK <-AMZRtot[AMZRtot$Treatment=="K",] 
pl1 <- plotregtrozosK(AMZRtotK,minimum,maximum,z='')

#Plot AMZR CO
AMZRtotCO2 <-AMZRtot[AMZRtot$Treatment=="CO2",] 
pl6 <- plotregtrozosCO2(AMZRtotCO2,minimum,maximum,z='')

#Plot SAM Control
OxygenSAM <-Metabolism[Metabolism$Stations=="SAM",]
SAMtot <- OxygenSAM[OxygenSAM$Time2<=48,]# Este paso se puede evitar pero se ha definido así la función
maximum <- max(SAMtot$Oxygen)
minimum <- min(SAMtot$Oxygen)
SAMtotK <-SAMtot[SAMtot$Treatment=="K",] 
pl2 <- plotregtrozosK(SAMtotK,minimum,maximum,z='')

#Plot SAM CO2
SAMtotCO2 <-SAMtot[SAMtot$Treatment=="CO2",] 
pl7 <- plotregtrozosCO2(SAMtotCO2,minimum,maximum,z='')

#Plot CAT Control
OxygenCAT <-Metabolism[Metabolism$Stations=="CAT",]
CATtot <- OxygenCAT[OxygenCAT$Time2<=48,]# Este paso se puede evitar pero se ha definido así la función
maximum <- max(CATtot$Oxygen)
minimum <- min(CATtot$Oxygen)
CATtotK <-CATtot[CATtot$Treatment=="K",] 
pl3 <- plotregtrozosK(CATtotK,minimum,maximum,z='')

#Plot CAT CO2
CATtotCO2 <-CATtot[CATtot$Treatment=="CO2",] 
pl8 <- plotregtrozosCO2(CATtotCO2,minimum,maximum,z='')

#Plot MAR Control
OxygenMAR <-Metabolism[Metabolism$Stations=="MAR",]
MARtot <- OxygenMAR[OxygenMAR$Time2<=48,]# Este paso se puede evitar pero se ha definido así la función
maximum <- max(MARtot$Oxygen)
minimum <- min(MARtot$Oxygen)
MARtotK <-MARtot[MARtot$Treatment=="K",] 
pl4 <- plotregtrozosK(MARtotK,minimum,maximum,z='')

#Plot MAR CO2
MARtotCO2 <-MARtot[MARtot$Treatment=="CO2",] 
pl9 <- plotregtrozosCO2(MARtotCO2,minimum,maximum,z='')

#Plot GD Control
OxygenGD <-Metabolism[Metabolism$Stations=="GD",]
GDtot <- OxygenGD[OxygenGD$Time2<=48,]# Este paso se puede evitar pero se ha definido así la función
GDtotK <-GDtot[GDtot$Treatment=="K",] 
maximum <- max(GDtot$Oxygen)
minimum <- min(GDtot$Oxygen)
pl5 <- plotregtrozosK(GDtotK,minimum,maximum,z='')

#Plot GD CO2
GDtotCO2 <-GDtot[GDtot$Treatment=="CO2",] 
pl10 <- plotregtrozosCO2(GDtotCO2,minimum,maximum,z='')#bold("GD"), para poner dentro del plot si se pone entre''

figure1 <- ggarrange(pl1,pl2,pl3,pl4,pl5, 
          labels = c("AMZR","SAM","CAT","MAR","GD"),
          ncol = 5, nrow = 1,common.legend=TRUE)
figure2 <- ggarrange(pl6,pl7,pl8,pl9,pl10, 
          labels = c("AMZR","SAM","CAT","MAR","GD"),
          ncol = 5, nrow = 1, common.legend=TRUE)
figure <- ggarrange(figure1,figure2,
          ncol = 1, nrow = 2)

#install.packages("grid")
library(grid)
###0
## Figure 4
png(file="/home/evaristo/Documentos/MOC2/REV/ScientiaMarina/Figure4.png",width=1500,height=1000)
annotate_figure(figure, left = textGrob(expression(bold("Oxygen ("*mu*mol*" "*O[2]*" "*L^-1*")")), rot = 90, vjust = 1, gp = gpar(cex = 1.3)),bottom = textGrob(expression(bold("Time (h)")), gp = gpar(cex = 1.3)))
dev.off()


####Modelos para hacer las pendientes de los timecourses, primero
Metabolism <- data.frame(Stations,Time2,Treatment,Light,Oxygen)# L = light, D = Dark
Metabolism <- na.omit(Metabolism)
tapply(Metabolism$Oxygen,list(Metabolism$Stations,Metabolism$Treatment,Metabolism$Light),mean)

attach(Metabolism)
y <-Metabolism$Oxygen
Time2 <-Metabolism$Time2
x <-as.factor(Metabolism$Stations)
w <-as.factor(Metabolism$Light)
z <-as.factor(Metabolism$Treatment)
f3 <- formula(y~Time2*x*w*z)

Mod0 <- gls(f3,weights=varIdent(form=~1|Time2*x*w*z),method="REML",data=Metabolism,correlation = corCompSymm(form =~Time2))
plot(Mod0)
anova(Mod0)
summary(Mod0)
Mod1 <- update(Mod0,.~.-w:z)##las concentraciones iniciales de oxígeno no dependian de si el tratamiento control o con adición de CO2 se ponia a la luz o a la oscuridad, dado que esa interacción no cosidera el tiempo y es, por tanto, el intercepto del modelo.
anova(Mod1)
summary(Mod1)
anova(Mod0,Mod1)

f3 <- formula(y~Time2*x*w*z)
Mod2 <- gls(f3,weights=varIdent(form=~1|Time2*x*w*z),method="REML",data=Metabolism,correlation = corCompSymm(form =~Time2|x))
anova(Mod2)
summary(Mod2)
Mod3 <- update(Mod2,.~.-w:z)
anova(Mod3)
anova(Mod0,Mod1,Mod2,Mod3)
Mod4 <- update(Mod3,.~.-Time2:x:z)
anova(Mod4)
Mod5 <- update(Mod4,.~.-Time2:w:z)
Mod6 <- update(Mod5,.~.-Time2:w:z)

library("mgcv")
library("itsadug")

model <- lm(y~x*w*z*Time2)
plot(model)
summary(model)

basic_model <- gam(y ~x+w+z+s(Time2,k=3), data = Metabolism, method = "REML")
summary(basic_model)
k.check(basic_model)
gam.check(basic_model)

par(mfrow = c(2, 2))
plot(basic_model, all.terms = TRUE)

basic_model2 <- gam(y ~x*w*z+s(Time2,k=3), data = Metabolism, method = "REML")
summary(basic_model2)
par(mfrow = c(2, 2))
plot(basic_model2, all.terms = TRUE)
plot_smooth(basic_model2, view = "Time2", rm.ranef = FALSE)

basic_model3 <- gam(y ~x*w*z+s(Time2,k=3,bs="re"), data = Metabolism, method = "REML")
summary(basic_model3)
par(mfrow = c(2, 2))
plot(basic_model3, all.terms = TRUE)









                                        #-------------pHs
pH_in_Situ <- c(8.04,8.05,8.03,8.13,7.97)
pH_Control <- c(8.04,7.97,8.065,8.05,7.95)
pH_Treatment <- c(7.7,7.82,7.75,7.5,7.65)
Tot_pH <- data.frame(cbind(c(rep(c(1,2,3),each=5)),rep(pH_in_Situ,3),c(pH_in_Situ,pH_Control,pH_Treatment)))
colnames(Tot_pH) <- c("Treatment","pHInsitu","pH")
Tot_pH$Treatment <- as.factor(Tot_pH$Treatment)

library(ggplot2)
# Basic box plot


colorPalette <- c("blue","skyblue","red")
p <- ggplot(Tot_pH, aes(x=Treatment, y=pH,col=Treatment))+geom_boxplot(notch=FALSE)+geom_jitter(shape=19, position=position_jitter(0),size=2.5)

plot1 <- p+labs(x="Treatment", y = "pH")+scale_x_discrete(labels=expression(italic("In situ"),Control,CO[2]))+theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+theme(plot.margin = margin(1.5,1.5,1.5,1.5, "cm"))+theme(axis.title.y = element_text(margin = margin(r = 0.5, unit="in")))+theme(axis.title.x = element_text(margin = margin(t = 0.5, unit="in")))+theme(panel.background = element_rect(fill = 'white', color = 'black'),panel.grid.major = element_line(color = 'black', linetype = 'dotted'))+scale_color_manual(name="",values=c(colorPalette),breaks=unique(levels(Tot_pH$Treatment)),labels=expression(italic("In situ"),"Control","CO\u2082"))+theme(legend.position="none")


shapiro.test(Tot_pH$pH)
model <- lm(Tot_pH$pH~as.factor(Tot_pH$Treat))
summary(model)
anova(model)

kruskal.test(Tot_pH$pH~as.factor(Tot_pH$Treat))
pairwise.wilcox.test(Tot_pH$pH,as.factor(Tot_pH$Treat),p.adjust.method = "BH")
## hay diferencias significativas entre los pH

df <- Tot_pH[,]
df1 <- Tot_pH[Tot_pH$Treatment==1,]
df2 <- Tot_pH[Tot_pH$Treatment==2,]
df3 <- Tot_pH[Tot_pH$Treatment==3,]

q <- ggplot(df,aes(df$pHInsitu,df$pH,col=df$Treatment))+geom_point(shape=19,size=3.5)+scale_color_manual(name="",values=c(colorPalette),breaks=unique(levels(df$Treatment)),labels=expression(italic("In situ"),"Control","CO\u2082"))+theme(axis.title.x = element_text(size = 16),axis.text.x = element_text(size = 14),axis.title.y = element_text(size = 16),axis.text.y = element_text(size = 14))+theme(plot.margin = margin(1.5,1.5,1.5,1.5, "cm"))+theme(axis.title.y = element_text(margin = margin(r = 0.5, unit="in")))+theme(axis.title.x = element_text(margin = margin(t = 0.5, unit="in")))+theme(panel.background = element_rect(fill = 'white', color = 'black'),panel.grid.major = element_line(color = 'black', linetype = 'dotted'))+theme(legend.text=element_text(size=14))+theme(panel.background = element_rect(fill = 'white', color = 'black'),panel.grid.major = element_line(color = 'black', linetype = 'dotted'))+theme(legend.position="none")

insitu <- expression(paste(italic("In situ")," pH"))
plot2 <- q+labs(x=insitu, y = "Treatment")+geom_point(shape=19,size=3.5)+geom_smooth(data=df1,aes(pHInsitu,pH),color="blue",method="glm")+geom_point(data=df2,aes(pHInsitu,pH),color="skyblue")+geom_point(shape=19,size=3.5)


library(ggpubr)
library(cowplot)
plotA <- plot1+geom_text(x=0.6, y=8.1, label="A", size=11.5, col="black")
plotB <- plot2+geom_text(x=7.98, y=8.1, label="B", size=11.5, col="black")

Figure2 <- ggarrange(plotA, plotB, 
          ncol = 2, nrow = 1)

ggexport(Figure2, filename = "figureNovapH.pdf",height=6,width=12)

model <- lm(df$pH~df$Treatment)
summary(model)

model <- lm(df1$pH~df1$pHInsitu)
model2 <- lm(df2$pH~df2$pHInsitu)
model3 <- lm(df2$pH~df2$pHInsitu+offset(df2$pHInsitu))
summary(model2)

mean(pH_in_Situ)
mean(pH_Control)
mean(pH_Treatment)
Diferencia_pH<- c(pH_in_Situ-pH_Control)
Diferencia_promedio<- mean(pH_in_Situ-pH_Control)## El botellón acidificaba ligeramente las muestras, o se tomaban en profundidades ligeramente más ácidas
mean(pH_Control-pH_Treatment)### Las muestras se acidificaron 10 veces más que la diferencia entre la muestra in situ y el botellón
max(pH_in_Situ-pH_Control)## diferencia máxima entre el pH del botellón o el segundo CTD y el CTD in situ
min(pH_in_Situ-pH_Control)## diferencia mínima entre el pH del botellón o el segundo CTD y el CTD in situ

max(pH_Control-pH_Treatment)## diferencia máxima entre el pH del botellón o el segundo CTD y el CTD in situ y la muestra acidificada
min(pH_Control-pH_Treatment)## diferencia máxima entre el pH del botellón o el segundo CTD y el CTD in situ y la muestra acidificada

plot(pH_Treatment-pH_Control~pH_in_Situ)## no hay una tendencia, la más básica no era la que más acidificamos)
model <- lm(pH_Treatment-pH_Control~pH_in_Situ)
abline(model)
summary(model)

plot(df2$pH~df2$pHInsitu,xlim=c(7.9,8.15),ylim=c(7.9,8.15))
model2 <- lm(df2$pH~df2$pHInsitu)
abline(model2)
abline(0,1, col="red")
summary(model3)#no hay diferencias entre e control y el pH i situ
