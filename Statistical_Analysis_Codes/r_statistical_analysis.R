######### Code for the anova statistical analysis of the statistical_analysis.txt file with fixation times of each fitness site scenario ##############

library(reshape2)
library(ggplot2)
library(fitdistrplus)
library(moments)
library(car)

setwd("/home/maria/Diplomatiki_parousiasi/statistical_analysis_new")

data<-read.csv(file="statistical_analysis.txt", sep='\t', header=TRUE)
head(data)

par(mfrow=c(4,1))
par(mfrow=c(1,4))
titles<-c("Linear_model", "Beta_model", "Constant_model","Sin_model")
for (col in 1:ncol(data)){
  hist(data[,col], main=substitute(paste('Histogram of ', a), list(a=titles[col])), xlab="Generations")
}
qqnorm(data$linear_model, pch = 1, frame = FALSE, main=" QQplot Linear_model")
qqline(data$linear_model, col = "steelblue", lwd = 2)

qqnorm(data$beta_model, pch = 1, frame = FALSE, main="QQplot Beta_model")
qqline(data$beta_model, col = "steelblue", lwd = 2)

qqnorm(data$constant_model, pch = 1, frame = FALSE, main="QQplot Constant_model")
qqline(data$constant_model, col = "steelblue", lwd = 2)

qqnorm(data$sin_model, pch = 1, frame = FALSE, main="QQplot Sin_model")
qqline(data$sin_model, col = "steelblue", lwd = 2)

linear_density<-density(data$linear_model)
plot(linear_density, xlab="Generations")

beta_density<-density(data$beta_model)
plot(beta_density, xlab="Generations")

constant_density<-density(data$constant_model)
plot(constant_density, xlab="Generations")

sin_density<-density(data$sin_model)
plot(sin_density, xlab="Generations")


shapiro.test(data$linear_model)
shapiro.test(data$beta_model)
shapiro.test(data$constant_model)
shapiro.test(data$sin_model)

descdist(data$linear_model, discrete = TRUE)
descdist(data$beta_model, discrete = TRUE) ### negative binomial distribution
descdist(data$constant_model, discrete = TRUE) ### negative binomial distribution

skewness(data$linear_model, na.rm = TRUE) ### positive skewnass so log transformation
skewness(data$beta_model, na.rm = TRUE)  ### positive skewnass so log transformation
skewness(data$beta_model, na.rm = TRUE)  ### positive skewnass so log transformation

data_anova<-read.csv(file='statistical_analysis_anova.txt', sep='\t', header=TRUE)
head(data_anova)

par(mfrow=c(1,1))
par(cex.axis=0.8)
boxplot(Time~Model,
        data=data_anova,
        main="Different boxplots of fixation times for each fitness model",
        ylab="Generations",
        col="orange",
        border="brown",
        las=2
)

levels(data_anova$Model)

fligner.test(Time ~ Model, data = data_anova) # The Fligner-Killeen test is one of the many tests for homogeneity of variances which is most robust against departures from normality. Null hypothesis  the variances in each of the groups (samples) are the same. For pvalue 0.01 the null hypothesis can not be rejected.

res.aov <- aov(Time ~ Model, data = data_anova) ### The one-way analysis of variance (ANOVA), also known as one-factor ANOVA, is an extension of independent two-samples t-test for comparing means in a situation where there are more than two groups. In one-way ANOVA, the data is organized into several groups base on one single grouping variable (also called factor variable).
### ANOVA test hypotheses:

###Null hypothesis: the means of the different groups are the same
###Alternative hypothesis: At least one sample mean is not equal to the others.

summary(res.aov) ### As the p-value is less than the significance level 0.05, we can conclude that there are significant differences between the  means of the groups highlighted with “*" in the model summary.
#levene.test(data_anova$Time, data_anova$Model, location="median", correction.method="zero.correction")

TukeyHSD(res.aov) ### It can be seen from the output, that only the difference between linear_model and constant_model is significant with an adjusted p-value of 0.042.

##### I used this site for anova guidance http://www.sthda.com/english/wiki/one-way-anova-test-in-r
