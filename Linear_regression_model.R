rm(list=objects())
graphics.off()
getwd()
setwd("C:/Users/ElAamrani/Desktop/Centrale_2�me_ann�e/S8/Statistiques_avanc�es/TP")


#----------------------------------------Lecture du jeu de donn�es-----------------------------------------
#Affichage des donn�es
df=read.table("pine.sup.data", header=TRUE)
dim(df)
head(df, 5)
str(df)
#Etude des grandeurs (min, max, moyenne, m�diane, quantile, �cart type) de chaque colonne
summary(df)
apply(df,2,mean)
apply(df,2,sd)
apply(df,2,function(x){c(mu=mean(x),sd=sd(x))})

#----------------------------------Etude des corr�lations entre variables----------------------------------
#Scatter plot des corr�lations
pairs(df)

#Scatter plot des corr�lations et distribution des covariables
library(GGally)
ggpairs(df)

#Calcul du coefficient de corr�lation
round(cor(df),3)

#Visulation graphique des coefficients de corr�lations
library(corrplot)
corrplot(cor(df))


#------------------------------R�gression lin�aire avec toutes les variables-------------------------------
#Extraction du plan d'exp�rience X et des valeurs de sortie Y
X=as.matrix(df)
X<-X[,c(11,1:10)]
colnames(X)[1]='x0'
colnames(X)
X[,1]=1 # intercept
Y=as.matrix(df[,11])
n=length(Y)
p=ncol(X)
X
Y

#Calcul de la matrice de projection
H=diag(X%*%solve(t(X)%*%X)%*% t(X))
#Calcul de l'estimateur de theta
theta.est=solve(t(X)%*%X)%*% t(X)%*%Y
round(theta.est,3)
#Calcul de l'estimateur non-biais� de sigma^2
sigma.est=sqrt(sum( (X%*%theta.est -Y)^2   )/(n-p))
round(sigma.est,3)
#Calcul de la matrice de covariance de la loi de theta
V=solve(t(X)%*%X) *sigma.est^2
#Calcul de l'�cart-type de chaque composante de theta
stddev=sqrt(diag(V))

#Calcul de la statistique de test de Student pour theta_j=0 contre theta_j#0
tobs=theta.est/stddev

#Calcul des valeurs estim�es dans le mod�le
Y.estim=X%*%theta.est

#V�rification de la coh�rence des r�sultats avec le r�sultat de la fonction 'lm' de R
res=lm(x11~.,data=df)
summary(res)
names(res)
res$residuals
Y-Y.estim
#Valeur estim�e de theta
theta.lm= res$coef
theta.lm-theta.est 
#Calcul de la p-value du test de Student pour l'intercept
round(2*pt(-abs(tobs),res$df),4)

#Calcul de la statistique de Fisher (F-Statistic)
A=diag(c(0,1,1,1,1,1,1,1,1,1,1))[2:11,]
fstat=1/((p-1)*sigma.est^2)*t(A%*%theta.est)%*%solve(A%*%solve(t(X)%*%X)%*%t(A))%*%A%*%theta.est
fstat
round(1-pf(fstat,p-1,n-p),4)

#------------------------------------Backward selection - Crit�re AIC--------------------------------------
library(MASS)
resf=stepAIC(res)
#V�rification : resf est de m�me classe que res
class(res)
class(resf)
n*log(deviance(res)/n)+2*11
#Explication des �tapes de resf
resf1=lm(x11~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=df)
summary(resf1)
deviance(resf1)
sum((Y-X%*%c(resf1$coefficients,0))^2)
n*log(deviance(resf1)/n)+2*10
#Explication (2) des �tapes de resf
resf2=lm(x11~x1+x2+x4+x5+x6+x7+x8+x9+x10,data=df)
n*log(deviance(resf2)/n)+2*10
deviance(res)
#Propri�t�s de resf
resf$coef
resf$anova
summary(resf)

#------------------------------------Backward selection - Crit�re BIC--------------------------------------
resg=stepAIC(res, k = log(n))
#V�rification : resg est de m�me classe que res
class(res)
class(resg)
#Propri�t�s de resg
resg$coef
resg$anova
#Explication des �tapes de resg
n*log(deviance(resf1)/n)+10*log(n)
-2*logLik(resf2)+10*log(n)


#------------------------------------Stepwise selction - Crit�re AIC---------------------------------------
resh=stepAIC(lm(x11~x3+x6, data=df), direction="both", scope=~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10)
#V�rification : resh est de m�me classe que res
class(res)
class(resh)
#Propri�t�s de resh
resh$coef
resh$anova


#------------------------------------------Recherche exhaustive-------------------------------------------
library(leaps)
recherche=regsubsets(x11~., int=TRUE, nbest=1, nvmax=10, method="exhaustive", data=df)
#Interpr�tation
plot(recherche,scale="bic", main="Recherche exhaustive par crit�re BIC")
resi=stepAIC(lm(x11~., data=df), direction="both", k = log(n))
#Interpr�tation
summary(recherche)
summary(recherche)$bic
plot(summary(recherche)$bic, main="Evolution du crit�re BIC")
#Recherche exhaustive pour BIC, Cp, adjr2, r2
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(recherche, scale="bic", main="Crit�re BIC")
plot(recherche, scale="Cp", main="Crit�re Cp")
plot(recherche, scale="adjr2", main="Crit�re adjr2")
plot(recherche, scale="r2", main="Crit�re r2")
title("Selection de mod�les selon diff�rents crit�res", outer=TRUE, font.main=4, font.lab=4, font.sub=4,
      cex.main=1.7, cex.lab=1.7, cex.sub=1.2)


#---------------------------------------------Cross-Validation--------------------------------------------
library(forecast)
resf1=lm(x11~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=df)
resf2=lm(x11~x1+x2+x4+x5+x6+x7+x8+x9+x10,data=df)
CV(resf1)
CV(resf2)
#Calcul des crit�res AIC et BIC
n*log(deviance(resf1)/n)+2*10
n*log(deviance(resf1)/n)+log(n)*10
n*log(deviance(resf2)/n)+2*10
n*log(deviance(resf2)/n)+log(n)*10
#Comparaison avec le crit�re AIC de CV
CV(resf1)[2]-(n*log(deviance(resf1)/n)+2*10)
CV(resf2)[2]-(n*log(deviance(resf2)/n)+2*10)
#Comparaison avec le crit�re BIC de CV
CV(resf1)[4]-(n*log(deviance(resf1)/n)+log(n)*10)
CV(resf2)[4]-(n*log(deviance(resf2)/n)+log(n)*10)
#Calcul et comparaison de PRESS dans le cas Leave One Out
resc1=lm(x11~.,data=df)
resc2=lm(x11~x2+x5+x6+x9,data=df)
CV(resc1)
(1/n)*sum(((Y-resc1$fitted)/(1-lm.influence(resc1)$hat))^2)
CV(resc2)
1/n*sum(((Y-resc2$fitted)/(1-lm.influence(resc2)$hat))^2)

#Calcul et comparaison de l'EQMP dans une cross-validation avec K=10
library(caret)
set.seed(42)

res10=train(x11~., data=df, method="lm",trControl=trainControl(method="cv",number=10, verboseIter=TRUE))
res10$results
res25=train(x11~., data=df, method="lm",trControl=trainControl(method="loocv", verboseIter=TRUE))
res25$resample
1/25*sum((res25$resample["RMSE"])^2)
CV(resc1)["CV"]

#Impl�mentation de la K-fold cross validation pour le mod�le complet
#Prendre K=n si on veut retrouver le PRESS
K=10
j=1
i=0
len = n%/%K
EQMP=0
for (i in 1:(K-1)){
  df_i=df[-(j:(j-1+len)),]
  model=lm(x11~., data=df_i)
  EQMP=EQMP+1/len*(sum((Y[j:(j-1+len)]-X[j:(j-1+len),]%*%model$coef)^2))
  print(1/(len)*(sum((Y[j:(j-1+len)]-X[j:(j-1+len),]%*%model$coef)^2)))
  j<-j+len
}
df_last=df[-(j:n),]
model=lm(x11~., data=df_last)
EQMP=EQMP+1/(n-(K-1)*len)*(sum((Y[j:n]-X[j:n,]%*%model$coef)^2))
EQMP=EQMP/K
EQMP

#Impl�mentation de la K-fold cross validation pour le mod�le retenu
K=10
j=1
i=0
len = n%/%K
EQMP=0
for (i in 1:(K-1)){
  df_i=df[-(j:(j-1+len)),]
  model=lm(x11~x2+x5+x6+x9, data=df_i)
  EQMP=EQMP+1/(len)*(sum((Y[j:(j-1+len)]-X[j:(j-1+len),c(1,3,6,7,10)]%*%model$coef)^2))
  print(1/(len)*(sum((Y[j:(j-1+len)]-X[j:(j-1+len),c(1,3,6,7,10)]%*%model$coef)^2)))
  j<-j+len
}
df_last=df[-(j:n),]
model=lm(x11~x2+x5+x6+x9, data=df_last)
EQMP=EQMP+1/(n-(K-1)*len)*(sum((Y[j:n]-X[j:n,c(1,3,6,7,10)]%*%model$coef)^2))
EQMP=EQMP/K
EQMP

