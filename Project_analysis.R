setwd("~/Documents/Measurements/")
#dat_full=read.csv("2017_Sample_log_ESS_212.csv", head=T)
dat_full=read.csv("data_cleaned.csv", head=T)
View(dat_full)

site=c("LT","SD","CF","LB","GC","SW","GC","GC","PL","SW","LT","GC","PL","SW","GC","SW","GC","GC","LB","SD","CF","GC","GC","PL","GC","SW")
date=c(1,1,1,1,2,2,3,4,4,4,5,5,5,5,6,6,7,7,8,8,8,8,9,10,10,10)
#Soil samples have been removed, data from 12/17/15 and 1/4/16 removed b/c didn't do elemental analysis

dat=dat_full[,7:37]
dat[is.na(dat)]=0 #convert NA's to 0's (last year's data converted negative measures to NA instead of 0)
dat[dat<0]=0 #convert negative entries to 0
#dat=data.frame(apply(dat, 2, as.numeric))


#dat=dat[,-28] #has NA's
#dat=dat[,-24] #remove Ag
#dat=dat[-## remove soil entries

### Average the triplicates ###
  #the first 42 rows are actually sextuplicates:
  b=rep(c(1:(42/6)), each=6)
  dat.p1=aggregate.data.frame(dat[1:42,], by=list(b), FUN='mean')[,-1]
  #last 57 rows are triplicates
  b=rep(c(1:(57/3)), each=3) #99-43+1=57
  dat.p2=aggregate.data.frame(dat[43:99,], by=list(b), FUN='mean')[,-1]
  dat=rbind(dat.p1, dat.p2) #26 x 31
  
### Remove elements with lots of NA's (or do something clever) ###
  remove_these=which(colnames(dat)=="delta15N" | colnames(dat)=="delta18O"
                     | colnames(dat)=="TN_umolperL" | colnames(dat)=="TP_umolperL" | colnames(dat)=="NO2_umolperL"
                     | colnames(dat)=="Ag" | colnames(dat)=="Cd" | colnames(dat)=="Cs" | colnames(dat)=="Tl")
  dat=dat[, -remove_these]


### PCA Regression ###
  library(pls)
  set.seed(30)
  pcr.fit=pcr(PO4~., data=dat, scale=T, validation="CV")
  summary(pcr.fit) #5 PC's seems to be the bet
  validationplot(pcr.fit, val.type="RMSEP")
  pca.pred=predict(pcr.fit, dat, ncomp=5)
  plot(dat$PO4, pca.pred, pch=16, col='steelblue', ylab="Predicted", xlab="Actual", main="PCA \n PO4 Predicted vs. Actual")
  r2=round(cor(pca.pred, dat$PO4)^2, 2)
  text(7,0.5, bquote(R^2 == .(r2)))
  abline(0,1, lty=2)

### Partial Least Squares ###
  set.seed(50)
  pls.fit=plsr(PO4~., data=dat, scale=T) #optoinal arguments to put in training and test set
  summary(pls.fit) #looks like 3 or 4 components best
  pls.pred=predict(pls.fit, dat, ncomp=4)
  plot(dat$PO4, pls.pred, pch=16, col='steelblue', ylab="Predicted", xlab="Actual", main="PLS \n PO4 Predicted vs. Actual")
  r2=round(cor(pls.pred, dat$PO4)^2, 2)
  text(5,1, bquote(R^2 == .(r2)))
  abline(0,1, lty=2)

##need to do a training/test subset to compare PCA regression and PLS

### LASSO ###
  library(glmnet)
  grid=10^seq(10, -2, length=100)
  lasso.mod=glmnet(x=as.matrix(dat[,c(1,2,4:22)]), y=dat[,3], alpha=1, lambda=grid)
  set.seed(9)
  cv.out=cv.glmnet(x=as.matrix(dat[,c(1,2,4:22)]), y=dat[,3], alpha=1)
  plot(cv.out)
  bestlam=cv.out$lambda.min
  lasso.pred=predict(lasso.mod, s=bestlam, newx=as.matrix(dat[,c(1,2,4:22)]))  
  plot(dat$PO4,lasso.pred, pch=16, col='steelblue', ylab="Predicted", xlab="Actual", main="LASSO \n PO4 Predicted vs. Actual")
  r2=round(cor(lasso.pred, dat$PO4)^2, 2)
  text(3.5,1, bquote(R^2 == .(r2)))
  abline(0,1, lty=2)
  
  lasso.coef=lasso.pred=predict(lasso.mod, s=bestlam, type="coefficients")
  #LASSO selects Cr, Mn, Ga, Se, Rb as significant predictors
  
  ## with outlier removed
  lasso.mod=glmnet(x=as.matrix(dat[-7,c(1,2,4:22)]), y=dat[-7,3], alpha=1, lambda=grid)
  set.seed(9)
  cv.out=cv.glmnet(x=as.matrix(dat[-7,c(1,2,4:22)]), y=dat[-7,3], alpha=1)
  plot(cv.out)
  bestlam=cv.out$lambda.min
  lasso.pred=predict(lasso.mod, s=bestlam, newx=as.matrix(dat[-7,c(1,2,4:22)]))  
  plot(dat$PO4[-7], lasso.pred, pch=16, col='steelblue', ylab="Predicted", xlab="Actual", main="LASSO \n PO4 Predicted vs. Actual")
  r2=round(cor(lasso.pred, dat$PO4[-7])^2, 2)
  text(3.5,1, bquote(R^2 == .(r2)))
  abline(0,1, lty=2)
  
  lasso.coef=lasso.pred=predict(lasso.mod, s=bestlam, type="coefficients")
  #LASSO selects NH4, Al, V, Cr, Mn, Ni, Ga, As, Se, Sr, Ba, Pb, and U
  #Cr, Ga, Rb if select a slightly larger lambda(lambda.1se from cv.out)
  
##### Unsupervised #####

  ### PCA ###
    pr.out=prcomp(dat, scale=T)
    biplot(pr.out, scale=0)
    biplot(pr.out, scale=0, xlabs=rep("", 26), xlim=c(-3,3))
    pr.var=pr.out$sdev^2
    pve=pr.var/sum(pr.var) #percentage variance explained by each component
    plot(pve, xlab="PC", ylab="Proportion of Variance Explained", ylim=c(0,1), pch=16, col='orange', bty='n') #5 principal components looks good
    plot(cumsum(pve), xlab="PC", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1), pch=16, col='orange', bty='n')
    
  
  ### K-means clustering ###
    #Cluster metals
    set.seed(63)
    km.out=kmeans(t(dat), centers=3, nstart=20) #centers = # of clusters
    km.out$cluster
    group.1=sort(colnames(dat)[km.out$cluster==1]); group.1
    group.2=sort(colnames(dat)[km.out$cluster==2]); group.2
    group.3=sort(colnames(dat)[km.out$cluster==3]); group.3
    #group.4=sort(colnames(dat)[km.out$cluster==4]); group.4
    #group.5=sort(colnames(dat)[km.out$cluster==5]); group.5

    #Cluster measurements
    set.seed(63)
    km.out=kmeans(dat, centers=3, nstart=20) #centers = # of clusters
    group.1=sort(site[km.out$cluster==1]); group.1
    group.2=sort(site[km.out$cluster==2]); group.2
    group.3=sort(site[km.out$cluster==3]); group.3
    group.4=sort(site[km.out$cluster==4]); group.4
    
    group.1=sort(date[km.out$cluster==1]); group.1 #no evidence they're clustering by date
    group.2=sort(date[km.out$cluster==2]); group.2
    group.3=sort(date[km.out$cluster==3]); group.3
    group.4=sort(date[km.out$cluster==4]); group.4


    
### K-Means But now with Soil! ###
soil=read.csv("data_cleaned_soil.csv", head=T)
soil.dat=soil[,7:37]
soil.dat[is.na(soil.dat)]=0 #convert NA's to 0's (last year's data converted negative measures to NA instead of 0)
soil.dat[soil.dat<0]=0 #convert negative entries to 0
remove_these=which(colnames(soil.dat)=="delta15N" | colnames(soil.dat)=="delta18O"
                     | colnames(soil.dat)=="TN_umolperL" | colnames(soil.dat)=="TP_umolperL" | colnames(soil.dat)=="NO2_umolperL"
                     | colnames(soil.dat)=="Ag" | colnames(soil.dat)=="Cd" | colnames(soil.dat)=="Cs" | colnames(soil.dat)=="Tl")
soil.dat=soil.dat[, -remove_these]
soil.dat=rbind(dat, soil.dat)
soil.ind=c(rep(0,26), rep(1, 6))
#Cluster measurements
    set.seed(50)
    km.out=kmeans(soil.dat, centers=2, nstart=20) #centers = # of clusters
    group.1=sort(soil.ind[km.out$cluster==1]); group.1 
    group.2=sort(soil.ind[km.out$cluster==2]); group.2 #soils 3 and 4
    
