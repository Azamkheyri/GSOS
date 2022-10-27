rm(list=ls())
library("mvtnorm")
library("glasso")
library("MASS")
##########################################
fun1 = function(L,CV){
Frobenius.dis = sqrt(sum((CV$Theta-Theta)^2))
dd=as.matrix(CV$Theta-Theta)
eigenvalues = eigen(t(dd)%*%dd)$values 
SP = sqrt(max(Re(eigenvalues[abs(Im(eigenvalues)) < 1e-6])))
TM = as.matrix(CV$Theta) %*% Sigma
QL = sum(diag((TM-diag(1,p))))^2
KL =  sum(diag(TM)) - determinant(TM, logarithm = TRUE)$modulus[1] -p
c(Frobenius.dis, SP, KL, QL)
}
#######################################Network#1
NetS20 <- c("CSSigma20.txt")
NetT20 <- c("CSTheta20.txt")
NetS50 <- c("CSSigma50.txt")
NetT50 <- c("CSTheta50.txt")
NetS100 <- c("CSSigma100.txt")
NetT100 <- c("CSTheta100.txt")
#####################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod1-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-d100.RData")
#############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod1-new-vI-d100.RData")
#######################################Network#2
NetS20 <- c("VerySparseSigma20.txt")
NetT20 <- c("VerySparseTheta20.txt")
NetS50 <- c("VerySparseSigma50.txt")
NetT50 <- c("VerySparseTheta50.txt")
NetS100 <- c("VerySparseSigma100.txt")
NetT100 <- c("VerySparseTheta100.txt")
############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,9))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL","f1score","MCC")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod2-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-d100.RData")
#############################################

cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}

#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod2-new-vI-d100.RData")
#######################################Network#3
NetS20 <- c("SampleSigma20.txt")
NetT20 <- c("SampleTheta20.txt")
NetS50 <- c("SampleSigma50.txt")
NetT50 <- c("SampleTheta50.txt")
NetS100 <- c("SampleSigma100.txt")
NetT100 <- c("SampleTheta100.txt")
#########################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod3-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-d100.RData")
#############################################

cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod3-new-vI-d100.RData")
#######################################Network#4
NetS20 <- c("StarSigma20.txt")
NetT20 <- c("StarTheta20.txt")
NetS50 <- c("StarSigma50.txt")
NetT50 <- c("StarTheta50.txt")
NetS100 <- c("StarSigma100.txt")
NetT100 <- c("StarTheta100.txt")
##########################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod4-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-d100.RData")
#############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod4-new-vI-d100.RData")
#######################################Network#5
NetS20 <- c("OtherMASigma20.txt")
NetT20 <- c("OtherMATheta20.txt")
NetS50 <- c("OtherMASigma50.txt")
NetT50 <- c("OtherMATheta50.txt")
NetS100 <- c("OtherMASigma100.txt")
NetT100 <- c("OtherMATheta100.txt")
########################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod5-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-d100.RData")
#############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod5-new-vI-d100.RData")
#######################################Network#6
NetS20 <- c("DiagonalDominantSigma20.txt")
NetT20 <- c("DiagonalDominantTheta20.txt")
NetS50 <- c("DiagonalDominantSigma50.txt")
NetT50 <- c("DiagonalDominantTheta50.txt")
NetS100 <- c("DiagonalDominantSigma100.txt")
NetT100 <- c("DiagonalDominantTheta100.txt")
#######################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
  alpha <- sort(alpha)

  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=NULL)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=NULL)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))
}
###########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)

save(A, file = "mod6-new-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam", "alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-d100.RData")
#############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  T <- diag(1,p)
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
   k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-I-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-I-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-I-d100.RData")
##############################################
cv <-function(nfold=5L,Y,lambda,alpha){
  
  if(is.null(Y)){ print("Y is required as input"); return(); }
  lambda <- sort(lambda)
   alpha <- sort(alpha)
  l = function(Theta,Sigma){
    res <- determinant(Theta, logarithm=TRUE)$modulus - sum(diag(Sigma%*%Theta))
    return(res)
  }
  n = nrow(Y)
  p = ncol(Y)
  
 ind=c(1:n)

  indnr <- rep(0L, nfold+1)
  
  indnr[1]=0
 rest = (n - (n %% nfold))/nfold
  for(i in 1:(nfold-1)){
    indnr[i+1]=i*rest
  }
  indnr[nfold+1]=n
  CV <- matrix(0,length(lambda),length(alpha))
  for(j in 1:length(lambda)){
  for(r in 1:length(alpha)){
        for(i in 1:nfold){
      indtest <- ind[(indnr[i]+1):indnr[i+1]]
      Strain <- crossprod(Y[-indtest,])/nrow(Y[-indtest,])
      T = 1/mean(diag(Strain))*diag(1,p)
      Thetatrain <- gsos(Strain,lambda=lambda[j],alpha=alpha[r],Target=T)$Theta
      CV[j,r]=CV[j,r]+l(Thetatrain, cov(Y[indtest,]))
    }
      }}
  k=which(CV==max(CV), arr.ind = TRUE)
  if( (k[1]==1) || (k[1]==length(lambda)) ){ warning("The optimal lambda value is a boundary value. The range of possible values should be increased") }
  if( (k[2]==1) || (k[2]==length(alpha)) ){ warning("The optimal alpha value is a boundary value. The range of possible values should be increased") }
  optimalla=lambda[k[1]]
  optimalal=alpha[k[2]]

    S=crossprod(Y)/n
    T = 1/mean(diag(S))*diag(1,p)
    object <- gsos(S=S,lambda=optimalla,alpha=optimalal,Target=T)$Theta
    return(list(optimal=c(optimalla,optimalal),lambda=lambda,CV=CV/nfold,Theta=object))

}
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-vI-d20.RData")
#########################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-vI-d50.RData")
#########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = cv( Y = X , lambda = lambda_grid, alpha=alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))

}
colMeans(A)
save(A, file = "mod6-new-vI-d100.RData")

