rm(list=ls())
library("GLassoElnetFast")
library("mvtnorm")
library("glasso")
library("MASS")
############################################

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
#######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod1-rope-vI-d100.RData")
#######################################Network#2
NetS20 <- c("VerySparseSigma20.txt")
NetT20 <- c("VerySparseTheta20.txt")
NetS50 <- c("VerySparseSigma50.txt")
NetT50 <- c("VerySparseTheta50.txt")
NetS100 <- c("VerySparseSigma100.txt")
NetT100 <- c("VerySparseTheta100.txt")
############################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod2-rope-vI-d100.RData")
#######################################Network#3
NetS20 <- c("SampleSigma20.txt")
NetT20 <- c("SampleTheta20.txt")
NetS50 <- c("SampleSigma50.txt")
NetT50 <- c("SampleTheta50.txt")
NetS100 <- c("SampleSigma100.txt")
NetT100 <- c("SampleTheta100.txt")
#########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod3-rope-vI-d100.RData")
#######################################Network#4
NetS20 <- c("StarSigma20.txt")
NetT20 <- c("StarTheta20.txt")
NetS50 <- c("StarSigma50.txt")
NetT50 <- c("StarTheta50.txt")
NetS100 <- c("StarSigma100.txt")
NetT100 <- c("StarTheta100.txt")
##########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod4-rope-vI-d100.RData")
#######################################Network#5
NetS20 <- c("OtherMASigma20.txt")
NetT20 <- c("OtherMATheta20.txt")
NetS50 <- c("OtherMASigma50.txt")
NetT50 <- c("OtherMATheta50.txt")
NetS100 <- c("OtherMASigma100.txt")
NetT100 <- c("OtherMATheta100.txt")
########################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod5-rope-vI-d100.RData")
#######################################Network#6
NetS20 <- c("DiagonalDominantSigma20.txt")
NetT20 <- c("DiagonalDominantTheta20.txt")
NetS50 <- c("DiagonalDominantSigma50.txt")
NetT50 <- c("DiagonalDominantTheta50.txt")
NetS100 <- c("DiagonalDominantSigma100.txt")
NetT100 <- c("DiagonalDominantTheta100.txt")
#######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-d100.RData")
######################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-I-d100.RData")
##################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-vI-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-vI-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,6))
colnames(A)= c("lam", "loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0,rope = TRUE,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))],fun1(L,CV))
}

save(A, file = "mod6-rope-vI-d100.RData")

