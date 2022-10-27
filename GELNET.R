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
##############################################
mycrossvalidation <- function (nfold = 5L, Y, lambda, alpha, ind = NULL, type = NULL, 
    Target = NULL, outer.maxit = 1000, outer.thr = 1e-05, inner.maxit = 1000, 
    inner.thr = outer.thr/10, penalize.diagonal = TRUE, cor = FALSE, 
    rope = FALSE) 
{
    if (is.null(Y)) {
        print("Y is required as input")
        return()
    }
    if (is.vector(lambda) == FALSE) {
        print("lambda should be a vector")
        return()
    }
    else {
        lambda <- sort(lambda)
           }
    l = function(Theta, Sigma) {
        res <- determinant(Theta, logarithm = TRUE)$modulus - 
            sum(diag(Sigma %*% Theta))
        return(res)
    }
    n = nrow(Y)
    p = ncol(Y)
    if (cor == TRUE) {
        f = function(x) {
            return(cor(x))
        }
    }
    else {
        f = function(x) {
            return(var(x))
        }
    }
    if (is.null(type)) {
        if (is.null(Target)) {
            Target <- matrix(0, p, p)
        }
        fTarget <- function(Y) {
            Target
        }
    }
    else {
        fTarget <- function(Y) {
            target(Y, type = type, cor = cor)
        }
    }
    if (is.null(ind)) {
        ind = c(1:n)
    }
    indnr <- rep(0L, nfold + 1)
    indnr[1] = 0
    rest = (n - (n%%nfold))/nfold
    for (i in 1:(nfold - 1)) {
        indnr[i + 1] = i * rest
    }
    indnr[nfold + 1] = n
    Targets <- array(NA, c(p, p, nfold))
    for (i in 1:nfold) {
        indtest <- ind[(indnr[i] + 1):indnr[i + 1]]
        Targets[, , i] <- fTarget(Y[-indtest, ])
    }
    CV <- matrix(0, length(lambda), length(alpha))
    for (j in 1:length(lambda)) {
     for (r in 1:length(alpha)) {
        for (i in 1:nfold) {
            indtest <- ind[(indnr[i] + 1):indnr[i + 1]]
            Strain <- f(Y[-indtest, ])
            if (rope == FALSE) {
                Thetatrain <- gelnet(S = Strain, lambda = lambda[j], 
                  alpha = alpha[r], outer.thr = outer.thr, penalize.diagonal = penalize.diagonal, 
                  Target = Targets[, , i])$Theta
            }
            else {
                Thetatrain <- rope(S = Strain, lambda = lambda[j], 
                  Target = Targets[, , i])
            }
            CV[j,r] = CV[j,r] + l(Thetatrain, f(Y[indtest, ]))
        }
    }}
    k = which(CV==max(CV), arr.ind =TRUE)
    if ((k[1] == 1) || (k[1] == length(lambda))) {
        warning("The optimal lambda value is a boundary value. The range of possible values should be increased")
    }
    if ((k[2] == 1) || (k[2] == length(alpha))) {
        warning("The optimal alpha value is a boundary value. The range of possible values should be increased")
    }
    optimalla = lambda[k[1]]
    optimalal = alpha[k[2]]

    if (rope == FALSE) {
        object <- gelnet(S = f(Y), lambda = optimalla, alpha = optimalal, 
            outer.thr = outer.thr, penalize.diagonal = penalize.diagonal, 
            Target = fTarget(Y))
        return(list(optimal = c(optimalla ,optimalal), lambda = lambda, CV = CV/nfold, 
            Theta = object$Theta, W = object$W, niter = object$niter, 
            del = object$del, conv = object$conv))
    }
    else {
        Theta <- rope(S = f(Y), lambda = optimal, Target = fTarget(Y))
        return(list(optimal = optimal, lambda = lambda, CV = CV/nfold, 
            Theta = Theta))
    }
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

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 57:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-d100.RData")
#################################################
n=50
p=20
t=100

Sigma = read.table(NetS20)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT20)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-I-d20.RData")
###############################################
n=50
p=50
t=100

Sigma = read.table(NetS50)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT50)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-I-d50.RData")
###########################################
n=50
p=100
t=100

Sigma = read.table(NetS100)
Sigma = as.matrix(Sigma)

Theta = read.table(NetT100)
Theta = as.matrix(Theta)

A = array(NA,c(t,7))
colnames(A)= c("lam","alpha","loglik","Frobenius.dis","SP","KL","QL")
lambda_grid <- seq(0.01,10,length.out=50)
alpha <- seq(0.01,0.99,length.out=20)

for (i in 1:t){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = mycrossvalidation( Y = X , lambda = lambda_grid , alpha =alpha,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod1-gelnet-vI-d100.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-d100.RData")
#################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-I-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-I-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod2-gelnet-vI-d100.RData")

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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-d100.RData")
#################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-I-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-I-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod3-gelnet-vI-d100.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal[1] ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-d100.RData")
#################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-I-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-I-d50.RData")
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
###############################75
for (i in 10:20){
set.seed(i)
X = mvrnorm(n,rep(0,p),Sigma)
X = scale(X,scale=F)

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal[1] ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
A
}

save(A, file = "mod4-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod4-gelnet-vI-d100.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-d100.RData")
#################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-I-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-I-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod5-gelnet-vI-d100.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5)
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-d100.RData")
#################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-I-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-I-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="Identity"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-I-d100.RData")
#######################################################
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-vI-d20.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-vI-d50.RData")
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

CV = crossvalidation( Y = X , lambda = lambda_grid , alpha =0.5,Target=target(Y=X,type="vI"))
A[i,] = c(CV$optimal ,CV$CV[which(CV$CV==max(CV$CV))][1],fun1(L,CV))
}

save(A, file = "mod6-gelnet-vI-d100.RData")



