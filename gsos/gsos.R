gsos <-function(S,lambda,alpha,zero=NULL,Theta=NULL,W=NULL,Target=NULL,outer.maxit=1000,outer.thr=1e-5,inner.maxit=1000,inner.thr=outer.thr,penalize.diagonal=TRUE,active=FALSE){
  
  BIG <- 10e9
  
  p <- nrow(S)
  #if (is.null(S)) { print("S is required as input"); return(); }else { p=nrow(S) }
  #check1<-(nrow(S)==ncol(S))&&(S==t(S));
  #if (check1==0) { print("S should be square and symmetric"); return(); }
  
  if(!is.matrix(lambda) & length(lambda) != 1 & length(lambda) != nrow(S))
  {stop("Wrong number of elements in lambda")}
  
  if(is.vector(lambda) && length(lambda) > 1){ lambda = matrix(sqrt(lambda))%*%sqrt(lambda)}
  if(length(lambda) == 1){lambda = matrix(lambda,ncol=p,nrow=p)}
  
  if(is.vector(Target) && length(Target) > 1){Target=diag(Target)}
  
  
  
  if(!is.null(zero)){
    if(!is.matrix(zero)){ zero=matrix(zero,nrow=TRUE)}
    for(k in 1:nrow(zero)){
      i=zero[k,1]
      j=zero[k,2]
      lambda[i,j]=BIG
      lambda[j,i]=BIG
    }
  }
  
  
  thr <- outer.thr; thr2 <- inner.thr
  lambda1 <- diag(lambda)*alpha; lambda2 <- diag(lambda)*(1-alpha)
  conv <- TRUE
  
  
  
  if(is.null(Target)) {Target <- matrix(0,p,p)}
  Theta <- matrix(0,p,p)
  for(i in 1:p){
          root1 <- (-S[i,i]-lambda1[i]+lambda2[i]*Target[i,i] + sqrt((S[i,i]+lambda1[i]-lambda2[i]*Target[i,i])^2 +4*lambda2[i]) ) / (2*lambda2[i])
          root2 <- (-S[i,i]-lambda1[i]+lambda2[i]*Target[i,i] - sqrt((S[i,i]+lambda1[i]-lambda2[i]*Target[i,i])^2 +4*lambda2[i]) ) / (2*lambda2[i])
          Theta[i,i] <- max(root1,root2)
        }
   W <- S; diag(W) <- diag(S) + lambda1 + lambda2*diag(Theta-Target)
  
  #p1 <- nrow(Theta); p11 <- ncol(Theta); p2 <- nrow(W); p22 <- ncol(W);
  #if (p1!=p11) { print("check dimensions of Theta"); return(); }
  #if (p2!=p22) { print("check dimensions of W"); return(); }
  
  #check2 <- (p1==p2)&&(p1==p);
  #if (check2==0) { print("check dimensions of Theta,W,S"); return(); }
  niter <- 0L
  dlz <- 0
  
  mode(p)="integer"
  mode(S)="double"
  mode(lambda)="double"
  mode(alpha)="double"
  mode(Theta)="double"
  mode(W)="double"
  mode(Target)="double"
  mode(outer.maxit)="integer"
  mode(thr)="double"
  mode(inner.maxit)="integer"
  mode(thr2)="double"
  mode(niter)="integer"
  mode(penalize.diagonal)="logical"
  mode(active)="logical"
  mode(dlz)="double"
  
  
  
  loop<-.Fortran("gsos",p,S,lambda,alpha,Theta,W,Target,outer.maxit,thr,inner.maxit,thr2,niter,penalize.diagonal,active,dlz)
  
  # Theta <- loop$TTh; W <- loop$Wm
  # del <- loop$dlz
  
  #if(anyNA(Theta) | anyNA(W)){ print("error"); return(); }
  
  if(loop[[12]] == outer.maxit){ conv <- FALSE }
  #return(list(Theta=Theta,W=W,niter=1,del=1,conv=1))
  return(list(Theta=loop[[5]],W=loop[[6]],niter=loop[[12]],del=loop[[15]],conv=conv))
}
