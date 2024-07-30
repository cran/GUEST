#' @export
LDA.boost=function(data, resp, theta,  sigma_e = 0.6,q = 0.8,lambda = 1, pi = 0.5){

  p=ncol(data)
  n=nrow(data)

  #record the distributions of the data
  model=NULL
  for (i in 1:p){
    if (sum(data[,i]%%1==0)==n & sum(data[,i]!=0 & data[,i]!=1)>0)
      model=c(model,"counts")
    else if (sum(data[,i]%%1==0)==n & sum(data[,i]!=0 & data[,i]!=1)==0)
      model=c(model,"binary")
    else
      model=c(model,"continuous")
  }

  #correction
  Xhat=data.frame()
  for (j in 1:p){
    if (model[j]=="continuous")
      Xhat[1:n,j]=mean(data[,j])+(sd(data[,j])-sigma_e)*(sd(data[,j]))^-1*(data[,j]-mean(data[,j]))
    else if (model[j]=="binary")
      Xhat[1:n,j]=(data[,j]-rep(1,n)+rbinom(n,1,q))/(2*rbinom(n,1,q)-rep(1,n))
    else
      Xhat[1:n,j]=((mean(data[,j])-lambda)/(1-pi))+((mean(data[,j])-lambda)/(1-pi))*(mean(data[,j])-lambda)/(lambda+(3*pi+1)/(1-pi)*(mean(data[,j])-lambda))*(data[,j]-mean(data[,j]))
  }
  Xhat=matrix(unlist(Xhat),nrow=n)

  x=Xhat
  y=as.factor(resp)

  classes <- length(unique(y)) #the number of y's levels

  #the within-group average
  mu <- NULL
  for (i in 1:classes){
    mu[[i]] <- colMeans(x[which(y==sort(unique(y))[i]),])
  }

  #the proportion of the group
  pi <- NULL
  for (i in 1:classes){
    pi[i]=summary(y)[i]/length(y)
  }

  #calculate the LDA scores
  delta <- matrix(rep(0,classes*length(y)),ncol=length(y))
  for(i in 1:classes){
    for(j in 1:length(y)){
      delta[i,j] <- log(pi[i])- 0.5 * (matrix(mu[[i]],nrow=1) %*% theta %*% matrix(mu[[i]],ncol=1))+ (x[j,]%*% theta %*% matrix(mu[[i]],ncol=1))
    }
  }

  #classification
  class_est <- NULL
  for(j in 1:length(y)){
    class_est[j] <- which(delta[,j]==max(delta[,j]))
  }
  class_est <- as.factor(class_est)
  levels(class_est)=levels(y)

  #output
  result_list <- list(
    score=delta,
    class=class_est
  )

  return(result_list)
}
