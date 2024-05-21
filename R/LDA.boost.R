#' @export
LDA.boost=function(data, resp, theta){

  x=data
  y=as.factor(resp)

  classes <- length(unique(y)) #the number of y's levels

  #the within-group average
  mu <- NULL
  for (i in 1:classes){
    mu[[i]] <- colMeans(x[which(y==unique(y)[i]),])
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
