#' @export
boost.graph<-function(data,ite1,ite2,ite3,thre,select = 0.9,inc = 10^(-3),
                      sigma_e = 0.6,q = 0.8,lambda = 1, pi = 0.5,rep = 100,cor=TRUE){

  p=ncol(data)
  n=nrow(data)

  #record every distribution of X_i
  model=NULL
  for (i in 1:p){
    if (sum(data[,i]%%1==0)==n & sum(data[,i]!=0 & data[,i]!=1)>0)
      model=c(model,"counts")
    else if (sum(data[,i]%%1==0)==n & sum(data[,i]!=0 & data[,i]!=1)==0)
      model=c(model,"binary")
    else
      model=c(model,"continuous")
  }

  if (rep!=1){
    Xhat_array=array(dim = c(n,p,rep))
    XI_array=array(dim = c(p,p,rep))
    XIstar_array=array(dim = c(p,p,rep))
    V_array=array(dim = c(p,p,rep))
    pair=NULL

    #bootstrap
    for (i in 1:rep){
      subject=sample(1:n,n,replace = T)
      X=data[subject,]

      #do correction according to the distribution
      if (cor==TRUE)
      {
        Xhat=data.frame()
        for (j in 1:p){
          if (model[j]=="continuous")
            Xhat[1:n,j]=mean(X[,j])+(sd(X[,j])-sigma_e)*(sd(X[,j]))^-1*(X[,j]-mean(X[,j]))
          else if (model[j]=="binary")
            Xhat[1:n,j]=(X[,j]-rep(1,n)+rbinom(n,1,q))/(2*rbinom(n,1,q)-rep(1,n))
          else
            Xhat[1:n,j]=((mean(X[,j])-lambda)/(1-pi))+((mean(X[,j])-lambda)/(1-pi))*(mean(X[,j])-lambda)/(lambda+(3*pi+1)/(1-pi)*(mean(X[,j])-lambda))*(X[,j]-mean(X[,j]))
        }
        Xhat=matrix(unlist(Xhat),nrow=n)
      }
      else
        Xhat=X

      Xhat_array[,,i]=Xhat

      #calculate XI
      XI = matrix(0,p,p)
      for(j in 1:p) {
        for(k in 1:j) {
          XI[j,k] = max(XICOR::xicor(Xhat[,j],Xhat[,k]),XICOR::xicor(Xhat[,k],Xhat[,j]))*1
          #XI[j,k] = max(abs(xicor(Xhat[,j],Xhat[,k])),abs(xicor(Xhat[,k],Xhat[,j])))*1
          XI[k,j]=XI[j,k]
        }
      }
      XI=abs(XI)

      #calculate XI_star without considering the diagonal entries
      temp=XI
      temp[which(col(temp)==row(temp))]=0
      XI_star=(XI-matrix(min(temp),p,p))/(max(temp)-min(temp))
      XI_star[which(col(XI_star)==row(XI_star))]=diag(XI)

      XI_array[,,i]=XI
      XIstar_array[,,i]=XI_star

      #record the chosen pairs:V[2,1]=1 means {x1,x2} is chosen
      V=matrix(0,p,p)
      for (j in 1:p){
        V[which(XI_star[,j]>thre),j]=1
      }
      V_array[,,i]=V

      #combine the chosen pairs (#total=rep)
      pair0 = NULL
      for(i in 1:p) {
        for(j in 1:p) {
          if(i<j && V[i,j] ==1) {
            pair0 = rbind(pair0,c(i,j))
          }
        }
      }
      pair=rbind(pair,pair0)
      pair=unique(pair)
      pair=pair[order(pair[,1]),]

      #record the XI values of the chosen pairs
      XI_new=NULL
      for (i in 1:length(pair[,1])){
        XI_new=c(XI_new,XI_star[pair[i,1],pair[i,2]])
      }
    }

    #boost.graph
    Theta_hat=rep(0,length(pair[,1]))
    for (i in 1:rep){
      theta=rep(0,length(pair[,1]))
      for (j in unique(pair[,1])){
        r=which(V_array[j,,i]==1)

        if (length(r)!=1){
          r=r[-which(r==j)]
          Xr=Xhat_array[,,i]

          for (k in 1:length(r)){
            Xr[which(col(Xr)==r[k])]=t(rep(0,n))
          }

          s=rep(0,p)
          theta0=rep(0,p) #starting point
          #iterate t times to calculate s

          if (sum(Xhat_array[j,i]%%1==0)==n & sum(Xhat_array[j,i]!=0 & Xhat_array[j,i]!=1)>0) #counts
            ite=ite3
          else if (sum(Xhat_array[j,i]%%1==0)==n & sum(Xhat_array[j,i]!=0 & Xhat_array[j,i]!=1)==0) #binary
            ite=ite2
          else
            ite=ite1

          for (l in 1:ite){
            if (sum(Xhat_array[j,i]%%1==0)==n)
              s=as.vector(t(Xhat_array[,,i]-Xr)%*%(Xhat_array[,j,i]-((Xhat_array[,,i]-Xr)%*%theta0)))/(2*n)
            else
              s=as.vector(t(Xhat_array[,,i]-Xr)%*%(Xhat_array[,j,i]-((Xhat_array[,,i]-Xr)%*%theta0)))/(-2*n)

            tau=matrix(0,length(r),p)
            for (k in 1:length(r)){
              temp=tau[k,]
              temp[which(abs(XIstar_array[,j,i]*s)>=select*max(XIstar_array[,j,i]*s))]=1
              temp[which(temp!=0)]=which(temp!=0)
              tau[k,]=temp
            }#one row for one r

            #update theta0
            for (k in 1:length(r)){
              if (sum(Xhat_array[j,i]%%1==0)==n)
                theta0[tau[k,]]=theta0[tau[k,]]+inc*abs(XIstar_array[tau[k,],r[k],i]*s[tau[k,]])
              else
                theta0[tau[k,]]=theta0[tau[k,]]+inc*sign(XIstar_array[tau[k,],r[k],i]*s[tau[k,]])
            }
          }
          theta[which(pair[,1]==j)]=theta[which(pair[,1]==j)]+theta0[pair[which(pair[,1]==j),2]]
        }
      }
      Theta_hat=Theta_hat+theta
    }
    Theta_hat=Theta_hat/rep
  }

  else{
    if (cor==TRUE){
      Xhat=data.frame()
      for (j in 1:p){
        if (model[j]=="continuous")
          Xhat[1:n,j]=mean(data[,j])+(sd(data[,j])-sigma_e)*(sd(data[,j]))^-1*(data[,j]-mean(data[,j]))
        else if (model[j]=="binary")
          Xhat[1:n,j]=(data[,j]-rep(1,n)+rbinom(n,1,q))/(2*rbinom(n,1,q)-rep(1,n))
        else
          Xhat[1:n,j]=((mean(data[,j])-lambda)/(1-pi))+((mean(data[,j])-lambda)/(1-pi))*(mean(data[,j])-lambda)/(lambda+(3*pi+1)/(1-pi)*(mean(data[,j])-lambda))*(data[,j]-mean(data[,j]))
      }
    }
    else
      Xhat=data

    Xhat=matrix(unlist(Xhat),nrow = n)

    XI = matrix(0,p,p)
    for(j in 1:p) {
      for(k in 1:j) {
        XI[j,k] = max(XICOR::xicor(Xhat[,j],Xhat[,k]),XICOR::xicor(Xhat[,k],Xhat[,j]))*1
        XI[k,j]=XI[j,k]
      }
    }
    XI=abs(XI)
    temp=XI
    temp[which(col(temp)==row(temp))]=0
    XI_star=(XI-matrix(min(temp),p,p))/(max(temp)-min(temp))
    XI_star[which(col(XI_star)==row(XI_star))]=diag(XI)

    #record the chosen pairs:V[2,1]=1 means {x1,x2} is chosen
    V=matrix(0,p,p)
    for (j in 1:p){
      V[which(XI_star[,j]>thre),j]=1
    }

    pair = NULL; XI_new=NULL
    for(i in 1:p) {
      for(j in 1:p) {
        if(i<j && V[i,j] ==1) {
          pair = rbind(pair,c(i,j))
          XI_new=c(XI_new,XI_star[i,j])
        }
      }
    }
    theta=rep(0,length(pair[,1]))
    Theta_hat=theta
    for (i in unique(pair[,1])){
      r=which(V[i,]==1)

      if (length(r)!=1){
        r=r[-which(r==i)]

        Xr=Xhat
        for (j in 1:length(r)){
          Xr[which(col(Xr)==r[j])]=t(rep(0,n))
        }

        s=rep(0,p)
        theta0=rep(0,p) #starting point
        #iterate t times to calculate s
        if (model[i]=="continuous")
          ite=ite1
        else if (model[i]=="binary")
          ite=ite2
        else
          ite=ite3

        for (l in 1:ite){
          if (model[i]=="continuous")
            s=as.vector(t(Xhat-Xr)%*%(Xhat[,i]-((Xhat-Xr)%*%theta0)))/(-2*n)
          else
            s=as.vector(t(Xhat-Xr)%*%(Xhat[,i]-((Xhat-Xr)%*%theta0)))/(2*n)

          tau=matrix(0,length(r),p)
          for (k in 1:length(r)){
            temp=tau[k,]
            temp[which(abs(XI_star[,i]*s)>=select*max(XI_star[,i]*s))]=1
            temp[which(temp!=0)]=which(temp!=0)
            tau[k,]=temp
          }#one row for one r

          #update theta
          for (k in 1:length(r)){
            if (model[i]=="continuous")
              theta0[tau[k,]]=theta0[tau[k,]]+inc*sign(XI_star[tau[k,],r[k]]*s[tau[k,]])
            else
              theta0[tau[k,]]=theta0[tau[k,]]+inc*abs(XI_star[tau[k,],r[k]]*s[tau[k,]])
          }
        }
        theta[which(pair[,1]==i)]=theta[which(pair[,1]==i)]+theta0[pair[which(pair[,1]==i),2]]
      }
    }
    Theta_hat=theta
  }

  #the final estimator
  theta.hat=diag(1,p,p)
  for (i in 1:length(pair[,1])){
    theta.hat[pair[i,1],pair[i,2]]=Theta_hat[i]
    theta.hat[pair[i,2],pair[i,1]]=theta.hat[pair[i,1],pair[i,2]]
  }

  #the network structure
  net = theta.hat
  net = network::network(net, directed = FALSE)
  network::network.vertex.names(net)=paste0("X",network::network.vertex.names(net))
  graph = GGally::ggnet2(net,size=3,node.color = "lightgray",label=T,label.size = 3,mode = "circle")

  #output
  result_list <- list(
    w = abs(theta.hat),
    p = pair,
    xi = XI_new,
    g = graph
  )

  return(result_list)
}
