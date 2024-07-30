\name{LDA.boost}

\alias{LDA.boost}

\title{Implementation of the linear discriminant function for multi-label classification.}

\description{This function applies the linear discriminant function to do classification for multi-label responses. The precision matrix, or the inverse of the covariance matrix, in the linear discriminant function can be estimated by \code{w} in the function \code{boost.graph}. In addition,  error-prone covariates in  the linear discriminant function are addressed by the regression calibration.
}
\usage{
LDA.boost(data, resp, theta, sigma_e = 0.6,q = 0.8,lambda = 1, pi = 0.5)
}

\arguments{
  \item{data}{
  An n (observations) times p (variables) matrix of random variables, whose distributions can be continuous, discrete, or mixed.
}
  \item{resp}{
  An n-dimensional vector of categorical random variables, which is the response in the data.
}
  \item{theta}{
  The estimator of the precision matrix.
}
\item{sigma_e}{
  The common value in the diagonal covariance matrix of the error for the classical measurement error model when \code{data} are continuous. The default value is 0.6.
}
  \item{q}{
  The common value used to characterize misclassification for binary random variables. The default value is 0.8.
}
  \item{lambda}{
  The parameter of the Poisson distribution, which is used to characterize error-prone count random variables. The default value is 1.
  }
  \item{pi}{
  The probability in the Binomial distribution, which is used to characterize error-prone count random variables. The default value is 0.5.
  }
}

\details{
The linear discriminant function used is as follow: \cr
  \deqn{
  \code{score}_{i,j} = \log (\pi _i) - 0.5\ \mu_{i}^\top\  \code{theta}\ \mu _{i} + \code{data}_{j}^\top\ \code{theta}\ \mu_{i},
  }\cr
for the class \eqn{i = 1, \cdots, I} with \eqn{I} being the number of classes in the dataset and subject \eqn{j = 1, \cdots, n}, where \eqn{\pi _i} is the proportion of subjects in the class \eqn{i}, \eqn{\code{data}_{j}} is the vector of covariates for the subject \eqn{j}, \eqn{\code{theta}} is the precision matrix of the covariates, and \eqn{\mu_{i}} is the empirical mean vector of the random variables in the class \eqn{i}.
}

\value{
  \item{score}{
  The value of the linear discriminant function (see details) with the estimator of the precision matrix accommodated.
}
  \item{class}{
  The result of predicted class for subjects.
}
}

\references{
Hui-Shan Tsao (2024). \emph{Estimation of Ultrahigh-Dimensional Graphical Models and Its Application to Dsicriminant Analysis.} Master Thesis supervised by Li-Pang Chen, National Chengchi University.
}

\author{
Hui-Shan Tsao and Li-Pang Chen\cr
  Maintainer: Hui-Shan Tsao \email{n410412@gmail.com}
}

\examples{
data(MedulloblastomaData)

X <- t(MedulloblastomaData[2:655,]) #covariates
Y <- MedulloblastomaData[1,] #response

X <- matrix(as.numeric(X),nrow=23)

p <- ncol(X)
n <- nrow(X)

#standarization
X_new=data.frame()
for (i in 1:p){
 X_new[1:n,i]=(X[,i]-rep(mean(X[,i]),n))/sd(X[,i])
}
X_new=matrix(unlist(X_new),nrow = n)

\donttest{
#estimate graphical model
result <- boost.graph(data = X_new, thre = 0.2, ite1 = 3, ite2 = 0, ite3 = 0, rep = 1)
theta.hat <- result$w

theta.hat[which(theta.hat<0.8)]=0 #keep the highly dependent pairs

#predict
pre <- LDA.boost(data = X_new, resp = Y, theta = theta.hat)
estimated_Y <- pre$class}
}

\keyword{classif}
