\name{boost.graph}

\alias{boost.graph}

\title{Estimation of precision matrix and detection of graphical structure}

\description{This function first applies the regression calibration to deal with measurement error effects. After that, the feature screening technique is employed to screen out independent pairs of random variables and reduce the dimension of random variables. Finally, we adopt the boosting method to detect informative pairs of random variables and estimate the precision matrix. This function can handle various distributions, such as normal, binomial, and Poisson distributions, as well as nonlinear effects among random variables.
}

\usage{
boost.graph(data,ite1,ite2,ite3,thre,select = 0.9,inc = 10^(-3),
sigma_e = 0.6,q = 0.8,lambda = 1,pi = 0.5,rep = 100,cor = TRUE)
}

\arguments{
  \item{data}{
  An n (observations) times p (variables) matrix of random variables, whose distributions can be continuous, discrete, or mixed.
}
  \item{ite1}{
  The number of iterations for continuous variables.
}
  \item{ite2}{
  The number of iterations for binary variables.
}
\item{ite3}{
  The number of iterations for count variables.
}
  \item{thre}{
  The treshold value for feature screening, whose value should be between 0 and 1.
}
\item{select}{
The treshold constant in the boosting algorithm, whose value should be between 0 and 1. The default value is 0.9.
}
  \item{inc}{
  The learning rate of the increment in the boosting algorithm, which shoud be a small value. The default value is 0.001.
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
  \item{rep}{
  The number of bootstrapping iterations. The default value is 100.
}
  \item{cor}{
  Measurement error correction when estimating the precision matrix. The default value is TRUE.
}
}

\value{
  \item{w}{
  The estimator of the precision matrix.
}
  \item{p}{
  The chosen pairs obtained by the feature screening.
}
  \item{xi}{
  The weights sorted with pairs in \code{p}.
}
  \item{g}{
  The visualization of the estimated network structure determined by \code{w}.
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

X <- t(MedulloblastomaData[2:656,]) #covariates
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
theta.hat <- result$w}
}

\keyword{multivariate}
\keyword{models}
\keyword{graphs}
