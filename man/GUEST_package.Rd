\name{GUEST_package}
\alias{GUEST_package}
\title{
Graphical Models in Ultrahigh-Dimensional and Error-Prone Data via Boosting Algorithm
}
\description{The package GUEST, referred to Graphical models in Ultrahigh-dimensional and Error-prone data via booSTing algorithm, is used to estimate the precision matrix and detect graphical structure for ultrahigh-dimensional, error-prone, and possibly nonlinear random variables. Given the estimated precision matrix, we further apply it to the linear discriminant function to deal with multi-classification. The precision matrix can be estimated by the function \code{boost.graph}, and the classification can be implemented by the function \code{LDA.boost}. Finally, we consider the medulloblastoma dataset to demonstrate the implementation of two functions. }
\details{To estimate the precision matrix and detect the graphical structure under our scenario, the function \code{boost.graph} first applies the regression calibration method to deal with measurement error in continuous, binary, or count random variables. After that, the feature screening technique is employed to reduce the dimension of random variable, and we then adopt the boosting algorithm to estimate the precision matrix. The estimated precision matrix also reflects the desired graphical structure. The function \code{LDA.boost} implements the linear discriminant function to do classification for multi-label classes, where the precision matrix, also known as the inverse of the covariance matrix, in the linear discriminant function can be estimated by the function \code{boost.graph}.
}
\value{GUEST_package
}
