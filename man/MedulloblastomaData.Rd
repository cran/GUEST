\name{MedulloblastomaData}
\alias{MedulloblastomaData}
\docType{data}
\title{The medulloblastoma dataset
}
\description{
The dataset, which is available on https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE468, contains 23 patients with medulloblastoma, and each patient has 2059 gene expression values. The response contains 2 classes: metastatic (M+) or non-metastatic (M0). After removing the missing and duplicate values, the dimension of remaining gene expressions is 655. The dataset is used to illustrate the usage of the \code{\link{boost.graph}} and \code{\link{LDA.boost}} functions.
}
\usage{data(MedulloblastomaData)}
\format{The dataset has 23 observations and 655 gene expression values.
}
\references{
MacDonald, T., Brown, K., LaFleur, B., Peterson K., Lawlor C., Chen Y., Packer RJ., Cogen P., Stephan DA.(2001). \emph{Expression profiling of medulloblastoma: PDGFRA and the RAS/MAPK pathway as therapeutic targets for metastatic disease}. Nat Genet, 29, 143–152.
}
\examples{
X <- t(MedulloblastomaData[2:655,]) #covariates
Y <- MedulloblastomaData[1,] #response

}
\keyword{datasets}
