\name{cv.lvmboost}
\alias{cv.lvmboost}
\title{Cross-validation}
\description{
Selection of the number of iterations of the boosting algorithm by cross-validation
}
\usage{
cv.lvmboost(obj.lvmboost, K = 5, trace = TRUE, seed = 1)
}
\arguments{
\item{obj.lvmboost}{an object of class \code{lvmboost} obtained from function \code{\link{lvmboost}}.}
\item{K}{number of folds.}
\item{trace}{logical; if \code{TRUE} tracing information is produced.}
\item{seed}{integer value used to specify the seed.}
}
\value{
List containing the following components:
\item{cv_err}{cross validation error.}
\item{sel}{iteration selected.}
\item{par_sel}{parameter estimates at the iteration selected.}
}
\seealso{
\code{\link{lvmboost}}
}

