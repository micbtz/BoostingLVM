\name{simdatalvm}
\alias{simdatalvm}
\title{Dataset simulation}
\description{
Simulation of a dataset from a latent variable model.
}
\usage{
simdatalvm(n, param, Phi = NULL, type = "binary")
}
\arguments{
 \item{n}{sample size.}
 \item{param}{matrix containing the parameters. See details.}
 \item{Phi}{Correlation matrix between latent variables. If NULL the dimensions are independent.}
 \item{type}{type of data. Currently only binary data are implemented}
}
\details{
The matrix param has one column for each item.
Each row corresponds to a different dimension.
In the last row there are the threshold parameters.
}
\value{
\code{data.frame} with dimension \eqn{n X p} containing the binary variables,
where $p$ is the number of manifest variables.
}
\examples{
# 10 items and 2 dimensions
param <- matrix(0, 3, 10)
param[1, 1:5] <- runif(5, 0.4, 0.7)
param[2, 6:10] <- runif(5, 0.4, 0.7)
param[3, ] <- runif(10, -0.4, 0.4)
dataset <- simdatalvm(n = 200,  param = param)
}

