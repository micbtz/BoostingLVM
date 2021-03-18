\name{lvmboost}
\alias{lvmboost}
\title{
Fit a latent variable model by likelihood-based boosting 
}
\description{
\code{lvmboost} is used to fit a latent variable model for binary responses by componentwise likelihood based boosting.
}
\usage{
lvmboost(data = NULL, trace = TRUE, nu = 0.05, tau = 2, Dmax = 5, kstop = 100, 
         obj.lvmboost = NULL, nfreepar = TRUE)
}
\arguments{
\item{data}{a data frame or matrix containing the binary item responses.}
\item{trace}{logical; if \code{TRUE} tracing information is produced. Set to \code{FALSE} to increase speadness.}
\item{nu}{step length.}
\item{tau}{value used to determine the choice between the Newton-type direction and the negative curvature direction.}
\item{Dmax}{maximum number of dimensions.}
\item{kstop}{number of maximum iterations.}
\item{obj.lvmboost}{output of function \code{lvmboost}. This is used to continue  iterating and no other argument is necessary.}
\item{nfreepar}{logical; if \code{TRUE} the iterations are stopped if the number 
of estimated parameters exceeds the number of free parameters.}
}
\details{
}
\value{
List containing the following components:
\item{parfin}{parameter estimates at the last iteration.}
\item{par}{list containing the parameter estimates at each iteration.}
\item{loglik}{vector of negative log-likelihood values at each iteration.}
\item{D}{number of dimensions at the last iteration.}
\item{nu}{step length.}
\item{Dmax}{maximum number of dimensions.}
\item{bifreq}{list containing the frequencies for pais of items.}
\item{data}{input data frame.}
}
\seealso{
\code{\link{cv.lvmboost}}, \code{\link{plot.lvmboost}}
}
\examples{
# generation of item parameters
nitems <- 20 # number of item parameters
a <- matrix(0, nitems, 2)
a[1:(nitems/2), 1] <- runif((nitems/2), 0.5, 0.8)
a[(nitems/2+1):nitems, 2] <- runif((nitems/2), 0.5, 0.8)
param <- rbind(t(a), 0.2)

# simulation of data
dataset <- simdatalvm(n = 2000, param = param)

# fit the model by boosting
res <- lvmboost(dataset, kstop = 20)
# continue iterating
res <- lvmboost(obj.lvmboost = res, kstop = 10)

# regularization path
plot(res)

# cross-validation
cv<-cv.lvmboost(res)
cv$par_sel
cv$sel # more iterations are needed
}


