\name{confint.cumres}
\alias{confint.cumres}
\title{Returns prediction bands for 'cumres' object}
\description{
  \code{confint} returns prediction bands for the cumulative residual
  process under the null.
}
\usage{
\method{confint}{cumres} (object, parm=1, level=0.95, cval=NULL, ...)
}

\arguments{
  \item{object}{Object produced by the function \code{cumres}}
  \item{parm}{vector of numbers indicating which processes from the \code{x} to
    calculate prediction bands for.} 
  \item{level}{The required prediction level.}
  \item{cval}{Overrules the level-parameter by calculating symmetric prediction
    bands defined by the standard error multiplied by \code{cval}.}
  \item{...}{Additional arguments.}
}
\value{
  list with the following members:

  \item{t}{Ordered values of variable that is used to cumulate
  residuals after}
  \item{yu}{Upper simultaneous confidence limit.}
}
\seealso{
  \code{\link[gof]{cumres}}
}
\examples{
n <- 500; x <- abs(rnorm(n,sd=0.2))+0.01; y <- sqrt(x) + rnorm(n,sd=0.2)
l <- lm(y ~ x)
g <- cumres(l, R=1000)
confint(g,1)
}
\author{Klaus K. Holst <kkho@biostat.ku.dk>}
\keyword{models}
\keyword{regression}
