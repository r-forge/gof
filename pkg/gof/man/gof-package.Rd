\name{gof-package}
\alias{gof-package}
\alias{gof}
\docType{package}
\title{Model-diagnostics based on cumulative residuals}
\description{
  Implementation of model-checking technique based on cumulative residuals
}
\details{
\tabular{ll}{
Package: \tab gof\cr
Type: \tab Package\cr
Version: \tab 0.8-2\cr
Date: \tab 2011-06-14\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
}
}
\author{
  Klaus. K. Holst <kkho@biostat.ku.dk>
}
\references{
D.Y. Lin and L.J. Wei and Z. Ying (2002) \emph{Model-Checking Techniques
  Based on Cumulative Residuals}. Biometrics, Volume 58, pp 1-12.

John Q. Su and L.J. Wei (1991) \emph{A lack-of-fit test for the mean function
  in a generalized linear model}. Journal. Amer. Statist. Assoc.,
Volume 86, Number 414, pp 420-426.
}
\keyword{ package }
\seealso{\code{\link[timereg]{cox.aalen}} in the \code{timereg}-package
  for similar GoF-methods for survival-data.}
\examples{
example(cumres)

}
