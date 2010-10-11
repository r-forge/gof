\name{cumres}
\alias{cumres}
\alias{cumres.lm}
\alias{cumres.glm}
\title{Calculates GOF measures for for generalized linear models}
\description{
  Currently linear regression models (lm) and logistic and poisson
  regression models are supported.  
}
\usage{
\method{cumres}{lm} (model, ...)
\method{cumres}{glm} (model,
         variable=c("predicted",colnames(model.matrix(model))),
         data=data.frame(model.matrix(model)),
         R=500, b=0, plots=min(R,50),
         breakties=1e-12,
         seed=round(runif(1,1,1e9)),...)
}
\arguments{
  \item{model}{Model object (\code{lm} or \code{glm})}
  \item{variable}{List of variable to order the residuals after}
  \item{data}{data.frame used to fit model (complete cases)}
  \item{R}{Number of samples used in simulation}
  \item{b}{Moving average bandwidth (0 corresponds to infinity =
    standard cumulated residuals)}
  \item{plots}{Number of realizations to save for use in the
    plot-routine}
  \item{breakties}{Add unif[0,breakties] to observations}
  \item{seed}{Random seed}
  \item{...}{additional arguments}
}
\value{
  Returns an object of class 'cumres'.
}
\author{Klaus K. Holst}

\seealso{\code{\link[timereg]{cox.aalen}} in the \code{timereg}-package
  for similar GOF-methods for survival-data.}
\examples{
sim1 <- function(n=100, f=function(x1,x2) {10+x1+x2^2}, sd=1, seed=1) {
  if (!is.null(seed))
    set.seed(seed)
  x1 <- rnorm(n);
  x2 <- rnorm(n)
  X <- cbind(1,x1,x2)
  y <- f(x1,x2) + rnorm(n,sd=sd)
  d <- data.frame(y,x1,x2)
  return(d)
}
d <- sim1(100); l <- lm(y ~ x1 + x2,d)
system.time(g <- cumres(l, R=100, plots=50))
g
\dontrun{plot(g)}
g1 <- cumres(l, c("y"), R=100, plots=50)
g1
g2 <- cumres(l, c("y"), R=100, plots=50, b=0.5)
g2

}
\keyword{models}
\keyword{regression}
