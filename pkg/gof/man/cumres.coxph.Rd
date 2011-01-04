\name{cumres.coxph}
\alias{cumres.coxph}
\title{Calculates GOF measures for Cox's propoportional hazard model for right
  censored survival times
}
\description{
  Calculates score processes and KS and Cvm tests for proportionaly of
  hazards via simulation (Martinussen and Scheike, 2006).
}
\usage{
\method{cumres}{coxph} (model,
         variable=c(colnames(model.matrix(model))),
         type=c("score","residual"),
         R=500, plots=min(R,50), seed=round(runif(1,1,1e9)), ...)
}
\arguments{
  \item{model}{Model object (\code{lm} or \code{glm})}
  \item{variable}{List of variable to order the residuals after}
  \item{R}{Number of samples used in simulation}
  \item{type}{Type of GOF-procedure},
  \item{plots}{Number of realizations to save for use in the
    plot-routine}
  \item{seed}{Random seed}
  \item{...}{additional arguments}
}
\value{
  Returns an object of class 'cumres'.
}
\references{
  Lin, D. Y. and Wei, L. J. and Ying, Z. (1993)
  \emph{Checking the Cox model with cumulative sums of martingale-based
    residuals} Biometrika, Volume 80, No 3, p. 557-572.

  Martinussen, Torben and Scheike, Thomas H.
  \emph{Dynamic regression models for survival data} (2006), 
  Springer, New York.

}
\seealso{
  \code{\link[gof]{cumres.glm}},
  \code{\link[survival]{coxph}}, and
  \code{\link[timereg]{cox.aalen}} in the \code{timereg} package
  for similar GOF-methods for survival-data. 
}
\examples{
library(survival)

simcox <- function(n=100, seed=1) {
  if (!is.null(seed))
    set.seed(seed)
  require(survival)
  time<-rexp(n); cen<-2*rexp(n); 
  status<-(time<cen); 
  time[status==0]<-cen[status==0]; 
  X<-matrix(rnorm(2*n),n,2)  
  return(data.frame(time=time, status=status, X))
}
n <- 100; d <- simcox(n); model <- coxph(Surv(time,status)~ X1 + X2, data=d, robust=TRUE)

dtimes <- sort(d$time[d$status==1])
system.time(a <- cumres(model))
summary(a)


## PBC example
data(pbc)
fit.cox <- coxph(Surv(time,status==2) ~ age + edema + bili + protime, data=pbc)
system.time(pbc.gof <- cumres(fit.cox,R=2000))
\donttest{
par(mfrow=c(2,2))
plot(pbc.gof, ci=TRUE, legend=NULL)
}

}
\author{Klaus K. Holst and Thomas Scheike}
\keyword{models}
\keyword{regression}

