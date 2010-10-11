confint.cumres <- function(object, parm=1, level=0.95, cval=NULL, ...) {
  t <- c(); yu <- c()
  for (idx in parm) {
    if (is.null(cval))
      cval <- quantile(object$cvalues[,idx], level)
    ##  y <- x$W[,i];
    t <- cbind(t,object$x[,idx])
    yu <- cbind(yu,cval*object$sd[,idx]);
  }
  return(list(t=t,yu=yu));
}
