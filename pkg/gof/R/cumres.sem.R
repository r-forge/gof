### `cumres` <-
###   function(model,...) UseMethod("cumres")

`Debug` <-
  function(msg, cond) {
    if (cond)
      print(paste(msg, collapse=" "))
  }

################################################################################

cumres.lvmfit <-
  function(model,y,x,cond,
           data=model.frame(model),p,
           R=100, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
           debug=FALSE, ...) {
    if (!require("lava")) stop("package 'lava' not available")
    if (!require("numDeriv")) stop("package 'numDeriv' not available")
    
    
    hatW.MC <- function(x) {
      ord <- order(x)
      output <- .C("W",
                   R=as.integer(R), ## Number of realizations
                   b=as.double(b), ## Moving average parameter
                   n=as.integer(n), ## Number of observations
                   npar=as.integer(nrow(Ii)), ## Number of parameters (columns in design)
                   xdata=as.double(x), ## Observations to cummulate after 
                   rdata=as.double(r), ## Residuals (one-dimensional)
                   betaiiddata=as.double(beta.iid), ## Score-process
                   etarawdata=as.double(grad), ## Eta (derivative of terms in cummulated process W)
                   plotnum=as.integer(plots), ## Number of realizations to save (for later plot)
                   seed=as.integer(seed), ## Seed (will probably rely on R's randgen in future version)
                   KS=as.double(0), ## Return: Kolmogorov Smirnov statistics for each realization
                   CvM=as.double(0), ## Return: Cramer von Mises statistics for each realization
                   Wsd=as.double(numeric(n)), ## Return: Pointwise variance of W(x)
                   cvalues=as.double(numeric(R)), ## Return: value for each realization s.t.  +/- cvalue * Wsd contains W*
                   Ws=as.double(numeric(plots*n)), ## Return: Saved realizations (for plotting function)
                   Wobs=as.double(numeric(n)) ## Observed process
                   , PACKAGE="gof")
      return(list(output=output,x=x[ord]))
    }

    x0 <- x
    n <- length(x0)
##    y <- u2~z1+z2+z3+y1+y2+y3
##    x <- ~u1
##    y <- "y2"
##    y0 <- predict(model,y=y,resid=FALSE)
    ##    y0 <- predict(model,y=y,resid=FALSE)
##    x0 <- predict(model,y=x,resid=FALSE)
##    x0 <- attributes(predict(model))$eta.x["u1",]
##    x0 <- predict(e,x=~y1+y2+y3,resid=TRUE)[,"u1"]
##    x0 <- predict(model,x=~y1+y3+y4+y5,resid=FALSE)[,"u1"]
##    x0 <- attributes(predict(model))$blup[,"u1"]
    
##    x0 <- attributes(predict(e))$epsilon.y[,"u1"]
##    x0 <- attributes(predict(e))$epsilon.y[,"u1"]
    if (missing(p))
      p <- lava::pars(model)

    myres <- function(p) {
      attributes(predict(model,p=p))$epsilon.y[,y]
##      predict(model,x=endogenous(model),resid=TRUE,p=p)[,"u2"]
##      predict(model,x=~z1+z2+z3+y1+y2+y3,resid=FALSE,p=p)[,"u2"]
    }
    grad <- -jacobian(myres,p,method=lava::lava.options()$Dmethod)
    r <- myres(p)
    ##    Ii <- solve(information(model,p), num=TRUE)
    Ii <- model$vcov
##    Debug(list("Ii=",Ii),debug)
##    Score <- score(model,data=data,p=p,indiv=TRUE)
##    beta.iid <- Ii%*%t(Score)    
    ord <- order(x0)
    x0 <- x0[ord]
    grad <- grad[ord,]
    r <- r[ord]
    Ii <- model$vcov
##    Debug(list("Ii=",Ii),debug)
    Score <- lava::score(model,data=data,p=p,indiv=TRUE,weight=lava::Weight(model))[ord,]
    beta.iid <- Ii%*%t(Score)

    onesim <- hatW.MC(x0)
    What <- matrix(onesim$output$Ws,nrow=n);
    ##    W <- cumsum(r[order(x0)]) 
    ##    matplot(x0,What,type="s", col="red", lty=1,pch=-1)
    ##    lines(onesim$output$Wobs ~ x0,type="s",lwd=2)

##    browser()
    
    res <- with(onesim$output,
                list(W=cbind(Wobs), What=list(What),
                     x=cbind(x0),
                     KS=KS, CvM=CvM,
                     R=R, n=n, sd=cbind(Wsd), 
                     cvalues=cbind(cvalues), variable="x",
                     type="sem",
                     model=class(model)[1]) ##, onesim$output$WW)
                )
    class(res) <- "cumres"
    
    return(res)    
  }



################################################################################

