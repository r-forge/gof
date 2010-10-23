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
    if (missing(p))
      p <- lava::pars(model)

    myres <- function(p) {
      attributes(predict(model,p=p))$epsilon.y[,y]
    }
    grad <- -numDeriv::jacobian(myres,p,method=lava::lava.options()$Dmethod)
    r <- myres(p)
    ##    Ii <- solve(information(model,p), num=TRUE)
    Ii <- model$vcov
    ord <- order(x0)

    x0 <- x0[ord]
    grad <- grad[ord,]
    r <- r[ord]
    Ii <- model$vcov
    Score <- lava::score(model,data=data,p=p,indiv=TRUE,weight=lava::Weight(model))[ord,]
    beta.iid <- Ii%*%t(Score)

    onesim <- hatW.MC(x0)
    What <- matrix(onesim$output$Ws,nrow=n);
    
    res <- with(onesim$output,
                list(W=cbind(Wobs), What=list(What),
                     x=cbind(x0),
                     KS=KS, CvM=CvM,
                     R=R, n=n, sd=cbind(Wsd), 
                     cvalues=cbind(cvalues), variable="x",
                     type="sem",
                     model=class(model)[1])
                )
    class(res) <- "cumres"
    
    return(res)    
  }



################################################################################

