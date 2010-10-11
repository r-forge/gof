`cumres.default` <-
  function(model,score,information,
           residual,variable,
           data=model.frame(model),p=pars(model),
           R=500, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
           debug=FALSE, ...) {

        
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
                   seed=as.integer(seed), ## Seed (will probably rely on R's rangen in future version)
                   KS=as.double(0), ## Return: Kolmogorov Smirnov statistics for each realization
                   CvM=as.double(0), ## Return: Cramer von Mises statistics for each realization
                   Wsd=as.double(numeric(n)), ## Return: Pointwise variance of W(x)
                   cvalues=as.double(numeric(R)), ## Return: value for each realization s.t.  +/- cvalue * Wsd contains W*
                   Ws=as.double(numeric(plots*n)), ## Return: Saved realizations (for plotting function)
                   Wobs=as.double(numeric(n)) ## Observed process
                   , PACKAGE="gof")
      return(list(output=output,x=x[ord]))
    }


  }
