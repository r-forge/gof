### `cumres` <-
###   function(model,...) UseMethod("cumres")

`Debug` <-
  function(msg, cond) {
    if (cond)
      print(paste(msg, collapse=" "))
  }

`cumres.lvmfit` <-
  function(model,variable=exogenous(model),index=with(model,c(endogenous(model),latent(model))),
           type="residual", data=model.frame(model),p=pars(model),
           R=500, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
           debug=FALSE, ...) {
    if (!require("gof")) stop("package 'gof' not available")
    if (!require("numDeriv")) stop("package 'numDeriv' not available")
    
    Debug(list("p=",p),debug)
  
##    Ii <- solve(information(model,p), num=TRUE)
    Ii <- model$vcov

    browser()
    Debug(list("Ii=",Ii),debug)
    Score <- score(model,data=data,p=p)
    beta.iid <- Ii%*%t(Score)
    
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

    n <- nrow(data)
    P <- predict(model)
    lat.var <- latent(model)
    lat.didx <- match(lat.var,names(data))
    y.var <- endogenous(model)
    if (!is.character(index)) 
      index <- names(data)[index]
    if (length(lat.var)<1)
      type <- "residuals"

    yp <- t(attributes(P)$Ey.x)
    etap <- t(attributes(P)$eta.x)
    Debug("Creating prediction data.frame (... conditional given x)", debug)
    pred <- cbind(data[,exogenous(model),drop=FALSE], yp, etap)
        
    W <- c()
    What <- c()
    Wsd <- c()  
    KS <- c()
    CvM <- c()
    mytype <- c()
    UsedVars <- c()
    UsedData <- c()
    allcvalues <- c();

    for (v in variable)
      for (ii in index)
##        if (!(tt=="residual" & ii%in%lat.var))
          {
            cat(ii, v, "\n")
###             if (tt=="residual") {
###               idx <- match(ii,y.var)
###               r <- residuals(model,data=data,p=p)[,ii]
###               grad <- -gradpredict(p,model,data)[((idx-1)*n+1):(idx*(n)),]
###             } else {
            {
              myres <- function(p) {
                attributes(predict(model,p=p))$epsilon.y[,ii]
              }
              r <- myres(p)
              grad <- -jacobian(myres,p)
            }
            x <- pred[,v]
            Debug("About to simulate...", debug)
            if (!is.null(x)) {
              UsedVars <- c(UsedVars, v)
              onesim <- hatW.MC(x)
              UsedData <- cbind(UsedData, onesim$x);  
              W <- cbind(W, onesim$output$Wobs)
              Wsd <- cbind(Wsd, onesim$output$Wsd)
              What <- c(What, list(matrix(onesim$output$Ws,nrow=n)));
              KS <- c(KS, onesim$output$KS);  CvM <- c(CvM, onesim$output$CvM)
              allcvalues <- cbind(allcvalues, onesim$output$cvalues)
              mytype <- c(mytype, paste(ii,sep=","))
            } else cat("Variable '", v , "' not found.\n",sep="")
          }
    
    if (length(UsedVars)<1)
      return(NULL)
    
    res <- list(W=W, What=What,
                x=UsedData,
                KS=KS, CvM=CvM,
                R=R, n=nrow(UsedData), sd=Wsd, 
                cvalues=allcvalues, variable=UsedVars,
                type=mytype,
                model=class(model)[1]) ##, onesim$output$WW)
    class(res) <- "cumres"
    
    return(res)    
  }


###{{{ R prototype

`Diag` <-
  function(model,variable,idx=1,data=model.frame(model),p=pars(model),
           R=500, b=0, plots=min(R,50), seed=round(runif(1,1,1e9)),
           debug=FALSE, ...) {

    n <- nrow(data)
    if (is.character(variable)) {
      x <- data[,variable]
    } else {
      x <- variable
      variable <- names(x)
    }

    r <- residuals(model,data=data,p=p)[idx,]
    grad <- -gradpredict(p,model,data)[((idx-1)*n+1):(idx*(n)), ]

    Debug(r,debug)
    
    ord <- order(x)
    ones <- matrix(0,n,n)
    for (i in 1:n)
      ones[i,] <- (x<=x[i])*1

    rs <- ones;
    for (i in 1:nrow(rs))
      rs[i,rs[i,]==1] <- r[rs[i,]==1]

    W <- 1/sqrt(n)*ones%*%r
    
    eta <- ones%*%grad
    Debug(list("p=",p),debug)
    
    Ii <- solve(information(model,p))
    Debug(list("Ii=",Ii),debug)
    S <- score(model,data=data,p=p)
    H0 <- eta%*%Ii%*%t(S)
    H1 <- (rs - H0)

    KS <- function(x) max(abs(x))
    CvM <- function(x) NA
    
    KS.obs <- KS(W)
    CvM.obs <- CvM(W)
    
    KSlist <- CvMlist <- cvalues <- c()
    last <- c(); for (i in 1:R) {
      N <- rnorm(n);
      W. <- 1/sqrt(n)*H1%*%N
      KSlist <- c(KSlist, ifelse( (KS(W.)>KS.obs), TRUE, FALSE))
      CvMlist <- c(CvMlist, ifelse( (CvM(W.)>CvM.obs), TRUE, FALSE))
      cvalues <- c(cvalues, 3)
      if (i<=plots)
        last <- cbind(last, W.)
    }

    W <- W[ord,]; What <- last[ord,]; x <- x[ord]; Wsd <- sd(t(H1))[ord]
    res <- list(W=cbind(W), What=list(What), x=cbind(x), KS=sum(KSlist)/R, CvM=sum(CvMlist)/R,
                R=R, n=n, sd=cbind(Wsd), cvalues=cbind(cvalues), variable=variable, type="residual")
    class(res) <- "cumres"
    return(res)
  }

###}}} R prototype
