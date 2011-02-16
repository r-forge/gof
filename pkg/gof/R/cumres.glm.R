`cumres` <-
function(model,...) UseMethod("cumres")

`cumres.lm` <- function(model,...) {
  cumres.glm(model,...)
}

`cumres.glm` <- function(model,
         variable=c("predicted",colnames(model.matrix(model))),
         data=data.frame(model.matrix(model)), 
         R=1000, b=0, plots=min(R,50),breakties=1e-12,
         seed=round(runif(1,1,1e9)),...) {
  dginv <- function(z) 1
  a.phi <- 1
  switch(class(model)[1],
         lm = {
           a.phi <- summary(model)$sigma^2
           h <- function(z) 1
         },
         glm = {
           f <- family(model)
           if (f$family=="gaussian") {
             a.phi <- summary(model)$dispersion
           }
           ## stop("Refit model using 'lm'")           

           g <- f$linkfun
           ginv <- f$linkinv
           dginv <- f$mu.eta ## D[linkinv]
           ##dg <- function(x) 1/dginv(g(x)) ## Dh^-1 = 1/(h'(h^-1(x)))
           canonf <- do.call(f$family,list())           
           caninvlink <- canonf$linkinv
           canlink <- canonf$linkfun
           Dcaninvlink <- canonf$mu.eta           
           Dcanlink <- function(x) 1/Dcaninvlink(canlink(x))
           ##gmu <- function(x) g(caninvlink(x))
           ##invgmu <- function(z) canlink(ginv(z))
           h <- function(z) Dcanlink(ginv(z))*dginv(z)                                
         },
         stop("Unsupported model!"))

  response <- all.vars(formula(model))[1]
  X <- model.matrix(model)
  n <- nrow(X)  
  r <- residuals(model, type="response") ## y-g^{-1}(Xb)
  yhat <- predict(model, type="response") ## g^{-1}(Xb)
  ##  Xbeta <- predict(model, type="link") ## X*b
  beta <- coef(model)
  if(any(is.na(beta))) stop("Over-parametrized model")
  Xbeta <- X%*%beta
 
  etaraw <- (as.numeric(dginv(Xbeta))*X)
  hatW.MC <- function(x) {
    myorder <- order(x)
    x <- x[myorder]
    Ii <- vcov(model)
    ##    S <- (X*(as.vector(h(Xbeta))*r))/a.phi
    A <- as.vector(h(Xbeta)*r)/a.phi 
    S <- apply(X,2,function(x) x*A)
    beta.iid <- Ii%*%t(S[myorder,,drop=FALSE])
    r0 <- r[myorder]
    D0 <- etaraw[myorder,,drop=FALSE]
    Wfun <- "W2"
    if (b!=0) { Wfun <- "W" }
    output <- .C(Wfun,
                 R=as.integer(R), ## Number of realizations
                 b=as.double(b), ## Moving average parameter
                 n=as.integer(n), ## Number of observations
                 npar=as.integer(nrow(Ii)),
                 xdata=as.double(x), ## Observations to cumulate after 
                 rdata=as.double(r0), ## Residuals (one-dimensional)
                 betaiiddata=as.double(beta.iid), ## Score-process
                 etarawdata=as.double(D0), ## Eta (derivative of terms in cummulated process W)
                 plotnum=as.integer(plots), ## Number of realizations to save (for later plot)
                 seed=as.integer(seed), ## Seed (will probably rely on R's rangen in future version)
                 KS=as.double(0), ## Return: Kolmogorov Smirnov statistics for each realization
                 CvM=as.double(0), ## Return: Cramer von Mises statistics for each realization
                 Wsd=as.double(numeric(n)), ## Return: Pointwise variance of W(x)
                 cvalues=as.double(numeric(R)), ## Return: value for each realization s.t.  +/- cvalue * Wsd contains W*
                 Ws=as.double(numeric(plots*n)), ## Return: Saved realizations (for plotting function)
                 Wobs=as.double(numeric(n)) ## Return: Observed process
                 )
    ##    W <- function() { 1/sqrt(n)*cumsum(r[myorder]) }
    return(list(output=output,x=x))
  }
  
  if (!is.na(match(response, variable))) variable[match(response, variable)] <- "predicted"
  variable <- unique(variable)
  UsedData <- X[,na.omit(match(variable, colnames(X))),drop=FALSE]
  myvars <- colnames(UsedData)[apply(UsedData,2,function(x) length(unique(x))>2)] ## Only consider variables with more than two levels
  if ("predicted"%in%variable) myvars <- c("predicted",myvars)

  
  untie <- runif(n,0,breakties)
  
  W <- c()
  What <- c()
  Wsd <- c()  
  KS <- c()
  CvM <- c()
  mytype <- c()
  UsedVars <- c()
  UsedData <- c()
  allcvalues <- c();
  for (v in myvars) {
    x <- NULL
    if (v=="predicted") {
      x <- yhat ## Xbeta
    } else if (v %in% colnames(X)) {
      x <- X[,v]       
    }
    if (!is.null(x)) {
      UsedVars <- c(UsedVars, v)
      onesim <- hatW.MC(x+untie) ## obs: break ties
      UsedData <- cbind(UsedData, onesim$x);  
      W <- cbind(W, onesim$output$Wobs)
      Wsd <- cbind(Wsd, onesim$output$Wsd)
      What <- c(What, list(matrix(onesim$output$Ws,nrow=n)));
      KS <- c(KS, onesim$output$KS);  CvM <- c(CvM, onesim$output$CvM)
      allcvalues <- cbind(allcvalues, onesim$output$cvalues)
      mytype <- c(mytype,"residual")
    } else cat("Variable '", v , "' not found.\n",sep="")
  }

  if (length(UsedVars)<1) 
    return(NULL)
  
  res <- list(W=W, What=What,
              x=UsedData,
              KS=KS, CvM=CvM,
              R=R, n=nrow(UsedData), sd=Wsd, 
              cvalues=allcvalues, variable=UsedVars,
              type=mytype, untie=untie,
              model=class(model)[1]) ##, onesim$output$WW)
  class(res) <- "cumres"
  res
}

