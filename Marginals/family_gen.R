################################################################################
### DESCRIPTION
### This is a helper function that creates for the different
### parameter sub-models, i.e. the parameters of the marginals and the 
### copula dependence parameter, the respective loss, risk and gradient functions.
### Finally, it creates and returns a parameter specific mboost Family objects .


### libraries 
library(numDeriv)



family_gen <- function(mu1 = NULL, sigma1 = NULL, nu1 = NULL, tau1 = NULL,
                       mu2 = NULL, sigma2 = NULL, nu2 = NULL, tau2 = NULL, 
                       rho = NULL,
                       stabilization,
                       loss_body, loss_args, 
                       risk_body, risk_args,
                       grad_body, grad_args,
                       offset_body, offset_args,
                       check_y_body, check_y_args,
                       response, name){
  
  
  
################################################################################
######################### Helpers from gamboostLSS #############################
################### To be deleted when integration into package ################
  
  ## function for weighted sd
  weighted.sd <- function(x, w, ...) {
    if (missing(w))
      w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
  }
  
  ## weighted median
  weighted.median <- function (x, w = 1, na.rm = FALSE) {
    if (length(w) == 1)
      w <- rep(w, length(x))
    
    ## remove observations with zero weights
    x <- x[w != 0]
    w <- w[w != 0]
    
    ## remove NAs if na.rm = TRUE
    if (na.rm) {
      keep <- !is.na(x) & !is.na(w)
      x <- x[keep]
      w <- w[keep]
    } else {
      if (any(is.na(x)) | any(is.na(w)))
        return(NA)
    }
    
    ## sort data and weights
    ind <- order(x)
    x <- x[ind]
    w <- w[ind]
    
    ## first time that fraction of weights is above 0.5
    ind1 <- min(which(cumsum(w)/sum(w) > 0.5))
    
    ## first time that fraction of weights is below 0.5
    ind2 <- ifelse(ind1 == 1, 1, max(which(cumsum(w)/sum(w) <= 0.5)))
    
    ## if sum of weights is an even integer
    if(sum(w) %% 1 == 0 && sum(w) %% 2 == 0)
      return(mean(c(x[ind1], x[ind2])))
    
    ## else return
    return(max(c(x[ind1], x[ind2])))
  }
  
  ## helper function that stabilizes the negative gradient if requested by the user
  stabilize_ngradient <- function(ngr, w = 1, stabilization) {
    ## set which to MAD if gamboostLSS_stab_ngrad = TRUE and which == "none"
    if (stabilization == "none" && getOption("gamboostLSS_stab_ngrad"))
      stabilization <- "MAD"
    ## stabilization using the mean absolute deviation (MAD)
    if (stabilization == "MAD") {
      div <- weighted.median(abs(ngr - weighted.median(ngr, w = w, na.rm = TRUE)),
                             w = w, na.rm = TRUE)
      div <- ifelse(div < 0.0001, 0.0001, div)
      ngr <- ngr / div
    }
    if (stabilization == "L2") {
      div <- sqrt(weighted.mean(ngr^2, w =w,  na.rm = TRUE))
      div <- ifelse(div < 1e-04, 1e-04, div)
      div <- ifelse(div > 1e+04, 1e+04, div)
      ngr <- ngr / div
    }
    ngr
  }
  
  
  ### pdf Singh-Maddala
  
  dSinghMaddala <- function (x, mu = 1, sigma = 1, tau = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", 
                 ""))
    
    z <- (x/mu)^sigma
    
    loglik <- log(z) + log(abs(sigma)) - log(x) + log(tau) - (1 + tau) * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    ft
  }
  
  ### cdf Singh-Maddala
  
  pSinghMaddala <- function (x, mu = 1, sigma = 1, tau = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(tau < 0)) 
      stop(paste("tau must be positive", "\n", 
                 ""))
    
    z <- (x/mu)^sigma
    
    log_cdf <- log(1 - (1 + z)^(-tau))
    
    if (log == FALSE) 
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  # pdf Dagum
  
  dDagum <- function (x, mu = 1, sigma = 1, nu = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(nu < 0)) 
      stop(paste("nu must be positive", "\n", ""))
    
    z <- (x/mu)^sigma
    loglik <- nu * log(z) + log(abs(sigma)) - log(x) - lgamma(nu) - 
      lgamma(1) + lgamma(nu + 1) - (nu + 1) * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    
    ft
  }
  
  ## cdf Dagum
  
  pDagum <- function (x, mu = 1, sigma = 1, nu = 1, log = FALSE){
    
    if (any(mu < 0))
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(nu < 0))
      stop(paste("nu must be positive", "\n",
                 ""))
    
    log_cdf <- log((1 + (x/mu)^-sigma)^-nu)
    
    if (log == FALSE)
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  
  ### pdf LogLogistic
  
  dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    z <- (x/mu)^(sigma)
    
    loglik <- log(sigma) + (sigma) * (log(x) - log(mu)) - log(x) - 2 * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    ft
  }
  
  ### cdf LogLogistic
  
  pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    log_cdf <- log(1/(1 + (x/mu)^(-sigma)))
    
    if (log == FALSE) 
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  weighted.sd <- function(x, w, ...) {
    if (missing(w))
      w <- rep(1, length(x))
    m <- weighted.mean(x, w, ...)
    var <- weighted.mean((x - m)^2, w, ...) * sum(w) / (sum(w) - 1)
    return(sqrt(var))
  }
  
  
  
################################################################################
######## Helper for numerical approximation of the gradient from gjrm ##########
################################################################################
  
  pdffz <- function(input){
    
    pdiff <- ifelse(input == 1, 0.9999999999999999, input) 
    pdiff <- ifelse(pdiff < 1e-16, 1e-16, pdiff) 
    
    return(pdiff)
    
  }
  
  dvffz <- function(derivative){
    
    deriv <- ifelse(is.na(derivative), .Machine$double.eps, derivative)
    deriv <- ifelse(deriv == Inf, 8.218407e+20, deriv)
    deriv <- ifelse(deriv == -Inf, -8.218407e+20, deriv)
    
    return(deriv)
  }
  
  # helper for numerical approximation from gjrm
  
  num_grad <- function(funcD, para){
    
    para <- c(para)
    
    if(length(para) == 1){
      
      fi <- jacobian(function(para) funcD(para), para)
      
    }
    
    if(length(para) > 1){
      
      fi <- grad(function(para) funcD(para), para)
      
    }
    
    return(fi)
    
    
  }
  
  check_eta_theta <- function(additivepredictor, coptype){
    
    if(coptype == "GAUSS"){
      
      checked_predictor <- ifelse(abs(additivepredictor) > 8.75, sign(additivepredictor) * 8.75, additivepredictor)
      
    }
    
    if(coptype == "FRANK"){
      
      epsilon <- c(sqrt(.Machine$double.eps))
      
      checked_predictor <- ifelse(abs(additivepredictor) < epsilon, epsilon, additivepredictor)
      
    }
    
    if(coptype %in% c("GUMBEL", "CLAYTON", "JOE")){
      
      checked_predictor <- ifelse(additivepredictor > 20, 20, additivepredictor)
      
      checked_predictor <- ifelse(additivepredictor < -17, -17, additivepredictor)
    }
    
    return(checked_predictor)
  }
  
  
  # GBS: Credit where it is due: Another incredibly useful piece of software from
  # gjrm:
  # 
  aux_matrix_disc <- function(response){
    
    ly1 <- length(response)
    y1m <- list()
    my1 <- max(response)
    
    for(i in 1:ly1){
      
      y1m[[i]] <- seq(0, response[i]);
      
      length(y1m[[i]]) <- my1 + 1
      
    }
    
    
    y1m <- do.call(rbind, y1m)
    
    
    return(y1m)
  }
  
  
################################################################################
######################### Function definition ##################################
################################################################################
  
  
  # defining the loss function
  loss <- function(){}
  body(loss) <- loss_body
  formals(loss) <- loss_args
  
  
  # defining the risk function
  risk <- function(){}   
  body(risk) <- risk_body
  formals(risk) <- risk_args
  
  
  # defining the gradient function
  ngradient <- function(){}
  body(ngradient) <- as.call(c(as.name("{"), grad_body)) # creating a multi-expression body (see: https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/body)
  formals(ngradient) <- grad_args
  
  # defining offset function
  offset <- function(){}
  body(offset) <- offset_body
  formals(offset) <- offset_args
  
  # defining the check_y function
  
  check_y <- function(){}
  body(check_y) <- check_y_body
  formals(check_y) <- check_y_args
  
  # setting the appropriate environment for offset, y_check and response function
  #environment(offset) <- current_env()
  environment(response) <- current_env()
  #environment(check_y) <- current_env()
  
  
  # creating the Family object
  mboost::Family(ngradient = ngradient,
                 risk = risk,
                 loss = loss,
                 response = response,
                 offset = offset,
                 name = name,
                 check_y = check_y)
  
  
}