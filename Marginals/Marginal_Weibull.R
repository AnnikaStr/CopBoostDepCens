################################################################################
### DESCRIPTION
### This file defines the Log-Logistic-Marginal function. It contains the 
### expressions of the pdf, cdf, its derivatives wrt the two parameters, names,
### response functions, offset functions, check_y functions required to create an mboost Family object.
### This function is handed over to a copula function and the expressions and functions
### are merged with the respective copula via the rlang package, in order to create 
### appropriate Families objects from gamboostLSS.


### libraries 
library(mboost)
library(gamboostLSS)
library(rlang)


### Weibull function

Weibull_Mar <- function(loc = NULL, offset_mu = NULL, offset_sigma = NULL, dep = FALSE){
  
  
  # check for appropriate offset values for parameter mu, sigma and tau
  
  if ((!is.null(offset_mu) && offset_mu <= 0))
    stop(sQuote("mu"), paste(" must be greater than zero in Marginal", loc))
  
  if ((!is.null(offset_sigma) && offset_sigma <= 0))                                          
    stop(sQuote("sigma"), paste(" must be greater than zero in Marginal", loc))
  
  
  
  
  # creating the location specific expressions
  y <- parse_expr(paste("y[,", loc, "]", sep = ""))                      
  mu <- parse_expr(paste("mu", loc, sep = ""))                           
  sigma <- parse_expr(paste("sigma", loc, sep = ""))                     
  
  if(dep)  y <- parse_expr(paste("y[,", 1, "]", sep = ""))
  
  ### generic functions
  
  # pdf
  
  pdf_gen <- expr(dweibull(x = !!y, scale = !!mu, shape = !!sigma))#expr(dLogLogistic(x = !!y, mu = !!mu, sigma = !!sigma))
  # cdf
  cdf_gen <- expr(pweibull(q = !!y, scale = !!mu, shape = !!sigma))#expr(pLogLogistic(x = !!y , mu = !!mu, sigma = !!sigma))
  
  
  l_generic <- list(pdf = pdf_gen,
                    cdf = cdf_gen)
  
  
  
  ### mu functions
  
  # pdf
  pdf_mu  <- expr(dweibull(x = !!y, scale = exp(f), shape = !!sigma))#expr(dLogLogistic(x = !!y, mu = exp(f), sigma = !!sigma))
  
  # derivativ logpdf
  derlpdf1.deretamu   <- expr( (- 1 / (exp(f)) + (!!sigma - 1) * (- 1 / (exp(f))) + (!!sigma) * (!!y)^(!!sigma)*(exp(f))^(- (!!sigma) - 1) ) * exp(f) )
  #expr(( -(!!sigma)/exp(f) + 2 * (!!sigma) * (!!y / exp(f))^(!!sigma) / (exp(f) * (1 + (!!y / exp(f))^(!!sigma))) ) * exp(f) )
  
  # cdf
  cdf_mu  <- expr(pweibull(q = !!y, scale = exp(f), shape = !!sigma))#expr(pLogLogistic(x = !!y, mu = exp(f), sigma = !!sigma))
  
  # derivative cdf
  dercdf.deretamu  <- expr(-(!!sigma * !!y * exp(- (!!y / (exp(f)) )^(!!sigma) ) * ( !!y/ (exp(f)) )^(!!sigma - 1) * (1 / (exp(f))^2 ) ) * exp(f) )
  #expr(- ( (!!sigma)/exp(f) * ((!!y)/exp(f))^(-(!!sigma)) / (1 + ((!!y)/exp(f))^(-(!!sigma)))^2 ) * exp(f)) 
  
  
  l_mu <- list(pdf = pdf_mu,
               derlpdf = derlpdf1.deretamu,
               cdf = cdf_mu,
               dercdf = dercdf.deretamu)
  
  
  
  ### sigma functions
  
  # pdf
  pdf_sig <- expr(dweibull(x = !!y, scale = !!mu, shape = exp(f)))#expr(dLogLogistic(x = !!y, mu = !!mu, sigma = exp(f)))
  # derivative pdf
  derlpdf1.deretasig  <- expr( ( (1/(exp(f))) + log(!!y) - log(!!mu) - (!!y/ (!!mu))^(exp(f)) * log(!!y / (!!mu) )  ) * exp(f) )
  #expr(( (1/exp(f)) + log(!!y) - log(!!mu) - 2/(1 + ((!!y) / (!!mu))^(exp(f))) * ((!!y) / (!!mu))^(exp(f)) * log((!!y)/(!!mu)) ) * exp(f))
  
  # cdf 
  cdf_sig <- expr(pweibull(q = !!y, scale = !!mu, shape = exp(f)))#expr(pLogLogistic( x = !!y, mu = !!mu, sigma = exp(f)))
  
  # derivative cdf 
  dercdf.deretasig <- expr( ( exp(- (!!y / (!!mu) )^(exp(f)) ) * ( (!!y / (!!mu) )^(exp(f)) * log(!!y / (!!mu) ) ) ) * exp(f) )
  #expr(( ((!!y) / (!!mu))^(-(exp(f))) * log((!!y)/(!!mu)) / (1 + ((!!y) / (!!mu))^(-(exp(f))))^2  ) * exp(f))
  
  
  l_sigma <- list(pdf = pdf_sig,
                  derlpdf = derlpdf1.deretasig,
                  cdf = cdf_sig,
                  dercdf = dercdf.deretasig)
  
  
  
  
  
  ### response functions 
  
  response_mu <- function(f) exp(f)
  response_sigma <- function(f) exp(f)
  
  
  l_response <- list(mu = response_mu,
                     sigma = response_sigma)
  
  ### offset functions     
  
  offset_mu <- expr({
    if (!is.null(!!mu)) {
      RET <- log(!!mu)
    }
    else {
      RET <- log(weighted.mean((!!y + weighted.mean(!!y, w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) # weighted version
    }
    return(RET)
  })
  
  offset_sigma <- expr({            # taken from gamlss
    if (!is.null(!!sigma)) {
      RET <- log(!!sigma)
    }
    else {
      sigma <- rep(0.1, length(!!y))             
      RET <- log(mean(sigma))
    }
    return(RET)
  })
  
  
  
  
  l_offset <- list(mu = offset_mu,
                   sigma = offset_sigma)
  
  ### names 
  
  name_mu <- "Weibull distribution: mu(log link)"
  name_sigma <- "Weibull distribution: sigma(log link)"
  
  l_names <- list(mu = name_mu,
                  sigma = name_sigma)
  
  ### check y function   
  
  
  check_y_mu <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
    y
  })
  
  check_y_sigma <- expr({
    if (!is.numeric(!!y))
      stop("response is not numeric but ", sQuote("WeibullMar()"))
    if (any(!!y < 0))
      stop("response is not positive but ", sQuote("WeibullMar()"))
    y
  })
  
  
  
  l_check_y <- list(mu = check_y_mu,
                    sigma = check_y_sigma)   
  
  
  
  
  # return list
  
  l = list(parameters = c(mu, sigma),
           parameter_names = c("mu", "sigma"),
           
           generic = l_generic,
           mu = l_mu,
           sigma = l_sigma,
           
           response = l_response,
           offset = l_offset,
           name = l_names,
           check_y = l_check_y,
           marg_name = "WeibullMarg")
  
  
  return(l)
  
}


