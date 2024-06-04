################################################################
################################################################
################################################################
# ------------ Application: colon cancer data set  ----------- #

library(copula)
library(gamboostLSS)
library(ggplot2)

# load copulas
source("Copulas/Clayton_dependentC.R")
source("Copulas/Gumbel_dependentC.R")
source("Copulas/Gauss_dependentC.R")


# load data
load("colcancer.RData")


ggplot(colcancer, aes(x = followup)) +
  geom_histogram(data = subset(colcancer, death == 0), aes(fill = "No"), binwidth = 5, alpha = 0.7) +
  geom_histogram(data = subset(colcancer, death == 1), aes(fill = "Yes"), binwidth = 5, alpha = 0.7) +
  scale_fill_manual(name = "Event", values = c("No" = "#238443", "Yes" = "#D9F0A3")) +
  theme_bw() + 
  ylab("Count") + 
  xlab("Follow-up time (months)") +
  theme(plot.title = element_text(size=14, hjust = 0.5),
        axis.title.x =  element_text(size = 15),
        axis.title.y =  element_text(size = 15),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 18),
        legend.text = element_text(size = 13), legend.title = element_text(size = 14)) 




# Calculation of the negative log-likelihood based on the left-out folds of
# 10-fold cross-validation for different combinations of copulas and margins
seed = 1234



set.seed(seed)
form <- as.formula(cbind(followup, death) ~ CTX + ASA.score + pUICC + age + LNE + LNR + R.status + preexisting.cancer + sex)


#### ---- Gauss Copula
mod <- glmboostLSS(form, data = colcancer, families = Gauss_Cop_depCens(marg1 = "WEI", marg2 = "WEI", dep = T, mu1=NULL, mu2 = NULL, sigma1 = NULL, sigma2 =  NULL, rho = NULL), 
                   control = boost_control(nu=0.01, mstop=1000, trace=T), method = "noncyclic")


folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds, grid = 1:20000)



#### ---- Gumbel Copula
mod <- glmboostLSS(form, data = colcancer, families = Gumbel_Cop_depCens(marg1 = "WEI", marg2 = "WEI", dep = T, mu1=NULL, mu2 = NULL, sigma1 = NULL, sigma2 =  NULL, rho = NULL), 
                   control = boost_control(nu=0.01, mstop=1000, trace=T), method = "noncyclic")


folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds, grid = 1:20000)



#### ---- Clayton Copula
mod <- glmboostLSS(form, data = colcancer, families = Clayton_Cop_depCens(marg1 = "WEI", marg2 = "WEI", dep = T, mu1=NULL, mu2 = NULL, sigma1 = NULL, sigma2 =  NULL, rho = NULL), 
                   control = boost_control(nu=0.01, mstop=10000, trace=T), method = "noncyclic")


folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds)

MSTOP <- mstop(CV_mod)




###################################################################################
############################### --- Final model --- ###############################

seed = 1234

set.seed(seed)
form <- as.formula(cbind(followup, death) ~ CTX + ASA.score + pUICC + age + LNE + LNR + R.status + preexisting.cancer + sex)


mod <- glmboostLSS(form, data = colcancer, families = Clayton_Cop_depCens(marg1 = "Wei", marg2 = "Wei", dep = T, mu1=NULL, mu2 = NULL, sigma1 = NULL, sigma2 =  NULL, rho = NULL), 
                   control = boost_control(nu=0.01, mstop=10000, trace=T), method = "noncyclic")


folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds)

MSTOP <- mstop(CV_mod)

mod <- glmboostLSS(form, data = colcancer, families = Clayton_Cop_depCens(marg1 = "Wei", marg2 = "Wei", dep = T, mu1=NULL, mu2 = NULL, sigma1 = NULL, sigma2 =  NULL, rho = NULL), 
                   control = boost_control(nu=0.01, mstop=MSTOP, trace=T), method = "noncyclic")


coef_last <- coef(mod, off2int=T)

res < - list(Coefficients = coef_last, mstop.opt = MSTOP, Model = mod)

save(res, file = "Results_colcancer.R")



###################################################################################
################################ --- Cox model --- ################################

set.seed(seed)

mod <- glmboost(form,data = colcancer, control = boost_control(mstop = 1000, nu =0.1), family = CoxPH())

folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds, grid = 1:20000)

MSTOP <- mstop(CV_mod)

mod <- glmboost(form, data = colcancer, family = CoxPH(), control = boost_control(nu=0.1, mstop=MSTOP))

coef_last <- coef(mod, off2int=T)


res_Cox <- list(Coefficients = coef_last, mstop.opt = MSTOP, Model = mod, CV_mod)



###################################################################################
################################ --- AFT model --- ################################


set.seed(seed)

mod <- glmboostLSS(form, data = colcancer, families = Custom_WeibullFamily(), 
                   control = boost_control(nu=0.01, mstop=10000, trace=T), method = "noncyclic")


folds <- cv(model.weights(mod), type = "kfold", B = 10)
CV_mod <- cvrisk(mod, folds = folds)

MSTOP <- mstop(CV_mod)

mod <- glmboostLSS(form, data = colcancer, families = Custom_WeibullFamily(), 
                   control = boost_control(nu=0.01, mstop=MSTOP, trace=T), method = "noncyclic")


coef_last <- coef(mod, off2int=T)

res_AFT <- list(Coefficients = coef_last, mstop.opt = MSTOP, Model = mod)






###################################################################################
pdffz <- function(input){
  
  pdiff <- ifelse(input > 0.999999, 0.999999, input) 
  pdiff <- ifelse(pdiff < 1e-16, 1e-16, pdiff) 
  
  return(pdiff)
  
}


dvffz <- function(derivative){
  
  deriv <- ifelse(is.na(derivative), .Machine$double.eps, derivative)
  deriv <- ifelse(deriv == Inf, 8.218407e+20, deriv)
  deriv <- ifelse(deriv == -Inf, -8.218407e+20, deriv)
  
  return(deriv)
}




### other stuff because these did not load for me... 
check_stabilization <- function(stabilization = c("none", "MAD", "L2")) {
  stabilization <- match.arg(stabilization)
  ## check if old stabilization interface is used and issue a warning
  if (getOption("gamboostLSS_stab_ngrad")) {
    warning("Usage of ", sQuote("options(gamboostLSS_stab_ngrad = TRUE)"),
            " is deprecated.\n", "Use argument ", sQuote("stabilization"),
            " in the fitting family. See ?Families for details.")
    if (stabilization == "none")
      warning(sQuote("stabilization"), " is set to ", dQuote("MAD"))
  }
  stabilization
}


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


### Create our own weibull and log-logistic margins for gamboostLSS: 
Custom_WeibullMu <- function(mu = NULL, sigma = NULL, stabilization){
  
  # y_check <- function(y){
  #   if ((is.matrix(y) && NCOL(y)!=2))
  #     stop("response should be a two-column matrix (y1 and y2) for this bivariate family")
  # }
  
  # neg. log-likelihood
  loss <- function(sigma, y, f=f){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(sigma = sigma, y = y, f = f))
  }
  
  # negative gradient
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = param, shape = sigma)
    
    SurvT <- pweibull(q = time, scale = param, shape = sigma, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dermu <- (- 1 / ( param ) + (sigma - 1) * (- 1 / (param)) + (sigma) * (time)^(sigma)*(param)^(- (sigma) - 1) ) 
    
    derS_dermu <- - ( -(sigma * time * exp(- (time / (param) )^( sigma ) ) * ( time/ (param) )^(sigma - 1) * (1 / (param)^2 ) ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dermu * exp(f) ) + (1-censind) * ( (1/SurvT) * derS_dermu * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  
  offset <- function(y, w){
    
    if (!is.null( mu )) {
      temp1 <- mu
      temp1 <- log( temp1 )
      
      RET <- temp1
    }else {
      
      
      RET <- log(weighted.mean((y[,1] + weighted.mean(y[,1], w = w, na.rm = TRUE))/2, w = w, na.rm = TRUE)) 
      
    }
    return(RET)
  }
  
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Weibull distribution: mu (log link)")
}

Custom_WeibullSigma <- function(mu = NULL, sigma = NULL, stabilization){
  
  
  # neg. log-likelihood
  loss <- function(mu, y, f=f){
    
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    lik <- censind * log(pdfT) + (1 - censind) * log(SurvT)
    
    neglik <- - lik
    
    return(neglik)
    
    
  }
  
  # risk is sum of loss
  risk <- function(y, f, w = 1) {
    sum(w * loss(mu = mu, y = y, f = f))
  }
  
  
  # ngradient 
  ngradient <- function(y, f, w = 1){
    
    time <- y[,1]
    censind <- y[,2]
    
    # in here this is MU or vartheta1
    param <- exp(f)
    
    pdfT <- dweibull(x = time, scale = mu, shape = param)
    
    SurvT <- pweibull(q = time, scale = mu, shape = param, lower.tail = FALSE)
    
    
    ## Small check
    pdfT <- dvffz(pdfT)
    
    SurvT <- pdffz(SurvT)
    
    
    derlpdf_dersigma <- ( (1/(param)) + log(time) - log(mu) - (time/ (mu))^(param) * log(time / (mu) )  )
    
    derS_dersigma <- - ( ( exp(- (time / (mu) )^(param) ) * ( (time / (mu) )^(param) * log(time / (mu) ) ) ) )
    
    
    #### negative gradient:
    ngr <- censind * ( derlpdf_dersigma * exp(f) ) + (1 - censind) * ( ( 1/SurvT ) * derS_dersigma * exp(f) )
    
    ngr <- stabilize_ngradient(ngr, w = w, stabilization)
    
    return(ngr)
    
    
  }
  
  offset <- function(y, w){
    
    if (!is.null( sigma )) {
      
      
      temp2 <- log( sigma )
      
      RET <- temp2
      
    }else{
      
      sigma_temp <- rep(0.1, length(y[,1]))
      
      RET <- log(mean(sigma_temp))
      
    }
    return(RET)
  }
  
  mboost:::Family(ngradient = ngradient, risk = risk, loss = loss,
                  response = function(f) exp(f), 
                  offset = offset,
                  name = "Weibull distribution: sigma (log link)")
}

#------------------ complete gamboostLSS families
Custom_WeibullFamily <- function (mu = NULL,  sigma = NULL, stabilization = c("none", "MAD", "L2")){
  
  stabilization <- check_stabilization(stabilization)
  
  Families(   mu = Custom_WeibullMu(mu = mu, sigma = sigma,  stabilization = stabilization ), 
              sigma = Custom_WeibullSigma(mu = mu, sigma = sigma, stabilization = stabilization ), 
              name = "Weibull distribution")
}


