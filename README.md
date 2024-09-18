# Boosting Dependent Censoring
A model-based boosting approach to deal with dependent censoring.


## Example
```
library(copula)
library(gamboostLSS)

source("Copulas/Gauss_dependentC.R")

p = 10
n = 2000

set.seed(seed+3)

# data generation
X_train <- matrix(NA, ncol = p, nrow = n)
for(i in 1:p)  X_train[,i] <- runif(n, -1, 1)
colnames(X_train) <- paste0("X", 1:p)

eta_mu1    <- 0.7 + 2 * X_train[,1] + 1 * X_train[,3]
eta_mu2    <-  rep(0.8, times = length( X_train[,2] ))  
mu1     <- exp(eta_mu1)
mu2     <- exp(eta_mu2)

eta_sigma1 <-  0.7  + 0.7 * X_train[,3]
eta_sigma2 <-   rep(0.5,times=length(X_train[,2]) ) 
sigma1  <- exp(eta_sigma1)
sigma2  <- exp(eta_sigma2)

eta_theta  <-2 + 1.5 * X_train[,5]
rho <- tanh(eta_theta)
  
y1 <- vector(length = length(mu1))
y2 <- vector(length = length(mu2))
for(i in 1:length(mu1)){
  cop = ellipCopula(family = "normal", dim = 2, param = rho[i])
  myMvd <- mvdc(copula = cop, margins = c("weibull", "weibull"), paramMargins = list(list( scale = (mu1[i]), shape = (sigma1[i])), 
                                                                                     list( scale = (mu2[i]), shape = (sigma2[i]))))
  vect <- rMvdc(1, myMvd)
  y1[i] <- vect[,1]
  y2[i] <- vect[,2]
}

Y <- pmin(y1,y2)
d <- as.numeric(Y == y1)
train <-  data.frame(Y,d, X_train)

# fit model via noncyclic gamboostLSS
mod <- glmboostLSS(cbind(Y,d)~., data= train, families = Gauss_Cop_depCens(marg1 = "WEI", marg2 = "WEI"), 
                   control = boost_control(nu=0.01, mstop=2000, risk = "oobag", trace = T
                   ), method = "noncyclic")


```
