library(VARgrowth)
library(splines)
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(nlme)
library(rstan)
library(bayesplot)
library(brms)

U <- 200 ### number of years
tpt <- 50 ### tpt per year
mr_obs_pct <- 0.25 ### percent last year observed
d <- 3
dat <- data.frame(myTime = rep(1:tpt, times = U),
                  myGroup = rep(factor(1:U), each = tpt))

trueObsVar <- 100^2
sigma_theta <- c(0.0001, 0.0001, 0.0001)^2
transThetaTrue <- thetaTrue <- transTheta <- theta <- matrix(nrow = U, ncol = d)

XTheta <- matrix(1, nrow = U)
beta <- apply(transTheta, 2, 
              function(col, XTheta){ qr.solve(XTheta, col)}, XTheta = XTheta)

# period <- 1.8*pi/U
# theta[,1] <- 1000*(sin(period*(1:U + 0)) + 2)
# theta[,2] <- 2*(sin(period*(1:U) + 3) + 4)
# theta[,3] <- 0.3*(sin(period*(1:U) + 2.8) + 2)

theta[,1] <- rep(20000, U)
theta[,2] <- rep(2, U)
theta[,3] <- rep(0.7, U)


transTheta[,1] <- log(theta[,1]) 
transTheta[,2] <- log(theta[,2])
transTheta[,3] <- log(theta[,3]/ (1 - theta[,3]))

group_dat <- data.frame(group = 1:U,
                        Asym = theta[,1],
                        offset = theta[,2],
                        growth = theta[,3])

eps <- mvtnorm::rmvnorm(U, sigma = diag(sigma_theta))
transThetaTrue <- XTheta %*% beta

set.seed(123123)

data <- VARGrowthData(dat, "myTime", "myGroup")
data <- data[1:(nrow(data) - (1 - mr_obs_pct)*tpt),] ### just remove 3/4 of most recent year

# ThetaTrend <- LinearModelTrend(data, ~ bs(as.numeric(group)))
ThetaTrend <- LinearModelTrend(data, ~ 1)
GrowthFunction <- GompertzFunction()

ThetaTransform <- ThetaTransformations(list(logTransform,
                                            logTransform,
                                            logitTransform),
                                       inverse = FALSE)

data$outcome <- SimObj$obs$obs

## Break into train/test set
train_data <- data %>% dplyr::filter(group %in% 1:(0.7*U))
test_data <- data %>% dplyr::filter(!(group %in% 1:(0.7*U)))