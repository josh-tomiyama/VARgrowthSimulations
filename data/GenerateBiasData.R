library(VARgrowth)
library(dplyr)

ntpt = c(5, 10, 20, 50)
U = c(10, 25, 50)
nreps = 100

generate_data <- function(seed, tstart, tend, beta, sigma_theta, obsVar, ntpt, U){
  dat <- data.frame(myTime = rep(seq(from = tstart,
                                     to = tend,
                                     by = (tend - tstart)/(ntpt - 1)),
                                     times = U),
                                 myGroup = rep(factor(1:U), each = ntpt)
                    )
  data <- VARGrowthData(dat, "myTime", "myGroup")
  ### FOR NOW: growth function and theta transform are fixed
  GrowthFunction <- GompertzFunction()
  ThetaTransform <- ThetaTransformations(list(logTransform,
                                              logTransform,
                                              logitTransform),
                                         inverse = FALSE)

  ### model 1: intercept model
  ThetaTrend <- LinearModelTrend(data, ~ 1)
  SimObj <- SimulateData(data,
                         param_list = list(ObsVar = obsVar,
                                           SigmaTheta = diag(sigma_theta),
                                           betaTheta = beta),
                         ThetaTrend,
                         GrowthFunction,
                         ThetaTransform,
                         simplify = FALSE,
                         seed = seed)

  data$outcome <- SimObj$obs$obs
  data
}

gompertz_half_point <- function(b2 = 2, b3 = 0.8){
  (log(log(2)) - log(b2))/log(b3)
}

deriv_gompertz <- function(t, Asym = 100, b2 = 2, b3 = 0.8, shift = 0.1*Asym){
  Asym*exp(-b2*b3^t)*(-b2*log(b3)*b3^t) - shift
}

set.seed(123123)
seeds <- sample.int(1000, nreps)

Asym <- 20000
offset <- c(0.5, 1, 2)
growth <- c(0.3, 0.5, 0.7)

mean_param_grid <- expand.grid(Asym, offset, growth)
var_param <- c(0.05)^2
var_param_grid <- expand.grid(var_param, var_param, var_param)
obs_vars <- c(100) ^ 2

Asym_trans <- log(Asym)
offset_trans <- log(offset)
growth_trans <- log(growth/(1-growth))

mean_param_grid <- expand.grid(Asym, offset, growth)

mean_param_grid_trans <- expand.grid(Asym_trans, offset_trans, growth_trans)
var_param_grid <- expand.grid(var_param, var_param, var_param)


t_half <- apply(mean_param_grid, 1, function(row){
  c(floor(gompertz_half_point(row[2], row[3])))
})

t_int_upr <- apply(mean_param_grid, 1, function(row){
  c(floor(gompertz_half_point(row[2], row[3])), floor(gompertz_half_point(row[2], row[3])) + 1000)
},
simplify = FALSE)

t_int_lwr <- apply(mean_param_grid, 1, function(row){
  c(floor(gompertz_half_point(row[2], row[3])) - 100, floor(gompertz_half_point(row[2], row[3])))
},
simplify = FALSE)

derivs <- apply(mean_param_grid, 1, function(row){
  function(t){
     row[1]*exp(-row[2]*row[3]^t)*(-row[2]*log(row[3])*row[3]^t) - 0.01*row[1]
  }
})

roots_lwr <- mapply(uniroot, derivs, t_int_lwr)
roots_lwr <- floor(unlist(roots_lwr[1,]))

roots_upr <- mapply(uniroot, derivs, t_int_upr)
roots_upr <- ceiling(unlist(roots_upr[1,]))

times <- cbind(mean_param_grid_trans, roots_lwr, roots_upr)
colnames(times) <- c("Asym", "offset", "growth", "tstart", "tend")

all_params <- expand.grid(ntpt, U, Asym_trans, offset_trans, growth_trans, var_param, var_param, var_param, obs_vars)
colnames(all_params) <- c("ntpt", "U", "Asym", "offset", "growth", "var_param1", "var_param2", "var_param3", "obs_vars")
nrow(all_params)
all_params <- left_join(all_params, times)

dir_name <- paste0("./InterceptBiasSimulation/")
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

for(i in 1:nrow(all_params)){
  curr_row <- all_params[i,]
  sims <- lapply(seeds,
                 generate_data,
                 tstart = curr_row$tstart,
                 tend = curr_row$tend,
                 beta = c(curr_row$Asym, curr_row$offset, curr_row$growth),
                 sigma_theta = c(curr_row$var_param1, curr_row$var_param2, curr_row$var_param3),
                 obsVar = curr_row$obs_vars,
                 ntpt = curr_row$ntpt,
                 U = curr_row$U)
  fname_prefix <- paste0("sim", i)
  sim_dir <- paste0(dir_name, "Sim", i, "/")
  if(!dir.exists(sim_dir)){
    dir.create(sim_dir)
  }
  true_param_file <- paste0(sim_dir, "Truth.RDS")
  data_file <- paste0(sim_dir, "SimData.RDS")

  saveRDS(curr_row, true_param_file)
  saveRDS(sims, data_file)
}





