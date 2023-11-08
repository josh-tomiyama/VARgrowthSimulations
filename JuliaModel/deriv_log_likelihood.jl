############################
### gompertz derivatives ###
############################

function gomp(;asym, offset, growth, time)
    asym*exp(-offset * exp(log(growth) * time))
end

## assuming the untransformed parameters: need additional terms for transformed parameters
## by chain rule
function gomp_deriv_asym(;asym, offset, growth, time)
    exp(-offset * growth ^ time)
end

function gomp_deriv_offset(;asym, offset, growth, time)
    -asym * exp(-offset * growth ^ time) * growth^time
end

function gomp_deriv_growth(;asym, offset, growth, time)
    if time == 0
        return 0
    end

    -asym * offset*exp(-offset * growth ^ time)*(time)*growth^(time-1)
end

function check_val(est, truth, msg::String)
    if abs(est - truth) > 10^-6
        print(msg*"\n")
    end
end

### checks
using ForwardDiff
using Distributions
asym = rand(Uniform(0, 100))
offset = rand(Uniform(0, 10))
growth = rand(Uniform(0, 1))
time = rand(Uniform(-1, 1)) ## allow negative times 

truth = ForwardDiff.derivative(x -> gomp(asym = x, offset = offset, growth = growth, time = time), asym) 
est = gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time)
check_val(est, truth, "error in asym derivative")

truth = ForwardDiff.derivative(x -> gomp(asym = asym, offset = x, growth = growth, time = time), offset) 
est = gomp_deriv_offset(asym = asym, offset = offset, growth = growth, time = time)
check_val(est, truth, "error in offset derivative")

truth = ForwardDiff.derivative(x -> gomp(asym = asym, offset = offset, growth = x, time = time), growth) 
est = gomp_deriv_growth(asym = asym, offset = offset, growth = growth, time = time)
check_val(est, truth, "error in growth derivative")

#########################################################################
### partial derivative of log likelihoods - non-centered param ##########
#########################################################################

### expect a single outcome for testing
function log_lik_raw2(;outcome, group, time, raw_asym, raw_offset, raw_growth, obs_sd, 
    beta_asym, X_asym, latent_sd_asym, 
    beta_offset, X_offset, latent_sd_offset, 
    beta_growth, X_growth, latent_sd_growth)

out = 0.0

if any(latent_sd_asym .<= 0)
return(-Inf)
end

if any(latent_sd_offset .<= 0)
return(-Inf)
end

if any(latent_sd_growth .<= 0)
return(-Inf)
end


asym = (X_asym*beta_asym)[group] + latent_sd_asym*raw_asym
println("asym: "*string(asym))
offset = (X_offset*beta_offset)[group] + latent_sd_offset*raw_offset
println("offset: "*string(offset))
growth = (X_growth*beta_growth)[group] + latent_sd_growth*raw_growth
println("growth: "*string(growth))
## observed process likelihood
mu = gomp(asym = asym, offset = offset, growth = growth, time = time)
println("mu: "* string(mu))
d = Normal(mu, obs_sd)
out += logpdf(d, outcome) 

d1 = Normal(0, 1)
out += logpdf(d1, raw_asym)
out += logpdf(d1, raw_offset)
out += logpdf(d1, raw_growth)
return(out) 


end

### non-centered versions 
### functions are written with respect to a scalar asym, offset, growth, not vectors
function deriv_asym_raw(;asym_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) * 
        latent_sd + ## term from XB + sd*raw = Asym
    ## asym_raw ~ N(0, 1)
    -asym_raw
end

function deriv_offset_raw(;offset_raw, asym, offset, growth, obs_sd, X, beta, latent_sd_asym, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_offset(asym = asym, offset = offset, growth = growth, time = time) * 
        latent_sd + ## term from transformation
    ## offset_raw ~ N(0, 1)
    -offset_raw
end

function deriv_growth_raw(;growth_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time)) / obs_sd^2 * 
        gomp_deriv_growth(asym = asym, offset = offset, growth = growth, time = time) * 
        latent_sd + ## term from transformation
    ## growth_raw ~ N(0, 1)
    -growth_raw
end

function deriv_beta_asym(;beta, asym, offset, growth, obs_sd, X, latent_sd, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) * X[group, :]
    ## ignore prior...for now
end

function deriv_beta_offset(;beta, asym, offset, growth, obs_sd, X, latent_sd, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_offset(asym = asym, offset = offset, growth = growth, time = time) * X[group, :]
    ## ignore prior...for now
end

function deriv_beta_growth(;beta, growth_raw, asym, offset, growth, obs_sd, X, latent_sd, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_growth(asym = asym, offset = offset, growth = growth, time = time) * X[group, :]
    ## ignore prior...for now
end

### NOTE FOR IMPLEMENTATION: proposals are on sd but priors are on the variance ###

function deriv_latent_sd_asym(;latent_sd, raw_asym, asym, offset, growth, obs_sd, X, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) * raw_asym
    ## ignore prior...for now
end

function deriv_latent_sd_offset(;latent_sd, raw_offset, asym, offset, growth, obs_sd, X, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) * raw_offset
    ## ignore prior...for now
end

function deriv_latent_sd_growth(;latent_sd, raw_growth, asym, offset, growth, obs_sd, X, outcome, group, time)
    ## obs process likelihood
    -(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) * raw_growth
    ## ignore prior...for now
end

function deriv_obs_sd(;obs_sd, asym, offset, growth, X, beta, latent_sd, outcome, group, time)
    -obs_sd^(-2) + (outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))^2 / obs_sd ^ 3
    ## ignore prior...for now
end

#################################################
### Test Non-Centered Derivatives scalar case ###
#################################################
asym = rand(Uniform(0, 100))
offset = rand(Uniform(0, 10))
growth = rand(Uniform(0, 1))
time = rand(Uniform(-1, 1))
outcome = gomp(asym = asym, offset = offset, growth = growth, time = time)
latent_sd_asym = rand(Uniform(0, 1))
latent_sd_offset = rand(Uniform(0, 1))
latent_sd_growth = rand(Normal(0, 0.001))
beta = [1, 2]
X = hcat([1, 1, 1], [1, 2, 3])
group = sample(X)
obs_sd = rand(Uniform(0, asym/10))

raw_asym = (asym - (X*beta)[group]) / latent_sd_asym
raw_offset = (offset - (X*beta)[group]) / latent_sd_offset
raw_growth = (growth - (X*beta)[group]) / latent_sd_growth

log_lik_raw2(outcome = outcome, group = group, time = time, 
            raw_asym = raw_asym, raw_offset = raw_offset, raw_growth = raw_growth, obs_sd = obs_sd, 
            beta_asym = beta, X_asym = X, latent_sd_asym = latent_sd_asym, 
            beta_offset = beta, X_offset = X, latent_sd_offset = latent_sd_offset, 
            beta_growth = beta, X_growth = X, latent_sd_growth = latent_sd_growth)

truth = ForwardDiff.derivative(p -> log_lik_raw2(outcome = outcome, group = group, time = time, 
                                    raw_asym = p, raw_offset = raw_offset, raw_growth = raw_growth, obs_sd = obs_sd, 
                                    beta_asym = beta, X_asym = X, latent_sd_asym = latent_sd_asym, 
                                    beta_offset = beta, X_offset = X, latent_sd_offset = latent_sd_offset, 
                                    beta_growth = beta, X_growth = X, latent_sd_growth = latent_sd_growth), raw_asym) 

est = deriv_asym_raw(asym_raw = raw_asym, asym = asym, offset = offset, growth = growth, obs_sd = obs_sd,
                    X = X, beta = beta, latent_sd = latent_sd_asym, outcome = outcome, group = group, time = time)
check_val(est, truth, "error in raw_asym derivative")

truth = ForwardDiff.gradient(p -> log_lik_raw2(outcome = outcome, group = group, time = time, 
                                    raw_asym = raw_asym, raw_offset = raw_offset, raw_growth = raw_growth, obs_sd = obs_sd, 
                                    beta_asym = p, X_asym = X, latent_sd_asym = latent_sd_asym, 
                                    beta_offset = beta, X_offset = X, latent_sd_offset = latent_sd_offset, 
                                    beta_growth = beta, X_growth = X, latent_sd_growth = latent_sd_growth), beta + [100000, -10]) 

est = deriv_beta_asym(beta = beta + [100000, -10], asym = asym, offset = offset, growth = growth, 
                     obs_sd = obs_sd, X = X, latent_sd = latent_sd_asym, outcome = outcome, 
                     group = group, time = time)


#####################################################################
### partial derivative of log likelihoods - centered param ##########
#####################################################################

### centered versions ###

### functions are written with respect to a scalar asym, offset, growth, not vectors
function deriv_asym(;asym, offset, growth, obs_sd, X, beta, latent_sd, outcome, group, time)
    ## obs process likelihood
    -2*(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_asym(asym = asym, offset = offset, growth = growth, time = time) +
    ## latent process likelihood
    -2*((X*beta)[group] - asym)/latent_sd^2
end

function deriv_offset(;asym, offset, growth, obs_sd, X, beta, latent_sd, outcome, group, time)
    ## obs process likelihood
    -2*(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_offset(asym = asym, offset = offset, growth = growth, time = time) +
    ## latent process likelihood
    -2*((X*beta)[group] - asym)/latent_sd^2
end

function deriv_growth(;asym, offset, growth, obs_sd, X, beta, latent_sd, outcome, group, time)
    ## obs process likelihood
    -2*(outcome - gomp(asym = asym, offset = offset, growth = growth, time = time))/obs_sd^2 * 
        gomp_deriv_growth(asym = asym, offset = offset, growth = growth, time = time) +
    ## latent process likelihood
    -2*((X*beta)[group] - asym)/latent_sd^2
end

