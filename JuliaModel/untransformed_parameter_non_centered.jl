using Distributions
using DataFrames
using BenchmarkTools
using LogExpFunctions
using XLSX
using StatsModels

using LinearAlgebra
# BLAS.set_num_threads(1)

# using ProfileView Profile.clear()
# @ProfileView.profview begin     
# for i in 1:100        
# a = ll2(zeros(size(X2)[2]))    
# end
# end

####################
## Format DataSet ##
####################

function gomp_data(dat, time_nm, year_nm, outcome_nm = nothing)
    out = copy(dat)
    out.time = dat[:, time_nm]
    out.year = dat[:, year_nm]
    select!(out, Not(time_nm, year_nm))
    if !isnothing(outcome_nm)
        out.outcome = out[:, outcome_nm]
        select!(out, Not(outcome_nm))
    end

    sort!(out, year, time)
end

function gomp_data!(dat, time_nm, year_nm, outcome_nm = nothing)
    dat.time = dat[:, time_nm]
    dat.year = dat[:, year_nm]
    select!(dat, Not(time_nm, year_nm))

    if !isnothing(outcome_nm)
        dat.outcome = dat[:, outcome_nm]
        select!(dat, Not(outcome_nm))
    end

    sort!(dat, year, time)
end

######################
## Prior Definition ##
######################
### TODO create a Dict of priors to be passed and called in fc functions
### aborted for now
### should be a Dict that maps priors to param names

priors = Dict(:asym =>,
            :offset =>,
            :growth =>,
            :asym_raw =>,
            :offset_raw =>,
            :growth_raw =>,
            :beta_asym => Normal(20000, 1),
            :beta_offset => Normal(2, 0.1),
            :beta_growth => Normal(9, 0.1),
            :latent_sd_asym => truncated(Normal(0, 1)),
            :latent_sd_offset => truncated(Normal(0, 1)), 
            :latent_sd_growth => truncated(Normal(0,1)),
            :obs_sd => InvGamma(1, 1))


#################
## Likelihoods ##
#################

### TODO fix up broadcasting data types using comprehensions
### TODO create a likelihood to be passed to other functions


### Likelihoods for observed process
### gompertz for scalar params asym offset growth
### use "trans_" for parameters on their native scale 
### otherwise use param name for real line
### for use with scalars
function gompertz_likelihood(dat, asym, offset, growth, obs_sd, log_lik = true)
    out = 0
    
    # asym log scale
    # growth logit scale
    # offset log scale

    # mu = exp.(asym[dat.group]).* exp.( -exp.(offset[dat.group]) .* 
    #     exp.( log.( logistic.(growth[dat.group]) ) .* dat.time) #change logistic_growth ^ dat.time to the log scale
    #     )

    ### check interaction of .=, but should be equivalent to above hopefully
    # mu = @. (exp(asym[dat.group]) * exp( -exp(offset[dat.group]) * 
    #     exp( log( logistic(growth[dat.group]) ) * dat.time) #change logistic_growth ^ dat.time to the log scale
    #     ))

    mu = @. (asym[dat.group] * exp( -offset[dat.group] * exp( log( growth[dat.group])  * dat.time) ) )

    

    d = Normal.(mu, obs_sd) ## sd assumed to be scalar
    if log_lik
        out += sum(logpdf.(d, dat.outcome))
    else
        out += prod(pdf.(d, dat.outcome))
    end


end

# prior <- dict(:beta_asym => Normal(0,1),
#               :beta_growth =>,
#               :beta_offset =>,
#               :latent_sd_asym =>,
#               :latent_sd_growth =>,
#               :latent_sd_offset,
#               :obs_sd => ,)

function log_posterior(dat, asym, offset, growth, obs_sd, beta_asym, X_asym, latent_sd_asym, beta_offset, X_offset, latent_sd_offset, beta_growth, X_growth, latent_sd_growth)
    out = 0.0
    ## observed process likelihood
    mu = @. (asym[dat.group] * exp( -offset[dat.group] * exp( log( growth[dat.group])  * dat.time) ) )
    out += sum(logpdf.(d, dat.outcome))

    ## latent process likelihood
    eta_asym = X_asym*beta_asym
    d_asym = Normal.(eta_asym, latent_sd_asym)
    out += sum(logpdf.(d_asym, asym))

    eta_offset = X_offset*beta_offset
    d_offset = Normal.(eta_offset, latent_sd_offset)
    out += sum(logpdf.(d_offset, offset))

    eta_growth = X_growth*beta_growth
    d_growth = Normal.(eta_growth, latent_sd_growth)
    out += sum(logpdf.(d_growth, growth))

    ## priors I will need to hard code this each time it changes...

    d2 = Uniform.(lwr, upr)
    out += sum(logpdf.(d2, gparam))

    # Prior on the VARIANCE
    d1 = truncated(Normal(prior_mean, prior_sd), lower = 0)
    out += sum(logpdf.(d1, latent_sd.^2))
    return(out) 

    # Prior on the VARIANCE
    d1 = InverseGamma(prior_shape, prior_scale)
    out += sum(logpdf.(d1, obs_sd.^2))
    return(out) 


end

### intended for scalars
# function gompertz_likelihood_scalar(outcome, asym, offset, growth, time, obs_sd, log_lik = true)
#     out = 0
#     mu = asym * exp( -offset * exp( log( growth)  * time) ) 

#     d = Normal(mu, obs_sd) ## sd assumed to be scalar
#     if log_lik
#         out += sum(logpdf(d, outcome))
#     else
#         out += prod(pdf(d, outcome))
#     end
# end

# function gompertz_likelihood(dat, params, log_lik = true)
#     out = 0.0
#     ## sd assumed to be scalar
#     # asym log scale
#     # growth logit scale
#     # offset log scale

#     mu = exp.(params.asym).*exp.( -exp.(params.offset) .* dat.time .* log.(logistic.(params.growth)))
#     d = Normal.(mu, params.obs_sd) 
#     if log_lik
#         out += sum(logpdf.(d, dat.outcome))
#     else
#         out += prod(pdf.(d, dat.outcome))
#     end


# end

### Latent likelihoods
### ideally: assign independently for each parameter of interest


# function latent_liklihood_linear(param, beta, model_matrix, sd, log_lik=true)
#     d = Normal(model_matrix * beta, sd)
#     if log_lik
#         out += sum(logpdf.(d, param))
#     else
#         out += sum(logpdf.(d, param))
#     end
# end

# function latent_likelihood_linear(param_name, dat, all_params::DataFrame, model_matrix, log_lik = true)
#     req_params = ["", "beta_", "sd_"] .* param_name
#     if !( issubset(req_params, names(all_params)) )
#         ErrorException("Required Parameters " * join(req_params, ",") * " not in param names")
#     end
# end



################################
## Full Conditional Functions ##
################################

## Note: all params are assumed independent a priori...I wonder if this is a problem
## read more: https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations

## Found something on priors for bspline regression models here
## https://link.springer.com/article/10.1007/s11222-012-9314-z


## Assume Asym, offset, growth are from the untransformed support (i.e. positive reals, interval 0 to 1)
## Assume these are vectors of length nYears
## TODO: Double check prior syntax relevant to my model

## General full conditional for the beta parameters for each gompertz parameter
## Assumes linear model for each parameter
function fc_beta(beta, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1)
    out = 0.0
    
    if (length(beta) != length(prior_mean)) & (length(prior_mean) > 1)
        error("Improper length of prior mean")
    end

    # Likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)  

    # Prior
    d1 = Normal(prior_mean[1], prior_intercept_sd)
    out += sum(logpdf.(d1, beta[1,:]))
    if((size(beta)[1]) > 1)
        d2 = Normal.(prior_mean[2:(size(beta)[1])], prior_sd)
        out += sum(logpdf.(d2, beta[2:(size(beta)[1]),:]))
    end
    
    return(out) 
end

function fc_beta(beta, dat, asym, offset, growth, obs_sd, prior)
    out = 0.0
    
    if (length(beta) != length(prior_mean)) & (length(prior_mean) > 1)
        error("Improper length of prior mean")
    end

    # Likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)  

    # Prior
    d1 = prior
    ## I would like beta to be a vector rather than a matrix here
    out += sum(logpdf.(d1, beta[1,:]))
    if((size(beta)[1]) > 1)
        d2 = Normal.(prior_mean[2:(size(beta)[1])], prior_sd)
        out += sum(logpdf.(d2, beta[2:(size(beta)[1]),:]))
    end
    
    return(out) 
end

## General full conditional for Gompertz parameters
## Change X, beta to corresponding param's X, beta, latent_sd
## no longer should be used in this file

function fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat, lwr, upr, lwr_lim, upr_lim)
    eta = X*beta
    out = 0.0

    if(length(lwr) > 1 & length(upr) > 1)
        if(length(lwr) != length(upr) | length(lwr) != length(gparam))
            stop("length for priors lwr and upr not same length as gparam")
        end
    end

    if(length(lwr) == 1 ⊻ length(upr) == 1)
        if(max(length(lwr), length(upr)) != length(gparam))
            error("length of either lwr or upper not same length as gparam")
        end
    end

    ## "prior"
    # d = truncated.(Normal.(eta, latent_sd), lwr_lim, upr_lim)
    ## remove truncation to avoid unintended consequences in reparameterization
    d = Normal.(eta, latent_sd)
    out += sum(logpdf.(d, gparam))
    ## uniform prior to constrain the parameters
    ## these should be on the real line scale
    # d2 = Uniform.(lwr, upr)
    # out += sum(logpdf.(d2, gparam))

    # Likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)                      

end

## replace asym/offset/growth with value of generate_gparam
function fc_gparam_raw(gparam_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, dat, lwr_lim, upr_lim)
    # gparam = gparam_raw*latent_sd + eta is non-centered param
    # generate that value before this function
    out = 0.0
    # gparam = Array{Float64}(undef, size(gparam_raw))
    # gparam = generate_gparam(gparam_raw, X, beta, latent_sd)

    # if any(gparam .< lwr_lim) | any(gparam .> upr_lim) 
    #     return(-Inf)
    # end

    ## "prior"
    d = Normal(0, 1)
    out += sum(logpdf.(d, gparam_raw))
    ## uniform prior to constrain the parameters
    ## these should be on the real line scale no priors for now...
    # d2 = Uniform.(lwr, upr)
    # out += sum(logpdf.(d2, gparam))

    # Likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)                      

end

## General full conditional for latent sd. Should be one for each gp param
## Test out truncated normal distributions
## gp_param = gp param of interest
function fc_latent_sd(latent_sd, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_sd = 1)
  if any(latent_sd .<= 0)
    return(-Inf)
  end
  out = 0.0
    # Likelihood
    ## basically all likelihoods replaced with the observed data likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd) 

    # Prior on the VARIANCE
    d1 = truncated(Normal(prior_mean, prior_sd), lower = 0)
    out += sum(logpdf.(d1, latent_sd.^2))
    return(out) 
end

## need to figure out how to cast arguments to simpler types
function fc_obs_sd(obs_sd, asym, offset, growth, dat, prior_shape = 0.1, prior_scale = 0.1)
    out = 0.0
    
    if any(obs_sd .<= 0)
        return(-Inf)
    end

    # Likelihood
    
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)  

    # Prior on the VARIANCE
    d1 = InverseGamma(prior_shape, prior_scale)
    out += sum(logpdf.(d1, obs_sd.^2))
    return(out) 
end

#############################
### Convenience Functions ###
#############################

## Convenience function to create random betas
function generate_init_betas(X, N, σ=0.5)
    return(rand(Normal(0, σ), (size(X)[2], N-1)))
end

## Big function to arrange data and store initial parameter values
## Adjust for X matrices for each gompertz parameter
## assumes a variable "fulldata" exists in parent environment
## will need to manually define X matrices here
## inits will contain the initial values for parameters
function generate_init_state(data, inits)

    ## don't need the missing data type but doesn't hurt
    X_asym = Matrix{Union{Float64, Missing}}(undef, length(unique(data.group)), 1)
    X_asym[:,1] .= 1

    X_offset = Matrix{Union{Float64, Missing}}(undef, length(unique(data.group)), 1)
    X_offset[:,1] .= 1

    X_growth = Matrix{Union{Float64, Missing}}(undef, length(unique(data.group)), 1)
    X_growth[:,1] .= 1
    
    inits = Dict(:outcome => data.outcome,
        :time => data.time,
        :group => data.group,
        :X_asym => X_asym,
        :X_offset => X_offset,
        :X_growth => X_growth,
        :asym => inits[:asym],
        :offset => inits[:offset],
        :growth => inits[:growth],
        :raw_asym => (inits[:asym] - X_asym*inits[:beta_asym]) ./ inits[:latent_sd_asym],
        :raw_offset => (inits[:offset] - X_offset*inits[:beta_offset]) ./ inits[:latent_sd_offset],
        :raw_growth => (inits[:growth] - X_growth*inits[:beta_growth]) ./ inits[:latent_sd_growth],
        :beta_asym => inits[:beta_asym],
        :beta_offset => inits[:beta_offset],
        :beta_growth => inits[:beta_growth],
        :latent_sd_asym => inits[:latent_sd_asym],
        :latent_sd_offset => inits[:latent_sd_offset],
        :latent_sd_growth => inits[:latent_sd_growth],
        :obs_sd => inits[:obs_sd],
        :accept_raw_asym => [0],
        :accept_raw_offset => [0],
        :accept_raw_growth => [0],
        :accept_beta_asym => [0],
        :accept_beta_offset => [0],
        :accept_beta_growth => [0],
        :accept_latent_sd_asym => [0],
        :accept_latent_sd_offset => [0],
        :accept_latent_sd_growth => [0],
        :accept_obs_sd => [0]
    )
    return(inits)
end


####################
## Draw Functions ##
####################

function perturb_normal!(vec, sd=0.5, proportion = 1)
    prop = Normal(0, sd) # proposal distribution
    if (proportion == 1)
        vec += rand(prop, size(vec))
    else
        idx = Int.(sort(sample(1:length(vec), Int(floor(length(vec)*proportion)), replace=false)))
        vec[idx] += rand(prop, length(idx))
    end
    return(vec)
end


function generate_gparam(gparam_raw, X, beta, latent_sd)
    # currently don't know how to use lwr_lim and upr_lim in simulation
    # gparam ~ N(eta, latent_sd)
    eta = X*beta
    
    gparam = gparam_raw .* latent_sd + eta
    return(gparam)
end

function generate_gparam_raw(beta, X, gparam, latent_sd)
    # currently don't know how to use lwr_lim and upr_lim in simulation
    eta = X*beta
    
    gparam_raw = (gparam - eta) ./ latent_sd
    return(gparam)
end

function draw_growth!(x0, dat, σ=0.1, track_accept = true)
    # proposed = perturb_normal!(copy(x0[:growth]), σ)
    proposed = perturb_normal!(copy(x0[:growth]), σ)
    ## fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat, lwr, upr, lwr_lim, upr_lim)
    f0 = fc_gparam(x0[:growth], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                   x0[:beta_growth], x0[:X_growth], x0[:latent_sd_growth], dat, 0.8, 0.95, 0, 1)

    f1 = fc_gparam(proposed, x0[:asym], x0[:offset], proposed, x0[:obs_sd], 
                    x0[:beta_growth], x0[:X_growth], x0[:latent_sd_growth], dat, 0.8, 0.95, 0, 1)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:growth] = proposed
        if track_accept
            x0[:accept_raw_growth] = x0[:accept_raw_growth] .+ 1
        end
    end
end

function draw_raw_asym!(x0, dat, σ=0.1, track_accept = false)
    ## generate_gparam(X, beta, latent_sd, lwr, upr, lwr_lim = -Inf, upr_lim = Inf)
    ## fc_gparam_raw(gparam_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, dat, lwr_lim, upr_lim)
    prop_gparam = Array{Float64}(undef, size(x0[:asym]))

    proposed = perturb_normal!(copy(x0[:raw_asym]), σ)
    prop_gparam = generate_gparam(proposed, x0[:X_asym], x0[:beta_asym], x0[:latent_sd_asym])

    f0 = fc_gparam_raw(x0[:raw_asym], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                        x0[:X_asym], x0[:beta_asym], x0[:latent_sd_asym], 
                      dat, 0, Inf)

    f1 = fc_gparam_raw(proposed, prop_gparam, x0[:offset], x0[:growth], x0[:obs_sd], 
                        x0[:X_asym], x0[:beta_asym], x0[:latent_sd_asym],
                        dat, 0, Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:raw_asym] = proposed
        x0[:asym] = prop_gparam
        if track_accept
            x0[:accept_raw_asym] = x0[:accept_raw_asym] .+ 1
        end
    end
end

function draw_raw_offset!(x0, dat, σ=0.1, track_accept = false)
    ## generate_gparam(X, beta, latent_sd, lwr, upr, lwr_lim = -Inf, upr_lim = Inf)
    ## fc_gparam_raw(gparam_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, dat, lwr_lim, upr_lim)
    prop_gparam = Array{Float64}(undef, size(x0[:offset]))

    proposed = perturb_normal!(copy(x0[:raw_offset]), σ)
    prop_gparam = generate_gparam(proposed, x0[:X_offset], x0[:beta_offset], x0[:latent_sd_offset])
    
    f0 = fc_gparam_raw(x0[:raw_offset], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
    x0[:X_offset], x0[:beta_offset], x0[:latent_sd_offset],
    dat, 0, Inf)

    f1 = fc_gparam_raw(proposed, x0[:asym], prop_gparam, x0[:growth], x0[:obs_sd], 
    x0[:X_offset], x0[:beta_offset], x0[:latent_sd_offset],
    dat, 0, Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:raw_offset] = proposed
        x0[:offset] = prop_gparam
        if track_accept
            x0[:accept_raw_offset] = x0[:accept_raw_offset] .+ 1
        end
    end
end


function draw_raw_growth!(x0, dat, σ=0.1, track_accept = false)
    ## generate_gparam(X, beta, latent_sd, lwr, upr, lwr_lim = -Inf, upr_lim = Inf)
    ## fc_gparam_raw(gparam_raw, asym, offset, growth, obs_sd, X, beta, latent_sd, dat, lwr_lim, upr_lim)
    # generate_gparam_raw(beta, X, gparam, latent_sd)
    prop_gparam = Array{Float64}(undef, size(x0[:growth]))

    proposed = perturb_normal!(copy(x0[:raw_growth]), σ)
    prop_gparam = generate_gparam(proposed, x0[:X_growth], x0[:beta_growth], x0[:latent_sd_growth])

    f0 = fc_gparam_raw(x0[:raw_growth], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                x0[:X_growth], x0[:beta_growth], x0[:latent_sd_growth],
                dat, 0, 1)
    f1 = fc_gparam_raw(proposed, x0[:asym], x0[:offset], prop_gparam, x0[:obs_sd], 
                    x0[:X_growth], x0[:beta_growth], x0[:latent_sd_growth],
                    dat, 0, 1)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:raw_growth] = proposed
        x0[:growth] = prop_gparam
        if track_accept
            x0[:accept_raw_growth] = x0[:accept_raw_growth] .+ 1
        end
    end
end

### Drawing 
### Be sure to change priors appropriately
function draw_beta_asym!(x0, dat, σ=0.1, track_accept = false)
    prop_gparam = Array{Float64}(undef, size(x0[:growth]))

    proposed = perturb_normal!(copy(x0[:beta_asym]), σ)
    ## need to update the gparam based on current proposal
    prop_gparam = generate_gparam(x0[:raw_asym], x0[:X_asym], proposed, x0[:latent_sd_asym])
    # fc_beta(beta, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1)
    ## TODO: CHECK IF GPARAM NEEDS TO BE UPDATED WITH PROPOSALS FOR BETA
    f0 = fc_beta(x0[:beta_asym], 
                dat, 
                prop_gparam, 
                x0[:offset], 
                x0[:growth], 
                x0[:obs_sd], 
                20000,
                10, 
                10)

    f1 = fc_beta(proposed, 
                dat, 
                prop_gparam, 
                x0[:offset], 
                x0[:growth], 
                x0[:obs_sd], 
                20000,
                10, 
                10)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_asym] = proposed
        x0[:asym] = prop_gparam
        if track_accept
            x0[:accept_beta_asym] = x0[:accept_beta_asym] .+ 1
        end
    end
end

function draw_beta_offset!(x0, dat, σ=0.1, track_accept = false)
    proposed = perturb_normal!(copy(x0[:beta_offset]), σ)
    ## need to update the gparam based on current proposal
    prop_gparam = generate_gparam(x0[:raw_offset], x0[:X_offset], proposed, x0[:latent_sd_offset])
    # fc_beta(beta, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1)
    f0 = fc_beta(x0[:beta_offset], 
                dat, 
                x0[:asym], 
                x0[:offset], 
                x0[:growth], 
                x0[:obs_sd], 
                2,
                0.2, 
                0.2)

    f1 = fc_beta(proposed, 
                dat, 
                x0[:asym], 
                prop_gparam, 
                x0[:growth], 
                x0[:obs_sd], 
                2,
                0.2, 
                0.2)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_offset] = proposed
        x0[:offset] = prop_gparam
        if track_accept
            x0[:accept_beta_offset] = x0[:accept_beta_offset] .+ 1
        end
    end
end

function draw_beta_growth!(x0, dat, σ=0.1, track_accept = false, noncenter = true)
    proposed = perturb_normal!(copy(x0[:beta_growth]), σ)
    prop_gparam = generate_gparam(x0[:raw_growth], x0[:X_growth], proposed, x0[:latent_sd_growth])
    # fc_beta(beta, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1)
    f0 = fc_beta(x0[:beta_growth], 
                dat, 
                x0[:asym], 
                x0[:offset], 
                x0[:growth], 
                x0[:obs_sd], 
                0.9,
                0.03, 
                0.03)
    if noncenter
        prop_gparam = generate_gparam(x0[:raw_growth], x0[:X_growth], proposed, x0[:latent_sd_growth])
        f1 = fc_beta(proposed, 
        dat, 
        x0[:asym], 
        x0[:offset], 
        prop_gparam, 
        x0[:obs_sd], 
        0.9,
        0.03, 
        0.03)
    else
        f1 = fc_beta(proposed, 
        dat, 
        x0[:asym], 
        x0[:offset], 
        x0[:growth], 
        x0[:obs_sd], 
        0.9,
        0.03, 
        0.03)
    end
    
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_growth] = proposed
        if noncenter
            x0[:growth] = prop_gparam
        end
        if track_accept
            x0[:accept_beta_growth] = x0[:accept_beta_growth] .+ 1
        end
    end
end


### latent standard deviations
function draw_latent_sd_asym!(x0, dat,  σ = 0.01, track_accept = false)
    prop = perturb_normal!(copy(x0[:latent_sd_asym]), σ)
    prop_gparam = generate_gparam(x0[:raw_asym], x0[:X_asym], x0[:beta_asym], prop)
    # function fc_latent_sd(latent_sd, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_sd = 1)
    f0 = fc_latent_sd(latent_sd_asym, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 1)
    f1 = fc_latent_sd(prop, dat, prop_gparam, x0[:offset], x0[:growth], x0[:obs_sd], 0, 1)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_asym] = prop
        x0[:asym] = prop_gparam
        if track_accept
            x0[:accept_latent_sd_asym] = x0[:accept_latent_sd_asym] .+ 1
        end
    end
end

function draw_latent_sd_offset!(x0, dat,  σ = 0.01, track_accept = false)
    prop = perturb_normal!(copy(x0[:latent_sd_offset]), σ)
    prop_gparam = generate_gparam(x0[:raw_offset], x0[:X_offset], x0[:beta_offset], prop)
    # function fc_latent_sd(latent_sd, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_sd = 1)
    f0 = fc_latent_sd(latent_sd_offset, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 0.1)
    f1 = fc_latent_sd(prop, dat, x0[:asym], prop_gparam, x0[:growth], x0[:obs_sd], 0, 0.1)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_offset] = prop
        x0[:offset] = prop_gparam
        if track_accept
            x0[:accept_latent_sd_offset] = x0[:accept_latent_sd_offset] .+ 1
        end
    end
end

function draw_latent_sd_growth!(x0, dat, σ = 0.01, track_accept = false, noncenter = true)
    prop = perturb_normal!(copy(x0[:latent_sd_growth]), σ)
    
    # function fc_latent_sd(latent_sd, dat, asym, offset, growth, obs_sd, prior_mean = 0, prior_sd = 1)
    f0 = fc_latent_sd(latent_sd_growth, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 0.1)

    if noncenter
        prop_gparam = generate_gparam(x0[:raw_growth], x0[:X_growth], x0[:beta_growth], prop)
        f1 = fc_latent_sd(prop, dat, x0[:asym], x0[:offset], prop_gparam, x0[:obs_sd], 0, 0.1)
    else
        f1 = fc_latent_sd(prop, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 0.1)
    end
    
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_growth] = prop
        if noncenter
            x0[:growth] = prop_gparam
        end
        if track_accept
            x0[:accept_latent_sd_growth] = x0[:accept_latent_sd_growth] .+ 1
        end
    end
end

## if i need to draw latent sd's jointly...
# function draw_latent_sd_all!(x0, σ = 0.1)
#     prop1 = perturb_normal!(copy(x0[:latent_sd_asym]), σ)
#     prop2 = perturb_normal!(copy(x0[:latent_sd_offset]), σ)
#     prop3 = perturb_normal!(copy(x0[:latent_sd_growth]), σ)
#     ##fc_latent_sd(latent_sd, beta, X, outcome, prior_sd = 1)
#     f0 = fc_latent_sd(x0[:latent_sd_asym], x0[:beta_asym], x0[:X_asym], x0[:asym], 1) + 
#         fc_latent_sd(x0[:latent_sd_offset], x0[:beta_offset], x0[:X_offset], x0[:offset], 1) +
#         fc_latent_sd(x0[:latent_sd_growth], x0[:beta_growth], x0[:X_growth], x0[:growth], 1)
#     f1 = fc_latent_sd(prop1, x0[:beta_asym], x0[:X_asym], x0[:asym], 1) + 
#         fc_latent_sd(prop2, x0[:beta_offset], x0[:X_offset], x0[:offset], 1) +
#         fc_latent_sd(prop3, x0[:beta_growth], x0[:X_growth], x0[:growth], 1)
#     a = (f1-f0)              
#     if (log(rand(Uniform(0,1),1)[1]) < a )
#         x0[:latent_sd_growth] = prop
#     end
# end

function draw_obs_sd!(x0, dat, σ = 0.1, track_accept = false)
    prop = perturb_normal!(copy(x0[:obs_sd]), σ)
    ##fc_obs_sd(obs_sd, asym, offset, growth, dat, prior_shape = 0.1, prior_scale = 0.1)
    f0 = fc_obs_sd(x0[:obs_sd], x0[:asym], x0[:offset], x0[:growth], dat, 1, 1)
    f1 = fc_obs_sd(prop, x0[:asym], x0[:offset], x0[:growth], dat, 1, 1)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:obs_sd] = prop
        if track_accept
            x0[:accept_obs_sd] = x0[:accept_obs_sd] .+ 1
        end
    end
end



function get_labels(obj, objname)
    if length(size(obj)) == 1
        return([objname * "["*string(i) * "]" for i in 1:size(obj)[1]])
    elseif length(size(obj)) == 2
        return([objname * "["*string(i) * "-" * string(j)*"]" for j in 1:size(obj)[2] for i in 1:size(obj)[1]])
    else
        error("Labels only implemented for 1 and 2 dimensional objects")
    end
end


function flatten(obj)
    if length(size(obj)) == 1
        return(obj)
    elseif length(size(obj)) == 2
        return(reshape(obj, length(obj)))
    else
        error("Flatten only implemented for 1 and 2 dimensional objects")
    end
end

function write_csv_line(x0; symbols, file, newfile=true)
    # Check to make sure we're opening a new file or appending to an existing one
    if (!newfile && typeof(file) != IOStream)
        error("When writing data, use an open IOStream rather than reopening the file")
    end
    labelvec = vcat([get_labels(x0[s], string(s)) for s in symbols]...)

    if (newfile)
        outfile = open(file, "w")
        write(outfile, join(labelvec, ",")*"\n")
        close(outfile)
    else
        valvec = vcat([flatten(x0[s]) for s in symbols]...)
        write(file, join([string(x) for x in valvec],",")*"\n")
    end
end




monitors = [:asym,
            :offset,
            :growth,
            :raw_asym,
            :raw_offset,
            :raw_growth,
            :beta_asym,
            :beta_offset,
            :beta_growth,
            :latent_sd_asym,
            :latent_sd_offset,
            :latent_sd_growth,
            :obs_sd]

accept_keys = [:accept_raw_asym ,
:accept_raw_offset ,
:accept_raw_growth ,
:accept_beta_asym ,
:accept_beta_offset ,
:accept_beta_growth ,
:accept_latent_sd_asym ,
:accept_latent_sd_offset ,
:accept_latent_sd_growth ,
:accept_obs_sd ]

function run_mcmc_chain(fname, x0, monitors, dat, warmup, run, thin=10, verbose = 1)
    ### probably should add checks for valid transformations raw -> gparam
    
    if (!isfinite(
            fc_gparam_raw(x0[:raw_asym], 
                          x0[:asym], 
                          x0[:offset], 
                          x0[:growth], 
                          x0[:obs_sd],
                          x0[:X_asym], x0[:beta_asym], x0[:latent_sd_asym],
                          dat, 0, Inf)
            )
        )
    error("raw_asym full conditional not fininte at init")
    end
    if (!isfinite(fc_gparam_raw(x0[:raw_offset], x0[:offset], x0[:offset], x0[:growth], x0[:obs_sd], 
        x0[:X_offset], x0[:beta_offset], x0[:latent_sd_offset],
        dat, 0, Inf)))
    error("raw_offset full conditional not fininte at init")
    end

    if (!isfinite(fc_gparam_raw(x0[:raw_growth], x0[:growth], x0[:offset], x0[:growth], x0[:obs_sd],
        x0[:X_growth], x0[:beta_growth], x0[:latent_sd_growth], dat, 0, 1)))
    error("raw_growth full conditional not fininte at init")
    end


    if (!isfinite(fc_beta(x0[:beta_asym], 
        dat, 
        x0[:asym], 
        x0[:offset], 
        x0[:growth], 
        x0[:obs_sd], 
        20000,
        10, 
        10)))
    error("Beta asym full conditional not finite at init")
    end

    if (!isfinite(fc_beta(x0[:beta_offset], 
        dat, 
        x0[:asym], 
        x0[:offset], 
        x0[:growth], 
        x0[:obs_sd], 
        2,
        0.2, 
        0.2)))
    error("Beta offset full conditional not finite at init")
    end

    if (!isfinite(
            fc_beta(x0[:beta_growth], 
            dat, 
            x0[:asym], 
            x0[:offset], 
            x0[:growth], 
            x0[:obs_sd], 
            0.85,
            0.03, 
            0.03)
            )
        )
        error("Beta growth full conditional not finite at init")
    end
    

    if (!isfinite(fc_latent_sd(latent_sd_asym, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 1)))
        error("asym latent sd full conditional not finite at init")
    end

    if (!isfinite(fc_latent_sd(latent_sd_offset, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 0.1)))
        error("offset latent sd full conditional not finite at init")
    end

    if (!isfinite(fc_latent_sd(latent_sd_growth, dat, x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 0, 0.1)))
        error("growth latent sd full conditional not finite at init")
    end

    if (!isfinite(fc_obs_sd(x0[:obs_sd], x0[:asym], x0[:offset], x0[:growth], dat, 100)))
        error("obs sd full conditional not finite at init")
    end


    if verbose > 0
        print("Initializing output file\n")
    end
    write_csv_line(x0,symbols=monitors, file=fname, newfile=true)
    outfile = open(fname, "a")
    write_csv_line(x0,symbols=monitors, file=outfile, newfile=false)
    
    if verbose > 0
        print("Beginning warmup\n")
    end
    
    # perturb_sd = Dict(:asym => 0.3, #asym
    #              :offset => 0.15, #offset
    #              :growth => 0.007, #growth
    #              :beta_asym => 30, #beta_asym
    #              :beta_offset => 0.0035,  #beta_offset
    #              :beta_growth => 0.00013, #beta_growth
    #              :latent_sd_asym => 1.5, #latent_sd_asym
    #              :latent_sd_offset => 0.003, #latent_sd_offset
    #              :latent_sd_growth => 0.00013, #latent_sd_growth
    #              :obs_sd => 5 #obs_sd
    # )
    perturb_sd = Dict(:asym => 0.3, #asym
                 :offset => 0.15, #offset
                 :growth => 0.007, #growth
                 :beta_asym => 30, #beta_asym
                 :beta_offset => 0.0035,  #beta_offset
                 :beta_growth => 0.005, #beta_growth
                 :latent_sd_asym => 1.5, #latent_sd_asym
                 :latent_sd_offset => 0.003, #latent_sd_offset
                 :latent_sd_growth => 0.005, #latent_sd_growth
                 :obs_sd => 5 #obs_sd
    )
    
    for i in 1:warmup
            
        if verbose > 0
            if i % 100 == 0
                print(string(i)*"/"*string(warmup)*"\n")
            end            
        end
        # draw_raw_asym!(x0, dat, perturb_sd[:asym])
        # draw_raw_offset!(x0, dat, perturb_sd[:offset])
        # draw_raw_growth!(x0, dat, perturb_sd[:growth])

        # draw_beta_asym!(x0, dat, perturb_sd[:beta_asym])
        # draw_beta_offset!(x0, dat, perturb_sd[:beta_offset])
        # draw_beta_growth!(x0, dat, perturb_sd[:beta_growth])

        # draw_latent_sd_asym!(x0, dat, perturb_sd[:latent_sd_asym])
        # draw_latent_sd_offset!(x0, dat, perturb_sd[:latent_sd_offset])
        # draw_latent_sd_growth!(x0, dat, perturb_sd[:latent_sd_growth])

        # draw_obs_sd!(x0, dat, perturb_sd[:obs_sd])
    end

    
    if verbose > 0
        print("Warmup complete, starting sample\n")
    end

    
    for i in 1:run
        
        if verbose > 0
        if i % 100 == 0
                print(string(i)*"/"*string(run)*"\n")
            end            
        end
        if i % thin == 0
            write_csv_line(x0,symbols=monitors, file=outfile, newfile=false)
        end
        # draw_raw_asym!(x0, dat, perturb_sd[:asym], true)
        # draw_raw_offset!(x0, dat, perturb_sd[:offset], true)
        # draw_raw_growth!(x0, dat, perturb_sd[:growth], true)
        draw_growth!(x0, dat, perturb_sd[:growth], true)
        # draw_beta_asym!(x0, dat, perturb_sd[:beta_asym], true)
        # draw_beta_offset!(x0, dat, perturb_sd[:beta_offset], true)
        draw_beta_growth!(x0, dat, perturb_sd[:beta_growth], true, false)

        # draw_latent_sd_asym!(x0, dat, perturb_sd[:latent_sd_asym], true)
        # draw_latent_sd_offset!(x0, dat, perturb_sd[:latent_sd_offset], true)
        draw_latent_sd_growth!(x0, dat, perturb_sd[:latent_sd_growth], true, false)

        # draw_obs_sd!(x0, dat, perturb_sd[:obs_sd], true)
    end
    close(outfile)

    ## output acceptance ratios
    print("\n")
    print("Acceptance Asym: ", string(x0[:accept_raw_asym]/run), "\n")
    print("Acceptance offset: ", string(x0[:accept_raw_offset]/run), "\n")
    print("Acceptance growth: ", string(x0[:accept_raw_growth]/run), "\n")
    # print("Acceptance growth: ", string(x0[:accept_growth]/run), "\n")
    print("Acceptance beta_asym: ", string(x0[:accept_beta_asym]/run), "\n")
    print("Acceptance beta_offset: ", string(x0[:accept_beta_offset]/run), "\n")
    print("Acceptance beta_growth: ", string(x0[:accept_beta_growth]/run), "\n")

    print("Acceptance latent_sd_asym: ", string(x0[:accept_latent_sd_asym]/run), "\n")
    print("Acceptance latent_sd_offset: ", string(x0[:accept_latent_sd_offset]/run), "\n")
    print("Acceptance latent_sd_growth: ", string(x0[:accept_latent_sd_growth]/run), "\n")
    print("Acceptance obs_sd: ", string(x0[:accept_obs_sd]/run), "\n")
end




    