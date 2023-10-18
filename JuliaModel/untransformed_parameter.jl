using Distributions
using DataFrames
using BenchmarkTools
using LogExpFunctions
using XLSX
using StatsModels
using PDMats
using LinearAlgebra

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
### aborted for now
### should be a names list that maps priors to param names
### It should also be expected that the params for evalutaiton are in a dataframe/named array


#################
## Likelihoods ##
#################

### TODO fix up broadcasting data types using iterators

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

function gompertz_likelihood(dat, params, log_lik = true)
    out = 0
    ## sd assumed to be scalar
    # asym log scale
    # growth logit scale
    # offset log scale

    mu = exp.(params.asym).*exp.( -exp.(params.offset) .* dat.time .* log.(logistic.(params.growth)))
    d = Normal.(mu, params.obs_sd) 
    if log_lik
        out += sum(logpdf.(d, dat.outcome))
    else
        out += prod(pdf.(d, dat.outcome))
    end


end

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


## Assume Asym, offset, growth are from the real line
## Assume these are vectors of length nYears
## TODO: Double check prior syntax relevant to my model
## TODO: May want more parameters in priors

## General full conditional for the beta parameters for each gompertz parameter
## Assumes linear model for each parameter
function fc_beta(beta, X, gparam, latent_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1, lwr_lim = 0, upr_lim = Inf)
    eta = X * beta
    out = 0.
    
    if (length(beta) != length(prior_mean)) & (length(prior_mean) > 1)
        error("Improper length of prior mean")
    end

    # Likelihood
    d = truncated.(Normal.(eta, latent_sd), lwr_lim, upr_lim)
    out += sum(logpdf.(d, gparam))

    # Prior
    d1 = Normal(prior_mean[1], prior_intercept_sd)
    out += sum(logpdf.(d1, beta[1,:]))
    if((size(beta)[1]) > 1)
        d2 = Normal(prior_mean[2:(size(beta)[1])], prior_sd)
        out += sum(logpdf.(d2, beta[2:(size(beta)[1]),:]))
    end
    
    return(out) 
end

## General full conditional for Gompertz parameters
## Change X, beta to corresponding param's X, beta, latent_sd
## uniform prior gomp params

function fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat, lwr, upr, lwr_lim, upr_lim)
    eta = X*beta
    out = 0.0

    if(length(lwr) > 1 & length(upr) > 1)
        if(length(lwr) != length(upr) | length(lwr) != length(gparam))
            stop("length for priors lwr and upr not same length as gparam")
        end
    end

    if(length(lwr) == 1 ⊻ length(upr) == 1)
        if(max(length(lqwr), length(upr)) != length(gparam))
            error("length of either lwr or upper not same length as gparam")
        end
    end

    ## "prior"
    d = truncated.(Normal.(eta, latent_sd), lwr_lim, upr_lim)
    out += sum(logpdf.(d, gparam))
    ## uniform prior to constrain the parameters
    ## these should be on the real line scale
    d2 = Uniform.(lwr, upr)
    out += sum(logpdf.(d2, gparam))

    # Likelihood
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)                      

end

## General full conditional for latent sd. Should be one for each gp param
## Test out truncated normal distributions
## gp_param = gp param of interest
function fc_latent_sd(latent_sd, beta, X, gp_param, prior_mean = 0, prior_sd = 1, lwr_lim = -Inf, upr_lim = Inf)
  if any(latent_sd .<= 0)
    return(-Inf)
  end

    eta = X * beta
    out = 0.0

    # Likelihood
    d = truncated.(Normal.(eta, latent_sd), lwr_lim, upr_lim)
    out += sum(logpdf.(d, gp_param))

    # Prior on the VARIANCE
    d1 = truncated(Normal(prior_mean, prior_sd), lower = 0)
    out += sum(logpdf.(d1, latent_sd.^2))
    return(out) 
end

## need to figure out how to cast arguments to simpler types
function fc_obs_sd(obs_sd, asym, offset, growth, dat, prior_mean = 0, prior_sd = 1)
    out = 0.0
    
    if any(obs_sd .<= 0)
        return(-Inf)
    end

    # Likelihood
    
    out += gompertz_likelihood(dat, asym, offset, growth, obs_sd)  

    # Prior on the VARIANCE
    # d1 = truncated(Normal(prior_mean, prior_sd), lower = 0)
    d1 = InverseGamma(prior_mean, prior_sd)
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
    ## Parameter model matrices

    ## Create example X matrix for condition, including missing values
    ## Note columns need to be manually changed here
    ## TODO change to Modelframe(form, data), modelmatrix(mf).m
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
        :beta_asym => inits[:beta_asym],
        :beta_offset => inits[:beta_offset],
        :beta_growth => inits[:beta_growth],
        :latent_sd_asym => inits[:latent_sd_asym],
        :latent_sd_offset => inits[:latent_sd_offset],
        :latent_sd_growth => inits[:latent_sd_growth],
        :obs_sd => inits[:obs_sd],
        :accept_asym => [0],
        :accept_offset => [0],
        :accept_growth => [0],
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
## assumes correlation same among variables and common variance
## hmmm need to know how to work for matrix inputs...
function perturb_multi_normal!(vec, sd::Vector{Float64}, corr::Float64)
    d = length(vec)
    sd_length = length(sd)

    if sd_length == 1
        sd_mat = Diagonal(repeat(sd, d))
    elseif sd_length == d
        sd_mat = Diagonal(sd)
    else
        error("improper sd length given")
    end

    corr_mat = fill(corr, d, d)
    corr_mat[diagind(corr_mat)] .= 1
    covar_mat = sd_mat * corr_mat * sd_mat 
    pd_mat = PDMat(covar_mat)

    prop = MvNormal(zeros(d), pd_mat) # proposal distribution
    vec += rand(prop, 1)
    return(vec)
end

function draw_asym!(x0, dat, σ=0.1, ρ::Float64 = 0.1, track_accept = true)
    # proposed = perturb_normal!(copy(x0[:asym]), σ)
    if ρ == 0
        proposed = perturb_normal!(copy(x0[:asym]), σ)
    else 
        proposed = perturb_multi_normal!(copy(x0[:asym]), σ, ρ)
    end
    ## fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat, lwr, upr, lwr_lim, upr_lim)
    f0 = fc_gparam(x0[:asym], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                   x0[:beta_asym], x0[:X_asym], x0[:latent_sd_asym], dat, 18000, 22000, 0, Inf)

    f1 = fc_gparam(proposed, proposed, x0[:offset], x0[:growth], x0[:obs_sd], 
                    x0[:beta_asym], x0[:X_asym], x0[:latent_sd_asym], dat,  18000, 22000, 0, Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:asym] = proposed
        if track_accept
            x0[:accept_asym] = x0[:accept_asym] .+ 1
        end
    end
end

function draw_offset!(x0, dat, σ=0.1,  ρ::Float64 = 0.1, track_accept = true)
    # proposed = perturb_normal!(copy(x0[:offset]), σ)
    if ρ == 0
        proposed = perturb_normal!(copy(x0[:offset]), σ)
    else 
        proposed = perturb_multi_normal!(copy(x0[:offset]), σ, ρ)
    end
    ## fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat)
    f0 = fc_gparam(x0[:offset], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                   x0[:beta_offset], x0[:X_offset], x0[:latent_sd_offset], dat, 1.5, 2.5, 0, Inf)

    f1 = fc_gparam(proposed, x0[:asym], proposed, x0[:growth], x0[:obs_sd], 
                    x0[:beta_offset], x0[:X_offset], x0[:latent_sd_offset], dat, 1.5, 2.5, 0, Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:offset] = proposed
        if track_accept
            x0[:accept_offset] = x0[:accept_offset] .+ 1
        end
    end
end

function draw_growth!(x0, dat, σ=0.1,  ρ::Float64 = 0, track_accept = true)
    # proposed = perturb_normal!(copy(x0[:growth]), σ)
    if ρ == 0
        proposed = perturb_normal!(copy(x0[:growth]), σ)
    else
        proposed = perturb_multi_normal!(copy(x0[:growth]), σ, ρ)
    end
    ## fc_gparam(gparam, asym, offset, growth, obs_sd, beta, X, latent_sd, dat)
    f0 = fc_gparam(x0[:growth], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
                   x0[:beta_growth], x0[:X_growth], x0[:latent_sd_growth], dat, 0.8, 0.95, 0, 1)

    f1 = fc_gparam(proposed, x0[:asym], x0[:offset], proposed, x0[:obs_sd], 
                    x0[:beta_growth], x0[:X_growth], x0[:latent_sd_growth], dat, 0.8, 0.95, 0, 1)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:growth] = proposed
        if track_accept
            x0[:accept_growth] = x0[:accept_growth] .+ 1
        end
    end
end

### Drawing 
### Be sure to change priors appropriately
function draw_beta_asym!(x0, σ=0.1, track_accept = true)
    proposed = perturb_normal!(copy(x0[:beta_asym]), σ)
    # fc_beta(beta, X, gparam, latent_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1, lwr_lim = 0, upr_lim = Inf)
    f0 = fc_beta(x0[:beta_asym], x0[:X_asym],
                x0[:asym],
                x0[:latent_sd_asym],
                20000,
                10, 
                10,
                0,
                Inf)

    f1 = fc_beta(proposed, x0[:X_asym], x0[:asym],
                        x0[:latent_sd_asym],
                        20000,
                        10, 
                        10,
                        0,
                        Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_asym] = proposed
        if track_accept
            x0[:accept_beta_asym] = x0[:accept_beta_asym] .+ 1
        end
    end
end

### Regression coefficients for imputation models

function draw_beta_offset!(x0, σ=0.1, track_accept = true)
    proposed = perturb_normal!(copy(x0[:beta_offset]), σ)
    # fc_beta(beta, X, gparam, latent_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1, lwr_lim = 0, upr_lim = Inf)
    f0 = fc_beta(x0[:beta_offset], x0[:X_offset],
                x0[:offset],
                x0[:latent_sd_offset],
                2,
                0.2, 
                0.2,
                0,
                Inf)

    f1 = fc_beta(proposed, x0[:X_offset],
                x0[:offset],
                x0[:latent_sd_offset],
                2,
                0.2, 
                0.2,
                0,
                Inf)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_offset] = proposed
        if track_accept
            x0[:accept_beta_offset] = x0[:accept_beta_offset] .+ 1
        end
    end
end

function draw_beta_growth!(x0, σ=0.1, track_accept = true)
    proposed = perturb_normal!(copy(x0[:beta_growth]), σ)
    # fc_beta(beta, X, gparam, latent_sd, prior_mean = 0, prior_intercept_sd = 1, prior_sd = 1, lwr_lim = 0, upr_lim = Inf)
    f0 = fc_beta(x0[:beta_growth], x0[:X_growth],
                x0[:growth],
                x0[:latent_sd_growth],
                0.85,
                0.03, 
                0.03,
                0,
                1)

    f1 = fc_beta(proposed, x0[:X_growth],
                x0[:growth],
                        x0[:latent_sd_growth],
                        0.85,
                        0.03, 
                        0.03,
                        0,
                        1)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_growth] = proposed
        if track_accept
            x0[:accept_beta_growth] = x0[:accept_beta_growth] .+ 1
        end
    end
end


### latent standard deviations
function draw_latent_sd_asym!(x0, σ = 0.01, track_accept = true)
    prop = perturb_normal!(copy(x0[:latent_sd_asym]), σ)
    # fc_latent_sd(latent_sd, beta, X, gp_param, prior_mean = 0, prior_sd = 1, lwr_lim = -Inf, upr_lim = Inf)
    f0 = fc_latent_sd(x0[:latent_sd_asym], x0[:beta_asym], x0[:X_asym], x0[:asym], 0, 1, 0, Inf)
    f1 = fc_latent_sd(prop, x0[:beta_asym], x0[:X_asym], x0[:asym], 0, 1, 0, Inf)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_asym] = prop
        if track_accept
            x0[:accept_latent_sd_asym] = x0[:accept_latent_sd_asym] .+ 1
        end
    end
end

function draw_latent_sd_offset!(x0, σ = 0.01, track_accept = true)
    prop = perturb_normal!(copy(x0[:latent_sd_offset]), σ)
    ##fc_latent_sd(latent_sd, beta, X, gp_param, prior_mean = 0, prior_sd = 1, lwr_lim = -Inf, upr_lim = Inf)
    f0 = fc_latent_sd(x0[:latent_sd_offset], x0[:beta_offset], x0[:X_offset], x0[:offset],  0, 0.1, 0, Inf)
    f1 = fc_latent_sd(prop, x0[:beta_offset], x0[:X_offset], x0[:offset], 0, 0.1, 0, Inf)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_offset] = prop
        if track_accept
            x0[:accept_latent_sd_offset] = x0[:accept_latent_sd_offset] .+ 1
        end
    end
end

function draw_latent_sd_growth!(x0, σ = 0.01, track_accept = true)
    prop = perturb_normal!(copy(x0[:latent_sd_growth]), σ)
    ##fc_latent_sd(latent_sd, beta, X, outcome, prior_sd = 1)
    f0 = fc_latent_sd(x0[:latent_sd_growth], x0[:beta_growth], x0[:X_growth], x0[:growth], 0,  0.1, 0, 1)
    f1 = fc_latent_sd(prop, x0[:beta_growth], x0[:X_growth], x0[:growth], 0, 0.1, 0, 1)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:latent_sd_growth] = prop
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

function draw_obs_sd!(x0, dat, σ = 0.1, track_accept = true)
    prop = perturb_normal!(copy(x0[:obs_sd]), σ)
    ##fc_obs_sd(obs_sd, asym, offset, growth, dat, prior_sd = 1)
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
            :beta_asym,
            :beta_offset,
            :beta_growth,
            :latent_sd_asym,
            :latent_sd_offset,
            :latent_sd_growth,
            :obs_sd]

accept_keys = [:accept_asym ,
:accept_offset ,
:accept_growth ,
:accept_beta_asym ,
:accept_beta_offset ,
:accept_beta_growth ,
:accept_latent_sd_asym ,
:accept_latent_sd_offset ,
:accept_latent_sd_growth ,
:accept_obs_sd ]

function run_mcmc_chain(fname, x0, monitors, dat, warmup, run, thin=10, verbose = 1)
    if (!isfinite(fc_gparam(x0[:asym], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
        x0[:beta_asym], x0[:X_asym], x0[:latent_sd_asym], dat, 18000, 22000, 0, Inf)
))
    error("asym full conditional not fininte at init")
    end
    if (!isfinite(fc_gparam(x0[:offset], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
        x0[:beta_offset], x0[:X_offset], x0[:latent_sd_offset], dat, 1.5, 2.5, 0, Inf)))
    error("offset full conditional not fininte at init")
    end

    if (!isfinite(fc_gparam(x0[:growth], x0[:asym], x0[:offset], x0[:growth], x0[:obs_sd], 
        x0[:beta_growth], x0[:X_growth], x0[:latent_sd_growth], dat, 0.8, 0.95, 0, 1)))
    error("growth full conditional not fininte at init")
    end


    if (!isfinite(fc_beta(x0[:beta_asym], x0[:X_asym], 
        x0[:asym],    
        x0[:latent_sd_asym],
        20000,
        1, 
        1,
        0,
        Inf)))
    error("Beta asym full conditional not finite at init")
    end

    if (!isfinite(fc_beta(x0[:beta_offset], x0[:X_offset],
        x0[:offset],
        x0[:latent_sd_offset],
        log(0.8),
        1, 
        1,
        0,
        Inf)))
    error("Beta offset full conditional not finite at init")
    end

    if (!isfinite(fc_beta(x0[:beta_growth], x0[:X_growth],
        x0[:growth],
        x0[:latent_sd_growth],
        logit(0.8),
        1, 
        1,
        0,
        1)))
    error("Beta growth full conditional not finite at init")
    end
    

    if (!isfinite(fc_latent_sd(x0[:latent_sd_asym], x0[:beta_asym], x0[:X_asym], x0[:asym], 500, 300, 0, Inf)))
        error("asym latent sd full conditional not finite at init")
    end

    if (!isfinite(fc_latent_sd(x0[:latent_sd_offset], x0[:beta_offset], x0[:X_offset], x0[:offset], 0, 1, 0, Inf)))
        error("offset latent sd full conditional not finite at init")
    end

    if (!isfinite(fc_latent_sd(x0[:latent_sd_growth], x0[:beta_growth], x0[:X_growth], x0[:growth], 0, 1, 0, 1)))
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
    ### best params for untransformed with these priors
    # perturb_sd = Dict(:asym => [0.22], #asym
    #              :offset => [0.0014], #offset
    #              :growth => [0.00007], #growth
    #              :beta_asym => 0.35, #beta_asym
    #              :beta_offset => 0.004,  #beta_offset
    #              :beta_growth => 0.005, #beta_growth
    #              :latent_sd_asym => 0.25, #latent_sd_asym
    #              :latent_sd_offset => 0.003, #latent_sd_offset
    #              :latent_sd_growth => 0.003, #latent_sd_growth
    #              :obs_sd => 5 #obs_sd
    # )

    perturb_sd = Dict(:asym => [0.24], #asym
                 :offset => [0.003], #offset
                 :growth => [0.0004], #growth
                 :beta_asym => 0.35, #beta_asym
                 :beta_offset => 0.005,  #beta_offset
                 :beta_growth => 0.005, #beta_growth
                 :latent_sd_asym => 0.25, #latent_sd_asym
                 :latent_sd_offset => 0.005, #latent_sd_offset
                 :latent_sd_growth => 0.003, #latent_sd_growth
                 :obs_sd => 20 #obs_sd
    )
    
    perturb_corr = Dict(:asym => 0.1, #asym
                 :offset => 0.15, #offset
                 :growth => 0.1, #growth
    )

    for i in 1:warmup
            
        if verbose > 0
            if i % 100 == 0
                print(string(i)*"/"*string(warmup)*"\n")
            end            
        end
        # draw_asym!(x0, dat, perturb_sd[:asym]) ### 
        # draw_offset!(x0, dat, perturb_sd[:offset])
        # draw_growth!(x0, dat, perturb_sd[:growth])

        # draw_beta_asym!(x0, perturb_sd[:beta_asym])
        # draw_beta_offset!(x0, perturb_sd[:beta_offset])
        # draw_beta_growth!(x0, perturb_sd[:beta_growth])

        # draw_latent_sd_asym!(x0, perturb_sd[:latent_sd_asym])
        # draw_latent_sd_offset!(x0, perturb_sd[:latent_sd_offset])
        # draw_latent_sd_growth!(x0, perturb_sd[:latent_sd_growth])

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
        # draw_asym!(x0, dat, perturb_sd[:asym], perturb_corr[:asym], true)
        draw_offset!(x0, dat, perturb_sd[:offset], perturb_corr[:offset], true)
        # draw_growth!(x0, dat, perturb_sd[:growth], perturb_corr[:growth], true)

        # draw_beta_asym!(x0, perturb_sd[:beta_asym], true)
        # draw_beta_offset!(x0, perturb_sd[:beta_offset], true)
        # draw_beta_growth!(x0, perturb_sd[:beta_growth], true)

        # draw_latent_sd_asym!(x0, perturb_sd[:latent_sd_asym], true)
        # draw_latent_sd_offset!(x0, perturb_sd[:latent_sd_offset], true)
        # draw_latent_sd_growth!(x0, perturb_sd[:latent_sd_growth], true)

        # draw_obs_sd!(x0, dat, perturb_sd[:obs_sd], true)
    end

    ## output acceptance ratios
    print("\n")
    print("Acceptance Asym: ", string(x0[:accept_asym]/run), "\n")
    print("Acceptance offset: ", string(x0[:accept_offset]/run), "\n")
    print("Acceptance growth: ", string(x0[:accept_growth]/run), "\n")

    print("Acceptance beta_asym: ", string(x0[:accept_beta_asym]/run), "\n")
    print("Acceptance beta_offset: ", string(x0[:accept_beta_offset]/run), "\n")
    print("Acceptance beta_growth: ", string(x0[:accept_beta_growth]/run), "\n")

    print("Acceptance latent_sd_asym: ", string(x0[:accept_latent_sd_asym]/run), "\n")
    print("Acceptance latent_sd_offset: ", string(x0[:accept_latent_sd_offset]/run), "\n")
    print("Acceptance latent_sd_growth: ", string(x0[:accept_latent_sd_growth]/run), "\n")
    print("Acceptance obs_sd: ", string(x0[:accept_obs_sd]/run), "\n")
end




    