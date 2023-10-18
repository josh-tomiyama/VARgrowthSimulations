using Distributions
using DataFrames
using XLSX
using BenchmarkTools

###############################################################
## categorical logit and partial categorical logit functions ##
###############################################################

function rPartialCatLogit(eta, nonzero_idx)
    p_last = 1/(1+sum(exp.(eta)))
    p = vcat(exp.(eta).*p_last,p_last)
    psub = [p[x] for x in nonzero_idx]
    psub = psub ./ sum(psub)
    return(nonzero_idx[rand(Categorical(psub))])
end

function dPartialCatLogit(x, eta, nonzero_idx, log=true)
    if (!(x in nonzero_idx))
        return(log ? -Inf : 0.0)
    end
    p_last = 1/(1+sum(exp.(eta)))
    p = vcat(exp.(eta).*p_last,p_last)
    if (any(.!isfinite.(p)))
        return(log ? -Inf : 0.0)
    end
    psub = [p[x] for x in nonzero_idx]
    psub = psub ./ sum(psub)
    xsub = findall(x .== nonzero_idx)[1]
    if (log)
        return(logpdf(Categorical(psub), xsub))
    end
    return(pdf(Categorical(psub), xsub))
end

function rCatLogit(eta)
    p_last = 1/(1+sum(exp.(eta)))
    p = vcat(exp.(eta).*p_last, p_last)
    return(rand(Categorical(p)))
end

function dCatLogit(x, eta, log=true)
    if (x < 1 || x > length(eta) + 1)
        return(log ? -Inf : 0.0)
    end
    p_last = 1/(1+sum(exp.(eta)))
    p = vcat(exp.(eta).*p_last,p_last)
    if (any(.!isfinite.(p)))
        return(log ? -Inf : 0.0)
    end
    if (log)
        return(logpdf(Categorical(p), x))
    end
    return(pdf(Categorical(p), x))
end

###################################
## Functions to Apply Imputation ##
###################################
(sqrt \circ)
function apply_single_impute!(X, row_idx, col_idx, val)
    X[row_idx, col_idx] = val
end
# Function to apply an interaction to certain rows of a design matrix
function apply_single_interact!(X, row_idx, col_idx, out_idx)
    X[row_idx, out_idx] = reduce(.*,[X[row_idx, i] for i in col_idx])
end

####################
## Format DataSet ##
####################

function gomp_data(dat, time_nm, year_nm)
    out = copy(dat)
    out.time = dat[:, time_nm]
    out.year = dat[:, year_nm]
    select!(out, Not(time_nm, year_nm))
    sort!(out, year, time)
end

function gomp_data!(dat, time_nm, year_nm)
    dat.time = dat[:, time_nm]
    dat.year = dat[:, year_nm]
    select!(dat, Not(time_nm, year_nm))
    sort!(dat, year, time)
end

######################
## Prior Definition ##
######################
### should be a names list that maps priors to param names
### It should also be expected that the params for evalutaiton are in a dataframe/named array

function assign_priors(dat, lik_func, )

#################
## Likelihoods ##
#################

### Likelihoods for observed process
### gompertz for scalar params asym offset growth
function gompertz_likelihood(dat, asym, offset, growth, obs_sd, log_lik = true)
    out = 0
    d = Normal(0, obs_sd) ## sd assumed to be scalar
    mu = asym*exp(growth.*dat.time)

    if log_lik
        out += sum(logpdf.(d, mu))
    else
        out += sum(pdf.(d, mu))
    end


end

### Latent likelihoods
### ideally: assign independently for each parameter of interest


function latent_liklihood_linear(param, beta, model_matrix, sd, log_lik=true)
    d = Normal(model_matrix * beta, sd)
    if log_lik
        out += sum(logpdf.(d, mu))
    else
        out += sum(logpdf.(d, mu))
    end
end

function latent_likelihood_linear(param_name, dat, all_params::DataFrame, model_matrix, log_lik = true)
    req_params = ["", "beta_", "sd_"] .* param_name
    if !( issubset(req_params, names(all_params)) )
        ErrorException("Required Parameters " * join(req_params, ",") * " not in param names")
    end
end



################################
## Full Conditional Functions ##
################################

## General full conditional for the beta parameter sets linked to
## multinomial likelihoods
function fc_beta(beta, X_imp, y, idx, prior_intercept_sd = 1, prior_sd = 1)
    eta = X_imp * beta
    out = 0.0
    # Likelihood
    # Liklihood only applies to a subset
    for i in idx
        out += dCatLogit(y[i], eta[i,:])
    end
    # Prior
    d1 = Normal(0, prior_intercept_sd)
    d2 = Normal(0, prior_sd)
    out += sum(logpdf.(d1, beta[1,:]))
    out += sum(logpdf.(d2, beta[2:(size(beta)[1]),:]))
    return(out) 
end

## Full conditional for beta_condition
function fc_beta_condition(beta_condition, X_imp_condition, y_imp_condition, 
                           prior_intercept_sd = 10, prior_sd = 3)
    eta = X_imp_condition*beta_condition
    out = 0.0
    
    # Likelihood
    for i in 1:(size(eta)[1])
        out += dCatLogit(y_imp_condition[i], eta[i,:])                       
    end
    # Prior
    d1 = Normal(0, prior_intercept_sd)
    d2 = Normal(0, prior_sd)
    out += sum(logpdf.(d1, beta_condition[1,:]))
    out += sum(logpdf.(d2, beta_condition[2:(size(beta_condition)[2]),]))
end

## Full conditional for diagnostic category
## Note: condition defined only for idx indices, the rest are observed
function fc_condition(condition, y_mRS90, idx,  X_condition, beta_condition,
                       eta_mRS_LVO, 
                       eta_mRS_Non_LVO,
                       eta_mRS_Hem,
                       eta_mRS_Mim)
    if (any(condition[idx] .> 2))
        error("Impossible condition value observed.")
    end
    ## Condition likelihoods
    eta_cond = X_condition[idx,:] * beta_condition
    out = sum([dPartialCatLogit(condition[idx[i]],  
                               eta_cond[i,:], 
                               1:2) for i in eachindex(idx)]) # LVO or non-LVO ischm)
    
    ## mRS likelihood (depends on condition)
    eta_mRS = ((condition[idx] .== 1) .* eta_mRS_LVO[idx,:] + 
               (condition[idx] .== 2) .* eta_mRS_Non_LVO[idx,:] + 
               (condition[idx] .== 3) .* eta_mRS_Hem[idx,:] + 
               (condition[idx] .== 4) .* eta_mRS_Mim[idx,:])
    #out += sum(dCatLogit.(y_mRS90[idx], eta_mRS[idx]))

    out += sum([dCatLogit(y_mRS90[idx[i]],  
                   eta_mRS[i,:]) for i in eachindex(idx)])
    return(out)
end

## Full conditionals for variables needing imputation
function fc_imp_RACE!(RACE, RACE_full, idx, 
                      condition,
                      BNIH_full,
                      ischemic,
                      X_condition, beta_condition,
                      X_imp_BNIH, beta_imp_BNIH,
                      X_imp_RACE, beta_imp_RACE,
                      sigma_RACE,
                      sigma_BNIH)
    
    # Apply imputations                   
    RACE_full[idx] .= RACE
    apply_single_impute!(X_condition, idx, 6, RACE)
    apply_single_impute!(X_imp_BNIH, idx, 4, RACE)

    eta_cond = X_condition * beta_condition
    eta_RACE = X_imp_RACE * beta_imp_RACE
    eta_BNIH = X_imp_BNIH * beta_imp_BNIH
    
    # Likeihood component - condition
    out = sum([dPartialCatLogit(condition[i],  
                               eta_cond[i,:], 
                               ischemic[i] ? (1:2) : (1:4)) 
                    for i in axes(eta_cond)[1]])
    

    # Likelihood component - RACE imputation model
    out += sum(logpdf.(Normal.(eta_RACE, sigma_RACE), RACE_full))

    # Likelihood component - BNIH imputation model
    out += sum(logpdf.(Normal.(eta_BNIH, sigma_BNIH), BNIH_full))
    
    return(out)
end


#fc_imp_RACE!(RACE, RACE_full, idx, 
#condition,
#BNIH_full,
#ischemic,
#X_condition, beta_condition,
#X_imp_BNIH, beta_imp_BNIH,
#X_imp_RACE, beta_imp_RACE,
#sigma_RACE,
#sigma_BNIH)

function fc_imp_BNIH!(BNIH, BNIH_full, idx, 
                      condition, y_mRS90,
                      RACE_full,
                      X_mRS_LVO, beta_mRS_LVO,
                      X_mRS_Non_LVO, beta_mRS_Non_LVO,
                      X_mRS_Hem, beta_mRS_Hem,
                      X_mRS_Mim, beta_mRS_Mim,
                      X_imp_BNIH, beta_imp_BNIH,
                      X_imp_RACE, beta_imp_RACE,
                      sigma_RACE,
                      sigma_BNIH)

                      
    # Apply imputations                   
    BNIH_full[idx] .= BNIH
    apply_single_impute!(X_imp_RACE, idx, 4, BNIH)
    apply_single_interact!(X_imp_RACE, 1:(size(X_imp_RACE)[1]), [4,4], 5) # BNIH^2

    apply_single_impute!(X_mRS_LVO, idx, 4, BNIH)
    apply_single_impute!(X_mRS_Non_LVO, idx, 4, BNIH)
    apply_single_impute!(X_mRS_Hem, idx, 4, BNIH)
    apply_single_impute!(X_mRS_Mim, idx, 4, BNIH)

    eta_RACE = X_imp_RACE * beta_imp_RACE
    eta_BNIH = X_imp_BNIH * beta_imp_BNIH
    
    eta_mRS_LVO = X_mRS_LVO * beta_mRS_LVO
    eta_mRS_Non_LVO = X_mRS_Non_LVO * beta_mRS_Non_LVO
    eta_mRS_Hem = X_mRS_Hem * beta_mRS_Hem
    eta_mRS_Mim = X_mRS_Mim * beta_mRS_Mim

    ## mRS likelihood (depends on condition)
    eta_mRS = ((condition[idx] .== 1) .* eta_mRS_LVO[idx,:] + 
       (condition[idx] .== 2) .* eta_mRS_Non_LVO[idx,:] + 
       (condition[idx] .== 3) .* eta_mRS_Hem[idx,:] + 
       (condition[idx] .== 4) .* eta_mRS_Mim[idx,:])

    out = sum([dCatLogit(y_mRS90[idx[i]],  
           eta_mRS[i,:]) for i in eachindex(idx)])

    # Likelihood component - RACE imputation model
    out += sum(logpdf.(Normal.(eta_RACE, sigma_RACE), RACE_full))

    # Likelihood component - BNIH imputation model
    out += sum(logpdf.(Normal.(eta_BNIH, sigma_BNIH), BNIH_full))
    
    return(out)
    
end

function fc_imp_SBP_DBP!(sysBP, diasBP, sysBP_full, daisBP_full,
                        idx_sbp, idx_dbp, idx_bp,
                        X_condition, beta_condition, isch, 
                        condition,
                        X_imp_sysBP, beta_imp_sysBP,
                        X_imp_diasBP, beta_imp_diasBP,
                        sigma_SBP,
                        sigma_DBP)


    sysBP_full[idx_sbp] .= sysBP
    daisBP_full[idx_dbp] .= diasBP
    ### sysBP Model
    apply_single_impute!(X_imp_sysBP, idx_dbp, 4, diasBP)
    ### diasBP Model
    apply_single_impute!(X_imp_diasBP,  idx_sbp, 4, sysBP)
    ### Imputation for condition model
    apply_single_impute!(X_condition, idx_sbp, 2, sysBP)
    apply_single_impute!(X_condition, idx_dbp, 3, diasBP)
    apply_single_interact!(X_condition, idx_bp, [2,3], 7)

    eta_imp_sysBP = X_imp_sysBP*beta_imp_sysBP
    eta_imp_diasBP = X_imp_diasBP*beta_imp_diasBP
    eta_condition  = X_condition*beta_condition
    
    ## Likelihood - condition
    out = sum([dPartialCatLogit(condition[i],  
                                eta_condition[i,:], 
                                isch[i] ? (1:2) : (1:4))
                                for i in axes(X_condition)[1]])
    
    ## Likelihood component - SBP imputation model
    out += sum(logpdf.(Normal.(eta_imp_sysBP, sigma_SBP), sysBP_full))

    
    ## Likelihood component - DBP imputation model
    out += sum(logpdf.(Normal.(eta_imp_diasBP, sigma_DBP), daisBP_full))
    

    return(out)
end

## Generic full conditional for imputation model betas
function fc_beta_imp(beta, X, sigma_imp, y, sigma_beta)
    # Likelihood
    eta = X*beta
    out = sum(logpdf.(Normal.(eta, sigma_imp), y))

    # Prior:
    out += sum([logpdf(Normal(0, sigma_beta[i]), beta[i]) for i in eachindex(beta)])
    return(out)
end

## Generic full conditional for imputation model sigmas
function fc_sigma_imp(sigma_imp, X, beta, y, prα, prβ)
    if (any(sigma_imp .<= 0))
        return(-Inf)
    end
    # Likelihood
    eta = X*beta
    out = sum(logpdf.(Normal.(eta, sigma_imp), y))

    # Prior:
    out += logpdf(Gamma(prα, prβ), sigma_imp[1])
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
function generate_init_state(data)
    ## Define the number of conditions
    N_condition = 4
    y_mRS90_actual = parse.(Int64, fulldata[!,:mRS90])
    
    ## Create example X matrix for condition, including missing values
    X_condition_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 7)
    X_condition_obs[:,1] .= 1
    X_condition_obs[:,2] .= fulldata[:,:sysBP]
    X_condition_obs[:,3] .= fulldata[:,:diasBP]
    X_condition_obs[:,4] .= fulldata[:,:age]
    X_condition_obs[:,5] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
    X_condition_obs[:,6] .= fulldata[:,:RACE_scale]
    X_condition_obs[:,7] .= fulldata[:,:sysBP] .*fulldata[:,:diasBP]

    
    ## Create model matrices for mRS outcomes, separately for each condition
    ## To do this, we create a "full" matrix with all possible terms and then perform column-wise subsets
    X_mRS_Full_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 13)
    X_mRS_Full_obs[:,1] .= 1
    X_mRS_Full_obs[:,2] .= fulldata[:,:age]
    X_mRS_Full_obs[:,3] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
    X_mRS_Full_obs[:,4] .= fulldata[:,:BNIH]
    X_mRS_Full_obs[:,5] .= fulldata[:,:TPA]
    X_mRS_Full_obs[:,6] .= fulldata[:,:TPATime]
    X_mRS_Full_obs[:,7] .= fulldata[:,:MT]
    X_mRS_Full_obs[:,8] .= fulldata[:,:TPA] .* fulldata[:,:MT]
    X_mRS_Full_obs[:,9] .= fulldata[:,:TPA] .* fulldata[:,:MT] .* fulldata[:,:MTTime]
    X_mRS_Full_obs[:,10] .= fulldata[:,:TPA] .* fulldata[:,:MT] .* fulldata[:,:TPATime]
    X_mRS_Full_obs[:,11] .= fulldata[:,:TPA] .* (.!fulldata[:,:MT]) .* fulldata[:,:TPATime]
    X_mRS_Full_obs[:,12] .= (.!fulldata[:,:TPA]) .* (fulldata[:,:MT]) .* fulldata[:,:MTTime]
    X_mRS_Full_obs[:,13] .= fulldata[:,:TimeToHosp]

    X_mRS_LVO_obs = X_mRS_Full_obs[:, [1,2,3,4,5,7,8,9,10,11,12]]
    X_mRS_Non_LVO_obs = X_mRS_Full_obs[:, [1,2,3,4,5,6]]
    X_mRS_Hem_obs = X_mRS_Full_obs[:, [1,2,3,4,13]]
    X_mRS_Mim_obs = X_mRS_Full_obs[:, [1,2,3,4]]

    ## Create model matrices for imputation 

    # From R: RACEformula <- formula("RACE_scale ~ age + sex + BNIH + I(BNIH^2)")
    X_imp_RACE_obs = X_mRS_Full_obs[:,[1,2,3,4,4]]
    apply_single_interact!(X_imp_RACE_obs, 1:(size(X_imp_RACE_obs)[1]), [4,4], 5) # BNIH^2

    # From R: BNIHformula <- formula("BNIH ~ age + sex + RACE_scale")
    X_imp_BNIH_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 4)
    X_imp_BNIH_obs[:,1] .= 1
    X_imp_BNIH_obs[:,2] .= fulldata[:,:age]
    X_imp_BNIH_obs[:,3] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
    X_imp_BNIH_obs[:,4] .= fulldata[:,:RACE_scale]

    # From R: sysBPformula <- formula("sysBP ~ age + sex + diasBP")
    X_imp_sysBP_obs = X_mRS_Full_obs[:, [1,2,3,3]]
    X_imp_sysBP_obs[:,4] .= fulldata[:,:diasBP]

    # From R: diasBPformula <- formula("diasBP ~ age + sex + sysBP")
    X_imp_diasBP_obs = X_mRS_Full_obs[:, [1,2,3,3]]
    X_imp_diasBP_obs[:,4] .= fulldata[:,:sysBP]

    impute_idx_BNIH = findall(ismissing.(X_mRS_LVO_obs[:,:4]))
    impute_idx_sbp = findall(ismissing.(X_condition_obs[:,:2]))
    impute_idx_dbp = findall(ismissing.(X_condition_obs[:,:3]))
    impute_idx_bp = sort(union(impute_idx_sbp, impute_idx_dbp))
    impute_idx_race = findall(ismissing.(X_condition_obs[:,:6]))

    ## Beta inits
    beta_condition  = rand(Normal(0,0.1), (size(X_condition_obs)[2], N_condition-1))
    beta_mRS_LVO = generate_init_betas(X_mRS_LVO_obs, 7, 0.5)
    beta_mRS_Non_LVO = generate_init_betas(X_mRS_Non_LVO_obs, 7, 0.5)
    beta_mRS_Hem = generate_init_betas(X_mRS_Hem_obs, 7, 0.5)
    beta_mRS_Mim = generate_init_betas(X_mRS_Mim_obs, 7, 0.5)

    ## Imputation value inits
    BNIH_missing = rand(Normal(10,1), length(impute_idx_BNIH))
    sysBP_missing = rand(Normal(0.5,1), length(impute_idx_sbp))
    diasBP_missing = rand(Normal(0.5,1), length(impute_idx_dbp))
    RACE_missing = rand(Normal(0.5,1), length(impute_idx_race))
    
    ## Design matrix inits based on impute inits
    X_condition = deepcopy(X_condition_obs)
    X_mRS_LVO = deepcopy(X_mRS_LVO_obs)
    X_mRS_Non_LVO = deepcopy(X_mRS_Non_LVO_obs)
    X_mRS_Hem = deepcopy(X_mRS_Hem_obs)
    X_mRS_Mim = deepcopy(X_mRS_Mim_obs)

    apply_single_impute!(X_mRS_LVO, impute_idx_BNIH, 4, BNIH_missing)
    apply_single_impute!(X_mRS_Non_LVO, impute_idx_BNIH, 4, BNIH_missing)
    apply_single_impute!(X_mRS_Hem, impute_idx_BNIH, 4, BNIH_missing)
    apply_single_impute!(X_mRS_Mim, impute_idx_BNIH, 4, BNIH_missing)

    apply_single_impute!(X_condition, impute_idx_sbp, 2, sysBP_missing)
    apply_single_impute!(X_condition, impute_idx_dbp, 3, diasBP_missing)
    apply_single_impute!(X_condition, impute_idx_race, 6, RACE_missing)
    apply_single_interact!(X_condition, impute_idx_bp, [2,3], 7)

    ## Imputation model design matrices
    X_imp_RACE = deepcopy(X_imp_RACE_obs)
    X_imp_BNIH = deepcopy(X_imp_BNIH_obs)
    X_imp_sysBP = deepcopy(X_imp_sysBP_obs)
    X_imp_diasBP = deepcopy(X_imp_diasBP_obs)
    
    beta_imp_RACE = generate_init_betas(X_imp_RACE, 2, 0.5)
    beta_imp_BNIH = generate_init_betas(X_imp_BNIH, 2, 0.5)
    beta_imp_sysBP = generate_init_betas(X_imp_sysBP, 2, 0.5)
    beta_imp_diasBP = generate_init_betas(X_imp_diasBP, 2, 0.5)

    ## Fill in imputed values for imputation models where needed
    ### RACE Model
    apply_single_impute!(X_imp_RACE, impute_idx_BNIH, 4, BNIH_missing)
    apply_single_interact!(X_imp_RACE, 1:(size(X_imp_RACE)[1]), [4,4], 5) # BNIH^2

    ### BNIH Model
    apply_single_impute!(X_imp_BNIH, impute_idx_race, 4, RACE_missing)

    ### sysBP Model
    apply_single_impute!(X_imp_sysBP, impute_idx_dbp, 4, diasBP_missing)

    ### diasBP Model
    apply_single_impute!(X_imp_diasBP,  impute_idx_sbp, 4, sysBP_missing)


    condition = [ismissing(fulldata[i,:type]) ? 0 : 
            fulldata[i,:type] == "LVO" ? 1 : 
            fulldata[i,:type] == "Not LVO" ? 2 :
            fulldata[i,:type] == "Hemorrhagic" ? 3 :
            fulldata[i,:type] == "Stroke mimic" ? 4 : 5 for i in eachindex(y_mRS90_actual)]

    sysBP_full = X_condition[:,2]
    diasBP_full = X_condition[:,3]
    BNIH_full = X_mRS_LVO[:,4]
    RACE_full = X_condition[:,6]
    sigma_RACE = rand(Gamma(1,1),1)
    sigma_BNIH = rand(Gamma(1,1),1)
    sigma_SBP = rand(Gamma(1,1),1)
    sigma_DBP = rand(Gamma(1,1),1)

    inits = Dict(:y_mRS90 => y_mRS90_actual,
        :condition => condition,
        :impute_idx_BNIH => impute_idx_BNIH,
        :impute_idx_sbp => impute_idx_sbp,
        :impute_idx_dbp => impute_idx_dbp,
        :impute_idx_bp => impute_idx_bp, 
        :impute_idx_race => impute_idx_race, 
        :ischemic => data[:, :Ischemic],
        :beta_condition => beta_condition, 
        :beta_mRS_LVO => beta_mRS_LVO,
        :beta_mRS_Non_LVO => beta_mRS_Non_LVO, 
        :beta_mRS_Hem => beta_mRS_Hem, 
        :beta_mRS_Mim => beta_mRS_Mim,
        :beta_imp_RACE => beta_imp_RACE,
        :beta_imp_BNIH => beta_imp_BNIH,
        :beta_imp_sysBP => beta_imp_sysBP,
        :beta_imp_diasBP => beta_imp_diasBP,
        :BNIH_missing => BNIH_missing,
        :sysBP_missing => sysBP_missing,
        :diasBP_missing => diasBP_missing,
        :RACE_missing => RACE_missing,
        :X_condition => X_condition,
        :X_mRS_LVO => X_mRS_LVO, 
        :X_mRS_Non_LVO => X_mRS_Non_LVO,
        :X_mRS_Hem => X_mRS_Hem,
        :X_mRS_Mim => X_mRS_Mim,
        :X_imp_RACE => X_imp_RACE,
        :X_imp_BNIH => X_imp_BNIH,
        :X_imp_sysBP => X_imp_sysBP,
        :X_imp_diasBP => X_imp_diasBP,
        :sysBP_full => sysBP_full,
        :diasBP_full => diasBP_full,
        :BNIH_full => BNIH_full,
        :RACE_full => RACE_full,
        :sigma_RACE => sigma_RACE,
        :sigma_BNIH => sigma_BNIH,
        :sigma_SBP => sigma_SBP,
        :sigma_DBP => sigma_DBP
    )
    return(inits)
end


####################
## Draw Functions ##
####################

function perturb_normal!(vec, sd=0.5, proportion = 1)
    prop = Normal(0, sd)
    if (proportion == 1)
        vec += rand(prop, size(vec))
    else
        idx = Int.(sort(sample(1:length(vec), Int(floor(length(vec)*proportion)), replace=false)))
        vec[idx] += rand(prop, length(idx))
    end
    return(vec)
end

### Regression coefficients for imputation models

function draw_beta_imp_sysBP!(x0, σ=0.1)
    proposed = perturb_normal!(copy(x0[:beta_imp_sysBP]), σ)
    
    f0 = fc_beta_imp(x0[:beta_imp_sysBP], x0[:X_imp_sysBP], 
                     x0[:sigma_SBP], x0[:sysBP_full], 
                     [10 for i in eachindex(x0[:beta_imp_sysBP])])

    f1 = fc_beta_imp(proposed, x0[:X_imp_sysBP], 
                     x0[:sigma_SBP], x0[:sysBP_full], 
                     [10 for i in eachindex(x0[:beta_imp_sysBP])])
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_imp_sysBP] = proposed
    end
end

function draw_beta_imp_diasBP!(x0, σ = 0.1)
    proposed = perturb_normal!(copy(x0[:beta_imp_diasBP]), σ)
    
    f0 = fc_beta_imp(x0[:beta_imp_diasBP], x0[:X_imp_diasBP], 
                     x0[:sigma_DBP], x0[:diasBP_full], 
                     [10 for i in eachindex(x0[:beta_imp_diasBP])])

    f1 = fc_beta_imp(proposed, x0[:X_imp_diasBP], 
                     x0[:sigma_DBP], x0[:diasBP_full], 
                     [10 for i in eachindex(x0[:beta_imp_diasBP])])
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_imp_diasBP] = proposed
    end
end


function draw_beta_imp_RACE!(x0, σ=0.1)
    proposed = perturb_normal!(copy(x0[:beta_imp_RACE]), σ)
    
    f0 = fc_beta_imp(x0[:beta_imp_RACE], x0[:X_imp_RACE], 
                     x0[:sigma_RACE], x0[:RACE_full], 
                     [10 for i in eachindex(x0[:beta_imp_RACE])])

    f1 = fc_beta_imp(proposed, x0[:X_imp_RACE], 
                     x0[:sigma_RACE], x0[:RACE_full], 
                     [10 for i in eachindex(x0[:beta_imp_RACE])])
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_imp_RACE] = proposed
    end
end

function draw_beta_imp_BNIH!(x0, σ = 0.1)
    proposed = perturb_normal!(copy(x0[:beta_imp_BNIH]), σ)
    
    f0 = fc_beta_imp(x0[:beta_imp_BNIH], x0[:X_imp_BNIH], 
                     x0[:sigma_BNIH], x0[:BNIH_full], 
                     [10 for i in eachindex(x0[:beta_imp_BNIH])])

    f1 = fc_beta_imp(proposed, x0[:X_imp_BNIH], 
                     x0[:sigma_BNIH], x0[:BNIH_full], 
                     [10 for i in eachindex(x0[:beta_imp_BNIH])])
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_imp_BNIH] = proposed
    end
end

### Variance terms for imputation models
function draw_sigma_imp_terms!(x0, σ = 0.1)
    prop1 = perturb_normal!(copy(x0[:sigma_SBP]), σ)
    prop2 = perturb_normal!(copy(x0[:sigma_DBP]), σ)
    prop3 = perturb_normal!(copy(x0[:sigma_RACE]), σ)
    prop4 = perturb_normal!(copy(x0[:sigma_BNIH]), σ)
    
    f0 = fc_sigma_imp(x0[:sigma_SBP], x0[:X_imp_sysBP], x0[:beta_imp_sysBP], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(x0[:sigma_DBP], x0[:X_imp_diasBP], x0[:beta_imp_diasBP], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(x0[:sigma_RACE], x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(x0[:sigma_BNIH], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], x0[:sysBP_full], 1, 1)
    f1 = fc_sigma_imp(prop1, x0[:X_imp_sysBP], x0[:beta_imp_sysBP], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(prop2, x0[:X_imp_diasBP], x0[:beta_imp_diasBP], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(prop3, x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sysBP_full], 1, 1) +
        fc_sigma_imp(prop4, x0[:X_imp_BNIH], x0[:beta_imp_BNIH], x0[:sysBP_full], 1, 1)
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:sigma_SBP] = prop1
        x0[:sigma_DBP] = prop2
        x0[:sigma_RACE] = prop3
        x0[:sigma_BNIH] = prop4
    end
end

## Imputation - sysBP
function draw_SBP_DBP!(x0, σ = 0.1, prop = 0.25)
    current_SBP = copy(x0[:sysBP_missing])
    current_DBP = copy(x0[:diasBP_missing])

    proposed_SBP = perturb_normal!(copy(current_SBP), σ, prop)
    proposed_DBP = perturb_normal!(copy(current_DBP), σ, prop)

    f0 = fc_imp_SBP_DBP!(x0[:sysBP_missing], x0[:diasBP_missing], 
               x0[:sysBP_full], 
               x0[:diasBP_full], x0[:impute_idx_sbp],
               x0[:impute_idx_dbp], x0[:impute_idx_bp], 
               x0[:X_condition],x0[:beta_condition], 
               x0[:ischemic], x0[:condition], x0[:X_imp_sysBP], 
               x0[:beta_imp_sysBP], x0[:X_imp_diasBP], 
               x0[:beta_imp_diasBP], x0[:sigma_SBP], 
               x0[:sigma_DBP])

    f1 = fc_imp_SBP_DBP!(proposed_SBP, proposed_DBP, 
               x0[:sysBP_full], 
               x0[:diasBP_full], x0[:impute_idx_sbp],
               x0[:impute_idx_dbp], x0[:impute_idx_bp], 
               x0[:X_condition],x0[:beta_condition], 
               x0[:ischemic], x0[:condition], x0[:X_imp_sysBP], 
               x0[:beta_imp_sysBP], x0[:X_imp_diasBP], 
               x0[:beta_imp_diasBP], x0[:sigma_SBP], 
               x0[:sigma_DBP])
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:sysBP_missing] = proposed_SBP
        x0[:diasBP_missing] = proposed_DBP
    else 
        # reset - TODO find a way to optimize this out
        fc_imp_SBP_DBP!(current_SBP, current_DBP, 
               x0[:sysBP_full], 
               x0[:diasBP_full], x0[:impute_idx_sbp],
               x0[:impute_idx_dbp], x0[:impute_idx_bp], 
               x0[:X_condition],x0[:beta_condition], 
               x0[:ischemic], x0[:condition], x0[:X_imp_sysBP], 
               x0[:beta_imp_sysBP], x0[:X_imp_diasBP], 
               x0[:beta_imp_diasBP], x0[:sigma_SBP], 
               x0[:sigma_DBP])
    end
end

## Imputation - RACE
function draw_RACE!(x0, σ=0.1, prop = 0.25)
    current_RACE = copy(x0[:RACE_missing])
    proposed_RACE = perturb_normal!(copy(current_RACE), σ, prop) 
    
    f0 = fc_imp_RACE!(x0[:RACE_missing], 
                        x0[:RACE_full], x0[:impute_idx_race], 
                        x0[:condition],x0[:BNIH_full], x0[:ischemic], x0[:X_condition], 
                        x0[:beta_condition], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                        x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sigma_RACE], x0[:sigma_BNIH])

    f1 = fc_imp_RACE!(proposed_RACE, 
                        x0[:RACE_full], x0[:impute_idx_race], 
                        x0[:condition],x0[:BNIH_full], x0[:ischemic], x0[:X_condition], 
                        x0[:beta_condition], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                        x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sigma_RACE], x0[:sigma_BNIH])
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:RACE_missing] = proposed_RACE
    else 
        # reset - TODO find a way to optimize this out
        f0 = fc_imp_RACE!(x0[:RACE_missing], 
                        x0[:RACE_full], x0[:impute_idx_race], 
                        x0[:condition],x0[:BNIH_full], x0[:ischemic], x0[:X_condition], 
                        x0[:beta_condition], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                        x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sigma_RACE], x0[:sigma_BNIH])
    end
end

## Imputation - BNIH
function draw_BNIH!(x0, σ, prop)
    current_BNIH = copy(x0[:BNIH_missing])
    proposed_BNIH = perturb_normal!(copy(current_BNIH), σ, prop) ## TODO - double check proposal is sensible
    
    f0 = fc_imp_BNIH!(x0[:BNIH_missing], x0[:BNIH_full], x0[:impute_idx_BNIH], 
                    x0[:condition],x0[:y_mRS90],x0[:RACE_full],
                    x0[:X_mRS_LVO], x0[:beta_mRS_LVO], 
                    x0[:X_mRS_Non_LVO], x0[:beta_mRS_Non_LVO], 
                    x0[:X_mRS_Hem], x0[:beta_mRS_Hem], 
                    x0[:X_mRS_Mim], x0[:beta_mRS_Mim], 
                    x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                    x0[:X_imp_RACE],x0[:beta_imp_RACE], 
                    x0[:sigma_RACE], x0[:sigma_BNIH])

    f1 = fc_imp_BNIH!(proposed_BNIH, x0[:BNIH_full], x0[:impute_idx_BNIH], 
                    x0[:condition],x0[:y_mRS90],x0[:RACE_full],
                    x0[:X_mRS_LVO], x0[:beta_mRS_LVO], 
                    x0[:X_mRS_Non_LVO], x0[:beta_mRS_Non_LVO], 
                    x0[:X_mRS_Hem], x0[:beta_mRS_Hem], 
                    x0[:X_mRS_Mim], x0[:beta_mRS_Mim], 
                    x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                    x0[:X_imp_RACE],x0[:beta_imp_RACE], 
                    x0[:sigma_RACE], x0[:sigma_BNIH])
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:BNIH_missing] = proposed_BNIH
    else 
        # reset - TODO find a way to optimize this out
        fc_imp_BNIH!(x0[:BNIH_missing], x0[:BNIH_full], x0[:impute_idx_BNIH], 
                    x0[:condition],x0[:y_mRS90],x0[:RACE_full],
                    x0[:X_mRS_LVO], x0[:beta_mRS_LVO], 
                    x0[:X_mRS_Non_LVO], x0[:beta_mRS_Non_LVO], 
                    x0[:X_mRS_Hem], x0[:beta_mRS_Hem], 
                    x0[:X_mRS_Mim], x0[:beta_mRS_Mim], 
                    x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
                    x0[:X_imp_RACE],x0[:beta_imp_RACE], 
                    x0[:sigma_RACE], x0[:sigma_BNIH])
    end
end

## Imputation function for condtion
function draw_condition!(x0, proportion=0.05, batches=1)
    current_condition = copy(x0[:condition])    

    # Compute current linear predictors        
    eta_mRS_LVO = x0[:X_mRS_LVO] * x0[:beta_mRS_LVO]
    eta_mRS_Non_LVO = x0[:X_mRS_Non_LVO] * x0[:beta_mRS_Non_LVO]
    eta_mRS_Hem = x0[:X_mRS_Hem] * x0[:beta_mRS_Hem]
    eta_mRS_Mim = x0[:X_mRS_Mim] * x0[:beta_mRS_Mim]
    proposed_condition = copy(current_condition)
    accepts = 0
    for batch in 1:batches
        # Choose what indices to propose changes 
        idx_prop = sort(sample(findall(x0[:ischemic]), 
                        Int(floor(proportion*sum(x0[:ischemic])))))
        
        # Flip everything in our proposal
        for i in idx_prop
            proposed_condition[i] = proposed_condition[i] == 1 ? 2 : 1
        end

        f0 = fc_condition(x0[:condition], x0[:y_mRS90], idx_prop, 
                    x0[:X_condition], x0[:beta_condition], 
                        eta_mRS_LVO, eta_mRS_Non_LVO, 
                        eta_mRS_Hem, eta_mRS_Mim)


        f1 = fc_condition(proposed_condition, x0[:y_mRS90], idx_prop, 
                    x0[:X_condition], x0[:beta_condition], 
                        eta_mRS_LVO, eta_mRS_Non_LVO, 
                        eta_mRS_Hem, eta_mRS_Mim)

        a = (f1-f0)              
        if (log(rand(Uniform(0,1),1)[1]) < a )
            accepts += 1
            for i in idx_prop
                x0[:condition][i] = proposed_condition[i]
            end
        else 
            for i in idx_prop
                proposed_condition[i] = current_condition[i]
            end
        end
    end
    #print(string(accepts)*" out of "*string(batches)*" accepts\n")
end





 
function draw_beta_condition!(x0, σ = 0.1)
    proposed = perturb_normal!(copy(x0[:beta_condition]), σ)
    
    f0 = fc_beta_condition(x0[:beta_condition], 
                            x0[:X_condition], 
                            x0[:condition])
    
    f1 = fc_beta_condition(proposed, 
                            x0[:X_condition], 
                            x0[:condition])
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[:beta_condition] = proposed
    end
end

function draw_beta_mRS!(x0, cond, symbol_beta, symbol_x, σ = 0.1)
    proposed = perturb_normal!(copy(x0[symbol_beta]), σ)
    idx = findall(x0[:condition] .== cond)
    f0 = fc_beta(x0[symbol_beta], x0[symbol_x], x0[:y_mRS90], idx)
    f1 =  fc_beta(proposed, x0[symbol_x], x0[:y_mRS90], idx)
    # Symmetric 
    a = (f1-f0)              
    if (log(rand(Uniform(0,1),1)[1]) < a )
        x0[symbol_beta] = proposed
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




monitors = [:beta_mRS_LVO, :beta_mRS_Non_LVO, :beta_mRS_Hem, 
            :beta_mRS_Mim, :beta_condition, :beta_imp_sysBP, 
            :beta_imp_diasBP, :beta_imp_RACE, :beta_imp_BNIH, 
            :sigma_SBP, :sigma_DBP, :sigma_RACE, :sigma_BNIH, 
            :BNIH_missing, :RACE_missing, :sysBP_missing, 
            :diasBP_missing]

function run_mcmc_chain(fname, x0, monitors, warmup, run, thin=10, verbose = 1)
    if (!isfinite(fc_beta_condition(x0[:beta_condition], 
        x0[:X_condition], 
        x0[:condition])))
    error("Beta condition full conditional not fininte at init")
    end

    if (!isfinite(fc_beta(x0[:beta_mRS_LVO], x0[:X_mRS_LVO], x0[:y_mRS90], findall(x0[:condition] .== 1))))
    error("Beta mRS LVO full conditional not finite at init")
    end

    if (!isfinite(fc_beta(x0[:beta_mRS_Non_LVO], x0[:X_mRS_Non_LVO], x0[:y_mRS90], findall(x0[:condition] .== 2))))
    error("Beta mRS Non LVO full conditional not finite at init")
    end

    if (!isfinite(fc_beta(x0[:beta_mRS_Hem], x0[:X_mRS_Hem], x0[:y_mRS90], findall(x0[:condition] .== 3))))
    error("Beta mRS Hem full conditional not finite at init")
    end
    
    if (!isfinite(fc_beta(x0[:beta_mRS_Mim], x0[:X_mRS_Mim], x0[:y_mRS90], findall(x0[:condition] .== 4))))
        error("Beta mRS Mim full conditional not finite at init")
    end


    eta_mRS_LVO = x0[:X_mRS_LVO] * x0[:beta_mRS_LVO]
    eta_mRS_Non_LVO = x0[:X_mRS_Non_LVO] * x0[:beta_mRS_Non_LVO]
    eta_mRS_Hem = x0[:X_mRS_Hem] * x0[:beta_mRS_Hem]
    eta_mRS_Mim = x0[:X_mRS_Mim] * x0[:beta_mRS_Mim]

    if (!isfinite(fc_condition(x0[:condition], x0[:y_mRS90], findall(x0[:ischemic]), 
    x0[:X_condition], x0[:beta_condition], 
    eta_mRS_LVO, eta_mRS_Non_LVO, 
    eta_mRS_Hem, eta_mRS_Mim)))
    error("condition full conditional not fininte at init")
    end


    if (!isfinite(fc_imp_RACE!(x0[:RACE_missing], 
    x0[:RACE_full], x0[:impute_idx_race], 
    x0[:condition],x0[:BNIH_full], x0[:ischemic], x0[:X_condition], 
    x0[:beta_condition], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
    x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sigma_RACE], x0[:sigma_BNIH])))
    error("RACE full conditional not finite at init")
    end

    if (!isfinite(fc_imp_SBP_DBP!(x0[:sysBP_missing], x0[:diasBP_missing], 
    x0[:sysBP_full], 
    x0[:diasBP_full], x0[:impute_idx_sbp],
    x0[:impute_idx_dbp], x0[:impute_idx_bp], 
    x0[:X_condition],x0[:beta_condition], 
    x0[:ischemic], x0[:condition], x0[:X_imp_sysBP], 
    x0[:beta_imp_sysBP], x0[:X_imp_diasBP], 
    x0[:beta_imp_diasBP], x0[:sigma_SBP], 
    x0[:sigma_DBP])))
    error("SBP/DBP full conditional not finite at init")
    end

    if (!isfinite(fc_imp_BNIH!(x0[:BNIH_missing], x0[:BNIH_full], x0[:impute_idx_BNIH], 
    x0[:condition],x0[:y_mRS90],x0[:RACE_full],
    x0[:X_mRS_LVO], x0[:beta_mRS_LVO], 
    x0[:X_mRS_Non_LVO], x0[:beta_mRS_Non_LVO], 
    x0[:X_mRS_Hem], x0[:beta_mRS_Hem], 
    x0[:X_mRS_Mim], x0[:beta_mRS_Mim], 
    x0[:X_imp_BNIH], x0[:beta_imp_BNIH], 
    x0[:X_imp_RACE],x0[:beta_imp_RACE], 
    x0[:sigma_RACE], x0[:sigma_BNIH])))
    error("BNIH full conditional not finite at init")
    end

    if (!isfinite(fc_beta_imp(x0[:beta_imp_sysBP], x0[:X_imp_sysBP], x0[:sigma_SBP], x0[:sysBP_full], [10 for i in eachindex(x0[:beta_imp_sysBP])])+
    fc_beta_imp(x0[:beta_imp_diasBP], x0[:X_imp_diasBP], x0[:sigma_DBP], x0[:diasBP_full], [10 for i in eachindex(x0[:beta_imp_diasBP])])+
    fc_beta_imp(x0[:beta_imp_RACE], x0[:X_imp_RACE], x0[:sigma_RACE], x0[:RACE_full], [10 for i in eachindex(x0[:beta_imp_RACE])])+
    fc_beta_imp(x0[:beta_imp_BNIH], x0[:X_imp_BNIH], x0[:sigma_BNIH], x0[:BNIH_full], [10 for i in eachindex(x0[:beta_imp_BNIH])])+

    fc_sigma_imp(x0[:sigma_SBP], x0[:X_imp_sysBP], x0[:beta_imp_sysBP], x0[:sysBP_full], 1, 1)+
    fc_sigma_imp(x0[:sigma_DBP], x0[:X_imp_diasBP], x0[:beta_imp_diasBP], x0[:sysBP_full], 1, 1)+
    fc_sigma_imp(x0[:sigma_RACE], x0[:X_imp_RACE], x0[:beta_imp_RACE], x0[:sysBP_full], 1, 1)+
    fc_sigma_imp(x0[:sigma_BNIH], x0[:X_imp_BNIH], x0[:beta_imp_BNIH], x0[:sysBP_full], 1, 1)))
    error("Imputation beta/sigma FCs not finite at init")
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
    
    for i in 1:warmup
            
        if verbose > 0
            if i % 100 == 0
                print(string(i)*"/"*string(warmup)*"\n")
            end            
        end
        draw_beta_imp_sysBP!(x0, 0.05)
        draw_beta_imp_diasBP!(x0, 0.05)
        draw_beta_imp_RACE!(x0, 0.05)
        draw_beta_imp_BNIH!(x0, 0.05)
        draw_sigma_imp_terms!(x0, 0.1)
        draw_SBP_DBP!(x0, 0.1, 0.25)
        draw_RACE!(x0, 0.1, 0.25)
        draw_BNIH!(x0, 0.1, 0.25)
        draw_condition!(x0, 0.01, 50)
        draw_beta_condition!(x0, 0.01)
        draw_beta_mRS!(x0, 1, :beta_mRS_LVO, :X_mRS_LVO, 0.05)
        draw_beta_mRS!(x0, 2, :beta_mRS_Non_LVO, :X_mRS_Non_LVO, 0.05)
        draw_beta_mRS!(x0, 3, :beta_mRS_Hem, :X_mRS_Hem, 0.05)
        draw_beta_mRS!(x0, 4, :beta_mRS_Mim, :X_mRS_Mim, 0.05)
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
        draw_beta_imp_sysBP!(x0, 0.05)
        draw_beta_imp_diasBP!(x0, 0.05)
        draw_beta_imp_RACE!(x0, 0.05)
        draw_beta_imp_BNIH!(x0, 0.05)
        draw_sigma_imp_terms!(x0, 0.1)
        draw_SBP_DBP!(x0, 0.1, 0.25)
        draw_RACE!(x0, 0.1, 0.25)
        draw_BNIH!(x0, 0.1, 0.25)
        draw_condition!(x0, 0.01, 50)
        draw_beta_condition!(x0, 0.01)
        draw_beta_mRS!(x0, 1, :beta_mRS_LVO, :X_mRS_LVO, 0.05)
        draw_beta_mRS!(x0, 2, :beta_mRS_Non_LVO, :X_mRS_Non_LVO, 0.05)
        draw_beta_mRS!(x0, 3, :beta_mRS_Hem, :X_mRS_Hem, 0.05)
        draw_beta_mRS!(x0, 4, :beta_mRS_Mim, :X_mRS_Mim, 0.05)
    end
end




############################
### Setup and Simulation ###
############################

## Read in the actual data
fulldata = DataFrame(XLSX.readtable("fulldata.xlsx", "Sheet 1"))

## Rescale BNIH, since it's out of wack for the sampler
fulldata[:,:BNIH] ./= 40


## Define the number of conditions
N_condition = 4

## Define the actual observed outcome - note that we won't be using this in the actual model fit (yet)
y_mRS90_actual = parse.(Int64, fulldata[!,:mRS90])

## Create example X matrix for condition, including missing values
X_condition_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 7)
X_condition_obs[:,1] .= 1
X_condition_obs[:,2] .= fulldata[:,:sysBP]
X_condition_obs[:,3] .= fulldata[:,:diasBP]
X_condition_obs[:,4] .= fulldata[:,:age]
X_condition_obs[:,5] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
X_condition_obs[:,6] .= fulldata[:,:RACE_scale]
X_condition_obs[:,7] .= fulldata[:,:sysBP] .*fulldata[:,:diasBP]

## Create model matrices for mRS outcomes, separately for each condition
## To do this, we create a "full" matrix with all possible terms and then perform column-wise subsets
X_mRS_Full_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 13)
X_mRS_Full_obs[:,1] .= 1
X_mRS_Full_obs[:,2] .= fulldata[:,:age]
X_mRS_Full_obs[:,3] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
X_mRS_Full_obs[:,4] .= fulldata[:,:BNIH]
X_mRS_Full_obs[:,5] .= fulldata[:,:TPA]
X_mRS_Full_obs[:,6] .= fulldata[:,:TPATime]
X_mRS_Full_obs[:,7] .= fulldata[:,:MT]
X_mRS_Full_obs[:,8] .= fulldata[:,:TPA] .* fulldata[:,:MT]
X_mRS_Full_obs[:,9] .= fulldata[:,:TPA] .* fulldata[:,:MT] .* fulldata[:,:MTTime]
X_mRS_Full_obs[:,10] .= fulldata[:,:TPA] .* fulldata[:,:MT] .* fulldata[:,:TPATime]
X_mRS_Full_obs[:,11] .= fulldata[:,:TPA] .* (.!fulldata[:,:MT]) .* fulldata[:,:TPATime]
X_mRS_Full_obs[:,12] .= (.!fulldata[:,:TPA]) .* (fulldata[:,:MT]) .* fulldata[:,:MTTime]
X_mRS_Full_obs[:,13] .= fulldata[:,:TimeToHosp]

X_mRS_LVO_obs = X_mRS_Full_obs[:, [1,2,3,4,5,7,8,9,10,11,12]]
X_mRS_Non_LVO_obs = X_mRS_Full_obs[:, [1,2,3,4,5,6]]
X_mRS_Hem_obs = X_mRS_Full_obs[:, [1,2,3,4,13]]
X_mRS_Mim_obs = X_mRS_Full_obs[:, [1,2,3,4]]

## Create model matrices for imputation 

# From R: RACEformula <- formula("RACE_scale ~ age + sex + BNIH + I(BNIH^2)")
X_imp_RACE_obs = X_mRS_Full_obs[:,[1,2,3,4,4]]
apply_single_interact!(X_imp_RACE_obs, 1:(size(X_imp_RACE_obs)[1]), [4,4], 5) # BNIH^2

# From R: BNIHformula <- formula("BNIH ~ age + sex + RACE_scale")
X_imp_BNIH_obs = Matrix{Union{Float64, Missing}}(undef, length(y_mRS90_actual), 4)
X_imp_BNIH_obs[:,1] .= 1
X_imp_BNIH_obs[:,2] .= fulldata[:,:age]
X_imp_BNIH_obs[:,3] .= 1.0 .* (fulldata[:,:sex] .== "FEMALE")
X_imp_BNIH_obs[:,4] .= fulldata[:,:RACE_scale]

# From R: sysBPformula <- formula("sysBP ~ age + sex + diasBP")
X_imp_sysBP_obs = X_mRS_Full_obs[:, [1,2,3,3]]
X_imp_sysBP_obs[:,4] .= fulldata[:,:diasBP]

# From R: diasBPformula <- formula("diasBP ~ age + sex + sysBP")
X_imp_diasBP_obs = X_mRS_Full_obs[:, [1,2,3,3]]
X_imp_diasBP_obs[:,4] .= fulldata[:,:sysBP]

## Create copies of design matrices to use for "true" imputation in simuation
X_condition_true = deepcopy(X_condition_obs)
X_mRS_LVO_true = deepcopy(X_mRS_LVO_obs)
X_mRS_Non_LVO_true = deepcopy(X_mRS_Non_LVO_obs)
X_mRS_Hem_true = deepcopy(X_mRS_Hem_obs)
X_mRS_Mim_true = deepcopy(X_mRS_Mim_obs)

## Simulate missing variables, and apply to "true" design matrices

### BNIH
impute_idx_BNIH = findall(ismissing.(X_mRS_LVO_obs[:,:4]))
BNIH_missing_true = rand(Normal(10,1), length(impute_idx_BNIH))
apply_single_impute!(X_mRS_LVO_true, impute_idx_BNIH, 4, BNIH_missing_true)
apply_single_impute!(X_mRS_Non_LVO_true, impute_idx_BNIH, 4, BNIH_missing_true)
apply_single_impute!(X_mRS_Hem_true, impute_idx_BNIH, 4, BNIH_missing_true)
apply_single_impute!(X_mRS_Mim_true, impute_idx_BNIH, 4, BNIH_missing_true)

# Imputation code for X_condition
impute_idx_sbp = findall(ismissing.(X_condition_obs[:,:2]))
impute_idx_dbp = findall(ismissing.(X_condition_obs[:,:3]))
impute_idx_bp = sort(union(impute_idx_sbp, impute_idx_dbp))
impute_idx_race = findall(ismissing.(X_condition_obs[:,:6]))

sysBP_missing_true = rand(Normal(0.5,1), length(impute_idx_sbp))
diasBP_missing_true = rand(Normal(0.5,1), length(impute_idx_dbp))
RACE_missing_true = rand(Normal(0.5,1), length(impute_idx_race))

apply_single_impute!(X_condition_true, impute_idx_sbp, 2, sysBP_missing_true)
apply_single_impute!(X_condition_true, impute_idx_dbp, 3, diasBP_missing_true)
apply_single_impute!(X_condition_true, impute_idx_race, 6, RACE_missing_true)
apply_single_interact!(X_condition_true, impute_idx_bp, [2,3], 7)

## Create example betas for condition, which we will treat as "true" values
beta_condition_true  = rand(Normal(0,0.5), (size(X_condition_true)[2], N_condition-1))

## Compute eta_condition
eta_condition = X_condition_true * beta_condition_true

## Determine some (simulated) ischemic indicators
isch_sim = [rand(Binomial(1, 0.25)) for i in 1:(size(fulldata)[1])] .> 0

## Simulate condition
condition_sim = [rPartialCatLogit(eta_condition[i,:], (isch_sim[i] ? (1:2) : (1:4))) 
                    for i in 1:(size(fulldata)[1])]


beta_mRS_LVO_true = generate_init_betas(X_mRS_LVO_obs, 7, 0.5)
beta_mRS_Non_LVO_true = generate_init_betas(X_mRS_Non_LVO_obs, 7, 0.5)
beta_mRS_Hem_true = generate_init_betas(X_mRS_Hem_obs, 7, 0.5)
beta_mRS_Mim_true = generate_init_betas(X_mRS_Mim_obs, 7, 0.5)

eta_mRS_LVO = X_mRS_LVO_true * beta_mRS_LVO_true
eta_mRS_Non_LVO = X_mRS_Non_LVO_true * beta_mRS_Non_LVO_true
eta_mRS_Hem = X_mRS_Hem_true * beta_mRS_Hem_true
eta_mRS_Mim = X_mRS_Mim_true * beta_mRS_Mim_true

# combine etas based on condition
eta_mRS = ((condition_sim .== 1) .* eta_mRS_LVO + 
(condition_sim .== 2) .* eta_mRS_Non_LVO + 
(condition_sim .== 3) .* eta_mRS_Hem + 
(condition_sim .== 4) .* eta_mRS_Mim)

## Simulate mRS outcome for analysis
y_mRS90 = [rCatLogit(eta_mRS[i,:]) for i in 1:(size(eta_mRS)[1])]


#########################
## MCMC Initialization ##
#########################



## Initialize state objects 
state0 = generate_init_state(copy(fulldata))
# update inits data with Simulation
state0[:y_mRS90] = copy(y_mRS90)
state0[:ischemic] = copy(isch_sim)
state0[:condition] = copy(condition_sim)

state1 = generate_init_state(copy(fulldata))
# update inits data with Simulation
state1[:y_mRS90] = copy(y_mRS90)
state1[:ischemic] = copy(isch_sim)
state1[:condition] = copy(condition_sim)

state2 = generate_init_state(copy(fulldata))
# update inits data with Simulation
state2[:y_mRS90] = copy(y_mRS90)
state2[:ischemic] = copy(isch_sim)
state2[:condition] = copy(condition_sim)


truevals = Dict(:beta_condition_true => beta_condition_true,
:beta_mRS_LVO_true => beta_mRS_LVO_true,
:beta_mRS_Non_LVO_true => beta_mRS_Non_LVO_true,
:beta_mRS_Hem_true => beta_mRS_Hem_true,
:beta_mRS_Mim_true => beta_mRS_Mim_true)
write_csv_line(truevals, symbols=[x for x in keys(truevals)], file="truevals.csv", newfile=true)
truefile = open("truevals.csv", "a")
write_csv_line(truevals, symbols=[x for x in keys(truevals)], file=truefile, newfile=false)
close(truefile)


# 2.73 hours
#using Threads
@time begin
Threads.@threads for chain in 1:3
    if chain == 1
        run_mcmc_chain("chain1.csv", state0, monitors, 5000, 100000, 20)
    elseif chain == 2
        run_mcmc_chain("chain2.csv", state1, monitors, 5000, 100000, 20)
    else
        run_mcmc_chain("chain3.csv", state2, monitors, 5000, 100000, 20)
    end
end
end
    