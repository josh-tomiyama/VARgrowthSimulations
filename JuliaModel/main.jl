## don't need this if open this folder as workspace in vscode
cd("H:/Thesis/VARgrowthSimulations/JuliaModel")
## include should have necessary packages
include("straight_forward.jl")
using CSV
using Random
############################
### Setup and Simulation ###
############################
data_dir = "../data/JuliaSim/Sim2"
## Read in the actual data
fulldata = CSV.read(data_dir * "/sim.csv", DataFrame)

## May replace with the last line of current output file
params = CSV.read(data_dir * "/true_params.csv", DataFrame)
names(params)

## Hopefully matrices work here
asym = params[:, startswith.(names(params), "asym_")]
asym = Matrix(asym)'

offset = params[:, startswith.(names(params), "offset_")]
offset = Matrix(offset)'

growth = params[:, startswith.(names(params), "growth_")]
growth = Matrix(growth)'

beta_asym = params[:, startswith.(names(params), "beta_asym")]
beta_asym = Matrix(beta_asym)'
beta_offset = params[:, startswith.(names(params), "beta_offset")]
beta_offset = Matrix(beta_offset)'
beta_growth = params[:, startswith.(names(params), "beta_growth")]
beta_growth = Matrix(beta_growth)'

latent_sd_asym = params[:, startswith.(names(params), "latent_sd_asym")]
latent_sd_asym = Matrix(latent_sd_asym)'
latent_sd_offset = params[:, startswith.(names(params), "latent_sd_offset")]
latent_sd_offset = Matrix(latent_sd_offset)'
latent_sd_growth = params[:, startswith.(names(params), "latent_sd_growth")]
latent_sd_growth = Matrix(latent_sd_growth)'
obs_sd = params[:, startswith.(names(params), "obs_sd")]
obs_sd = Matrix(obs_sd)'

inits_truth = Dict(:asym => asym,
:offset => offset,
:growth => growth,
:beta_asym => beta_asym,
:beta_offset => beta_offset,
:beta_growth => beta_growth,
:latent_sd_asym => latent_sd_asym,
:latent_sd_offset => latent_sd_offset,
:latent_sd_growth => latent_sd_growth,
:obs_sd => obs_sd)

Random.seed!(1234)

inits2 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in beta_asym],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in beta_offset],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in beta_growth ],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in latent_sd_asym],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in latent_sd_offset],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in latent_sd_growth],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in obs_sd])

inits3 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in beta_asym],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in beta_offset],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in beta_growth ],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in latent_sd_asym],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in latent_sd_offset],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in latent_sd_growth],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in obs_sd])

inits4 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in beta_asym],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in beta_offset],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in beta_growth ],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in latent_sd_asym],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in latent_sd_offset],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in latent_sd_growth],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in obs_sd])

#########################
## MCMC Initialization ##
#########################

## Initialize state objects 
# state0 = generate_init_state(copy(fulldata), inits)
# @time run_mcmc_chain(data_dir * "/chains1.csv", state0, monitors, fulldata, 5000, 100000, 20)

state0 = generate_init_state(copy(fulldata), inits_truth)
state1 = generate_init_state(copy(fulldata), inits2)
state2 = generate_init_state(copy(fulldata), inits3)
state3 = generate_init_state(copy(fulldata), inits4)
# function run_mcmc_chain(fname, x0, monitors, dat, warmup, run, thin=10, verbose = 1)
@time run_mcmc_chain(data_dir * "/reparam.csv", state0, monitors, fulldata, 1, 10000, 20)

# write_csv_line(inits_truth, symbols=[x for x in keys(inits_truth)], file="truevals_model.csv", newfile=true)
# truefile = open("truevals_model.csv", "a")
# write_csv_line(inits_truth, symbols=[x for x in keys(inits_truth)], file=truefile, newfile=false)
# close(truefile)


# 2.73 hours
# using Base.Threads
# # define number of threads in settings or set computer environment variable
if nthread < 4
    error("Set number of threads before running parallel")
end
@time begin
Threads.@threads for chain in 1:4
    if chain == 1
        run_mcmc_chain(data_dir * "/chains1.csv", state1, monitors, fulldata, 10000, 100000, 20)
    elseif chain == 2
        run_mcmc_chain(data_dir * "/chains2.csv", state2, monitors, fulldata, 10000, 100000, 20)
    elseif chain == 3
        run_mcmc_chain(data_dir * "/chains3.csv", state3, monitors, fulldata, 10000, 100000, 20)
    else 
        run_mcmc_chain(data_dir * "/chains0.csv", state0, monitors, fulldata, 10000, 100000, 20)
    end
end
end