## don't need this if open this folder as workspace in vscode
cd("H:/Thesis/VARgrowthSimulations/JuliaModel")
## include should have necessary packages
include("untransformed_parameter.jl")
using CSV
using Random
############################
### Setup and Simulation ###
############################
## change directory based on sim
data_dir = "../data/JuliaSim/Sim2"
## Read in the actual data
fulldata = CSV.read(data_dir * "/sim.csv", DataFrame)

########################################
### Sim Different gp parameterization ##
########################################
Random.seed!(1234)

t = 50
U = 50

beta_asym = 20000.0
beta_offset = 2.0
beta_growth = 0.83
X = repeat([1], outer = U) 

latent_sd_asym = 1000
latent_sd_offset = 0.01
latent_sd_growth = 0.01

obs_sd = 10

asym = [rand(d, 1)[1] for d in Normal.(X*beta_asym, latent_sd_asym)] 
offset = [rand(d, 1)[1] for d in Normal.(X*beta_offset, latent_sd_offset)] 
growth = [rand(d, 1)[1] for d in Normal.(X*beta_growth, latent_sd_growth)] 

df = DataFrame(group = repeat(collect(1:U), inner = t),
               time = repeat(collect(1:t), outer = U)
)

mu = Vector{Float64}(undef, nrow(df))

# for i in 1:nrow(df)
#     mu[i] = asym[df.group[i]] * exp(-offset[df.group[i]] * exp(df.time[i]*log(growth[df.group[i]]) ) )
# end

# why doesn't log work here
mu = [asym[i][1] for i in df.group] .* 
        exp.(-[offset[i][1] for i in df.group] .* exp.(df.time .* [log(growth[i][1]) for i in df.group]) ) 
df.mu .= mu 
[rand(Normal(temp_mu, obs_sd), 1) for temp_mu in df.mu]
outcome = [rand(Normal(temp_mu, obs_sd), 1)[1] for temp_mu in df.mu]

df.outcome .= outcome



## the [1] is to ensure we don't get a vector of vectors
inits2 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in [latent_sd_offset]],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_growth]],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in [obs_sd]])

inits3 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in [latent_sd_offset]],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_growth]],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in [obs_sd]])

inits4 = Dict(:asym => [rand(Normal(a, 100), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.1), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 100), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.1), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 100), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.1), lower = 0), 1)[1] for a in [latent_sd_offset]],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_growth]],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in [obs_sd]])

init_truth = Dict(:asym => asym,
:offset => offset,
:growth => growth,
:beta_asym => [beta_asym ;;],
:beta_offset => [beta_offset ;;],
:beta_growth => [beta_growth ;;],
:latent_sd_asym => [latent_sd_asym],
:latent_sd_offset => [latent_sd_offset],
:latent_sd_growth => [latent_sd_growth],
:obs_sd => [obs_sd])



#########################
## MCMC Initialization ##
#########################

## Initialize state objects 
# state0 = generate_init_state(copy(fulldata), inits)
# @time run_mcmc_chain(data_dir * "/chains1.csv", state0, monitors, fulldata, 5000, 100000, 20)

state0 = generate_init_state(copy(df), init_truth)
state1 = generate_init_state(copy(df), inits2)
state2 = generate_init_state(copy(df), inits3)
state3 = generate_init_state(copy(df), inits4)
# function run_mcmc_chain(fname, x0, monitors, dat, warmup, run, thin=10, verbose = 1)
# @time run_mcmc_chain(data_dir * "/reparam.csv", state1, monitors, fulldata, 10000, 100000, 20)
# truevals = Dict(:beta_condition_true => beta_condition_true,
# :beta_mRS_LVO_true => beta_mRS_LVO_true,
# :beta_mRS_Non_LVO_true => beta_mRS_Non_LVO_true,
# :beta_mRS_Hem_true => beta_mRS_Hem_true,
# :beta_mRS_Mim_true => beta_mRS_Mim_true)
true_file_name = data_dir * "/truevals_reparam.csv"
write_csv_line(init_truth, symbols=[x for x in keys(init_truth)], file= true_file_name, newfile=true)
truefile = open(true_file_name, "a")
write_csv_line(init_truth, symbols=[x for x in keys(init_truth)], file=truefile, newfile=false)
close(truefile)


# 2.73 hours
using Base.Threads
# define number of threads in settings or set computer environment variable
if nthreads() < 4
    error("Set number of threads before running parallel")
end


@time begin
@threads for chain in 1:4
    if chain == 1
        run_mcmc_chain(data_dir * "/reparam1.csv", state1, monitors, df, 10000, 100000, 20)
    elseif chain == 2
        run_mcmc_chain(data_dir * "/reparam2.csv", state2, monitors, df, 10000, 100000, 20)
    elseif chain == 3
        run_mcmc_chain(data_dir * "/reparam3.csv", state3, monitors, df, 10000, 100000, 20)
    else 
        run_mcmc_chain(data_dir * "/reparam0.csv", state0, monitors, df, 10000, 100000, 20)
    end
end
end