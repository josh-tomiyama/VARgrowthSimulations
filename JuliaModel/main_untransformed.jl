## don't need this if open this folder as workspace in vscode
cd("H:/Thesis/VARgrowthSimulations/JuliaModel")
## include should have necessary packages
include("untransformed_parameter.jl")
# include("untransformed_parameter.jl")
using CSV
using Random
############################
### Setup and Simulation ###
############################
## change directory based on sim
data_dir = "../data/JuliaSim/Sim2"
## Read in the actual data
# fulldata = CSV.read(data_dir * "/sim.csv", DataFrame)

########################################
### Sim Different gp parameterization ##
########################################
Random.seed!(1234)

t = 50
U = 50

beta_asym = 20000.0
beta_offset = 2.0
beta_growth = 0.9
X = repeat([1], outer = U) 

latent_sd_asym = 1
latent_sd_offset = 0.01
latent_sd_growth = 0.01

obs_sd = 100

# asym = [rand(d, 1)[1] for d in truncated.(Normal.(X*beta_asym, latent_sd_asym), lower = 0, upper = Inf)] 
# offset = [rand(d, 1)[1] for d in truncated.(Normal.(X*beta_offset, latent_sd_offset), lower = 0, upper = Inf)] 
# growth = [rand(d, 1)[1] for d in truncated.(Normal.(X*beta_growth, latent_sd_growth), lower = 0, upper = 1)] 

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
## function is exactly the same
# mu_test = @. (asym[df.group] * exp( -offset[df.group] * exp( log( growth[df.group])  * df.time) ) )
# why doesn't log work here
mu = [asym[i][1] for i in df.group] .* 
        exp.(-[offset[i][1] for i in df.group] .* exp.(df.time .* [log(growth[i][1]) for i in df.group]) ) 


df.mu .= mu 

outcome = [rand(Normal(temp_mu, obs_sd), 1)[1] for temp_mu in df.mu]

df.outcome .= outcome

CSV.write(data_dir*"/untransform_data_sim.csv", df)



## the [1] is to ensure we don't get a vector of vectors
inits2 = Dict(:asym => [rand(Normal(a, 0.5), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.01), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 0.5), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.01), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 0.5), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_offset]],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_growth]],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in [obs_sd]])

inits3 = Dict(:asym => [rand(Normal(a, 0.5), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.01), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 0.5), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.01), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 0.5), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_offset]],
:latent_sd_growth =>  [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_growth]],
:obs_sd =>  [rand(truncated(Normal(a, 10), lower = 0), 1)[1]  for a in [obs_sd]])

inits4 = Dict(:asym => [rand(Normal(a, 0.5), 1)[1] for a in asym],
:offset => [rand(Normal(a, 0.01), 1)[1] for a in offset],
:growth => [rand(Normal(a, 0.01), 1)[1] for a in growth],
:beta_asym => [rand(Normal(a, 0.5), 1)[1] for a in [beta_asym ;;]],
:beta_offset => [rand(Normal(a, 0.01), 1)[1] for a in [beta_offset ;;]],
:beta_growth =>  [rand(Normal(a, 0.01), 1)[1] for a in [beta_growth ;;]],
:latent_sd_asym => [rand(truncated(Normal(a, 0.5), lower = 0), 1)[1] for a in [latent_sd_asym]],
:latent_sd_offset => [rand(truncated(Normal(a, 0.01), lower = 0), 1)[1] for a in [latent_sd_offset]],
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
run_mcmc_chain(data_dir * "/reparam.csv", state0, monitors, df, 1, 20000, 20)
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

# using Base.Threads
# # define number of threads in settings or set computer environment variable
# if nthreads() < 4
#     error("Set number of threads before running parallel")
# end
# function write_accepts(state, accept_file_name, accept_keys)
#     write_csv_line(state, symbols=accept_keys, file= accept_file_name, newfile=true)
#     accept_file = open(accept_file_name, "a")
#     write_csv_line(state, symbols=accept_keys, file=accept_file, newfile=false)
#     close(accept_file)
# end

# accept_keys = [:accept_raw_asym ,
# :accept_raw_offset ,
# :accept_raw_growth ,
# :accept_beta_asym ,
# :accept_beta_offset ,
# :accept_beta_growth ,
# :accept_latent_sd_asym ,
# :accept_latent_sd_offset ,
# :accept_latent_sd_growth ,
# :accept_obs_sd ]

# warmup = 10000
# run = 50000

# @time begin
# @threads for chain in 1:4
#     if chain == 1
#         run_mcmc_chain(data_dir * "/reparam1.csv", state1, monitors, df, warmup, run, 20)
#         write_accepts(state1, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     elseif chain == 2
#         run_mcmc_chain(data_dir * "/reparam2.csv", state2, monitors, df, warmup, run, 20)
#         write_accepts(state2, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     elseif chain == 3
#         run_mcmc_chain(data_dir * "/reparam3.csv", state3, monitors, df, warmup, run, 20)
#         write_accepts(state3, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     else 
#         run_mcmc_chain(data_dir * "/reparam0.csv", state0, monitors, df, warmup, run, 20)
#         write_accepts(state0, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     end
# end
# end



# [state0[s] = [state0[s]] for s in accept_keys]
# [state1[s] = [state1[s]] for s in accept_keys]
# [state2[s] = [state2[s]] for s in accept_keys]
# [state3[s] = [state3[s]] for s in accept_keys]

# function write_accepts(state, accept_file_name, accept_keys)
#     write_csv_line(state, symbols=accept_keys, file= accept_file_name, newfile=true)
#     accept_file = open(accept_file_name, "a")
#     write_csv_line(state, symbols=accept_keys, file=accept_file, newfile=false)
#     close(accept_file)
# end
# accept_file_name = data_dir * "/accept_reparam0.csv"


# accept_file_name = data_dir * "/accept_reparam1.csv"
# write_csv_line(state1, symbols=accept_keys, file= accept_file_name, newfile=true)
# accept_file = open(accept_file_name, "a")
# write_csv_line(state1, symbols=accept_keys, file=accept_file, newfile=false)
# close(accept_file)

# accept_file_name = data_dir * "/accept_reparam2.csv"
# write_csv_line(state2, symbols=accept_keys, file= accept_file_name, newfile=true)
# accept_file = open(accept_file_name, "a")
# write_csv_line(state2, symbols=accept_keys, file=accept_file, newfile=false)
# close(accept_file)

# accept_file_name = data_dir * "/accept_reparam3.csv"
# write_csv_line(state3, symbols=accept_keys, file= accept_file_name, newfile=true)
# accept_file = open(accept_file_name, "a")
# write_csv_line(state3, symbols=accept_keys, file=accept_file, newfile=false)
# close(accept_file)