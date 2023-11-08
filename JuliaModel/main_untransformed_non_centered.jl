## don't need this if open this folder as workspace in vscode
cd("H:/Thesis/VARgrowthSimulations/JuliaModel")
## include should have necessary packages
include("untransformed_parameter_non_centered.jl")
# include("untransformed_parameter.jl")
using CSV
using Random
using Plots
using Optim
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

# beta_asym = log(20000.0)
# beta_offset = log(2.0)
# beta_growth = log(0.9/(1-0.9))

beta_asym = 20000.0
beta_offset = 2.0
beta_growth = 0.9
X = repeat([1], outer = U) 

# latent_sd_asym = 0.1
# latent_sd_offset = 0.001
# latent_sd_growth = 0.001

latent_sd_asym = 1
latent_sd_offset = 0.001
latent_sd_growth = 0.001

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

mu = [asym[i][1] for i in df.group] .* 
        exp.(-[offset[i][1] for i in df.group] .* exp.(df.time .* [log(growth[i][1]) for i in df.group]) ) 

# mu = [exp(asym[i][1]) for i in df.group] .* 
#         exp.(-[exp(offset[i][1]) for i in df.group] .* 
#         exp.(df.time .* [log( logistic(growth[i][1])) for i in df.group]) ) 
df.mu .= mu 

outcome = [rand(Normal(temp_mu, obs_sd), 1)[1] for temp_mu in df.mu]

df.outcome .= outcome

CSV.write(data_dir*"/untransform_data_sim.csv", df)

function perturb_init(init, perturb_sd::Dict)
        Dict(:asym => [rand(Normal(a, perturb_sd[:asym]), 1)[1] for a in init[:asym] ],
        :offset => [rand(Normal(a, perturb_sd[:offset]), 1)[1] for a in init[:offset] ],
        :growth => [rand(Normal(a, perturb_sd[:growth]), 1)[1] for a in init[:growth]],
        :beta_asym => [rand(Normal(a, perturb_sd[:beta_asym]), 1)[1] for a in init[:beta_asym] ],
        :beta_offset => [rand(Normal(a, perturb_sd[:beta_offset]), 1)[1] for a in init[:beta_offset] ],
        :beta_growth =>  [rand(Normal(a, perturb_sd[:beta_growth]), 1)[1] for a in init[:beta_growth] ],
        :latent_sd_asym => [rand(truncated(Normal(a, perturb_sd[:latent_sd_asym]), lower = 0), 1)[1] for a in init[:latent_sd_asym]],
        :latent_sd_offset => [rand(truncated(Normal(a, perturb_sd[:latent_sd_offset]), lower = 0), 1)[1] for a in init[:latent_sd_offset]],
        :latent_sd_growth =>  [rand(truncated(Normal(a, perturb_sd[:latent_sd_growth]), lower = 0), 1)[1] for a in init[:latent_sd_growth]],
        :obs_sd =>  [rand(truncated(Normal(a, perturb_sd[:obs_sd]), lower = 0), 1)[1]  for a in init[:obs_sd]])

end

init_sd_trans = Dict(:asym => 0.1,
:offset => 0.01,
:growth => 0.01,
:beta_asym => 0.05,
:beta_offset => 0.01,
:beta_growth =>  0.01,
:latent_sd_asym => 0.01,
:latent_sd_offset => 0.01,
:latent_sd_growth =>  0.01,
:obs_sd =>  10)

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

inits2 = perturb_init(init_truth, init_sd_trans)
inits3 = perturb_init(init_truth, init_sd_trans)
inits4 = perturb_init(init_truth, init_sd_trans)


#########################
## MCMC Initialization ##
#########################

## Initialize state objects 
state0 = generate_init_state(copy(df), init_truth)

# function run_mcmc_chain(fname, x0, monitors, dat, warmup, run, thin=10, verbose = 1)
# latent_sd_offset not looking so good
perturb_sd = Dict(:asym => 0.25, #asym
:offset => 0.15, #offset
:growth => 0.07, #growth
:beta_asym => 6.5, #beta_asym
:beta_offset => 0.0035,  #beta_offset
:beta_growth => 0.0001, #beta_growth
:latent_sd_asym => 2.3, #latent_sd_asym
:latent_sd_offset => 0.003, #latent_sd_offset
:latent_sd_growth => 0.0001, #latent_sd_growth
:obs_sd => 5 #obs_sd
)

# perturb_sd_trans = Dict(:asym => 0.002, #asym
# :offset => 0.03, #offset
# :growth => 0.032, #growth
# :beta_asym => 0.00035, #beta_asym
# :beta_offset => 0.0015,  #beta_offset
# :beta_growth => 0.0013, #beta_growth
# :latent_sd_asym => 0.0004, #latent_sd_asym
# :latent_sd_offset => 0.0015, #latent_sd_offset
# :latent_sd_growth => 0.0017, #latent_sd_growth
# :obs_sd => 5 #obs_sd
# )

# function run_mcmc_chain(fname, x0, monitors, dat, priors, perturb_sd, warmup, run, thin=10, verbose = 1)

@time run_mcmc_chain(data_dir * "/reparam.csv", state0, monitors, df, priors, 
                        gompertz_likelihood, perturb_sd, 1,  10000, 20)


##################################
######## Parallel MCMC ###########
##################################

state0 = generate_init_state(copy(df), init_truth)
state1 = generate_init_state(copy(df), inits2)
state2 = generate_init_state(copy(df), inits3)
state3 = generate_init_state(copy(df), inits4)

# true_file_name = data_dir * "/truevals_reparam.csv"
# write_csv_line(init_truth, symbols=[x for x in keys(init_truth)], file= true_file_name, newfile=true)
# truefile = open(true_file_name, "a")
# write_csv_line(init_truth, symbols=[x for x in keys(init_truth)], file=truefile, newfile=false)
# close(truefile)

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

# warmup = 1000
# run = 100000

# @time begin
# @threads for chain in 1:4
#     if chain == 1
#         run_mcmc_chain(data_dir * "/reparam1.csv", state1, monitors, df, perturb_sd, warmup, run, 20)
#         write_accepts(state1, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     elseif chain == 2
#         run_mcmc_chain(data_dir * "/reparam2.csv", state2, monitors, df, perturb_sd, warmup, run, 20)
#         write_accepts(state2, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     elseif chain == 3
#         run_mcmc_chain(data_dir * "/reparam3.csv", state3, monitors, df, perturb_sd, warmup, run, 20)
#         write_accepts(state3, data_dir * "/accept_reparam" * string(chain) * ".csv", accept_keys)
#     else 
#         run_mcmc_chain(data_dir * "/reparam0.csv", state0, monitors, df, perturb_sd, warmup, run, 20)
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