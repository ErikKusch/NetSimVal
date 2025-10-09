Network_igraph <- Sim.Network(n_spec = 10, NetworkType = "Association", Sparcity = 0.5, MaxStrength = 1, seed = 42)
use_data(Network_igraph)

CarryingK_vec <- Sim.CarryingK(n_spec = 10, k_range = c(200,200), seed = 42)
use_data(CarryingK_vec)

Niches_vec <- Sim.Niche(n_spec = 10, Env_range = c(0, 10), seed = 42)
use_data(Niches_vec)

Initialise_df <- Sim.Initialise(n_spec = 10, n_individuals = 4e2,n_mode = "total", Env_range = c(0, 10), 
                                Trait_means = Niches_vec, Trait_sd = 1, seed = 42)
use_data(Initialise_df)

Env_mat <- Sim.Space(
    x_range = c(0, 10), y_range = c(0, 10),
    ncol = 1e3, nrow = 1e3,
    x_gradient = function(x) x,
    y_gradient = function(y) y
)
use_data(Env_mat)

SimulationOutput <- Sim.Compute(d0 = 0.4,
                                b0 = 0.6,
                                env.xy = Env_mat,
                                t_max = 10,
                                t_inter = 0.1,
                                sd = 2.5,
                                migration = 0.2,
                                Effect_Dis = 0.5,
                                Network_igraph,
                                k_vec = CarryingK_vec,
                                ID_df = Initialise_df,
                                seed = 42,
                                verbose = TRUE, # whether to print progress in time as current time
                                RunName = "TestRun"
)
use_data(SimulationOutput)

Inferred_igraph <- Sim.Network(n_spec = 10, NetworkType = "Association", Sparcity = 0.5, MaxStrength = 1, seed = 43)
use_data(Inferred_igraph)
