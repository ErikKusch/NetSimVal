#' Update Dynamic Death Rates
#'
#' This is a helper function for updating of dynamic death rates in the simulation run.
#'
#' @param ID_df Data frame of initialising individuals with columns ID, Trait, X, Y, and Species. Output of Sim.Initialise()$ID_df.
#' @param which Either Character or Numeric. If Character and == "Initial" then dynbamic death rate is computed for all individuals. If numeric, dynamic death rate is only computed for indviduals affected by the event applied to the individual in the indexed row.
#' @param event Character. Either "Birth" or "Death".
#' @param d0 Numeric. background death rate (gets altered by the environment and interactions).
#' @param b0 Numeric. background birth rate (remains constant).
#' @param env.xy Environmental matrix as produced by Sim.Space().
#' @param sd Numeric. Habitat suitability in death rate function. Higher values allow individuals to persist in areas of greater environmental maladaptation.
#' @param Effect_Mat Weighted adjacency matrix of association/interaction network.
#' @param Effect_Dis Distance to which association/interaction effects are computed around individuals.
#' @param k_vec Named vector containing carrying capacity for each species. Output of Sim.CarryingK().
#' @param seed Numeric. Seed for random processes.
#'
#' @return An updated data frame of individuals and their dynamic death rates.
#'
Sim.d0Update <- function(
    ID_df, which = "Initial", event = NULL,
    env.xy,
    d0 = 0.4, b0 = 0.6, sd = 2.5,
    Effect_Mat,
    k_vec,
    Effect_Dis = 0.5,
    seed) {
  # set.seed(seed)

  ## dynamic death rate functions
  d0P <- function(N, b0, d0, k) {
    N * (b0 - d0) / k
  }
  d0E <- function(Tr, Env, sd) {
    exp((abs(Tr - Env) / sd)^2)
  }
  d0Omega <- function(Effect_Mat, Effect_Dis, ID_df, i) {
    Abundances_i <- rep(0, ncol(Effect_Mat))
    names(Abundances_i) <- colnames(Effect_Mat)
    Abundances_obs <- table(ID_df[-i, ][
      abs(ID_df[i, "X"] - ID_df[-i, "X"]) <= Effect_Dis &
        abs(ID_df[i, "Y"] - ID_df[-i, "Y"]) <= Effect_Dis,
      "Species"
    ])
    Abundances_i[match(names(Abundances_obs), names(Abundances_i))] <- Abundances_obs
    ### Extract all effects that the focal species is subject to in the interaction matrix
    Effects_i <- Effect_Mat[which(rownames(Effect_Mat) == ID_df[i, "Species"]), ]
    ### Weigh effects by abundances
    WeightedEffects_i <- Abundances_i * Effects_i
    ### Calculate final effect
    if (sum(Abundances_i) != 0) {
      FinalEffect_i <- sum(WeightedEffects_i) / sum(Abundances_i)
    } else {
      FinalEffect_i <- 0
    }
    return(FinalEffect_i)
  }
  dt <- function(d0, d0P, d0E, d0Omega) {
    d0 + d0P * d0E - d0Omega
  }
  Get.Environment <- function(x, y, env_mat) {
    col <- which.min(abs(x - as.numeric(colnames(env_mat))))[1]
    row <- which.min(abs(y - as.numeric(rownames(env_mat))))[1]
    env_mat[row, col]
  }
  ## queried d0 update
  if (which[1] == "Initial") {
    ## dynamic death rate components
    ### population dynamics
    N_vec <- table(ID_df$Species)
    d0P_vec <- d0P(
      N = N_vec, b0 = b0, d0 = d0,
      k = k_vec[match(names(k_vec), names(N_vec))]
    )
    ID_df$d0P <- as.numeric(d0P_vec[match(ID_df$Species, names(d0P_vec))])
    ### environment
    ID_df$Env <- apply(ID_df, 1, FUN = function(row){
      # print(row)
      Get.Environment(
        x = as.numeric(row["X"]),
        y = as.numeric(row["Y"]),
        env_mat = env.xy
      )
    })
    ID_df$d0E <- d0E(
      Tr = ID_df$Trait, Env = ID_df$Env,
      sd = sd
    )
    ### Interactions
    ID_df$d0Omega <- sapply(1:nrow(ID_df), FUN = function(i) {
      d0Omega(
        Effect_Mat = Effect_Mat, Effect_Dis = Effect_Dis,
        ID_df = ID_df, i = i
      )
    })
  } else {
    ## updating population dynamics
    if (sum(ID_df$Species == which$Species) > 0) {
      d0P <- d0P(
        N = sum(ID_df$Species == which$Species),
        b0 = b0, d0 = d0,
        k = k_vec[names(k_vec) == which$Species]
      )
      ID_df[ID_df$Species == which$Species, "d0P"] <- d0P
    }
    ## environment
    if (event == "Birth") {
      which$Env <- Get.Environment(
          x = which$X,
          y = which$Y,
          env_mat = env.xy)
      ID_df$d0E[ID_df$ID == which$ID] <- d0E(
        Tr = which$Trait,
        Env = which$Env,
        sd = sd
      )
    }
    ## interactions
    ### individuals affected by addition or removal of other indiivudal
    newinterac <- ID_df[abs(ID_df$X - which$X) <= Effect_Dis &
      abs(ID_df$Y - which$Y) <= Effect_Dis, "ID"]
    if (length(newinterac) > 0) {
      ID_df$d0Omega[ID_df$ID %in% newinterac] <- sapply(newinterac, FUN = function(i) {
        d0Omega(
          Effect_Mat = Effect_Mat, Effect_Dis = Effect_Dis,
          ID_df = ID_df, i = i
        )
      })
    }
  }
  ID_df$dt <- dt(
    d0 = d0,
    d0P = ID_df$d0P,
    d0E = ID_df$d0E,
    d0Omega = ID_df$d0Omega
  )
  ID_df$dt[ID_df$dt < 0] <- 0 # make sure probabilities are never < 0
  return(ID_df)
}
