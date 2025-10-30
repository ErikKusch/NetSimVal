#' Carry out Demographic Simulation
#'
#' This function carries out spatially explicit, individual-based demographic simulation for the generation of data products ready for network inference.
#'
#' @param d0 Numeric. background death rate (gets altered by the environment and interactions).
#' @param b0 Numeric. background birth rate (remains constant).
#' @param env.xy Environmental matrix as produced by Sim.Space().
#' @param t_max Numeric. Maximum simulation time.
#' @param t_inter Numeric. Interval length at which to record simulation outputs (measured in simulation time).
#' @param sd Numeric. Habitat suitability in death rate function. Higher values allow individuals to persist in areas of greater environmental maladaptation.
#' @param migration Numeric. Standard deviation of 0-centred normal distribution from which natal dispersal is drawn.
#' @param range Numeric. Distance around a birth event within which the chance of a birth to occur are high (i.e. the flat-top portion of a bivariate flat-top normal distribution).
#' @param Effect_Dis Numeric. Distance within which neighbouring individuals interact with a focal individual.
#' @param Network_igraph An igraph object with association/interaction strength stored as "weight" attribute of edges. Output of Sim.Network().
#' @param k_vec Named vector containing carrying capacity for each species. Output of Sim.CarryingK().
#' @param ID_df Data frame of initialising individuals with columns ID, Trait, X, Y, and Species. Output of Sim.Initialise()$ID_df.
#' @param seed Numeric. Seed for random processes.
#' @param verbose Logical. Whether to print simulation time and sampling interval.
#' @param RunName Character. Name for temporary .RData object written to disk.
#'
#' @return A list containing data frame object with the same columns as ID_df at each sampling interval defined via n_inter until t_max is reached.
#'
#' @importFrom lubridate seconds_to_period
#' @importFrom igraph as_adjacency_matrix
#'
#' @examples
#' data("Initialise_df")
#' data("CarryingK_vec")
#' data("Network_igraph")
#' data("Env_mat")
#' 
#' SimResult <- Sim.Compute(
#'     d0 = 0.4,
#'     b0 = 0.6,
#'     env.xy = Env_mat,
#'     t_max = 5,
#'     t_inter = 0.1,
#'     sd = 2.5,
#'     migration = 0.2,
#'     range = 0.05,
#'     Effect_Dis = 0.5,
#'     Network_igraph = Network_igraph,
#'     k_vec = CarryingK_vec,
#'     ID_df = Initialise_df,
#'     seed = 42,
#'     verbose = TRUE, # whether to print progress in time as current time
#'     RunName = "Trial"
#' )
#'
#' @export
Sim.Compute <- function(
    d0 = 0.4,
    b0 = 0.6,
    env.xy, # space matrix
    t_max = 10,
    t_inter = 0.1,
    sd = 2.5,
    migration = 0.2,
    range = 0.05,
    Effect_Dis = 0.5,
    Network_igraph,
    k_vec,
    ID_df,
    seed = 42,
    verbose = TRUE, # whether to print progress in time as current time
    RunName = "") {
  call_info <- match.call()
  set.seed(seed)

  ## Network as adjacency matrix
  Effect_Mat <- as_adjacency_matrix(Network_igraph, attr = "weight") # columns affect rows
  rownames(Effect_Mat) <- colnames(Effect_Mat) <- names(k_vec)

  ## environmental range
  Env_range_x <- range(as.numeric(colnames(Env_mat)))
  Env_range_y <- range(as.numeric(rownames(Env_mat)))

  ## dynamic death rate initialisation
  ID_df <- Sim.d0Update(
    ID_df = ID_df, which = "Initial",
    env.xy = env.xy, d0 = d0, b0 = b0, sd = sd,
    Effect_Mat = Effect_Mat, k_vec = k_vec,
    Effect_Dis = Effect_Dis, seed = seed
  )

  ## list object to store individuals at each time step
  ID_ls <- list(ID_df)

  ## setting start times
  t <- 0 # start time at 1
  names(ID_ls)[length(ID_ls)] <- t
  ## progress bar
  if (!verbose) {
    pb <- txtProgressBar(min = 0, max = t_max, style = 3)
  }
  TimeStart <- Sys.time()
  ## simulation loop over time steps
  while (t < t_max) {
    # if(verbose){print(t)}
    ## vectors for storing birth and death probabilities for each individual
    birth_prob <- rep(b0, nrow(ID_df))
    death_prob <- ID_df$dt
    names(birth_prob) <- names(death_prob) <- ID_df$ID

    ## event identification
    EventSample_vec <- paste(rep(c("Birth", "Death"), each = nrow(ID_df)), names(birth_prob), sep = "_")
    ProbSample_vec <- c(birth_prob, death_prob)
    if (any(ProbSample_vec == Inf)) {
      event <- EventSample_vec[which(ProbSample_vec == Inf)[1]]
    } else {
      event <- sample(
        EventSample_vec,
        size = 1,
        prob = ProbSample_vec
      )
    }
    ## event evaluation
    event_eval <- strsplit(event, split = "_")
    event_ID <- event_eval[[1]][2]
    event_EV <- event_eval[[1]][1]
    if (event_EV == "Birth") {
      append_df <- ID_df[ID_df$ID == event_ID, ]
      append_df$ID <- max(ID_df$ID) + 1
      newloc.xy <- rFlatTopNorm(1, 
                                mean = c(append_df$X, append_df$Y),
                                sd = migration,
                                range = range,
                                truncDist = Effect_Dis)
      newloc.x <- newloc.xy[1]
      newloc.y <- newloc.xy[2]

      ## ensuring species don't disperse beyond the environmental limit
      newloc.x <- ifelse(newloc.x < Env_range_x[1], Env_range_x[1], newloc.x)
      newloc.x <- ifelse(newloc.x > Env_range_x[2], Env_range_x[2], newloc.x)
      newloc.y <- ifelse(newloc.y < Env_range_y[1], Env_range_y[1], newloc.y)
      newloc.y <- ifelse(newloc.y > Env_range_y[2], Env_range_y[2], newloc.y)
      append_df$X <- newloc.x
      append_df$Y <- newloc.y
      ID_df <- rbind(ID_df, append_df)
      affected_row <- ID_df[nrow(ID_df), ]
    }
    if (event_EV == "Death") {
      affected_row <- ID_df[ID_df$ID == as.numeric(event_ID), ]
      ID_df <- ID_df[ID_df$ID != event_ID, ]
    }
    ## dynamic death rate recalculation
    ID_df <- Sim.d0Update(
      ID_df = ID_df, which = affected_row, event = event_EV,
      env.xy = env.xy, d0 = d0, b0 = b0, sd = sd,
      Effect_Mat = Effect_Mat, k_vec = k_vec,
      Effect_Dis = Effect_Dis, seed = round(t, 3) * 1e4 + seed
    )
    ## Gillespie time
    ### identify by how much time advances
    tadvance <- rexp(1, rate = sum(c(birth_prob, death_prob)))
    t <- t + tadvance
    ### record data only if interval is met
    if (t - as.numeric(names(ID_ls)[length(ID_ls)]) >= t_inter) {
      # if(verbose){message(t)}
      ID_ls <- c(ID_ls, list(ID_df))
      names(ID_ls)[length(ID_ls)] <- t
      saveobj <- list(Call = call_info, Network = Network_igraph, K = k_vec, Simulation = ID_ls)
      save(saveobj,
        file = paste0("TEMP_SIM_", RunName, "-", seed, ".RData")
      )
    }
    ## update progress
    if (!verbose) {
      setTxtProgressBar(pb, t)
    }
    # estimator of finish
    if (verbose) {
      TimeEnd <- Sys.time()
      PercRem <- ((t_max - t) / t_max) * 100
      PercDone <- 100 - PercRem
      if (PercDone > 100) {
        PercDone <- 100
      }
      TimeElapsed <- difftime(TimeEnd, TimeStart, units = "secs")
      TimeRemaining <- (TimeElapsed / PercDone) * (100 - PercDone)
      # EstimatedFinish <- TimeStart+TimeRemaining
      cat("\r", sprintf("%05s", format(round(PercDone, 2), nsmall = 2)), "%  - ", paste("Estimated time remaining =", round(seconds_to_period(TimeRemaining), 0), "         "))
      flush.console()
    }
    if (nrow(ID_df) == 0) {
      warning("All species went extinct")
      break
    }
  }
  ID_ls <- c(ID_ls, list(ID_df))
  names(ID_ls)[length(ID_ls)] <- t
  unlink(paste0("TEMP_SIM_", RunName, "-", seed, ".RData"))
  return(ID_ls)
}
