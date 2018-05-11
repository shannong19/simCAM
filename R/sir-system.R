
#' Wrapper for deSolve ODE solver for the deterministic SIR-CM model
#'
#' @param T total time
#' @param init_vals vector of length 3, (S0, I0, R0)
#' @param beta between 0 and 1 - infection rate
#' @param gamma between 0 and 1 - recovery rate
#' @param step what times to return.  Default is 1
#' @param inner_fxn default is SIR_inner which does the basic SIR model
#' @param do_plot logical. Default is TRUE.  Makes a plot fo the SIR curves
#' @return data frame of values
SIR <- function(T, init_vals, beta,
                 gamma,
                 step = 1, inner_fxn = SIR_inner,
                 do_plot = TRUE){

  params <- c(beta=beta,
              gamma=gamma)
  state <- c(X = init_vals[1],
             Y = init_vals[2],
             Z = init_vals[3])
  times <- seq(0, T, by = step)


  ## use ODE to integrate and get the simulation
  ode_results <- deSolve::ode(state, times, func = inner_fxn, parms = params,
                              maxsteps = 50000, method = "rk4")
  ## Plot results
  if(do_plot){
    plot(1, 1, type = "n", xlim = range(ode_results[, 1]), ylim = range(ode_results[, 2:4]),
         xlab = "Time", ylab = "# of Individuals",
         main = paste0("SIR curve \n beta = ", beta,
                       "; gamma = ", gamma
         ))
    lines(ode_results[, 1], ode_results[, 2], col = "blue", lwd = 2) ## S
    lines(ode_results[, 1], ode_results[, 3], col = "red", lwd = 2) ## I
    lines(ode_results[,1], ode_results[,4], col = "darkgreen", lwd=2) ## R
    legend("topright", c("# Susceptible",
                         "# Infectious",
                         "# Recovered "),
           col = c(4, 2, "darkgreen"), lwd = 2)
  }
  colnames(ode_results) <- c("t", "S", "I", "R")
  results <- as.data.frame(ode_results)
  return(results)
}

#' Making the SIR function for deSolve
#'
#' @param t a vector of times.  Should be strictly monotonically increasing.
#' @param state the initial conditions, a vector. (X, Y, Z)= (S, I, R)
#' @param params the parameters, here beta and gamma
#' @return the function used in ode()
SIR_inner <- function(t, state, params){
  N <- sum(state)
  with(as.list(c(state, params)), {
    dX <- - params[1] * X * Y / N
    dZ <- params[2] * Y
    dY <- -dX - dZ
    list(c(dX, dY, dZ))
  })
}


#' Extract the probabilities from the SIR for updating individuals/agents
#'
#' @param ode_results object returned from SIR. Specifically the columns are t, S, I, R
#' @return probabilities at each time step
extract_probs_sir <- function(ode_results, beta, gamma){
  p <- ode_results[, -1]
  p$S <- beta * ode_results$I / (N)
  p$R <- gamma
  return(p)

}



#' Run the CM of the SIR model
#'
#' @param T - total time
#' @param init_vals vectof of initial values for S, I, R
#' @param p matrix of probabilities
#' @param L L total number of runs
#' @param sim_list list of results
#' @return updated sim_list
run_cm_sir <- function(T, init_vals, p, L=1, sim_list){
  for(run in 1:L){
    for(tt in 0:(T-1)){

      gamma <- p$R[1]
      S_cur <- sim_list$S[tt+1, run]
      I_cur <- sim_list$I[tt+1, run]
      Z_cur <- sim_list$R[tt+1, run]
      ## S1
      sim_list$S[tt+2, run] <-  S_cur - rbinom(1, S_cur, p$S[tt+1])
      sim_list$R[tt+2, run] <- sim_list$R[tt+1, run] + rbinom(1, I_cur, gamma)
      sim_list$I[tt+2, run] <- N - sim_list$S[tt+2, run] - sim_list$R[tt+2, run]

    }
  }
  return(sim_list)

}

#' Set initial values in the simulation list for the SIR model
#'
#' @param init_vals  vectof of initial values for S, I, R
#' @param sim_list list of results
#' @return sim_list with init_vals in the right places
fill_inits_sir <- function(init_vals, sim_list){
    sim_list$S[1, ] <- init_vals[1]
    sim_list$I[1, ] <- init_vals[2]
    sim_list$R[1, ] <- init_vals[3]
    return(sim_list)
}


#' Compute the row variances
#'
#' @param mat matrix
#' @return vector of row variances
rowVar <- function(mat){
    t(apply(mat, 1, var))
}

#' Compute the row average
#'
#' @param mat matrix
#' @return vector of row averages
rowAvg <- function(mat){
    t(apply(mat, 1, mean))
}

#' Summarize the simulations from CM or AM
#'
#' @param sim_list - output from run_cm_* or run_am_*
#' @param fxn - summary function for each compartment.  Default is the rowAvg, which computes the mean number in each compartment over all the simulations
#' @param var_names - What we call the resulting columns in the new df
#' @return data frame for plotting of the summaries of the compartments
summary_sims <- function(sim_list, fxn=rowAvg,
                         var_names = c("S1", "S2", "I", "R1", "R2")){
    out <- t(plyr::laply(sim_list,
                   .fun=fxn))
    out <- as.data.frame(out)
    colnames(out) <- var_names
    return(out)
}


#' Run the agent-based model for the SIR simulatio
#'
#' @param T - total time
#' @param N - total number of agents
#' @param p matrix of probabilities
#' @param L L total number of runs
#' @param agents list of L x T+1 x N
#' @return updated agents
run_am_sir <- function(T, p, L=1, agents){
    gamma <- p$R[1]
    for(run in 1:L){
        for(tt in 0:(T-1)){
            agent_current <- agents[[run]][tt+1,]
            StoI <- rbinom(N, 1, p$S[tt+1])
            ItoR <- rbinom(N, 1, gamma)
            agents[[run]][tt+2, ] <- ifelse(agent_current == 1,
                                             agent_current + StoI,
                                      ifelse(agent_current == 2,
                                             agent_current + ItoR,
                                             agent_current)) 
        }
        
    }
    return(agents)


}

#' Initialize agents for the AM
#'
#' @param T - total time
#' @param N - total number of agents
#' @param init_vals vectof of initial values for S1, S2, I, R1, R2
#' @param p matrix of probabilities
#' @param L L total number of runs
#' @param agents list of L x T+1 x N
#' @return agents initialized so agents are in right place
initialize_agents <- function(T, N, init_vals, agents){
    init_vec <- rep(1:length(init_vals), times = init_vals)
    agents <- lapply(agents, function(x){
        mat <- matrix(0, nrow = T+1, ncol = N)
        mat[1, ] <- init_vec
        return(mat)
    })
    return(agents)
}


#' Put the agents in a format to look at summary stats
#'
#' @param agents list of L x T+1 x N
#' @param n_states number of states (e.g. 3 for SIR)
summarize_agents <- function(agents, n_states=3){
    sims <- lapply(1:n_states, function(x){
        counts <- sapply(agents, function(mat){
            apply(mat, 1, function(row){
                length(which(row == x))
            })
        })
        return(counts)
    })

}
