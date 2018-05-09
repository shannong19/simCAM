
#' Wrapper for deSolve ODE solver for the deterministic SIR-CM model
#'
#' @param T total time
#' @param init_vals vector of length 3, (S0, I0, R0)
#' @param beta between 0 and 1 - infection rate
#' @param gamma between 0 and 1 - recovery rate
#' @param step what times to return.  Default is 1
#' @param inner_fxn default is SIR_inner which does the basic SIR model
#' @return data frame of values
SIR<- function(T, init_vals, beta1, beta2,
                 gamma1, gamma2,
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
    dY <- -dX - dY
    list(c(dX, dY, dZ))
  })
}


#' Extract the probabilities from the SIR for updating individuals/agents
#'
#' @param ode_results object returned from SIR. Specifically the columns are t, S, I, R
#' @return probabilities at each time step
extract_probs <- function(ode_results, beta, gamma){
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
      S_cur <- sim_list$S[tt+1]
      I_cur <- sim_list$I[tt+1]
      Z_cur <- sim_list$R[tt+1]
      ## S1
      sim_list$S[tt+2, run] <-  S_cur - rbinom(S_cur, 1, p$S[tt+1])
      sim_list$R1[tt+2, run] <- sim_list$R1[tt+1, run] + rbinom(I_cur, gamma)
      sim_list$I[tt+2, run] <- N - sim_list$S[tt+2, run] - sim_list$R[tt+2, run]

    }
  }
  return(sim_list)

}

