

#' Find the time points corresponding to differences in  matrix of simulations
#'
#' @param val the value of the change
#' @param x matrix where entry ij is the jth simulation's compartment value at time i
#' @return data frame of the time and frequency with which the change is equal to val
change_val <- function(val, x, L){
    mat <- apply(x, 2, function(y){
        ifelse(diff(y) == val, 1, 0)
    })
    tab <- rowSums(mat) / L
    df <- data.frame(Time=1:(length(tab)), Percent=tab, Change = val) 
    return(df)

}




#' Set initial values in the simulation list
#'
#' @param init_vals  vectof of initial values for S1, S2, I, R1, R2
#' @param sim_list list of results
#' @return sim_list with init_vals in the right places
fill_in_init_vals <- function(init_vals, sim_list){
    sim_list$S1[1, ] <- init_vals[1]
    sim_list$S2[1, ] <- init_vals[2]
    sim_list$I[1, ] <- init_vals[3]
    sim_list$R1[1, ] <- init_vals[4]
    sim_list$R2[1, ] <- init_vals[5]
    
    return(sim_list)
}




#' @param T - total time
#' @param N - total number of agents
#' @param p matrix of probabilities
#' @param L L total number of runs
#' @param agents list of L x T+1 x N
#' @return updated agents
run_am <- function(T, p, L=1, agents){
    
    for(run in 1:L){
        for(tt in 0:(T-1)){
            bern1 <- rbinom(1,1, p$S1[tt+1])
            bern2 <- rbinom(1,1, p$S2[tt+1])
            Z <- rmultinom(1, 1, c(gamma1, gamma2, 1 - gamma1 - gamma2))
            agents[[run]][tt+2, ] <- ifelse(agents[[run]][tt+1,] == 1,
                                            agents[[run]][tt+1,] + bern1 * 2,
                                     ifelse(agents[[run]][tt+1,] == 2,
                                            agents[[run]][tt+1,] + bern2 * 1,
                                     ifelse(agents[[run]][tt+1,] == 3,
                                            agents[[run]][tt+1,] + sum(Z[1:2] * 1:2),
                                            agents[[run]][tt+1,])))
            ## for(n in 1:N){
            ##     agent_current <- agents[[run]][tt+1, n]
            ##     ## Update the agent based on its current state
            ##     ## 1 = S1, 2=s2, 3=i, r=r1, 5=r2
            ##     if(agent_current== 1){
            ##         agents[[run]][tt+2, n] <- agent_current + bern1 * 2
            ##     } else if(agent_current == 2){
            ##         agents[[run]][tt+2, n] <- agent_current + bern2 * 1
            ##     } else if (agent_current == 3){
            ##         agents[[run]][tt+2, n] <- agent_current + sum(Z[1:2] * 1:2)
            ##     } else {
            ##         agents[[run]][tt+2, n] <- agent_current
            ##     }

            ## }
        }
    }
    return(agents)


}

#' Run the CM of the S^2IR^2 model
#'
#' @param T - total time
#' @param init_vals vectof of initial values for S1, S2, I, R1, R2
#' @param p matrix of probabilities
#' @param L L total number of runs
#' @param sim_list list of results
#' @return updated sim_list
run_cm <- function(T, init_vals, p, L=1, sim_list){
    for(run in 1:L){
        for(tt in 0:(T-1)){

            gamma1 <- p$R1[1]
            gamma2 <- p$R2[1]
            ## S1
            sim_list$S1[tt+2, run] <- sim_list$S1[tt+1, run] * (1 - rbinom(1, 1, p$S1[tt+1]))
            sim_list$S2[tt+2, run] <- sim_list$S2[tt+1, run] * (1 - rbinom(1, 1, p$S2[tt+1]))
            Z <- rmultinom(1, 1, c(gamma1, gamma2, 1 - gamma1 - gamma2))
            sim_list$R1[tt+2, run] <- sim_list$R1[tt+1, run] + sim_list$I[tt+1, run] * Z[1]
            sim_list$R2[tt+2, run] <- sim_list$R2[tt+1, run] + sim_list$I[tt+1, run] * Z[2]
            sim_list$I[tt+2, run] <- N - sim_list$S1[tt+2, run] -
                sim_list$S2[tt+2, run] -
                sim_list$R1[tt+2, run] -
                sim_list$R2[tt+2, run] 
        }
    }
    return(sim_list)
    
}

#' Extract the probabilities from the ODEs for updating
#'
#' @param ode_results object returned from SIR2. Specifically the columns are t, S1, S2, I, R1, R2.
#' @return probabilities
extract_probs <- function(ode_results){
    p <- ode_results[, -1]
    p$S1 = beta1 * ode_results$I / (N )
    p$S2 = beta2 * ode_results$I / (N )
    p$R1 <- gamma1
    p$R2 <- gamma2
    return(p)
    
}


#' Wrapper for DeSolve ODE solver
#'
#' @param T total time
#' @param init_vals vector of length 5, (S10, S20, I0, R10, R20)
#' @param beta1 between 0 and 1
#' @param beta2 between 0 and 1
#' @param gamma1 between 0 and 1
#' @param gamma2 between 0 and 1 
#' @param step what times to return.  Default is 1
#' @param inner_fxn default is SIR_inner which does the basic SIR model
#' @param do_plot Default is FALSE. will make a pretty plot for testing
#' @return data frame of values
SIR2 <- function(T, init_vals, beta1, beta2,
                 gamma1, gamma2,
                step = 1, inner_fxn = SIR2_inner,
                do_plot = TRUE){

    params <- c(beta1=beta1, beta2=beta2,
                gamma1=gamma1, gamma2=gamma2)
    state <- c(X1 = init_vals[1],
               X2 = init_vals[2],
               Y = init_vals[3],
               Z1 = init_vals[4],
               Z2 = init_vals[5])
    times <- seq(0, T, by = step)

    
    ## use ODE to integrate and get the simulation
    ode_results <- deSolve::ode(state, times, func = inner_fxn, parms = params,
                                maxsteps = 50000, method = "rk4")
    ## Plot results
    if(do_plot){
        plot(1, 1, type = "n", xlim = range(ode_results[, 1]), ylim = range(ode_results[, 2:6]),
             xlab = "Time", ylab = "# of Individuals",
             main = paste0("S2IR2 curve \n beta1 = ", beta1,
                          "; beta2 = ", beta2,
                          "; gamma1 = ", gamma1,
                          "; gamma2 = ", gamma2
                          ))
        lines(ode_results[, 1], ode_results[, 2], col = "blue", lwd = 2) ## S1
        lines(ode_results[, 1], ode_results[, 3], col = "orange", lwd = 2) ## S2
        lines(ode_results[, 1], ode_results[, 4], col = "red", lwd = 2) ## I
        lines(ode_results[,1], ode_results[,5], col = "darkgreen", lwd=2) ## R1
        lines(ode_results[,1], ode_results[,6], col = "purple", lwd=2) ## R2
        legend("topright", c("# Susceptible 1", "# Susceptible 2",
                             "# Infectious",
                             "# Recovered 1", "# Recovered 2"),
               col = c(4, "orange", 2, "darkgreen", "purple"), lwd = 2)
    }
    colnames(ode_results) <- c("t", "S1", "S2", "I", "R1", "R2")
    results <- as.data.frame(ode_results)
    return(results) 
}

#' Making the ODE function for deSolve
#'
#' @param t a vector of times.  Should be strictly monotonically increasing.
#' @param state the initial conditions, a vector. (X1, X2, Y, Z1, Z2)
#' @param params the parameters, here beta1, beta2, gamma1, gamma2
#' @return the function used in ode
SIR2_inner <- function(t, state, params){
    N <- sum(state)
    with(as.list(c(state, params)), {
        dX1 <- - params[1] * X1 * Y / N
        dX2 <- - params[2] * X2 * Y / N
        dZ1 <- params[3] * Y
        dZ2 <- params[4] * Y
        dY <- -dX1 - dX2 - dZ1 - dZ2
        list(c(dX1, dX2, dY, dZ1, dZ2))
    })
}

