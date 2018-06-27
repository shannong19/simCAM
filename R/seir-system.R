
#' Wrapper for deSolve ODE solver for the deterministic SEIR-CM model
#'
#' @param T total time
#' @param init_vals vector of length 3, (S0, I0, R0)
#' @param beta between 0 and 1 - infection rate
#' @param gamma between 0 and 1 - recovery rate
#' @param alpha between 0 and 1 - rate of exposed to infectious
#' @param step what times to return.  Default is 1
#' @param inner_fxn default is SEIR_inner which does the basic SEIR model
#' @param do_plot logical. Default is TRUE.  Makes a plot fo the SEIR curves
#' @return data frame of values
SEIR <- function(T, init_vals, beta,
                 gamma, alpha,
                 step = 1, inner_fxn = SEIR_inner,
                 do_plot = TRUE){

    params <- c(beta=beta,
                gamma=gamma,
                alpha = alpha)
    
    state <- c(X = init_vals[1],
               E = init_vals[2],
               Y = init_vals[3],
               Z = init_vals[4])
    times <- seq(0, T, by = step)


    ## use ODE to integrate and get the simulation
    ode_results <- deSolve::ode(state, times, func = inner_fxn, parms = params,
                                maxsteps = 50000, method = "rk4")
    ## Plot results
    if(do_plot){
        plot(1, 1, type = "n", xlim = range(ode_results[, 1]), ylim = range(ode_results[, 2:5]),
             xlab = "Time", ylab = "# of Individuals",
             main = paste0("SEIR curve \n beta = ", beta,
                           "; gamma = ", gamma, "; alpha = ", alpha
                           ))
        lines(ode_results[, 1], ode_results[, 2], col = "blue", lwd = 2) ## S
        lines(ode_results[, 1], ode_results[, 3], col = "orange", lwd = 2) ## E
        lines(ode_results[,1], ode_results[,4], col = "red", lwd=2) ## I
        lines(ode_results[,1], ode_results[,5], col = "darkgreen", lwd=2) ## R
        legend("topright", c("# Susceptible",
                             "# Exposed",
                             "# Infectious",
                             "# Recovered "),
               col = c(4, "orange",  2, "darkgreen"), lwd = 2)
    }
    colnames(ode_results) <- c("t", "S", "E", "I", "R")
    results <- as.data.frame(ode_results)
    return(results)
}

#' Making the SIR function for deSolve
#'
#' @param t a vector of times.  Should be strictly monotonically increasing.
#' @param state the initial conditions, a vector. (X, Y, Z)= (S, I, R)
#' @param params the parameters, here beta and gamma
#' @return the function used in ode()
SEIR_inner <- function(t, state, params){
    N <- sum(state)
    with(as.list(c(state, params)), {
        dX <- - params[1] * X * Y / N
        dE <- -dX - params[3] * E
        dZ <- params[2] * Y
        dY <- -dX - dZ - dE
        list(c(dX, dE, dY, dZ))
    })
}
