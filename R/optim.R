#' Compute sum of squared error for the SEIR model
#'
#' @param params c(beta, gamma, alpha)
#' @param data output from the SEIR() function
#' @param fxn function to be optimized. Default is SEIR
#' @param init_vals c(X, E, Y, Z)
#' @param step time step.  Default is 1
#' @return sse across X, E, and Y
SSE_SEIR <- function(params, data,
                     fxn=SEIR, init_vals=c(9950,0,50,0),
                     step=1){
    data_sub <- data[, -c(1, ncol(data))] # take out time column and last column
    T <- max(data[, 1])
    ## Fit model with params
    beta <- params[1]
    gamma <- params[2]
    alpha <- params[3]
    new_data <- fxn(T, init_vals, beta, gamma, alpha,
                    step=1, SEIR_inner, do_plot = FALSE)
    new_sub <- new_data[, -c(1, ncol(new_data))]
    sse <- sum(sapply(1:ncol(new_sub), function(col){
        sum((new_data[, col] - data[, col])^2)
    }))
    return(sse)
    
    

}


#' Compute sum of squared error for the SIR model
#'
#' @param params c(beta, gamma, alpha)
#' @param data output from the SIR() function
#' @param fxn function to be optimized. Default is SIR
#' @param init_vals c(X, Y, Z)
#' @param step time step.  Default is 1
#' @return sse across X, and Y
SSE_SIR <- function(params, data,
                     fxn=SIR, init_vals=c(9950,50,0),
                     step=1){
    data_sub <- data[, -c(1, ncol(data))] # take out time column and last column
    T <- max(data[, 1])
    ## Fit model with params
    beta <- params[1]
    gamma <- params[2]
    new_data <- fxn(T, init_vals, beta, gamma, 
                    step=1, SIR_inner, do_plot = FALSE)
    new_sub <- new_data[, -c(1, ncol(new_data))]
    sse <- sum(sapply(1:ncol(new_sub), function(col){
        sum((new_data[, col] - data[, col])^2)
    }))
    return(sse)
    
    

}
