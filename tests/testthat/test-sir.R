context("SIR functions")
library(gridExtra)

test_that("SIR from deSolve is correct", {

    T <- 100

    init_vals <- c(950, 50, 0)
    beta <- .1
    gamma <- .03
    N <- sum(init_vals)
    step <- 1

    ## Plot the SIR
    do_plot = TRUE
    results <- SIR(T, init_vals, beta,
                   gamma, step, SIR_inner, do_plot)

    ## extract the probabilities of transition
    probs <- extract_probs_sir(results, beta, gamma)
})


test_that("the CM simulation for the SIR", {
    L <- 5000 # number of runs
    S <- matrix(0, nrow = T+1, ncol = L)
    I <- matrix(0, nrow = T+1, ncol = L)
    R <- matrix(0, nrow = T+1, ncol = L)
    cm_sim_list <- list(S=S, I=I, R=R)
    cm_sim_list <- fill_inits_sir(init_vals, cm_sim_list)




    ## Run the simulation
    cm_sim_list <- run_cm_sir(T, init_vals,
                              probs, L, cm_sim_list)
    counts <- cm_sim_list$S


    ## Plotting the average
    cm_mean <- summary_sims(cm_sim_list, var_names = c("S", "I", "R")) * 100 / N
    cm_mean$time <- 0:(nrow(cm_mean)-1)
    df_melt <- reshape2::melt(cm_mean, id = "time", variable.name = "sim")
    ggplot2::ggplot(df_melt, ggplot2::aes(x=time, y=value * 100 / N, col=sim)) +
        ggplot2::geom_line(size=3)

})

test_that("AM functions", {
    L <- 5000
    T <- 100
    agents <- vector("list", L)
    agents <- initialize_agents(T, N, init_vals, agents)

    t <- proc.time()[3]
    out_agents <- run_am_sir(T, probs, L, agents)
    proc.time()[3] - t

    ## saveRDS(out_agents, "../../sims/out_agents.RDS")
    out_agents <- readRDS("../../sims/out_agents.RDS")
    agents_sims <- summarize_agents(out_agents)
    am_sim_list <- agents_sims

    ## Plotting the average
    am_mean <- summary_sims(am_sim_list, var_names = c("S", "I", "R")) * 100 / N
    am_mean$time <- 0:(nrow(am_mean) - 1)
    df_melt <- reshape2::melt(am_mean, id = "time", variable.name = "sim")
    ggplot2::ggplot(df_melt, ggplot2::aes(x=time, y=value * 100 / N, col=sim)) +
        ggplot2::geom_line(size=3)



})

test_that("overlays", {

    ## Average
    
    sum1 <- cm_mean
    sum2 <- am_mean
    plot_overlap(sum1, sum2, plot_dash = TRUE,
                 size = 2)
   ## ggplot2::ggsave("../../images/sir-mean.pdf",
   ##                 width=10,height=8)

    ## Variance
    ## Plotting the variance
    am_var <- summary_sims(am_sim_list, fxn=rowVar,
                           var_names = c("S", "I", "R"))
    cm_var <- summary_sims(cm_sim_list, fxn=rowVar,
                           var_names = c("S", "I", "R"))
    plot_overlap(cm_var, am_var,
                 summary_name = "Variance in states",
                 ylim=c(0,max(max(cm_var), max(am_var))),
                 plot_dash = FALSE, size = 2)
    ggplot2::ggsave("../../images/sir-var.pdf",
                    ## width=10,height=8)
    
})

test_that("plotting sir simulations", {

    ## CM
    ## An individual state
    g_s_draws <- plot_draws_sir(cm_sim_list$S,
                                beta, gamma,
                                N, L)
    print(g_s_draws)

    ## All at once
    cols <- c("Blues", "Reds", "Greens")
    cm_titles <- c("Susceptible",
                   "Infectious",
                   "Recovered")
    cm_symbols <-c("\\hat{S}(t)",
                   "\\hat{I}(t)",
                   "\\hat{R}")
    g_list_cm <- lapply(1:length(cm_sim_list),
                        function(ind){
                            plot_draws_sir(cm_sim_list[[ind]],
                                           beta,
                                           gamma,
                                           N, L,
                                           col = cols[ind],
                                           cat_title = cm_titles[ind],
                                           tex_symbol = cm_symbols[ind],
                                           approach = "CM")
                        })


    g_list_am <- lapply(1:length(am_sim_list),
                        function(ind){
                            plot_draws_sir(cm_sim_list[[ind]],
                                           beta,
                                           gamma,
                                           N, L,
                                           col = cols[ind],
                                           cat_title = cm_titles[ind],
                                           tex_symbol = cm_symbols[ind],
                                           approach = "AM")
                        })


    ## Plot png
    ## png("../../images/draws-sir.png", height=1000, width =800)
    ## do.call("grid.arrange", c(g_list_cm, g_list_am, ncol=2, as.table = FALSE))
    ## dev.off()

                 

})

