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
         ##          width=10,height=8)

    ## Variance
    ## Plotting the variance
    am_var <- summary_sims(am_sim_list, fxn=rowVar,
                           var_names = c("S", "I", "R"))
    cm_var <- summary_sims(cm_sim_list, fxn=rowVar,
                           var_names = c("S", "I", "R"))
    plot_overlap(cm_var, am_var,
                 summary_name = "Variance in states",
                 ylim=c(0,max(max(cm_var), max(am_var))),
                 plot_dash = FALSE, size = 2,
                 ylab = "Variance of state totals")
    ## ggplot2::ggsave("../../images/sir-var.pdf",
       ##          width=10,height=8)
    
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

test_that("Looking at s2ir2", {
    T <- 50
    init_vals <- c(250,500, 250, 0, 0)
    N <- sum(init_vals)
    beta1 <- .25
    beta2 <- .5
    gamma1 <- .05
    gamma2 <- .1
    step <- 1
    inner_fxn <- SIR2_inner
    results <- SIR2(T, init_vals, beta1,
                    beta2, gamma1, gamma2)
    head(results)
    ## Extract the probabilities of transition
    p <- extract_probs(results)


    ## Set up the CM list
    L <- 5000 # number of runs
    S1 <- matrix(0, nrow= T+1, ncol = L)
    S2 <- matrix(0, nrow= T+1, ncol = L)
    I <- matrix(0, nrow= T+1, ncol = L)
    R1 <- matrix(0, nrow= T+1, ncol = L)
    R2 <- matrix(0, nrow= T+1, ncol = L)
    cm_sim_list <- list(S1=S1, S2=S2, I=I, R1=R1, R2=R2)
    cm_sim_list <- fill_in_init_vals(init_vals, cm_sim_list)
    cm_sim_list <- run_cm(T, init_vals, p, L=L, cm_sim_list)


    ## Save it
    sim_list_all <- list(sim_list=cm_sim_list, beta1=beta1, beta2=beta2,
                     gamma1=gamma1, gamma2=gamma2, N=N, L=L)
    ## saveRDS(sim_list_all, "../../sims/cm-sim-list-sir2.RDS")

    ## read_cm <- readRDS("../../sims/cm-sim-list-sir2.RDS")
    cm_sim_list <- read_cm$sim_list
    ## Plotting all the draws
    cols <- c("Blues", "Oranges", "Reds", "Greens", "Purples")
    cm_titles <- c("Susceptible 1", "Susceptible 2",
                   "Infectious",
                   "Recovered 1", "Recovered 2")
    cm_symbols <-c("\\hat{S}_1(t)", "\\hat{S}_2(t)",
                   "\\hat{I}(t)",
                   "\\hat{R}_1(t)", "\\hat{R}_2(t)")
    g_list_cm <- lapply(1:length(cm_sim_list),
                        function(ind){
                            plot_draws_s2ir2(cm_sim_list[[ind]],
                                       beta1, beta2,
                                       gamma1, gamma2,
                                       N, L,
                                       col = cols[ind],
                                       cat_title = cm_titles[ind],
                                       tex_symbol = cm_symbols[ind])
                        })

    ## png("../../images/cm-sir2.png", height=1000)
    ## do.call("grid.arrange", c(g_list_cm, ncol=1))
    ## dev.off()

    
    ## Averages
    cm_mean <- summary_sims(cm_sim_list)


    g_avg_cm <- plot_summary(cm_mean, L =L, beta1=beta1, beta2=beta2,
                             gamma1=gamma1, gamma2=gamma2)
    g_avg_cm
    ## ggplot2::ggsave("../../images/cm-sir2-avg.pdf", height=6, width = 8)

    ## Variance

    cm_var <- summary_sims(cm_sim_list, fxn=rowVar)
    g_var_cm <- plot_summary(cm_var, sum_name = "Variance",
                             ylab="Variance within Compartment",
                             L=L, beta1=beta1, gamma1=gamma1,
                             beta2=beta2, gamma2=gamma2)
    g_var_cm
    ## ggplot2::ggsave("../../images/cm-sir2-var.pdf", height=6, width=8)

    
    ## Time distribution

    cols <- c("blue", "orange", "darkred", "darkgreen", "purple")
    cm_titles <- c("Susceptible 1", "Susceptible 2",
                   "Infectious",
                   "Recovered 1", "Recovered 2")
    cm_symbols <-c("\\hat{S}_1(t)", "\\hat{S}_2(t)",
                   "\\hat{I}(t)",
                   "\\hat{R}_1(t)", "\\hat{R}_2(t)")
    g_list_cm_times<- lapply(1:length(sim_list),
                             function(ind){
                                 plot_time_dist(sim_list[[ind]],
                                                beta1, beta2,
                                                gamma1, gamma2,
                                                N, L,
                                                col = cols[ind],
                                                cat_title = cm_titles[ind],
                                                tex_symbol = cm_symbols[ind])
                             })


    ###################################3
    ## The AM whee
    L <- 5000
    T <- 50
    agents <- vector("list", L)
    agents <- initialize_agents(T, N, init_vals, agents)

    t <- proc.time()[3]
    out_agents <- run_am(T, p, L, agents)
    proc.time()[3] - t

    ## saveRDS(out_agents, "../../sims/agents-sir2.RDS")

    ## out_agents <- readRDS("../../sims/agents-sir2.RDS")
    agents_sims <- summarize_agents(out_agents, n_states = 5)
    am_sim_list <- agents_sims



    ## Plotting all the draws
    cols <- c("Blues", "Oranges", "Reds", "Greens", "Purples")
    cm_titles <- c("Susceptible 1", "Susceptible 2",
                   "Infectious",
                   "Recovered 1", "Recovered 2")
    cm_symbols <-c("\\hat{S}_1(t)", "\\hat{S}_2(t)",
                   "\\hat{I}(t)",
                   "\\hat{R}_1(t)", "\\hat{R}_2(t)")
    g_list_am <- lapply(1:length(am_sim_list),
                        function(ind){
                            plot_draws_s2ir2(am_sim_list[[ind]],
                                       beta1, beta2,
                                       gamma1, gamma2,
                                       N, L,
                                       col = cols[ind],
                                       cat_title = cm_titles[ind],
                                       tex_symbol = cm_symbols[ind],
                                       approach = "AM")
                        })

    ## png("../../images/am-sir2.png", height=1000)
    ## do.call("grid.arrange", c(g_list_am, ncol=1))
    ## dev.off()

    ## Side by side
    ## png("../../images/cm-am-s2ir2.png", height=1000, width=1000)
    ## do.call("grid.arrange", c(g_list_cm, g_list_am, ncol = 2, as.table = FALSE))
    ## dev.off()

    ## Averages
    am_mean <- summary_sims(am_sim_list)


    g_avg_am <- plot_summary(am_mean, approach = "AM",
                             L=L, beta1=beta1, gamma1=gamma1,
                             beta2=beta2, gamma2=gamma2)
    g_avg_am
    ## ggplot2::ggsave("../../images/am-sir2-avg.pdf", height=6, width = 8)

    ## Variance

    am_var <- summary_sims(am_sim_list, fxn=rowVar)
    g_var_am <- plot_summary(am_var, sum_name = "Variance",
                             ylab="Variance within Compartment",
                             approach = "AM",
                             L=L, beta1=beta1, gamma1=gamma1,
                             beta2=beta2, gamma2=gamma2)
    g_var_am
    ## ggplot2::ggsave("../../images/am-sir2-var.pdf", height=6, width = 8)


    ## Time distribution

    cols <- c("blue", "orange", "darkred", "darkgreen", "purple")
    cm_titles <- c("Susceptible 1", "Susceptible 2",
                   "Infectious",
                   "Recovered 1", "Recovered 2")
    cm_symbols <-c("\\hat{S}_1(t)", "\\hat{S}_2(t)",
                   "\\hat{I}(t)",
                   "\\hat{R}_1(t)", "\\hat{R}_2(t)")
    g_list_am_times <- lapply(1:length(am_sim_list),
                             function(ind){
                                 plot_time_dist(am_sim_list[[ind]],
                                                beta1, beta2,
                                                gamma1, gamma2,
                                                N, L,
                                                col = cols[ind],
                                                cat_title = cm_titles[ind],
                                                tex_symbol = cm_symbols[ind])
                             })

    png("../../images/am-sir2-time.png", height=1000)
    do.call("grid.arrange", c(g_list_am_times, ncol=1))
    dev.off()


    ## major rearrangement of time steps

    ## s1 and s2
    pdf("../../images/am-cm-s1s2.pdf", height=10, width = 12)
    grid.arrange(g_list_cm_times[[1]], g_list_am_times[[1]],
                 g_list_cm_times[[2]], g_list_am_times[[2]],
                 ncol = 2)
    dev.off()


    ## I
    pdf("../../images/am-cm-i.pdf", height=5, width = 12)
    grid.arrange(g_list_cm_times[[3]], g_list_am_times[[3]],
                 ncol = 2)
    dev.off()


    ## r1
    pdf("../../images/am-cm-r1.pdf", height=5, width = 12)
    grid.arrange(g_list_cm_times[[4]], g_list_am_times[[4]],
                 ncol = 2)
    dev.off()


    ## r2
    pdf("../../images/am-cm-r2.pdf", height=5, width = 12)
    grid.arrange(g_list_cm_times[[5]], g_list_am_times[[5]],
                 ncol = 2)
    dev.off()


})
