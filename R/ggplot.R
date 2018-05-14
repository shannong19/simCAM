
#' Plot the functions to overlap
#'
#' @param sum1 summary df of the first set of  simulations
#' @param sum2 summary df of the second set of simulations
#' @param col1 color palette for first set of curves
#' @param col2 color palette for second set of curves
#' @param cat_names1 curve names for first set
#' @param cat_names2 curve names for second set
#' @param beta infection parameter
#' @param gamma second parameter
#' @param total number of individuals
#' @param L total number of runs
#' @param approach1 name of first approach
#' @param approach2 name of second approach
#' @param summary_name title of graph
#' @param xlab xlab title
#' @param ylab ylab title
#' @param ylim range of values for y to plot
#' @param size line size
#' @param plot_dash logical. if true, plot dashed line for second set of lines, otherwise use solid
plot_overlap <- function(sum1, sum2, col1=c("blue", "darkred", "darkgreen"),
                         col2 = c("lightblue", "darksalmon", "olivedrab1"),
                         cat_names1 = c("S", "I", "R"),
                         cat_names2 = c("S", "I", "R"),
                         beta=.1, gamma=.03,
                         N=1000, L=5000,
                         tex_symbol = c("\\hat{S}(t)", 
                                        "\\hat{I}(t)",
                                        "\\hat{R}(t)"),
                         approach1="CM",
                         approach2="AM",
                         summary_name = "Mean proportion in states",
                         xlab = "Time",
                         ylab = "% of population",
                         ylim =c(0,100),
                         size = 3,
                         plot_dash = TRUE){
                             
    df <- data.frame(time=0:(nrow(sum1)-1), sum1[, 1:3], sum2[, 1:3])
    cat_names <- c(paste0(cat_names1, "-", approach1), paste0(cat_names2, "-", approach2))
    colnames(df) <- c("time", cat_names)
    df_melt <- reshape2::melt(df, id="time", variable.name="Type")
    if(plot_dash){
        df_melt$type <- ifelse(grepl("CM", df_melt$Type), 1, 2)
    } else {
        df_melt$type <- 1  # no dash
    }
    cols <- c(col1, col2)

    g <- ggplot2::ggplot(data = df_melt, ggplot2::aes(x = time, y = value, colour=Type,
                                                      linetype=factor(type), group=Type))+
        ggplot2::geom_line(size=size) +
        ggplot2::scale_color_manual(values = cols, labels = cat_names) +
        ggplot2::guides(linetype=FALSE) +
        my_theme() +
        ggplot2::labs(title = summary_name,
                      subtitle = latex2exp::TeX(sprintf("%d agents; %d runs; $\\beta = %.2f$; $\\gamma = %.2f$",
                                                        N, L,  beta, gamma)),
                      x=xlab, y=ylab) +
        ggplot2::ylim(ylim)
                      
                      

    print(g)
    return(g)

}



#' My custom ggplot background
#'
#' @return ggplot theme
my_theme <- function(){
    ggplot2::theme_bw() + ggplot2::theme(
                     axis.text.x = ggplot2::element_text(size = 16, family = "Palatino"),
                     axis.text.y= ggplot2::element_text(size = 16, family = "Palatino"),
                     axis.title.x= ggplot2::element_text(size = 18, family = "Palatino"),
                     axis.title.y= ggplot2::element_text(size = 18, family = "Palatino"),
                     plot.title = ggplot2::element_text(size = 24, family = "Palatino"),
                     legend.title = ggplot2::element_text(size = 20, family = "Palatino"),
                     legend.text = ggplot2::element_text(family = "Palatino", size = 16),
                     legend.key.size = ggplot2::unit(3, "line"),
                     plot.subtitle = ggplot2::element_text(size=16, family = "Palatino")
                     )
}



#' Plot individual draws for a single state for the SIR
#'
#' @param counts data frame of simulations and counts for a state
#' @param beta infection parameter between 0 and 1
#' @param gamma recovery parameter between 0 and 1
#' @param N number of total individuals
#' @param L number of total simulations
#' @param col color pallette used in brewer.pal from RColorBrewer
#' @param cat_title title of the state (e.g. Infectious)
#' @param tex_symbol latex symbol associated with cat_title
#' @param approach either CM or AM
plot_draws_sir <- function(counts,
                       beta,
                       gamma,
                       N,
                       L,
                       col="Blues",
                       cat_title = "Susceptible",
                       tex_symbol = "\\hat{S}(t)",
                       approach = "CM" ){
    
    df <- data.frame(time = 0:(nrow(counts)-1), counts)
    df_melt <- reshape2::melt(df, id = "time", variable.name="sim")
    my_cols <- RColorBrewer::brewer.pal(9, col)

    ## The GGplot
    g <- ggplot2::ggplot(data = df_melt, ggplot2::aes(x = time,
                                                     y= value * 100/N, group = sim)) +
        ggplot2::geom_line(size = 1,
                           colour = rep(my_cols, each = nrow(counts),
                                        length.out = nrow(df_melt))) +
    ggplot2::ylim(0,100) +
    my_theme() +
    ggplot2::xlab("Time") + ggplot2::ylab(paste0("% ", cat_title)) +
    ggplot2::labs(title = latex2exp::TeX(sprintf("Simulations of %s",
                                                 tex_symbol)),
          subtitle = latex2exp::TeX(sprintf("$\\beta = %.2f$; $\\gamma = %.2f$; %s approach",
                                 beta, gamma, approach))) +
    ggplot2::theme(legend.position = "none")
    return(g)


}
