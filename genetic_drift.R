#==============================================================================
#    Genetic-Drift.R : Genetic Drift Simulator
#    Copyright (C) 2023  Bruno Toupance <bruno.toupance@mnhn.fr>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by    
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#==============================================================================


library("ggplot2")
library("colorspace")
# library("scales")
# library("MASS")


#==============================================================================
# check_integer
#==============================================================================
check_integer <- function(
        x = 0,
        label = "x", 
        min_x = -Inf, 
        max_x = +Inf, 
        min_inc = FALSE,
        max_inc = FALSE,
        rtn = list(check_flag = TRUE, check_msg = ""))
{
    if (is.numeric(x)) {
        if (x == round(x)) {
            out_of_bounds_flag <- FALSE
            if (min_inc) {
                if (x <= min_x) {
                    out_of_bounds_flag <- TRUE
                }
            } else {
                if (x < min_x) {
                    out_of_bounds_flag <- TRUE
                }
            }
            if (max_inc) {
                if (x >= max_x) {
                    out_of_bounds_flag <- TRUE
                }
            } else {
                if (x > max_x) {
                    out_of_bounds_flag <- TRUE
                }
            }
            if (out_of_bounds_flag) {
                rtn$check_flag <- FALSE
                rtn$check_msg <- sprintf("%s\nFAIL: [%s] out of bounds", 
                                         rtn$check_msg, label)
            }
        } else {
            rtn$check_flag <- FALSE
            rtn$check_msg <- sprintf("%s\nFAIL: [%s] not integer", 
                                     rtn$check_msg, label)
        }
    } else {
        rtn$check_flag <- FALSE
        rtn$check_msg <- sprintf("%s\nFAIL: [%s] not numeric", 
                                 rtn$check_msg, label)
    }

    return(rtn)
}


#==============================================================================
# check_real
#==============================================================================
check_real <- function(
        x = 0.0,
        label = "x", 
        min_x = -Inf, 
        max_x = +Inf, 
        min_inc = FALSE,
        max_inc = FALSE,
        rtn = list(check_flag = TRUE, check_msg = "")) 
{
    if (is.numeric(x)) {
        out_of_bounds_flag <- FALSE
        if (min_inc) {
            if (x <= min_x) {
                out_of_bounds_flag <- TRUE
            }
        } else {
            if (x < min_x) {
                out_of_bounds_flag <- TRUE
            }
        }
        if (max_inc) {
            if (x >= max_x) {
                out_of_bounds_flag <- TRUE
            }
        } else {
            if (x > max_x) {
                out_of_bounds_flag <- TRUE
            }
        }
        if (out_of_bounds_flag) {
            rtn$check_flag <- FALSE
            rtn$check_msg <- sprintf("%s\nFAIL: [%s] out of bounds", 
                                     rtn$check_msg, label)
        }
    } else {
        rtn$check_flag <- FALSE
        rtn$check_msg <- sprintf("%s\nFAIL: [%s] not numeric", 
                                 rtn$check_msg, label)
    }

    return(rtn)
}



#==============================================================================
# check_parameters
#==============================================================================
check_parameters <- function(
        nb_ind = 10, 
        ini_p = 0.5, 
        nb_gen = 15, 
        nb_rep = 10)
{
    check_list <- list(check_flag = TRUE, check_msg = "")

    check_list <- check_integer(nb_ind, label = "population size", 
                                min_x = 1, rtn = check_list)
    
    check_list <- check_integer(nb_gen, label = "number of generations", 
                                min_x = 1, rtn = check_list)
    
    check_list <- check_real(ini_p, label = "initial frequency", 
                             min_x = 0, max_x = 1, rtn = check_list)
    
    check_list <- check_integer(nb_rep, label = "number of repetitions", 
                                min_x = 1, rtn = check_list)

    return(list(msg = check_list$check_msg, 
                flag = check_list$check_flag, 
                nb_ind = nb_ind, 
                ini_p = ini_p, 
                nb_gen = nb_gen, 
                nb_rep = nb_rep))
}





#==============================================================================
# simulate_data
#  Simulate frequency evolution under drift for allele A
#==============================================================================
simulate_data <- function(
        nb_ind = 10, 
        ini_p = 0.5, 
        nb_gen = 15, 
        nb_rep = 10, 
        diploid_flag = TRUE) 
{
    param_checking <- check_parameters(nb_ind, ini_p, nb_gen, nb_rep)
    # print(param_checking)

    if (param_checking$flag) {
        if (diploid_flag) {
            nb_copies <- 2 * nb_ind
        } else {
            nb_copies <- nb_ind
        }
        freq_mat <- matrix(NA, nrow = nb_rep, ncol = nb_gen + 1)
        freq_mat[, 1] <- ini_p

        count_mat <- matrix(0, nrow = nb_copies + 1, ncol = nb_gen + 1)
        count_mat[round(ini_p * nb_copies) + 1, 1] <- nb_rep

        mean_p <- rep(NA, times = nb_gen + 1)
        var_p <- rep(NA, times = nb_gen + 1)
        nb_A_fixed <- rep(NA, times = nb_gen + 1)
        nb_A_lost <- rep(NA, times = nb_gen + 1)
        nb_polym <- rep(NA, times = nb_gen + 1)
        time_A_fixed <- c()
        time_A_lost <- c()

        j_max <- nb_gen + 1
        for (i_rep in 1:nb_rep) {
            i <- i_rep
            p_A <- ini_p
            gen <- 1
            while (gen <= nb_gen) {
                nb_A <- rbinom(1, nb_copies, p_A)
                p_A <- nb_A / nb_copies
                j <- gen + 1
                k <- nb_A + 1
                if (nb_A == nb_copies) {
                    time_A_fixed <- append(time_A_fixed, gen)
                    freq_mat[i, j:j_max] <- p_A
                    count_mat[k, j:j_max] <- count_mat[k, j:j_max] + 1
                    gen <- nb_gen
                } else {
                    if (nb_A == 0) {
                        time_A_lost <- append(time_A_lost, gen)
                        freq_mat[i, j:j_max] <- p_A
                        count_mat[k, j:j_max] <- count_mat[k, j:j_max] + 1
                        gen <- nb_gen
                    } else {
                        freq_mat[i, j] <- p_A
                        count_mat[k, j] <- count_mat[k, j] + 1
                    }
                }
                gen <- gen + 1
            }
        }
        
        rownames(count_mat) <- 0:nb_copies
        colnames(count_mat) <- 0:nb_gen

        nb_A_fixed <- apply(
            freq_mat, 
            MARGIN = 2, 
            FUN = function(x) { 
                return(sum(x == 1)) 
            })
        
        nb_A_lost <- apply(
            freq_mat, 
            MARGIN = 2, 
            FUN = function(x) {
                return(sum(x == 0))
            })
        
        nb_polym <- apply(
            freq_mat, 
            MARGIN = 2, 
            FUN = function(x) { 
                return(sum(x != 0 & x != 1)) 
            })
        
        mean_p <- apply(
            freq_mat, 
            MARGIN = 2, 
            FUN = function(x) { 
                return(mean(x)) 
            })
        
        if (nb_rep > 1) {
            var_p <- apply(
                freq_mat, 
                MARGIN = 2, 
                FUN = function(x) { 
                    return(var(x)) 
                })
        }

        sim_data <- list(freq_mat = freq_mat, 
                         count_mat = count_mat, 
                         nb_copies = nb_copies, 
                         ini_p = ini_p, 
                         nb_gen = nb_gen, 
                         nb_rep = nb_rep, 
                         nb_A_fixed = nb_A_fixed, 
                         nb_A_lost = nb_A_lost, 
                         nb_polym = nb_polym, 
                         mean_p = mean_p, 
                         var_p = var_p, 
                         time_A_fixed = time_A_fixed, 
                         time_A_lost = time_A_lost, 
                         param = param_checking)
    } else {
        sim_data <- list(param = param_checking)
    }

    return(sim_data)
}



#==============================================================================
# get_count_df
#==============================================================================
get_count_df <- function(sim_data)
{
    if (sim_data$param$flag) {
        count_mat <- sim_data$count_mat
        count_df <- as.data.frame(count_mat)
        nb_gen <- sim_data$nb_gen
        nb_copies <- sim_data$nb_copies
        xi <- 0:nb_copies
        colnames(count_df) <- paste("ni_gen_", 0:nb_gen, sep = "")
        count_df <- cbind(xi, count_df)
        count_df <- count_df[nrow(count_df):1, ]
        
        return(count_df)
    } else {
        return(NULL)
    }
}



#==============================================================================
# plot_error
#==============================================================================
plot_error <- function(msg)
{
    plot(c(0, 1), c(0, 1), type = "n", 
         xlab = "", ylab = "", main = "", 
         xaxt = "n", yaxt = "n", bty = "n")
    
    err_msg <- sprintf("ERROR: check parameter values%s", msg)
    text(0.5, 0.5, err_msg, col = "red", adj = 0.5)
}



#==============================================================================
# plot_freq
#==============================================================================
plot_freq <- function(sim_data, fix_flag = FALSE)
{
    if (sim_data$param$flag) {
        freq_mat <- sim_data$freq_mat
        nb_gen <- sim_data$nb_gen
        nb_rep <- sim_data$nb_rep

        if (nb_rep > 1000) {
            msg <- sprintf("\ntoo many replicates to display")
            plot_error(msg)
        } else {
            par_bak <- par(no.readonly = TRUE)

            layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

            par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, 
                        top = 0.5 + 0, right = 0.5 + 5))
            box_y <- c(0, 1)
            box_x <- c(0, nb_gen)
            plot(box_x, box_y, type = "n",
                 xlim = c(0, nb_gen), 
                 main = "", 
                 xlab = "Time (generations)", 
                 ylab = "Frequency P", 
                 bty = "n")
            
            abline(h = 0, col = "gray")
            abline(h = 1, col = "gray")
        
            # color_rep <- hue_pal()(nb_rep)  # 'scales' package
            # color_rep <- rainbow(nb_rep)
            color_rep <- rainbow_hcl(nb_rep)   # 'colorspace' package
            color_2 <- rainbow_hcl(2)   # 'colorspace' package
            color_A_fixed <- color_2[1]
            color_A_lost <- color_2[2]
            
            val_x <- 0:nb_gen
            pos_A_fixed <- freq_mat[, nb_gen + 1] == 1
            pos_A_lost <- freq_mat[, nb_gen + 1] == 0
            pos_polym <- !(pos_A_fixed | pos_A_lost) 
            for (i in 1:nb_rep) {
                val_y <- freq_mat[i, ]
                if (fix_flag) {
                    pos <- which((val_y != 0) & (val_y != 1))
                    pos_0 <- which(val_y == 0)
                    pos_1 <- which(val_y == 1)
                    if (length(pos_0) > 1) {
                        pos_0 <- pos_0[1]
                        pos <- append(pos, pos_0)
                        lines(val_x[pos], val_y[pos], col = color_A_lost)
                        points(val_x[pos_0], val_y[pos_0], pch = 20, 
                               col = color_A_lost)
                    } else {
                        if (length(pos_1) > 1) {
                            pos_1 <- pos_1[1]
                            pos <- append(pos, pos_1)
                            lines(val_x[pos], val_y[pos], col = color_A_fixed)
                            points(val_x[pos_1], val_y[pos_1], pch = 20, 
                                   col = color_A_fixed)
                        } else {
                            lines(val_x, val_y, col = "black")
                        }
                    }
                } else {
                    lines(val_x, val_y, col = color_rep[i])
                }
            }
            val_x <- 0
            val_y <- sim_data$ini_p
            points(val_x, val_y, pch = 15, col = "black", cex = 1)

            if (fix_flag) {

                A_lost_label <- sprintf("A lost\n%d", sum(pos_A_lost))
                polym_label <- sprintf("Polym\n%d", sum(pos_polym))
                A_fixed_label <- sprintf("A fixed\n%d", sum(pos_A_fixed))
                axis(4, at = c(0.0), label = A_lost_label, las = 2, 
                     lwd = 0, col.axis = color_A_lost, 
                     cex.axis = 1)
                axis(4, at = c(0.5), label = polym_label, las = 2, 
                     lwd = 0, col.axis = "black", 
                     cex.axis = 1)
                axis(4, at = c(1.0), label = A_fixed_label, las = 2, 
                     lwd = 0, col.axis = color_A_fixed, 
                     cex.axis = 1)

                if (sum(sim_data$time_A_fixed) > 0) {
                    rug(sim_data$time_A_fixed, side = 3, col = color_A_fixed)
                    rug(mean(sim_data$time_A_fixed), side = 3, col = "black", 
                        lwd = 2)
                }
                if (sum(sim_data$time_A_lost) > 0) {
                    rug(sim_data$time_A_lost, side = 1, col = color_A_lost)
                    rug(mean(sim_data$time_A_lost), side = 1, col = "black", 
                        lwd = 2)
                }
            }

            par(par_bak)
        }
    } else {
        plot_error(sim_data$param$msg)
    }
}




#==============================================================================
# plot_freq_density
#==============================================================================
plot_freq_density <- function(sim_data, count_flag = FALSE) 
{
    if (sim_data$param$flag) {
        count_mat <- sim_data$count_mat
        nb_copies <- sim_data$nb_copies
        ini_p <- sim_data$ini_p
        nb_gen <- sim_data$nb_gen
        nb_rep <- sim_data$nb_rep

        par_bak <- par(no.readonly = TRUE)

        layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, 
                    top = 0.5 + 0, right = 0.5 + 5))
        box_y <- c(0 - 0.5, nb_copies + 0.5)
        box_x <- c(0 - 0.5, nb_gen + 0.5)
        plot(box_x, box_y, type = "n", 
             main = "", 
             xlab = "Time (generations)", 
             ylab = "Allele A count", 
             bty = "n")
        
        if (nb_gen <= 40) {
            axis(1, at = 0:nb_gen, labels = FALSE)
        }
        if (nb_copies <= 40) {
            axis(2, at = 0:nb_copies, labels = FALSE)
        }
        scale_step <- 1 / (nb_rep + 1)
        scale_array <- (1 - (log10(seq(from = scale_step, to = 1, 
                                       by = scale_step)) / log10(scale_step))) 
        color_array <- rgb(1, 165 / 255, 0, alpha = scale_array) 
        for (nb_A in 0:nb_copies) {
            for (k in 1:nb_gen) {
                i <- nb_A + 1
                j <- k + 1
                val <- count_mat[i, j]
                rect_color <- color_array[val + 1]
                x1 <- k - 0.5
                y1 <- nb_A - 0.5
                x2 <- k + 0.5
                y2 <- nb_A + 0.5
                rect(x1, y1, x2, y2, col = rect_color, border = NA)
                if (count_flag) {
                    text(k, nb_A, val)
                }
            }
        }
        k <- 0
        nb_A <- nb_copies * ini_p
        x1 <- k - 0.5
        y1 <- nb_A - 0.5
        x2 <- k + 0.5
        y2 <- nb_A + 0.5
        rect(x1, y1, x2, y2, col = "gray", border = NA)
        val <- nb_rep
        if (count_flag) {
            text(k, nb_A, val)
        }

        par(par_bak)
    } else {
        plot_error(sim_data$param$msg)
    }
}




#==============================================================================
# plot_mean_P
#==============================================================================
plot_mean_P <- function(sim_data, expected_mean_flag = FALSE) 
{
    if (sim_data$param$flag) {
        ini_p <- sim_data$ini_p
        nb_gen <- sim_data$nb_gen
        mean_p <- sim_data$mean_p

        par_bak <- par(no.readonly = TRUE)

        layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, 
                    top = 0.5 + 0, right = 0.5 + 5))

        box_y <- c(0, 1)
        box_x <- c(0, nb_gen)

        plot(box_x, box_y, type = "n", 
             main = "", 
             xlab = "Time (generations)", 
             ylab = "Mean of P", 
             bty = "n")
        
        abline(h = seq(from = 0, to = 1, by = 0.1), col = "gray", lty = 2)

        if (expected_mean_flag) {
            axis(4, at = ini_p, label = expression(p[0]), las = 2, lwd = 0)
            abline(h = ini_p, lty = 2, lwd = 2, col = "red")
        }

        abline(h = 0, col = "gray")
        abline(h = 1, col = "gray")

        val_x <- 0:nb_gen
        val_y <- mean_p
        lines(val_x, val_y)

        if (nb_gen < 50) {
            points(val_x, val_y, pch = 15, col = "gray", cex = 1)
        }
        
        par(par_bak)
    } else {
        plot_error(sim_data$param$msg)
    }
}




#==============================================================================
# plot_variance_P
#==============================================================================
plot_variance_P <- function(sim_data, expected_var_flag = FALSE) 
{
    if (sim_data$param$flag) {
        nb_copies <- sim_data$nb_copies
        ini_p <- sim_data$ini_p
        nb_gen <- sim_data$nb_gen
        nb_rep <- sim_data$nb_rep
        var_p <- sim_data$var_p

        par_bak <- par(no.readonly = TRUE)

        layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, 
                    top = 0.5 + 0, right = 0.5 + 5))

        max_y <- ini_p * (1 - ini_p)
        if (nb_rep > 1) {
            max_y <- max(c(var_p, max_y))
        }
        max_y <- ceiling(max_y * 20) / 20
        box_y <- c(0, max_y)
        box_x <- c(0, nb_gen)

        plot(box_x, box_y, type = "n", 
             main = "", 
             xlab = "Time (generations)", 
             ylab = "Variance of P", 
             bty = "n")
        
        if (expected_var_flag) {
            axis(4, at = ini_p * (1 - ini_p), 
                 label = expression(p[0](1 - p[0])), las = 2, lwd = 0)
            abline(h = ini_p * (1 - ini_p), lty = 2, lwd = 2, col = "gray")
            val_x <- 0:nb_gen
            val_y <- ini_p * (1 - ini_p) * (1 - (1 - 1 / nb_copies)^val_x)
            lines(val_x, val_y, lty = 2, lwd = 2, col = "red")
        }

        abline(h = 0, col = "gray")

        if (nb_rep > 1) {
            val_x <- 0:nb_gen
            val_y <- var_p
            lines(val_x, val_y)

            if (nb_gen < 50) {
                points(val_x, val_y, pch = 15, col = "gray", cex = 1)
            }
        
        } else {
            text(mean(box_x), mean(box_y), "NOT AVAILABLE", col = "red")
        }

        par(par_bak)
    } else {
        plot_error(sim_data$param$msg)
    }
}



#==============================================================================
# plot_fixation_prob
#==============================================================================
plot_fixation_prob <- function(sim_data, expected_prob_flag = FALSE) 
{
    if (sim_data$param$flag) {
        ini_p <- sim_data$ini_p
        nb_gen <- sim_data$nb_gen
        nb_rep <- sim_data$nb_rep
        nb_A_fixed <- sim_data$nb_A_fixed
        nb_A_lost <- sim_data$nb_A_lost
        nb_polym <- sim_data$nb_polym

        par_bak <- par(no.readonly = TRUE)

        color_2 <- rainbow_hcl(2)  # 'colorspace' package
        color_A_fixed <- color_2[1]
        color_A_lost <- color_2[2]
        color_polym <- "black"

        layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, 
                    top = 0.5 + 0, right = 0.5 + 5))

        box_y <- c(0, 1)
        box_x <- c(0, nb_gen)

        plot(box_x, box_y, type = "n", 
             main = "", 
             xlab = "Time (generations)", 
             ylab = "Proportion", 
             bty = "n")
        
        abline(h = seq(from = 0, to = 1, by = 0.1), col = "gray", lty = 2)
        
        if (expected_prob_flag) {
            abline(h = ini_p, lty = 2, lwd = 2, col = color_A_fixed)
            abline(h = 1 - ini_p, lty = 2, lwd = 2, col = color_A_lost)
            if (ini_p == 0.5) {
                axis(4, at = ini_p, label = expression(p[0]), las = 2, 
                     lwd = 0, cex.axis = 1)
            } else {
                axis(4, at = ini_p, label = expression(p[0]), las = 2, lwd = 0, 
                     cex.axis = 1)
                axis(4, at = 1 - ini_p, label = expression(1 - p[0]), las = 2, 
                     lwd = 0, cex.axis = 1)
            }
        }

        abline(h = 0, col = "gray")
        abline(h = 1, col = "gray")

        gen_max <- 50
        
        val_x <- 0:nb_gen
        val_y <- nb_A_fixed / nb_rep
        lines(val_x, val_y, col = color_A_fixed)
        if (nb_gen < gen_max) {
            points(val_x, val_y, pch = 15, col = color_A_fixed, cex = 1)
        }

        val_x <- 0:nb_gen
        val_y <- nb_A_lost / nb_rep
        lines(val_x, val_y, col = color_A_lost)
        if (nb_gen < gen_max) {
            points(val_x, val_y, pch = 15, col = color_A_lost, cex = 1)
        }

        val_x <- 0:nb_gen
        val_y <- nb_polym / nb_rep
        lines(val_x, val_y, col = color_polym)
        if (nb_gen < gen_max) {
            points(val_x, val_y, pch = 15, col = color_polym, cex = 1)
        }

        legend("topright", 
               legend = c("Polym", "A fixed", "A lost"), 
               col = c(color_polym, color_A_fixed, color_A_lost), 
               lty = 1, bg = "white")
        
        par(par_bak)
    } else {
        plot_error(sim_data$param$msg)
    }
}



#==============================================================================
# plot_fixation_time
#==============================================================================
plot_fixation_time <- function(sim_data, expected_time_flag = FALSE)
{
    if (sim_data$param$flag) {
        nb_copies <- sim_data$nb_copies
        ini_p <- sim_data$ini_p
        nb_gen <- sim_data$nb_gen
        nb_rep <- sim_data$nb_rep
        time_A_fixed <- sim_data$time_A_fixed
        time_A_lost <- sim_data$time_A_lost

        mean_time_A_fixed <- mean(time_A_fixed)
        mean_time_A_lost <- mean(time_A_lost)

        # color_2 <- hue_pal()(2)  # 'scales' package
        color_2 <- rainbow_hcl(2)  # 'colorspace' package

        par_bak <- par(no.readonly = TRUE)

        layout(matrix(c(1), nrow = 1, ncol = 1, byrow = TRUE))

        par(mar = c(bottom = 0.5 + 4, left = 0.5 + 4, top = 0.5 + 0, 
                    right = 0.5 + 5))

        expected_time_A_fixed <- -2 * nb_copies * (1 - ini_p) / ini_p * 
            log(1 - ini_p)
        expected_time_A_lost <- -2 * nb_copies * (ini_p) / (1 - ini_p) * 
            log(ini_p)


        A_label <- sprintf("A fixed\nObserved mean = %5.1f", mean_time_A_fixed)
        B_label <- sprintf("A lost\nObserved mean = %5.1f", mean_time_A_lost)
        if (expected_time_flag) {
            A_label <- sprintf("%s\nExpected mean = %5.1f", A_label, 
                               expected_time_A_fixed)
            B_label <- sprintf("%s\nExpected mean = %5.1f", B_label, 
                               expected_time_A_lost)
        }
    
        if (length(time_A_fixed) > 0 | length(time_A_lost) > 0) {
            plot_df <- data.frame(
                Time = c(time_A_fixed, time_A_lost), 
                Allele = as.factor(c(rep(A_label, 
                                         times = length(time_A_fixed)), 
                                     rep(B_label, 
                                         times = length(time_A_lost)))))
            fix_plot <- ggplot(plot_df, 
                               aes(x = Time, colour = Allele, fill = Allele))
            fix_plot <- fix_plot + geom_histogram(
                aes(y = ..density..), 
                alpha = 0.5, 
                position = "identity")
            fix_plot <- fix_plot + xlim(c(0, nb_gen))
            fix_plot <- fix_plot + geom_density(alpha = 0.2)
            fix_plot <- fix_plot + labs(x = "Time (generations)", 
                                        y = "Density")
            fix_plot <- fix_plot + geom_vline(
                xintercept = c(mean_time_A_fixed, mean_time_A_lost), 
                colour = color_2, size = 1)
            if (expected_time_flag) {
                fix_plot <- fix_plot + geom_vline(
                    xintercept = c(expected_time_A_fixed, expected_time_A_lost), 
                    colour = color_2, linetype = "dotted", size = 1)
            }
            fix_plot <- fix_plot + theme(legend.justification = c(1, 1), 
                                         legend.position = c(1, 1))
            
            print(fix_plot)
        } else {
            plot(c(0, 1), c(0, 1), type = "n", xlab = "", ylab = "", main = "", 
                 xaxt = "n", yaxt = "n", bty = "n")
            err_msg <- "NO FIXATION"
            text(0.5, 0.5, err_msg, col = "red", adj = 0.5)
        }

        par(par_bak)
    } else {
        plot_error(sim_data$param$msg)
    }
}


