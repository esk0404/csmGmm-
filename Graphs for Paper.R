library(ggplot2)

#------------------------------Run first before simulating-------------------------------------

# Run symm_fit_ind_fixed() and calc_dens_ind_hl() in symm_fit_ind_fixed.R

pi <- 0.5 #mixture weight for simulating +/- mu under the alt

#Function to generate alt and null stats
generate_alt_stats <- function(n, mu, pi, hl) {
  K <- length(mu)
  Z <- matrix(NA, nrow = n, ncol = K)
  
  for (j in 1:K) {
    if (hl[j] == 0) {
      Z[, j] <- rnorm(n, 0, 1)
    } else {
      sign_vec <- sample(c(-1, 1), n, replace = TRUE, prob = c(1 - pi, pi))
      Z[, j] <- rnorm(n, mean = sign_vec * mu[j], sd = 1)
    }
  }
  return(Z)
}
#-----------------------------------Run simulations-----------------------------------------------------
fixedPi_3 <- c(0.9, 0.04, 0.03, 0.03)
fixedPi_40 <- c(0.4, 0.1, 0.1, 0.4)
fixedPi_90 <- c(0.03, 0.03, 0.04, 0.9)

# 3% Prior
set.seed(123)
run_simulation_3 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(30, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(900, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 900), rep(0, 40), rep(0, 30), rep(1, 30))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_3)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 30  # 30 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

set.seed(123)
run_simulation_40_pi3 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(100, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(100, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(400, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(400, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 400), rep(0, 100), rep(0, 100), rep(1, 400))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_3)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 400  # 400 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

set.seed(123)
run_simulation_90_pi3 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(900, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 30), rep(0, 30), rep(0, 40), rep(1, 900))  # 0 = null, 1 = alt, -1 = partial alt (ignore)
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_3)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 900  # 900 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

# 40% Prior
set.seed(123)
run_simulation_40 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(100, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(100, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(400, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(400, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 400), rep(0, 100), rep(0, 100), rep(1, 400))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_40)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 400  # 400 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

set.seed(123)
run_simulation_3_pi40 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(30, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(900, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 900), rep(0, 40), rep(0, 30), rep(1, 30))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_40)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 30  # 30 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}


set.seed(123)
run_simulation_90_pi40 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(900, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 30), rep(0, 30), rep(0, 40), rep(1, 900))  # 0 = null, 1 = alt, -1 = partial alt (ignore)
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_40)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 900  # 900 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

# 90% Prior
set.seed(123)
run_simulation_90 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(900, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 30), rep(0, 30), rep(0, 40), rep(1, 900))  # 0 = null, 1 = alt, -1 = partial alt (ignore)
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_90)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 900  # 900 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

set.seed(123)
run_simulation_3_pi90 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(40, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(30, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(30, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(900, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 900), rep(0, 40), rep(0, 30), rep(1, 30))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_90)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 30  # 30 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}

set.seed(123)
run_simulation_40_pi90 <- function(fixedMu, n_sim = 1000) {
  fdp_vec <- numeric(n_sim)
  power_vec <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    alt_1 <- generate_alt_stats(100, fixedMu, pi, hl = c(1, 0))
    alt_2 <- generate_alt_stats(100, fixedMu, pi, hl = c(0, 1))
    alt_3 <- generate_alt_stats(400, fixedMu, pi, hl = c(1, 1))
    null_stats <- generate_alt_stats(400, fixedMu, pi, hl = c(0, 0))
    testStats_new <- rbind(null_stats, alt_1, alt_2, alt_3)
    true_labels <- c(rep(0, 400), rep(0, 100), rep(0, 100), rep(1, 400))  # 0 = null, 1 = alt
    
    fit_result_new <- symm_fit_ind_fixed(testStats_new, fixedMu, fixedPi_90)
    lfdr <- fit_result_new$lfdrResults
    
    sorted_lfdr <- sort(lfdr)
    cum_avg <- cumsum(sorted_lfdr) / seq_along(sorted_lfdr)
    r <- max(which(cum_avg <= 0.1))
    
    if (r > 0) {
      rejected_indices <- order(lfdr)[1:r]
      num_fp <- sum(true_labels[rejected_indices] == 0) # false positives
      num_tp <- sum(true_labels[rejected_indices] == 1) 
      fdp_vec[i] <- num_fp / r
      power_vec[i] <- num_tp / 400  # 400 alternatives
    } else {
      fdp_vec[i] <- 0
      power_vec[i] <- 0
    }
  }
  
  list(mean_fdp = mean(fdp_vec),
       power = mean(power_vec))
}


#-------------------------------90% Prior FDP Curve----------------------------------------------------------------------
mu_values <- seq(0, 8, by = 0.25)
# 90% alt
fixedPi_90 <- c(0.03, 0.03, 0.04, 0.9)
fdp_curve_90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90(fixedMu)$mean_fdp
})

# 3% alt, 90% alt prior
fdp_curve_3_pi90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3_pi90(fixedMu)$mean_fdp
})

# 40% alt, 90% alt prior
fdp_curve_40_pi90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40_pi90(fixedMu)$mean_fdp
})


# Put data into a long format first
df <- data.frame(
  mu = rep(mu_values, 3),
  fdp = c(fdp_curve_3_pi90, fdp_curve_40_pi90, fdp_curve_90),
  group = rep(c("3% Alternative", "40% Alternative", "True Alternative"),
              each = length(mu_values))
)

ggplot(df, aes(x = mu, y = fdp, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "FDP") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotdash", "dashed", "solid"))

#----------------------------90% Prior Power Curve-------------------------------------------
#Plot all power curves in one graph
mu_values <- seq(0, 8, by = 0.25)

# 90% Alt, 90% Prior
fixedPi_90 <- c(0.03, 0.03, 0.04, 0.9)
power_curve_90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90(fixedMu)$power
})

# 3% Alt, 90% Alt Prior
power_curve_3_pi90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3_pi90(fixedMu)$power
})

# 40% alt, 90% Alt Prior
fixedPi_40_pi90 <- c(0.4, 0.1, 0.1, 0.4)
power_curve_40_pi90 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40_pi90(fixedMu)$power
})

# Put data into long format
df_power_90 <- data.frame(
  mu = rep(mu_values, 3),
  power = c(power_curve_3_pi90, power_curve_40_pi90, power_curve_90),
  group = rep(c("3% Alternative", "40% Alternative", "True Alternative"),
              each = length(mu_values))
)

# Plot with JSS style
ggplot(df_power_90, aes(x = mu, y = power, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "Power") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotdash", "dashed", "solid"))


#------------------------------40% Prior FDP Curve--------------------------------------------------
# 40% alt
fixedPi_40 <- c(0.4, 0.1, 0.1, 0.4)
fdp_curve_40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40(fixedMu)$mean_fdp
})

# 3% alt, 40% alt prior
fdp_curve_3_pi40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3_pi40(fixedMu)$mean_fdp
})

# 90% alt, 40% alt prior
fdp_curve_90_pi40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90_pi40(fixedMu)$mean_fdp
})
# Put data into long format
df2 <- data.frame(
  mu = rep(mu_values, 3),
  fdp = c(fdp_curve_3_pi40, fdp_curve_40, fdp_curve_90_pi40),
  group = rep(c("3% Alternative", "True Alternative", "90% Alternative"),
              each = length(mu_values))
)

# Plot with JSS style
ggplot(df2, aes(x = mu, y = fdp, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "FDP") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dashed", "dotdash", "solid"))


#--------------------------40% Prior Power Curve---------------------------------------------------
#Plot all power curves in one graph
mu_values <- seq(0, 8, by = 0.25)
fixedPi_40 <- c(0.4, 0.1, 0.1, 0.4)

# 40% Alt, 40% Prior
power_curve_40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40(fixedMu)$power
})

# 3% Alt, 40% Alt Prior
power_curve_3_pi40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3_pi40(fixedMu)$power
})

# 90% alt, 40% Alt Prior
power_curve_90_pi40 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90_pi40(fixedMu)$power
})

# Put data into long format
df_power_40 <- data.frame(
  mu = rep(mu_values, 3),
  power = c(power_curve_3_pi40, power_curve_40, power_curve_90_pi40),
  group = rep(c("3% Alternative", "True Alternative", "90% Alternative"),
              each = length(mu_values))
)

# Plot with JSS style
ggplot(df_power_40, aes(x = mu, y = power, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "Power") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dashed", "dotdash", "solid"))


#--------------------------3% Prior FDP Curve-------------------------------------------------------
# 3% alt
fixedPi_3 <- c(0.9, 0.04, 0.03, 0.03)
fdp_curve_3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3(fixedMu)$mean_fdp
})

# 90% alt, 3% alt prior
fdp_curve_90_pi3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90_pi3(fixedMu)$mean_fdp
})

# 40% alt, 3% alt prior
fdp_curve_40_pi3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40_pi3(fixedMu)$mean_fdp
})

# Put data into long format
df3 <- data.frame(
  mu = rep(mu_values, 3),
  fdp = c(fdp_curve_3, fdp_curve_40_pi3, fdp_curve_90_pi3),
  group = rep(c("True Alternative", "40% Alternative", "90% Alternative"),
              each = length(mu_values))
)

# Plot with JSS style
ggplot(df3, aes(x = mu, y = fdp, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "FDP") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotdash", "dashed", "solid"))

#------------------------3% Prior Power Curve-------------------------------------------------
# 3% Alt, 3% Prior
fixedPi_3 <- c(0.9, 0.04, 0.03, 0.03)
power_curve_3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_3(fixedMu)$power
})

# 90% Alt, 3% Alt Prior
power_curve_90_pi3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_90_pi3(fixedMu)$power
})

# 40% alt, 3% Alt Prior
fixedPi_40_pi3 <- c(0.4, 0.1, 0.1, 0.4)
power_curve_40_pi3 <- sapply(mu_values, function(mu_val) {
  fixedMu <- rep(mu_val, 2)
  run_simulation_40_pi3(fixedMu)$power
})
# Put data into long format
df_power <- data.frame(
  mu = rep(mu_values, 3),
  power = c(power_curve_3, power_curve_40_pi3, power_curve_90_pi3),
  group = rep(c("True Alternative", "40% Alternative", "90% Alternative"),
              each = length(mu_values))
)

# Plot with JSS style
ggplot(df_power, aes(x = mu, y = power, linetype = group, color = group)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0.8, linetype = "dashed") +
  labs(x = "Effect Size (mu)", y = "Power") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title = element_blank()
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_linetype_manual(values = c("dotdash", "dashed", "solid"))


