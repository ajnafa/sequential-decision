#----------------------AR(1) Thompson Sampling-----------------------#
#--------------------------------------------------------------------#
# License: MIT--------------------------------Author: A. Jordan Nafa #

# Set session options
options(
    scipen = 999,
    digits = 4,
    mc.cores = 8
)

# Load libraries
pacman::p_load(
    data.table,
    cmdstanr,
    ggplot2,
    ggdist,
    extraDistr
)

# Load helper functions
source("functions/helpers.R")

# Set the ggplot2 theme
theme_set(plot_theme(
    xaxis_size = 22,
    yaxis_size = 22,
    title_size = 25,
    caption_size = 18,
    axis_text_size = 20,
    strip_size = 20,
    legend_text_size = 20,
    legend.position = "top",
    strip_face = "bold",
    caption.hjust = 0,
    caption.vjust = -1,
    subtitle_size = 20,
    base_size = 20,
    transparent = FALSE,
    plot.margin = margin(2, 8, 4, 5, "mm")
))

#--------------------------------------------------------------------#
#--------------------------Data Simulation----------------------------
#--------------------------------------------------------------------#

# Set the cmdstan path
set_cmdstan_path("E:/Users/Documents/.cmdstan/cmdstan-2.33.1")

# Compile the model
model = cmdstan_model("models/binomial_bandit_ar1.stan")

# Define the dimensions of the experiment
n = 10000              # Average number of new players per day
T = 90                 # Maximum Number of days to run the experiment
het = FALSE            # Whether to include temporal heterogeneity
k = 0.01               # Decay Rate for the number of new players
decay = FALSE

# Define the treatment effect and baseline conversion rate
tau = 0.8              # Treatment Effect
alpha = -4.5           # Baseline Conversion Rate
prob = 0.5             # Allocation Probability at period 1
rho = 0.00             # Correlation between the potential outcomes

# Define the data frame to store the results
data = data.table(
    period = integer(),
    treat = integer(),
    arm_id = integer(),
    purchases = integer(),
    players = integer()
)

# Define the data frame to store the results
track_data = data.table(
    next_period = integer(),
    arm_id = integer(),
    current_period = integer(),
    allocation_probs = numeric()
)

for (t in 1:T) {

    if (decay) {
        # Simulate the number of new players
        n_t = rpois(1, n * exp(-k * t))
    } else {
        n_t = rpois(1, n)
    }

    # Simulate the allocation of new players to each arm
    X = rbinom(n_t, 1, prob)

    # Treatment Assignment Mechanism
    muA = alpha + 0*tau; sdA = 1  # Mean and Std. Dev for the Control
    muB = alpha + 1*tau; sdB = 1  # Mean and Std. Dev for Arm 1

    # True Difference in Means, True Long-Run ATE is 0.01314
    delta = inv_logit(muB) - inv_logit(muA)

    # Simulate the potential outcomes for each arm
    mu = c(muA, muB)

    # Define the covariance matrix
    cov = matrix(
        c(sdA^2, rho*sdA*sdB, 
        rho*sdA*sdB, sdB^2), 
        nrow = 2, ncol = 2)

    # Simulate the Potential Outcomes for each arm
    mu = MASS::mvrnorm(n_t, mu, cov)

    # Simulate the observed outcomes for each arm
    YA = mu[, 1]        # Potential Outcome if A = 1, A' = 0
    YB = mu[, 2]        # Potential Outcome if A = 0, B = 1

    # Realization of the Potential Outcomes on the Latent Scale
    Y_obs = YA * (1 - X) + YB * X
    Y_mis = YA * X + YB * (1 - X)

    # Simulate the observed outcomes at time t
    Y = rbinom(n_t, 1, inv_logit(Y_obs))

    # Store the results in a data frame
    data_t = data.table(
        period = rep(t, length(Y)),
        treat = X,
        arm_id = X + 1,
        purchase = Y,
        players = rep(1, length(Y))
    )[
        , .(purchases = sum(purchase), players = sum(players)), 
        by = .(period, treat, arm_id)
    ][sort(arm_id, decreasing = TRUE)]

    # Append the data to the main data frame
    data = rbind(data, data_t)

    # Need at leat 3 periods to estimate the AR(1) model
    if (t > 3) {

        # Define the data list 
        stan_data <- list(
            N = nrow(data),
            D = max(data$arm_id),
            T = max(data$period),
            M = 3,
            K = data$players,
            Y = data$purchases,
            dd = data$arm_id,
            tt = data$period,
            truth = delta,
            debug = 0
        )

        # Set the initial number of draws to 0
        min_draws <- 0
        
        # Repeat the run if any of the chains stop unexpectedly due 
        # to a weird cmdstanr issue where the chains stop unexpectedly 
        # when running models repeatedly
        while (min_draws < (1000 * 4)) {
            
            # Fit the Stan Model
            fit <- model$sample(
                data = stan_data,
                chains = 4,
                parallel_chains = 4,
                iter_warmup = 2000,
                iter_sampling = 2000,
                refresh = 0
            )
            
            # Update the check
            min_draws <- posterior::ndraws(fit$draws())
        }

        # Extract the posterior draws
        posterior <- fit$draws('prob_best', format = 'draws_df')

        # Extract the posterior draws for prob_best
        prob_best_draws = melt(
            as.data.table(posterior),
            id.vars = c('.iteration', '.draw'),
            measure.vars = patterns('prob_best')
        )[
            , `:=` (
                'arm_id' = stringr::str_extract(variable, ",[0-9]+") |>
                    stringr::str_remove(",") |> as.numeric(),
                'time' = stringr::str_extract(variable, "[0-9]+") |> 
                    as.numeric(),
                'post_prob' = value
            )
        ][
            , .(post_prob = mean(post_prob)), by = .(time, arm_id)
        ]
        
        # Extract the draws for the next period
        prob_best_draws <- prob_best_draws[prob_best_draws$time > t][
            , .(
            allocation_probs = mean(post_prob),
            next_period = min(time),
            current_period = t
            ), by = .(arm_id)
        ]

        track_data <- rbind(track_data, prob_best_draws)

        # Update the allocation probability
        prob = prob_best_draws[arm_id == 2, allocation_probs]

        # Sleep for 2 seconds to prevent chains from failing
        Sys.sleep(2)
    }
}

# Save the final model's output files
fit$save_output_files(
    dir = "output/fit/", 
    basename = "binomial_bandit_ar1_fit_t90"
    )

# Save the final model object
fit$save_object(file = "output/fit/binomial_bandit_ar1_fit_t90.Rds")

# Save the model data
arrow::write_parquet(
    data, 
    "output/data/binomial_bandit_ar1_data.parquet"
)

# Save the tracked estimates
arrow::write_parquet(
    track_data, 
    "output/data/binomial_bandit_ar1_track_data.parquet"
)

#--------------------------------------------------------------------#
#-----------------------Graphing the Results--------------------------
#--------------------------------------------------------------------#

# Read in the model object
fit <- readr::read_rds("output/fit/binomial_bandit_ar1_fit_t90.Rds")

# Read in the tracked allocation probabilities
prob_data <- arrow::read_parquet(
    "output/data/binomial_bandit_ar1_track_data.parquet"
)

# Extract the results of the final model
posterior <- fit$draws(c('prob_best', 'theta_diff', 'bias'), format = 'draws_df')

# Extract the posterior draws for theta
ate_draws = melt(
    as.data.table(posterior),
    id.vars = c('.iteration', '.draw'),
    measure.vars = patterns('theta_diff'),
    value.name = 'theta_diff'
)[
    , `:=` (
        'time' = stringr::str_extract(variable, "[0-9]+") |> 
            as.numeric()
    )
]

# Plot the results
ate_plot <- ggplot(ate_draws, aes(x = time, y = theta_diff, group = time)) +
    stat_gradientinterval(
        aes(
            slab_fill = after_stat(y > 0),
            point_fill = after_stat(y > 0),
            slab_alpha = after_stat(pdf / ave(pdf, group, datatype, FUN = "max"))
        ),
        .width = c(0.50, 0.80),
        shape = 22,
        fill_type = "gradient",
        point_interval = mean_qi,
        scale = 1
    ) +
    scale_fill_manual(
        "Direction of Effect",
        values = c("#F8766D", "#00BFC4"),
        aesthetics = "slab_fill",
        guide = guide_legend(
            override.aes = list(
                shape = 22,
                size = 4,
                fill = c("#F8766D", "#00BFC4"),
                alpha = 1
            )
        ),
        labels = c("Negative", "Positive")
    ) +
    scale_fill_manual(
        values = c("#F8766D", "#00BFC4"),
        aesthetics = "point_fill",
        guide = "none"
    ) +
    scale_slab_alpha_continuous(guide = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(
        breaks = scales::pretty_breaks(n = 10),
        limits = c(-0.2, 0.2)
    ) +
    labs(
        x = "Time",
        y = "Daily Average Treatment Effect",
        title = "AR(1) Thompson Sampling for Bernoulli Bandit Problem",
        subtitle = "Posterior Predictive Distribution of the Estimated Treatment Effect Over Time"
    )

# Save the plot
ggsave(
    filename = "figures/ate_plot_t90.png",
    plot = ate_plot,
    width = 15.5,
    height = 10,
    units = "in",
    dpi = "retina",
    bg = 'white',
    type = "cairo"
)

# Plot the posterior probabiltiies
prob_plot <- ggplot(
    prob_data, 
    aes(
        x = next_period, 
        y = allocation_probs,
        group = arm_id,
        fill = factor(arm_id)
    )) +
    geom_point(size = 2, shape = 21) +
    scale_fill_manual(
        "Condition",
        values = c("#F8766D", "#00BFC4"),
        labels = c("Control", "Treatment"),
        aesthetics = "fill"
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    labs(
        x = "Time",
        y = "Allocation Probability for Time t + 1",
        title = "AR(1) Thompson Sampling for Bernoulli Bandit Problem",
        subtitle = "Evolution of the Allocation Probabilities Over the Course of the Experiment"
    )

# Save the plot
ggsave(
    filename = "figures/allocation_probs_t90.png",
    plot = prob_plot,
    width = 15.5,
    height = 10,
    units = "in",
    dpi = "retina",
    bg = 'white',
    type = "cairo"
)


#--------------------------------------------------------------------#
#-------------------------Basic Diagnostics---------------------------
#--------------------------------------------------------------------#

# Basic Diagnostics
fit$diagnostic_summary()

# Check the Rhat values
summ <- posterior::summarise_draws(fit$draws())

# Check the Rhat values
bayesplot::mcmc_rhat_hist(summ$rhat)
