/* 
    Model: Binomial Bandit for Sequential Bayesian A/B Testing 
    with an AR(1) Formulation for the Latent Mean
    Author: A. Jordan Nafa
    License: MIT
    Date: 2023-08-07
*/


data {
    // Data Dimensions
    int<lower=1> N;                   // Number of observations
    int<lower=1> D;                   // Number of Policies
    int<lower=1> T;                   // Number of Observed Time Periods
    int<lower=1> M;                   // Number of Periods to Forecast

    // Input Data
    array[N] int<lower=1> K;          // Number of Trials
    array[N] int<lower=0> Y;          // Number of Successes
    array[N] int<lower=1,upper=D> dd; // Policy to Observation Mapping
    array[N] int<lower=1,upper=T> tt; // Time to Observation Mapping
    real<lower=-1,upper=1> truth;     // True Long-Run Treatment Effect

    // Debugging Flag
    int<lower=0> debug;               // Prior Predictive Checks
}

transformed data {
    int<lower=T + 1> J;               // Number of Time Periods 
    J = T + M;
}

parameters {
   matrix[D, T] mu;
   vector<lower=0>[D] sigma;
   // @TODO: Might be worth extending the model to 
   // accomodate non-linearities via a gaussian process
}

transformed parameters {
   // Probability of Success
   matrix[D, T] theta;
   for (d in 1:D) {
        theta[d, 1:T] = inv_logit(mu[d, 1:T]);
   } 
   // @TODO: Investigate using a non-linear state-space approach here
}

model {
    // Priors on the Parameters
    target += normal_lpdf(sigma | 0, 1) - 2 * normal_lccdf(0 | 0, 1);
    for (n in 1:N) {
        if (tt[n] == 1) {
            target += normal_lpdf(mu[dd[n], tt[n]] | 0, sigma[dd[n]]);
        } else {
            // @TODO: Consider reparameterizing this
            target += normal_lpdf(mu[dd[n], tt[n]] | mu[dd[n], tt[n] - 1], sigma[dd[n]]);
        }
    }

    // Binomial Likelihood
    if (!debug) {
        for (n in 1:N) {
            target += binomial_lpmf(Y[n] | K[n], theta[dd[n], tt[n]]);
        }
    }
}

generated quantities {
    // Posterior Predictive Check
    array[N] int y_rep;
    for (n in 1:N) {
        y_rep[n] = binomial_rng(K[n], theta[dd[n], tt[n]]);
    }

    // Posterior Predictive Distribution for the Full Series
    array[D, J] real mu_pred;             // Posterior Predictive Distribution for mu
    array[D*J] real mu_pred_vec;
    array[D, J] real theta_pred;          // Posterior Predictive Distribution for Theta
    array[D*J] real theta_pred_vec;

    for (d in 1:D) {
        // @TODO: Maybe this could be recast as MVN
        mu_pred[d, 1:T] = normal_rng(mu[d, 1:T], sigma[d]);
        for (t in (T+1):J) {
            mu_pred[d, t] = normal_rng(mu_pred[d, t-1], sigma[d]);
        }
        // @TODO: Consider making this a function to make this less ugly
        mu_pred_vec[((d-1)*J + 1):(d*J)] = mu_pred[d, 1:J];
        theta_pred[d, 1:J] = inv_logit(mu_pred[d, 1:J]);
        theta_pred_vec[((d-1)*J + 1):(d*J)] = theta_pred[d, 1:J];
    }

    // Difference Between the Series (Daily Average Treatment Effect)
    vector[J] theta_diff;
    vector[J] bias;
    for (t in 1:J) {
        theta_diff[t] = theta_pred[2, t] - theta_pred[1, t];
        bias[t] = theta_diff[t] - truth;
    }

    // Probability of Best at Each Period
    array[J] simplex[D] prob_best;
    {
        for (t in 1:J) {
            real best_arm = max(theta_pred[1:D, t]);
            for (d in 1:D) {
                prob_best[t, d] = (theta_pred[d, t] >= best_arm);
            }
            prob_best[t, 1:D] = prob_best[t, 1:D] / sum(prob_best[t, 1:D]);
        }
    }
}

