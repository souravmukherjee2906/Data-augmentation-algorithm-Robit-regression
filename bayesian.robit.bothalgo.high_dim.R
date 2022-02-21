## Analysis of convergence of DA and sandwich DA algorithms in low-dimensional data.


#------------------------------------------------------------------------------------
## Clear the environment ##
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(TruncatedNormal)) install.packages("TruncatedNormal")
if(!require(truncnorm)) install.packages("truncnorm")
if(!require(crch)) install.packages("crch")
if(!require(coda)) install.packages("coda")
if(!require(tidyverse)) install.packages("tidyverse")
if(!require(mvtnorm)) install.packages("mvtnorm")


#------------------------------------------------------------------------------------
## Load the packages into the R environment ##
#------------------------------------------------------------------------------------
library(TruncatedNormal)
library(truncnorm)
library(crch)   # For generating samples from truncated student-t distribution.
library(coda)
library(tidyverse)
library(mvtnorm) # For generating from the multivariate normal distribution.


#------------------------------------------------------------------------------------
## Source all the functions from "Functions_robit.R" ##
#------------------------------------------------------------------------------------
source("Functions_robit.R")


#------------------------------------------------------------------------------------
## Dataset (High dimensional case) ##
#------------------------------------------------------------------------------------
data(prostate, package = "spls")
X_high.dim <- cbind(1, prostate$x[, 1:150])
y_high.dim <- prostate$y


#------------------------------------------------------------------------------------

# User provided values
nsim <- 3e5
burn.in <- 2e5
after.burnin.vec <- 1:(nsim +1 - burn.in)
nu.small <- 3 
nu.large <- 1000

#------------------------------------------------------------------------------------

# Call the appropriate Markov chain algorithm.
run.chain.high.dim <- function(X, y, nu = 3, beta.start = rep(0, ncol(X)),
                               Sigma_a = diag(1, nrow = ncol(X)),
                               nsim = 1000, burn.in = 100,
                               model = "probit", method = "DA", ...){
  
  if ((model == "probit") && (method == "DA")){
    final.fn <- run_probit_beta_MCMC
  }else if ((model == "probit") && (method == "sandwich")){
    final.fn <- run_probit_sandwich_beta_MCMC
  }else if ((model == "robit") && (method == "DA")){
    final.fn <- run_robit_beta_MCMC
  }else{
    final.fn <- run_robit_sandwich_beta_MCMC
  }
  out <- final.fn(X = X, y = y, nu = nu, beta.start = beta.start,
                  Sigma_a = Sigma_a, nsim = nsim,
                  burn.in = burn.in, model = model, method = method, ...)
  return(out)
}


all.params.high.dim <- expand_grid(algo = c("probit", "robit with small \u03bd",
                                            "robit with large \u03bd"),
                                   method = c("DA", "sandwich") )


all.runs.high.dim <- all.params.high.dim %>%
  mutate(all.chains.high.dim = pmap(
    list(algo, method),
    function(this.algo, this.method){
      msg <- glue::glue("\nNow running, algo = {this.algo}, method = {this.method}")
      message(msg)
      this.model <- ifelse(grepl("probit", this.algo), "probit", "robit")
      this.nu <- ifelse(grepl("large", this.algo), nu.large, nu.small)
      
      run.chain.high.dim(X = X_high.dim, y = y_high.dim, nu = this.nu,
                         nsim = nsim, burn.in = burn.in,
                         model = this.model, method = this.method)
    }
  )
  )

#------------------------------------------------------------------------------------

# Function to calculate log-likelihood and log-posterior for both probit and robit
# regression models.
log.lik.post.both <- function(X, y, nu, beta, beta_a = rep(0, ncol(X)),
                              Sigma_a = diag(1, ncol(X)), model, ...){
  if (model == "probit"){
    final.fn <- log.lik.post.probit
  }else{
    final.fn <- log.lik.post.robit
  }
  
  out <- final.fn(X = X, y = y, nu = nu, beta = beta, beta_a = beta_a,
                  Sigma_a = Sigma_a)
  return(out)
}

#------------------------------------------------------------------------------------

## Getting the log-likelihood chains and log-posterior chains calculated from
## from the generated Markov chains of beta values for all possible model and
## markov chain types.

daf <- all.runs.high.dim %>% 
  mutate(
    log.lik.post.chain = pmap(
      list(all.chains.high.dim, algo),
      function(beta.mat, this.algo){
        this.model <- ifelse(grepl("probit", this.algo), "probit", "robit")
        this.nu <- ifelse(grepl("large", this.algo), nu.large, nu.small)
        
        log.lik.chain <- apply(beta.mat, 1,
                               FUN = function(beta){
                                 log.lik.post.both(X = X_high.dim,
                                                   y = y_high.dim,
                                                   nu = this.nu,
                                                   beta = beta,
                                                   model = this.model
                                 )$logLikelihood
                               }
        )
        log.post.chain <- apply(beta.mat, 1,
                                FUN = function(beta){
                                  log.lik.post.both(X = X_high.dim,
                                                    y = y_high.dim,
                                                    nu = this.nu,
                                                    beta = beta,
                                                    model = this.model
                                  )$logPosterior
                                }
        )
        output <- list(log.lik.chain = log.lik.chain,
                       log.post.chain = log.post.chain)
        return(output)
      }
    )
  )

# save the R environment
save.image(file = "probit_robit_all_cases_nu=3_and_nu=1000_high_dim.RData")

#------------------------------------------------------------------------------------
## Autocorrelation plots of log-Likelihood and log-Posterior. ##
#------------------------------------------------------------------------------------

daf1.acf <- daf %>% 
  mutate(
    loglik.acf = lapply(log.lik.post.chain,
                        function(which.list){
                          acf(which.list$log.lik.chain, lag.max = 50, plot = F)$acf[ ,1,1]
                        }
    ),
    logpost.acf = lapply(log.lik.post.chain,
                         function(which.list){
                           acf(which.list$log.post.chain, lag.max = 50, plot = F)$acf[ ,1,1]
                         }
    )
  ) %>% 
  mutate(
    model = algo
  )



daf1.acf.formatted <- daf1.acf %>% 
  select(algo, model, method, loglik.acf, logpost.acf) %>%
  pivot_longer(c(loglik.acf, logpost.acf), names_to = "lik_Or_post", values_to = "ACF") %>% 
  mutate(index = lapply(1:n(), function(x){0:50} )) %>% 
  unnest(c(ACF, index)) %>% 
  mutate(
    lik_Or_post_acf = ifelse(
      grepl(pattern = "loglik", lik_Or_post),
      "Autocorrelation plot for log-likelihood",
      "Autocorrelation plot for log-posterior"
    )
  )


daf1.acfplot <- daf1.acf.formatted %>% 
  ggplot(aes(x = index, y = ACF, color = model, linetype = method)) +
  geom_line() +
  facet_wrap(~ lik_Or_post_acf, scales = "free") +
  labs(x = "Lag", y = "Autocorrelation", color = "Model:",
       linetype = "Markov Chain Type:") + 
  scale_color_manual(values = c("chocolate", "grey", "blue")) + 
  theme_bw() + 
  theme(legend.position = "top")


daf1.acfplot

# Save the plot
ggsave(filename = "lik_Or_post.acfplot.pdf", plot = daf1.acfplot, height = 4,
       width = 8, device = cairo_pdf)


#------------------------------------------------------------------------------------
## Running mean plots of log-Likelihood and log-Posterior. ##
#------------------------------------------------------------------------------------

daf2.cummean <- daf %>% 
  mutate(
    loglik.cummean = lapply(log.lik.post.chain,
                            function(which.list){
                              cummean(which.list$log.lik.chain)
                            }
    ),
    logpost.cummean = lapply(log.lik.post.chain,
                             function(which.list){
                               cummean(which.list$log.post.chain)
                             }
    )
  ) %>% 
  mutate(
    model = algo
  )


daf2.cummean.formatted <- daf2.cummean %>% 
  select(algo, model, method, loglik.cummean, logpost.cummean) %>%
  pivot_longer(c(loglik.cummean, logpost.cummean),
               names_to = "lik_Or_post", values_to = "cum_mean") %>% 
  mutate(index = lapply(1:n(), function(x){0:(nsim - burn.in)} )) %>% 
  unnest(c(cum_mean, index)) %>% 
  mutate(
    lik_Or_post_mean = ifelse(
      grepl(pattern = "loglik", lik_Or_post),
      "Running mean plot for log-likelihood",
      "Running mean plot for log-posterior"
    )
  )


daf2.running.mean.plot <- daf2.cummean.formatted %>% 
  ggplot(aes(x = index, y = cum_mean, color = model, linetype = method)) +
  scale_x_continuous(limits = c(0, (nsim - burn.in)) ) +
  geom_line() +
  facet_wrap(~ lik_Or_post_mean, scales = "free") +
  labs(x = "Iteration", y = "Running mean", color = "Model:",
       linetype = "Markov Chain Type:") +
  scale_color_manual(values = c("chocolate", "grey", "blue")) + 
  theme_bw() + 
  theme(legend.position = "top")


daf2.running.mean.plot

# Save the plot
ggsave(filename = "lik_Or_post.meanplot.png", plot = daf2.running.mean.plot,
       height = 4, width = 8)

