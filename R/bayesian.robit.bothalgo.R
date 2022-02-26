## Analysis of convergence of DA and sandwich DA algorithms in low-dimensional
## lupus dataset.


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
## Dataset (Low dimensional case) ##
#------------------------------------------------------------------------------------
data("lupus", package = "TruncatedNormal")
x1 <- lupus[,"x1"]
x2 <- lupus[,"x2"]
const <- lupus[,"const"]
X <- lupus[,-1]
y <- lupus[,"response"]

#------------------------------------------------------------------------------------

# Prior choice
I <- diag(1, ncol(X), ncol(X))
g.small <- 3.49  # Choosing a small value of g in the Zellner's g prior
                 # (see Chakraborty and Khare (2017)). g = 3.49 ensures
                 # the trace-class property of the probit DA algorithm
                 # discussed in Chakraborty and Khare (2017).

g.large <- 1000  # Choosing a large value of g in the Zellner's g prior
                 # (see Chakraborty and Khare (2017)). This high value of
                 # g induces a diffuse prior on the coefficient vector \beta.


# User provided values
nsim <- 3e6
burn.in <- 2e6
after.burnin.vec <- 1:(nsim +1 - burn.in)
beta.start <- c(-1.778,4.374,2.428)
nu.small <- 3 
nu.large <- 1000

#------------------------------------------------------------------------------------

# Call the appropriate Markov chain algorithm.
run.chain <- function(X, y, nu = 3, beta.start = rnorm(ncol(X)),
                      nsim = 1000, burn.in = 100, g,
                      model = "probit", method = "DA", ...){
  Sigma_a <-  (t(X)%*%X)/g
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
                  burn.in = burn.in, g = g, model = model, method = method, ...)
  return(out)
}



all.params <- expand_grid(algo = c("probit", "robit with small \u03bd",
                                   "robit with large \u03bd" ),
                          method = c("DA", "sandwich"),
                          g = c(g.small, g.large) )



all.runs <- all.params %>%
  mutate(all.chains = pmap(
    list(algo, method, g),
    function(this.algo, this.method, this.g){
      msg <- glue::glue("\nNow running, algo = {this.algo}, method = {this.method}, g = {this.g}")
      message(msg)
      this.model <- ifelse(grepl("probit", this.algo), "probit", "robit")
      this.nu <- ifelse(grepl("large", this.algo), nu.large, nu.small)
      
      run.chain(X = X, y = y, nu = this.nu, beta.start = beta.start,
                nsim = nsim, burn.in = burn.in, g = this.g,
                model = this.model, method = this.method)
    }
  )
  )



daf <- all.runs

# save the R environment
save.image(file = "probit_robit_all_cases_nu=3_and_nu=1000.RData")


#------------------------------------------------------------------------------------
## Autocorrelation plots ##
#------------------------------------------------------------------------------------
daf1 <- daf %>% 
  mutate(beta.acf = lapply(all.chains, function(x){acf(x, lag.max = 50, plot = F) })) %>% 
  mutate(
    beta_1.acf = lapply(beta.acf, function(x){x$acf[ ,2,2] }),
    beta_2.acf = lapply(beta.acf, function(x){x$acf[ ,3,3] })
  ) %>% 
  mutate(
    model = algo
  )


daf1.beta <- daf1 %>% 
  select(algo, model, method, beta_1.acf, beta_2.acf, g) %>%
  pivot_longer(c(beta_1.acf, beta_2.acf), names_to = "beta_type", values_to = "ACF") %>% 
  mutate(index = lapply(1:n(), function(x){0:50} )) %>% 
  unnest(c(ACF, index)) %>% 
  mutate(
    beta_type_formatted = ifelse(
      grepl(pattern = "beta_1", beta_type),
      "'Autocorrelation plot for '*beta[1]",
      "'Autocorrelation plot for '*beta[2]"
      ),
    g.formatted = glue::glue("'g = {g}'")
    )


daf1.acfplot <- daf1.beta %>% 
  ggplot(aes(x = index, y = ACF, color = model, linetype = method)) +
  geom_line() +
  facet_grid(g.formatted ~ beta_type_formatted, labeller = label_parsed) +
  labs(x = "Lag", y = "Autocorrelation", color = "Model:",
       linetype = "Markov Chain Type:") + 
  scale_color_manual(values = c("chocolate", "grey", "blue")) + 
  theme_bw() + 
  theme(legend.position = "top")


daf1.acfplot

# Save the plot
ggsave(filename = "beta.acfplot.pdf", plot = daf1.acfplot, height = 6, width = 8,
       device = cairo_pdf)


#------------------------------------------------------------------------------------
## Running mean plots ##
#------------------------------------------------------------------------------------
daf2 <- daf %>% 
  mutate(
    beta.running.sum = lapply(all.chains, function(beta_mat){apply(beta_mat, 2, cumsum)})
  ) %>% 
  mutate(
    beta1.running.mean = lapply(beta.running.sum,
                              function(beta_mat_sum){beta_mat_sum[,2]/after.burnin.vec}),
    beta2.running.mean = lapply(beta.running.sum,
                             function(beta_mat_sum){beta_mat_sum[,3]/after.burnin.vec})
  ) %>% 
  mutate(
    model = algo
  )


daf2.beta <- daf2 %>%
  select(algo, model, method, beta1.running.mean, beta2.running.mean, g) %>% 
  pivot_longer(c(beta1.running.mean, beta2.running.mean),
               names_to = "beta_type", values_to = "cum_mean" ) %>% 
  mutate(index = lapply(1:n(), function(x){0:(nsim - burn.in)} )) %>% 
  unnest(c(cum_mean, index)) %>% 
  mutate(
    beta_type_formatted = ifelse(
      grepl(pattern = "beta1", beta_type),
      "'Running mean plot for '*beta[1]",
      "'Running mean plot for '*beta[2]"
    ),
    g.formatted = glue::glue("'g = {g}'")
  )


daf2.running.mean.plot <- daf2.beta %>% 
  ggplot(aes(x = index, y = cum_mean, color = model, linetype = method)) +
  scale_x_continuous(limits = c(0, (nsim - burn.in)) ) +
  geom_line(size = 0.5) +
  facet_grid(g.formatted ~ beta_type_formatted, labeller = label_parsed) +
  labs(x = "Iteration", y = "Running mean", color = "Model:",
       linetype = "Markov Chain Type:") +
  scale_color_manual(values = c("chocolate", "grey", "blue")) +
  theme_bw() +
  theme(legend.position = "top")


daf2.running.mean.plot

# Save the plot
ggsave(filename = "beta.meanplot.png", plot = daf2.running.mean.plot,
       height = 6, width = 8)

