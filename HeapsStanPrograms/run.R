## Functions for reading in transformed, differenced data
my_path <- "./HeapsStanPrograms" # path to directory containing R code and Stan programs,
# starting in current working directory
source(paste(my_path, "/read.R", sep = ""))

## Construct data set
my_m <- 3 # m = 3 example
my_omit <- c(1:2, 197:200) # constructs quarterly time-series described in paper
my_Nahead <- 40 # hold-back last 40 time-points, as described in paper
yraw_trim <- process_data(my_m, my_omit, my_Nahead, paste(my_path, "/data", sep = "")) # construct data
yraw_trim_hb <- yraw_trim$y

## Run HMC using rstan
# All the Stan programs except statrmlprior.stan can be run using rstan;
# the latter requires version 2.26 of Stan which is, at the time of writing, available
# in cmdstanr (instructions below) but not rstan
library(rstan) # load the rstan package
options(mc.cores = parallel::detectCores()) # allow use of multiple cores
rstan_options(auto_write = TRUE)

# As an example, run the code for the stationary, exchangeable prior
my_p <- 4 # fit model of order 4
my_data <- list(
    m = my_m, p = my_p, N = nrow(yraw_trim_hb), y = yraw_trim_hb, es = c(0, 0),
    fs = sqrt(c(0.455, 0.455)), gs = c(1.365, 1.365), hs = c(0.071175, 0.071175),
    scale_diag = 1, scale_offdiag = 0,
    df = my_m + 4
) # create named list containing everything from the data block
my_iter <- 2000

output <- stan(paste(my_path, "/statprior.stan", sep = ""),
    data = my_data, chains = 4,
    warmup = my_iter / 2, iter = my_iter
) # Run 4 HMC chains, each for 1000 warm-up
# iterations, then 2000 post warm-up iterations

## Extract MCMC output as 3-dimensional array of form: [iterations, chains, parameters]
draws_arr <- as.array(output)

saveRDS(draws_arr, file = "./HeapsStanPrograms/rstanStatPriorSamps.rds")

## Alternatively, run HMC using cmdstanr
# All the Stan programs can be run using cmdstanr
library(cmdstanr) # load the cmdstanr package

# As an example, run the code for the stationary, exchangeable prior
my_p <- 4 # fit model of order 4
my_data <- list(
    m = my_m, p = my_p, N = nrow(yraw_trim_hb), y = yraw_trim_hb, es = c(0, 0),
    fs = sqrt(c(0.455, 0.455)), gs = c(1.365, 1.365), hs = c(0.071175, 0.071175),
    scale_diag = 1, scale_offdiag = 0,
    df = my_m + 4
) # create named list containing everything from the data block
my_iter <- 2000
mod <- cmdstan_model(paste(my_path, "/statprior.stan", sep = ""))
output <- mod$sample(
    data = my_data, chains = 4, parallel_chains = 4, iter_warmup = my_iter / 2,
    iter_sampling = my_iter / 2
) # Run 4 HMC chains, each for 1000 warm-up
# iterations, then 2000 post warm-up iterations

## Extract MCMC output as 3-dimensional array of form: [iterations, chains, parameters]
draws_arr <- output$draws()
draws_arr <- unclass(draws_arr)
dimnames(draws_arr) <- list(
    iterations = NULL, chains = paste("chain:", 1:dim(draws_arr)[2], sep = ""),
    parameters = dimnames(draws_arr)[[3]]
)
saveRDS(draws_arr, file = "./HeapsStanPrograms/cmdstanrStatPriorSamps.rds")