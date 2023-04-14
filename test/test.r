# library(BMIMIC)
library(rstan)

for(i in fs::dir_ls("R")) {source(i)}

generated_data <-
  mkData(N1         = 100,
         N2         = 100,
         nitem      = 10,
         mu_2       = 0.2,
         var_2      = 2,
         uni_dif    = c(0,0,0,0,0.2,0.2,0.2,0,0,0),
         nonuni_dif = c(0,0,0,0,0,0,0,0.2,0.2,0.2),
         anchor_n   = 4)

# stan_fit <-
#   rstan::stan(
#     file = file.path(system.file("stan", package = "BMIMIC"),
#                      "BMIMIC_example.stan"),
#     data   = generated_data$stan_data
#   )

# summary(stan_fit)

# https://github.com/stan-dev/rstan/wiki/Configuring-C---Toolchain-for-Windows#r40
# Sys.getenv("BINPREF")
# readLines("~/.Rprofile")
#
# writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
# Sys.which("make")
#
#
# cat("CXX14FLAGS += -O3 -mmmx -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2",file = "~/.R/Makevars.win", sep = "\n", append = FALSE)

# install.packages("StanHeaders", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
# install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

generated_data$stan_data$resp <- generated_data$response

stan_fit <-
  rstan::stan(
    file = file.path("test/BMIMIC_example_1.stan"),
    data   = generated_data$stan_data
  )
