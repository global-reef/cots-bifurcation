# R/functions_bifurcation.R
# Helper functions for CoTS bifurcation models


# Basic *cubic* bifurcation-style population function
# N = CoTS density
# r = growth rate
# K = upper carrying capacity / outbreak equilibrium
# A = Allee threshold / lower unstable threshold
# m = mortality/removal term

cots_growth <- function(N, r, K, A, m = 0) {
  
  dN_dt <- r * N * (1 - N / K) * (N / A - 1) - m * N
  
  return(dN_dt)
}
# This means: 
# N < A        negative growth
# A < N < K    positive growth
# N > K        negative growth
# and has three equilibria: 
# N = 0  stable 
# N = A  unstable threshold 
# N = K stable high-density equilibrium

# Generate model output across a sequence of CoTS densities
generate_bifurcation_curve <- function(
    density_seq,
    r,
    K,
    A,
    m = 0
) {
  
  output <- data.frame(
    density = density_seq,
    growth_rate = cots_growth(
      N = density_seq,
      r = r,
      K = K,
      A = A,
      m = m
    )
  )
  
  return(output)
}

##### ashna's framework 
# core model 
dN/dt = r(T) * (N - State1) * (1 - N / State3) * (N / State2 - 1)
# where: 
# State 1 = endemic equilibrium, e.g. 3 CoTS/ha
# State 2 = unstable tipping threshold, e.g. 5–6 CoTS/ha = K_t = proportion of K_c defined by % live coral 
# State 3 = outbreak equilibrium/carrying capacity, e.g. 16 CoTS/ha = K_c

# and then our simulations should test different combinations of parameteres for state 1, 2 and 3 

# for example 
state_scenarios <- list(
  ashna_original = list(S1 = 3,  S2 = 6,  S3 = 16),
  koh_tao_low    = list(S1 = 5,  S2 = 15, S3 = 40),
  koh_tao_mid    = list(S1 = 10, S2 = 30, S3 = 75),
  koh_tao_high   = list(S1 = 15, S2 = 50, S3 = 100)
)
