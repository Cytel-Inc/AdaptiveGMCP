# File for testing survival simulations for the p-value combination method.

meth <- "CombPValue"
alp <- 0.025
N <- 200
D <- 88
TSCont <- NA
TSBin <- NA
CC <- NA
CER_Z2 <- NA

# num_arms <- 2 # 1 treatment arm versus one control arm
num_arms <- 3 # 2 treatment arms versus one control arm

num_eps <- 1 # one endpoint
ep_type <- list(EP1 = "Survival")
arm_means <- NA
arm_stdevs <- NA
arm_props <- NA
# arm_lmbdas <- list(EP1 = c(0.03466, 0.03466)) # hazard rates
arm_lmbdas <- list(EP1 = c(0.03466, 0.03466, 0.03466)) # hazard rates
accr_start <- 0
accr_rates <- 8

# alloc_ratio <- c(1, 1) # balanced design
alloc_ratio <- c(1, 1, 1) # balanced design
corr <- NA
# weights <- 1 # Only 1 hypothesis to test
weights <- rep(1 / 2, 2)
# trans_mat <- matrix(0, nrow = 1, ncol = 1) # Only 1 hypothesis to test
trans_mat <- rbind(H1 = c(0, 1), H2 = c(1, 0))

test <- "Dunnett"
t <- 1 # Fixed sample design
sel <- FALSE
nSim <- 10
seed <- 1234

# Calling the simulation function
out1 <- simMAMSMEP(
  Method = meth,
  alpha = alp,
  SampleSize = N,
  totalEvents = D,
  TestStatCont = TSCont,
  TestStatBin = TSBin,
  UseCC = CC,
  FWERControl = CER_Z2,
  nArms = num_arms,
  nEps = num_eps,
  lEpType = ep_type,
  Arms.Mean = arm_means,
  Arms.std.dev = arm_stdevs,
  Arms.Prop = arm_props,
  Arms.Haz.Rates = arm_lmbdas,
  Accr.Start.Times = accr_start,
  Accr.Rates = accr_rates,
  Arms.alloc.ratio = alloc_ratio,
  EP.Corr = corr,
  WI = weights,
  G = trans_mat,
  test.type = test,
  info_frac = t,
  Selection = sel,
  nSimulation = nSim,
  Seed = seed,
  Parallel = FALSE,
  Verbose = FALSE
)

out1
