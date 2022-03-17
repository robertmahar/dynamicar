# Simulates and analyses a degenerate two-stage SMART.

# ============================================================================
# Cluster specific setup.

# Create indicator of whether this is a cluster job.
pbsScript <- ifelse(Sys.getenv("PBS_ARRAYID") == "", F, T) 

# Get any existing PBS index.
tIndex    <- as.numeric(Sys.getenv("PBS_ARRAYID"))

# Set the working directory to the current PBS directory.
if(pbsScript == T) 
{
  # Builds package on server.
  library(devtools)
  devtools::load_all(Sys.getenv("PATH_TO_PKG_ON_CLUSTER"))
  
  # Anchors the current working directory to the that of the pbs script.
  setwd(Sys.getenv("PBS_O_WORKDIR"))
  
  # Define file paths for batch job.
  dirOut     <- file.path("..", "..", "output")
  pathToStan <- file.path("..", "stan", "qModel.stan")
  fileOut <- file.path(dirOut, "runSimBatch", 
                       paste0("runSimBatch", tIndex, ".Rdata"))
  
} else {
  
  # Define file paths for local job.
  dirOut     <- file.path("analysis", "output")
  pathToStan <- file.path("analysis", "scripts", "stan", "qModel.stan")
  fileOut <- file.path(dirOut, "runSim.Rdata")
}

# ============================================================================
# Simulation setup.

# Number of simulations.
nSims <- 10

# Preference ordering of utilities. 
utils <- c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1)

# List of randomisation scheme.
scheme <-  data.frame(
  stage = c(   1,    1,   2,   2,   2,   2),
  int   = c(   1,    1,   1,   1,   1,   1),
  a1    = c(   0,    1,   0,   0,   1,   1),
  y1    = c(   0,    0,   1,   1,   1,   1),
  a2    = c(  NA,   NA,   0,   1,   0,   1),
  a1a2  = c(  NA,   NA,   0,   0,   0,   1),
  rho   = c( 0.5,  0.5, 0.5, 0.5, 0.5, 0.5))

# Create the first-stage event probability pairings, drop dups.
pi_1 <- expand.grid(pbo_1 = seq(0, 100, 5), trt_1 = seq(0, 100, 5))

# Create the second-stage marginal event probability pairings.
pi_2 <- expand.grid(pbo_21 = c(5, 10, 20, 40, 60, 80, 90, 95), 
                    pbo_22 = c(5, 10, 20, 40, 60, 80, 90, 95))
pi_2 <- cbind(pi_2, "index" = seq(nrow(pi_2)))

# Nest the first stage pairings within the second stage pairings.
pi <- do.call(rbind, lapply(pi_2[["index"]], function(x) cbind(pi_1, pi_2[x,])))
pi[["trt_21"]] <- pi[["pbo_21"]]
pi[["trt_22"]] <- pi[["pbo_22"]]

# Remove the `building blocks'.
rm(pi_1, pi_2) 

# Add in adaptive information to make the scripts work.
grid <- rbind(cbind(pi, "interimSize" = 500,  "interimNumber" = 4),
              cbind(pi, "interimSize" = 2000, "interimNumber" = 1))
grid[["tuning"]] <- with(grid, ifelse(interimNumber == 1, 0, 1))

# Add myopic trigger.
grid <- rbind(cbind(grid, "myopic" = FALSE),
              cbind(grid, "myopic" = TRUE))

# Now take the fractions.
piNames       <- grep(colnames(grid), pattern = "_", value = T)
grid[piNames] <- grid[piNames]/100 

if(pbsScript == T) 
{
  # If batch job, use the grid defined by batch index.
  grid <- grid[grid[["index"]] == tIndex,]  
} else{
  # If local job, just use the first two examples of the entire grid.
  grid <- grid[sample(seq(nrow(grid)), 2),]
}
# ============================================================================
# Compile the stan model.
qModel <- stan_model(file = pathToStan)

# ============================================================================
# Run the simulation and place output in list.
simulations <- list()

# The simulation loop in all of its glory. 
# Note that if this is batch script each embedded simulation 'truth'
# Is sent to compute on an individual HPC core (or local machine). 
for(i in seq(nrow(grid)))
{
  scenario <- grid[i,]
  
  scheme[["pi"]] <- t(scenario[c("pbo_1", 
                                 "trt_1", 
                                 "pbo_21", 
                                 "trt_21",  
                                 "pbo_22", 
                                 "trt_22")])

  details <- paste(colnames(scenario), "=", scenario, collapse = ", ")
  
  simulations[[details]] <- lapply(seq(nSims),
    function(x) simulateTrial(
    scheme        = scheme,
    u             = utils,
    tuning        = scenario[["tuning"]],
    myopic        = scenario[["myopic"]],
    interimSize   = scenario[["interimSize"]],
    interimNumber = scenario[["interimNumber"]]))
}

# ============================================================================
# Save the output. 

save(simulations, file = fileOut)

# End of script.

