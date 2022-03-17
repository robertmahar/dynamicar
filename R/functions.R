#' Simulate SMART patient
#'
#' \code{simulatePatients} simulates a single patient for a two-stage, two-treatment, degenerate SMART.
#' 
#' @param scheme A dataframe that defines the structure of the SMART in terms of stages (`stage`),
#'  treatments and interactions (`a1`, `a2`, `a1a2`), the intermediate outcome (`y1`), and initial randomisation
#'  probabilities (`rho`), and event rates for the specified treatments (`pi`). 
#'
#' @return A dataframe containing simulated data for a single trial participant.
#' @export
#'
#' @examples 
#'
#' For a two-stage, two-treatment, degenerate (i.e. participants exit the trial if event) we simulate a participant
#' with equal randomisation and event rates as follows:
#'
#' scheme <- data.frame(
#' stage = c(   1,    1,   2,   2,   2,   2),
#' int   = c(   1,    1,   1,   1,   1,   1),
#' a1    = c(   0,    1,   0,   0,   1,   1),
#' y1    = c(   0,    0,   1,   1,   1,   1),
#' a2    = c(  NA,   NA,   0,   1,   0,   1),
#' a1a2  = c(  NA,   NA,   0,   0,   0,   1),
#' rho   = c( 0.5,  0.5, 0.5, 0.5, 0.5, 0.5),
#' pi    = c( 0.8,  0.7, 0.5, 0.4, 0.5, 0.4))
#' 
#' simulatePatients(scheme)
#' 
simulatePatients <- function(scheme)
{
  rho <- scheme[["rho"]][which(scheme[["a1"]] == 1 & scheme[["stage"]] == 1)]
  a1  <- rbinom(1, 1, rho) 
  
  pi <- scheme[["pi"]][which(scheme[["a1"]] == a1 & scheme[["stage"]] == 1)]
  y1 <- rbinom(1, 1, pi)
  if(y1 == 0)
  {
    a2 <- NA
    y2 <- NA
  }
  if(y1 == 1)
  {
    rho  <- scheme[["rho"]][which(scheme[["a1"]] == a1 & scheme[["y1"]] == 1 & scheme[["a2"]] == 1)]
    a2   <- rbinom(1, 1, rho)
    if(a2 == 0)
    {
      pi <- scheme[["pi"]][which(scheme[["a1"]] == a1 & scheme[["a2"]] == 0)]
      y2 <- rbinom(1, 1, pi)
    }
    if(a2 == 1)
    {
      pi <- scheme[["pi"]][which(scheme[["a1"]] == a1 & scheme[["a2"]] == 1)]
      y2 <- rbinom(1, 1, pi)
    }
  }
  df <- data.frame("int" = 1, a1, y1, a2, "a1a2" = a1 * a2, y2)
  df
}

#' Analyse SMART data
#'
#' \code{doAnalysis} analyses simulated SMART data using a Q-learning decision-theoretic model implemented in Rstan. Requires the Rstan model (`qModel`) to be compiled before running. 
#' 
#' @param data   A dataframe containing simulated SMART data. 
#' @param scheme A dataframe that defines the structure of the SMART in terms of stages (`stage`),
#'  treatments and interactions (`a1`, `a2`, `a1a2`), the intermediate outcome (`y1`), and initial randomisation
#'  probabilities (`rho`), and event rates for the specified treatments (`pi`). 
#' @param u A vector of utilities.
#' @param myopic A logical. Does this analysis use 'myopic' response adaptive randomisation (as opposed to 'dynamic' response adaptive randomisation)? Defaults to FALSE.
#' @param tuning A real values scalar used to define whether any randomisation is equal (`tuning` equals zero) or adaptive (`tuning` greater than zero, with adaptive algorithm being 'greedier' for higher values). Defaults to zero.
#'
#' @return
#' @export
#' @import rstan 
#' 
#' @examples
#' 
#' For a two-stage, two-treatment, degenerate (i.e. participants exit the trial if event) we simulate a participant
#' with equal randomisation and event rates as follows:
#'
#' scheme <- data.frame(
#' stage = c(   1,    1,   2,   2,   2,   2),
#' int   = c(   1,    1,   1,   1,   1,   1),
#' a1    = c(   0,    1,   0,   0,   1,   1),
#' y1    = c(   0,    0,   1,   1,   1,   1),
#' a2    = c(  NA,   NA,   0,   1,   0,   1),
#' a1a2  = c(  NA,   NA,   0,   0,   0,   1),
#' rho   = c( 0.5,  0.5, 0.5, 0.5, 0.5, 0.5),
#' pi    = c( 0.8,  0.7, 0.5, 0.4, 0.5, 0.4))
#' 
#' data <- lapply(1:1000, function(x) simulatePatients(scheme))
#' data <- do.call(rbind, data)
#' 
#' doAnalysis(data, scheme)
#' 
doAnalysis <- function(data, 
                       scheme, 
                       u = c(0, 0, 0, 0, 1,  1,  1,  1,  1,  1), 
                       myopic = F, 
                       tuning = 0)
{
  fits <- list()

  # ===========================
  # Mung the second stage data.
  # Filter the second stage data
  q <- data[complete.cases(data),]
  
  # Get unique predictor matrix
  predictor <- scheme[c("int", "a1", "a2", "a1a2")] 
  predictor <- predictor[complete.cases(predictor),]
  
  # Myopic thinking in stage 2 (i.e. don't use information from the first stage).
  if(myopic == T)
  {
    predictor <- predictor[c("int", "a2")]
  }

  # Put everything in a list.
  stanData <- list("X"          = q[colnames(predictor)],
                   "Y"          = as.array(q[[ncol(q)]]),
                   "iDesign"    = nrow(q),
                   "jDesign"    = ncol(predictor),
                   "iPredictor" = nrow(predictor),
                   "jPredictor" = ncol(predictor),
                   "predictor"  = predictor,
                   "utils"      = matrix(ncol = 2, byrow = T, c(u[8], u[4],
                                                                u[7], u[3],
                                                                u[6], u[2],
                                                                u[5], u[1])))
  
  # Function to keep things `simple'.
  f1 <- function()
  {
    stanFit <- sampling(qModel, 
                        stanData, 
                        chains = 4, 
                        iter = 1000, 
                        refresh = 0)
    
    # Extract the expected/maximum/total/relative utilities and optimal decisions.
    expectedUtility      <- summary(stanFit, "u")[["summary"]][,"mean"]
    utilities            <- data.frame(expectedUtility, "action" = c(0, 1))
    utilities[["point"]] <- rep(1:(nrow(utilities)/2), each = 2)
    utilities            <- ddply(utilities, "point", mutate, expectedUtilityTotal = sum(expectedUtility))
    utilities[["rho"]]   <- with(utilities, (expectedUtility ^ tuning)/
                                   (expectedUtility ^ tuning + (expectedUtilityTotal - expectedUtility) ^ tuning))
    utilities            <- ddply(utilities, "point", mutate, expectedUtilityMax = max(expectedUtility))
    
    utilities <- utilities[, c("point", "action", "expectedUtility", "expectedUtilityMax", "expectedUtilityTotal", "rho")]
    
    return(list("utilities" = utilities, 
                "stanFit"   = stanFit))
  }
  
  f <- f1()
  
  # Save the results.
  results <- cbind("stage" = 2, f[["utilities"]])
  fits[["q2"]] <- f[["stanFit"]]
  
  # ===========================
  # Subset the second stage data
  q <- data[-grep(colnames(data), pattern = "2")]
  
  predictor <- scheme[c("int", "a1", "a2", "a1a2")]
  predictor <- predictor[!complete.cases(predictor),]
  predictor <- predictor[c("int", "a1")]
  
  stanData <- list("X"          = q[colnames(predictor)],
                   "Y"          = as.array(q[[ncol(q)]]),
                   "iDesign"    = nrow(q),
                   "jDesign"    = ncol(predictor),
                   "iPredictor" = nrow(predictor),
                   "jPredictor" = ncol(predictor),
                   "predictor"  = predictor,
                   "utils"      = cbind(c(u[10], u[9]), unique(f[["utilities"]][["expectedUtilityMax"]])))
  
  # Myopic thinking in stage 1.
  if(myopic == T)
  {
    foo <- 0
    bar <- 0
    stanData[["utils"]][,2] <- c(foo, bar)
  }
  
  f <- f1()
  
  results <- rbind(cbind("stage" = 1, f[["utilities"]]), results)
  
  fits[["q1"]] <- f[["stanFit"]]
  return(list("inference" = fits, 
              "decisions" = results))
}


#' Get known utilities.
#' 
#' \code{getKnownUtilities} obtains the true utilities for a given utility set and trial scheme.
#'
#' @param scheme A dataframe that defines the structure of the SMART in terms of stages (`stage`),
#'  treatments and interactions (`a1`, `a2`, `a1a2`), the intermediate outcome (`y1`), and initial randomisation
#'  probabilities (`rho`), and event rates for the specified treatments (`pi`). 
#' @param u A vector of utilities.
#'
#' @return
#' @export
#' @import plyr
#' 
#' @examples
#'
#' For a two-stage, two-treatment, degenerate (i.e. participants exit the trial if event) we simulate a participant
#' with equal randomisation and event rates as follows:
#'
#' scheme <- data.frame(
#' stage = c(   1,    1,   2,   2,   2,   2),
#' int   = c(   1,    1,   1,   1,   1,   1),
#' a1    = c(   0,    1,   0,   0,   1,   1),
#' y1    = c(   0,    0,   1,   1,   1,   1),
#' a2    = c(  NA,   NA,   0,   1,   0,   1),
#' a1a2  = c(  NA,   NA,   0,   0,   0,   1),
#' rho   = c( 0.5,  0.5, 0.5, 0.5, 0.5, 0.5),
#' pi    = c( 0.8,  0.7, 0.5, 0.4, 0.5, 0.4))
#' 
#' getKnownUtility(scheme)
#' 
getKnownUtility <- function(scheme, 
                            u = c(0, 0, 0, 0, 1,  1,  1,  1,  1,  1))
{ 
  scheme[["eventNeg"]] <- u[10:5]
  scheme[["eventPos"]] <- c(NA, NA, u[4:1])
  scheme[["expectedUtility"]] <- (1 - scheme[["pi"]]) * scheme[["eventNeg"]] + scheme[["pi"]] * scheme[["eventPos"]]  
  scheme <- ddply(scheme, .(stage, a1), mutate, optimalUtility = max(expectedUtility))
  qUtility <- ddply(scheme, .(stage, a1), summarise, "eventPos" = mean(optimalUtility))
  qUtility <- qUtility[complete.cases(qUtility),]
  
  ind      <- which(scheme[["stage"]] == qUtility[["stage"]] - 1 & scheme[["a1"]] == qUtility[["a1"]])
  scheme[["eventPos"]][ind] <- qUtility[["eventPos"]]
  scheme[["expectedUtility"]] <- (1 - scheme[["pi"]]) * scheme[["eventNeg"]] + scheme[["pi"]] * scheme[["eventPos"]]  
  
  ind   <- which(is.na(scheme[["optimalUtility"]])) 
  scheme[["optimalUtility"]][ind] <- max(scheme[["expectedUtility"]][ind])
  scheme
}

#' Simulate a decision-theoretic SMART with response adaptive randomisation.
#'
#' \code{simulateTrial} simulates a decision-theoretic SMART with either fixed randomisation, or 'myopic' or 'dynamic' response adaptive randomisation. 
#' 
#' @param scheme A dataframe that defines the structure of the SMART in terms of stages (`stage`),
#'  treatments and interactions (`a1`, `a2`, `a1a2`), the intermediate outcome (`y1`), and initial randomisation
#'  probabilities (`rho`), and event rates for the specified treatments (`pi`). 
#' @param interimSize A positive integer indicating the size of (equal sized and complementary) analysis cohorts used in the analysis.
#' @param interimNumber A positive integer indicating the number of interim analyses to be performed (include the final analysis). 
#' @param u A vector of utilities.
#' @param myopic A logical. Does this analysis use 'myopic' response adaptive randomisation (as opposed to 'dynamic' response adaptive randomisation)? Defaults to FALSE.
#' @param tuning A real values scalar used to define whether any randomisation is equal (`tuning` equals zero) or adaptive (`tuning` greater than zero, with adaptive algorithm being 'greedier' for higher values). Defaults to zero.
#'
#' @return
#' @export
#' @import rstan
#'
#' @examples
#' 
#' For a two-stage, two-treatment, degenerate (i.e. participants exit the trial if event) we simulate a participant
#' with equal randomisation and event rates as follows:
#'
#' scheme <- data.frame(
#' stage = c(   1,    1,   2,   2,   2,   2),
#' int   = c(   1,    1,   1,   1,   1,   1),
#' a1    = c(   0,    1,   0,   0,   1,   1),
#' y1    = c(   0,    0,   1,   1,   1,   1),
#' a2    = c(  NA,   NA,   0,   1,   0,   1),
#' a1a2  = c(  NA,   NA,   0,   0,   0,   1),
#' rho   = c( 0.5,  0.5, 0.5, 0.5, 0.5, 0.5),
#' pi    = c( 0.8,  0.7, 0.5, 0.4, 0.5, 0.4))
#' 
#' utils <- c(0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
#' simulateTrial(scheme, interimSize = 500, interimNumber = 4, u = utils)
#' 
simulateTrial <- function(scheme, 
                          interimSize, 
                          interimNumber, 
                          tuning = 0, 
                          myopic = F, 
                          u) 
{
    
  # Possible histories.
  histories <- c("1 0 0 NA NA NA",
                 "1 0 1 0 0 0",
                 "1 0 1 0 0 1",
                 "1 0 1 1 0 0",
                 "1 0 1 1 0 1",
                 "1 1 0 NA NA NA", 
                 "1 1 1 0 0 0", 
                 "1 1 1 0 0 1",
                 "1 1 1 1 1 0",
                 "1 1 1 1 1 1")  
  
  util_order <- c(10, 8, 4, 7, 3, 9, 6, 2, 5, 1)
  
  # Set up holders for things we want to collect.
  results      <- list()
  data         <- c()
  expectations <- c()
  utilities    <- c()   
  
  if(tuning == 0)
  {
    interimSize   <- interimSize * interimNumber
    interimNumber <- 1
  }
  
  for(i in seq(interimNumber))
  {
    data           <- rbind(data, do.call(rbind, lapply(seq(interimSize), function(x) simulatePatients(scheme = scheme)))) 
    results[[i]]   <- doAnalysis(data = data, scheme = scheme, u = u, myopic = myopic, tuning = tuning)

    expectations[[i]] <- cbind(i, 
                            results[[i]][["decisions"]], 
                            "knownOptimalUtility"  = getKnownUtility(scheme = scheme, u = u)[["optimalUtility"]], 
                            "knownExpectedUtility" = getKnownUtility(scheme = scheme, u = u)[["expectedUtility"]])

    expectations[[i]][["rho"]] <- ifelse(expectations[[i]][["rho"]] < 0, 0, 
                                  ifelse(expectations[[i]][["rho"]] > 1, 1, expectations[[i]][["rho"]]))
      
    # Update the randomisation probabilities.
    if(tuning > 0) 
    {
      scheme[["rho"]] <-  expectations[[i]][["rho"]]
    } 

    # Get the utility of patients on the trial.
    utilities <- c(utilities, sapply(match(apply(data, 1, paste, collapse = " "), histories[order(util_order)]), function(x) u[x]))
  }
  
  expectations <- do.call(rbind, expectations)
  utilities    <- mean(utilities)

  return(list("expectations" = expectations, 
              "utilities"    = utilities))
}

# End of script.



