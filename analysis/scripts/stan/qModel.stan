// Define functions. 
functions 
{
  /*
  * Inverse logit transform.*
  * @param x Real. Log odds.  
  */
  real expit(real x) 
  {
  return exp(x)/(1 + exp(x));
  }
  
  /*
  * Utility function.*
  * @param p        Ratio. Event rate for treatment
  * @param utility0 Real. Utility of non-event
  * @param utility1 Real. Utility of event
  */
  real utility(real p, real utility0, real utility1) 
  {
  return utility0 * (1 - p) + utility1 * p;
  }
}

// Define data. 
data 
{
  int iDesign;
  int jDesign;
  int iPredictor;
  int jPredictor;
  matrix[iDesign, jDesign] X;
  int Y[iDesign];
  matrix[iPredictor, jPredictor] predictor;
  matrix[iPredictor, 2] utils;

}

// Define model parameters. 
parameters 
{
  vector[jPredictor] beta;
}

// Define transformed parameters.
transformed parameters
{
  vector[iPredictor] p;
  vector[iPredictor] u;
  
  // Generate posterior probabilities under each treatment sequence.
  p =  predictor * beta;
  for(i in 1:iPredictor) p[i] = expit(p[i]);
  
  // Weight the utilities of each decision by the probabilites.
  for(i in 1:iPredictor) u[i] = utility(p[i], utils[i, 1], utils[i, 2]);
}

// Define model likelihood and priors. 
model 
{
  beta ~ normal(0, 5);
  if(iDesign == 0)
  {
    # This is a hack to get around the noted bug of stan's underlying C++ libs not allowing zero length arrays to be multiplied.
    # Basically we just ignore the likelihood and sample from the priors. 
  } else{
    Y ~ bernoulli_logit(X * beta);  
  }
  
}
