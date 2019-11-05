// Program 2(b):  Poisson likelihood, log link, Fixed Effects (Dietary fat) 
// Poisson likelihood, log link 
// Fixed effects model for multi-arm trials 


functions {
  real to_real(int r){ return r; }    // for generating real values for pD statistics
}



data {
  int<lower=1> NO; // Total number of observations/data points conforms [i,k]; Tip: "n_ik" hier ungeeignet => nur so kann man "for i in 1:NO" bzw "[1:NO]" setzen
  int<lower=1> NT; // Total number of treatments; Restriction because of begin natural numbers without 0
  int<lower=1> NS; // Total number of studies
  
  real<lower=0> E[NO];   // Exposure times in years
  int<lower=0> r[NO]; // Number of events ("successes") in the control group (r_1) and under active treatment (r_2)  
  int<lower=1, upper=NT> t[NO]; // // Indicates, which treatment is just observed
  int<lower=1, upper=NS> s[NO]; // indicates seperate study   
  int<lower=1, upper=NT> base[NO]; // Indicates baseline (first or "base-case" treatment)
} 



parameters {  
  vector[NT] d;      // seperate treatment effect  
  vector[NS] mu; 
  vector[NT] A; 
}  



transformed parameters {
  vector[NT] d_II=d;                     // seperate treatment effect + d[1]=0 
  vector[NO] lambda;     
  vector[NO] theta;     
  
  d_II[1] = 0.0;                         // treatment effect is zero for reference treatment 
  

  for(i in 1:NO){ 
 	  lambda[i] = exp(mu[s[i]] + d_II[t[i]] - d_II[base[i]]);  // model for linear predictor 
  	theta[i] = lambda[i] * E[i];              // failure rate * exposure 
  }

}  



model {
  A ~ normal(-3, 1/1.77); 
  // target += normal_lpdf(A| -3,1/1.77);    

  mu[s[1:NO]] ~ normal(0, 1);
  // target += normal_lpdf(mu[s[1:NO]]| 0,1);    
  
  
  target += normal_lpdf(d[t[1:NO]]| 0,1);                        // vague priors for treatment effects 
  // durch "target" erreichen viel mehr Werte ein Rhat von 1

  r[1:NO] ~ poisson(theta[1:NO]); // Poisson likelihood 

}



generated quantities { // further diagnostics
  vector[NO] dev;
  vector[NT] T;        // for absolute effects
  real totalresdev;

  for(i in 1:NO) {     
		  dev[i] = 2*((theta[i]-to_real(r[i])) + to_real(r[i])*log(to_real(r[i])/theta[i])) ;     // Deviance contribution 
  }
  
  totalresdev   = sum(dev);        // Total Residual Deviance 
  
  //Provide estimates of treatment effects T[k] on the natural (probability) scale  
  // Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA
  T[1:NT] = exp(A + d_II[1:NT]); 
  
}




