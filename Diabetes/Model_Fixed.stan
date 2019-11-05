// Program 3(b):  Binomial likelihood, cloglog link, Fixed Effects (Diabetes example) 

// Binomial likelihood, cloglog link 
// Fixed effects model for multi-arm trials 




functions {
  real to_real(int r){ return r; }    // for generating real values for pD statistics
  real to_real_II(int n){ return n; } // for generating real values for leverage statistics
}




data {
  int<lower=1> NO; // Total number of observations/data points conforms [i,k]; Tip: "n_ik" hier ungeeignet => nur so kann man "for i in 1:NO" bzw "[1:NO]" setzen
  int<lower=1> NT; // Total number of treatments; Restriction because of begin natural numbers without 0
  int<lower=1> NS; // Total number of studies
  
  int<lower=1> r[NO]; // Number of events ("successes") in the control group (r_1) and under active treatment (r_2)  
  int<lower=1, upper=NT> t[NO]; // // Indicates, which treatment is just observed
  int<lower=1, upper=NS> s[NO]; // indicates seperate study   
  int<lower=1, upper=NT> base[NO]; // Indicates baseline (first or "base-case" treatment)
  
  int<lower=1> n[NO]; // Number of observations at each trial in each single study
  real<lower=1> time[NO];   // Exposure times in years
} 




parameters {  
  vector[NT] d;      // seperate treatment effect  
  vector[NS] mu; 
  vector[NT] A; 
}  




transformed parameters {
  vector[NT] d_II=d;                     // seperate treatment effect + d[1]=0 
  
  real p[NO];
  real<lower=0.00000001> p_cloglog[NO];
 
  d_II[1] = 0.0;                         // treatment effect is zero for reference treatment 
  
 
  for (i in 1:NO){ 
    p[i] = log(time[i]) +  mu[s[i]] + d_II[t[i]] - d_II[base[i]] ;
    p_cloglog[i] = inv_cloglog(p[i]);
  } 
}  




model {
  A ~ normal(-4.2, 0.9491579); 
  // target += normal_lpdf(A| -3,1/1.77);    

  mu[s[1:NO]] ~ normal(0.000001, 1);
  // target += normal_lpdf(mu[s[1:NO]]| 0,1);    
  
  d[2:NT] ~ normal(0, 1);
  //  target += normal_lpdf(d[t[1:NO]]| 0,1);                        // vague priors for treatment effects 
  // durch "target" erreichen viel mehr Werte ein Rhat von 1

  for (i in 1:NO)
    r[i] ~ binomial(n[i], p_cloglog[i]); // binomial likelihood 

}



generated quantities { // further diagnostics
  vector[NO] rhat;   
  vector[NO] dev;
  vector[NT] T;        // for absolute effects
  vector[NT] T_cloglog;
  real totalresdev;

  for(i in 1:NO) {     
     rhat[i] = p_cloglog[i] * n[i] ;
       dev[i]   =  2 *(               // Deviance contribution 
     			        to_real(r[i]) *  ( log( to_real(r[i]))  - log(rhat[i]) )
      		 	        + ( to_real_II(n[i]) - to_real(r[i])) * ( log (to_real_II(n[i]) - to_real(r[i]))  - log (to_real_II(n[i]) - rhat[i]) )); 
  }
  
  totalresdev   = sum(dev);        // Total Residual Deviance 
  
  //Provide estimates of treatment effects T[k] on the natural (probability) scale  
  // Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA
  T[1:NT] = log(3) + A + d_II[1:NT] ; 
  T_cloglog [1:NT] = inv_cloglog(T[1:NT]);
}



