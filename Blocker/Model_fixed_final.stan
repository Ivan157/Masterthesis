// Program 1(d): Binomial likelihood, logit link, Fixed Effects (Blocker example) 

// Binomial likelihood, logit link 
// Fixed effects model for multi-arm trials 

functions {
  real to_real(int r){ return r; }    // for generating real values for leverage statistics
  real to_real_II(int n){ return n; } // for generating real values for leverage statistics
}



data {
  int<lower=1> NO; // Total number of observations/data points conforms [i,k]; Tip: "n_ik" hier ungeeignet => nur so kann man "for i in 1:NO" bzw "[1:NO]" setzen
  int<lower=1> NT; // Total number of treatments; Restriction because of begin natural numbers without 0
  int<lower=1> NS; // Total number of studies
  
  int<lower=0> n[NO]; // Number of observations at each trial in each single study
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
  d_II[1] = 0.0;                         // treatment effect is zero for reference treatment 

}  

model {
  A ~ normal(-2.2, 1/3.3);
  mu[s[1:NO]] ~ normal(0, 10);
  d[t[1:NO]] ~ normal(0, 10);
  r[1:NO] ~ binomial_logit(n[1:NO], mu[s[1:NO]] + d_II[t[1:NO]] - d_II[base[1:NO]]); // binomial likelihood 

}



generated quantities { // further diagnostics
  vector[NO] rhat;     // for leverage plot
  vector[NO] p_II;     // for leverage plot
  vector[NO] dev;
  vector[NT] T;        // for absolute effects
  real totalresdev;

  for(i in 1:NO) {     // generate single values for leverage plot
    p_II[i] = exp(mu[s[i]] + d_II[t[i]] - d_II[base[i]]) / (1 + exp(mu[s[i]] + d_II[t[i]] - d_II[base[i]])); // model for linear predictor 
    
    rhat[i] = p_II[i] * n[i];      // expected value of the numerators
    
     dev[i]   =  2 *(               // Deviance contribution 
			        to_real(r[i]) * log( to_real(r[i]) /rhat[i]) 
		 	        + ( to_real_II(n[i]) - to_real(r[i])) * log( (to_real_II(n[i]) - to_real(r[i])) /(to_real_II(n[i]) - rhat[i])   ) ); 
  }
  
  totalresdev   = sum(dev);        // Total Residual Deviance 

  T[1:NT] = exp(A + d_II[1:NT]) ./(1 + exp(A + d_II[1:NT])); 


}




