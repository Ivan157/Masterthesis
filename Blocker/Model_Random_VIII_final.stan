// Program 1(c): Binomial likelihood, logit link, Random Effects (Blocker example) 

// Binomial likelihood, logit link 
// Random effects model for multi-arm trials 

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
  vector[NO] delta;   
  real<lower=0> s_d; 
  vector[NT] A; 
}  



transformed parameters {
  vector[NT] d_II=d;                     // seperate treatment effect + d[1]=0 
  real<lower=0> tau=sqrt(s_d);           // between-trial precision
  real w[NS,NT]=rep_array(0.0, NS, NT);  // & ff.: for NMA
  vector[NO] sw_vec;                     // vector of sw; for NMA
  vector[NO] vard_vec;
  vector[NO] variance_d;
  vector[NO] md;
  
  vector[NO] delta_II=delta;
  d_II[1] = 0.0;                         // treatment effect is zero for reference treatment 
  
  
  for (i in 1:NO) {         // LOOP THROUGH STUDIES 
    if (t[i] == 1) {        //  LOOP THROUGH ARMS; adjustment for multi-arm trials and treatment effect is zero for control arm, so all further steps for [,1] are set to zero
      delta_II[i] = 0.0;    // treatment effect is zero for control arm 
      w[s[i],t[i]] = 0.0;   // adjustment for multi-arm trials is zero for control arm
      vard_vec[i] = 1E-2;
      sw_vec[i] = 0.0;      // cumulative adjustment for multi-arm trials 
    } else {
      w[s[i],t[i]] = delta_II[i] - d_II[t[i]] + d_II[base[i]]; // adjustment for multi-arm RCTs
      vard_vec[i] = (pow(tau,2) * t[i])/(2 * (t[i]-1));        // Transformation of precision to variance
      sw_vec[i] = sum(w[s[i],1:(t[i]-1)])/(t[i]-1); 
    }
    variance_d[i] = sqrt(vard_vec[i]);                         // Variance of LOR distributions (with multi-arm trial correction); "conforms" variance_d in TSD2
  }
  md[1:NO] = d_II[t[1:NO]] - d_II[base[1:NO]] + sw_vec[1:NO];  // mean of LOR distributions (with multi-arm trial correction) 
}  



model {
  A ~ normal(-2.2, 1/3.3);
 r[1:NO] ~ binomial_logit(n[1:NO], mu[s[1:NO]] + delta_II[1:NO]); // binomial likelihood 

  delta[1:NO] ~ normal(md[1:NO], variance_d[1:NO]);                // trial-specific LOR distributions 
  //mu[s[1:NO]] ~ normal(0, 10000);                          // vague priors for all trial baselines 
  //d[t[1:NO]] ~ normal(0, 10000);                           // vague priors for treatment effects 
  //s_d ~ uniform(0, 10);                                             // vague prior for between-trial SD. // beta=5 oder =10 machen nichts aus
  s_d ~ uniform(0, 5);  
  //mu[s[1:NO]] ~ normal(0, 10);                          // vague priors for all trial baselines 
  //d[t[1:NO]] ~ normal(0, 10);                           // vague priors for treatment effects 
  //s_d ~ uniform(0, 10);                                             // vague prior for between-trial SD.
  //mu[s[1:NO]] ~ normal(0, 1000);                          // vague priors for all trial baselines 
  mu[s[1:NO]] ~ normal(0, 10);
  // d[t[1:NO]] ~ normal(0, 100);                           // vague priors for treatment effects 
  // d[r[1:NT]] ~ normal(0, 10000);                           // vague priors for treatment effects 
  //d[t[1:NO]] ~ normal(0, 1000);                           // vague priors for treatment effects
  //d[t[1:NO]] ~ normal(0, 500);                           // vague priors for treatment effects
  //   ///////////////////////////
  //d[t[1:NO]] ~ normal(0, 10000);                           // vague priors for treatment effects 
  //d[t[1:NO]] ~ normal(0, 1000);                           // vague priors for treatment effects 
  //d[t[1:NO]] ~ normal(0, 100);                           // vague priors for treatment effects 
  d[t[1:NO]] ~ normal(0, 10);                           // vague priors for treatment effects 
  //d[t[1:NO]] ~ normal(0, 200);                           // vague priors for treatment effects 
  //d[t[1:NO]] ~ normal(0, 150);                           // vague priors for treatment effects   
}



generated quantities { // further diagnostics
  vector[NO] rhat;     // for leverage plot
  vector[NO] p_II;     // for leverage plot
  vector[NO] dev;
  vector[NT] T;        // for absolute effects
  real totalresdev;

  for(i in 1:NO) {     // generate single values for leverage plot
    p_II[i] = exp(mu[s[i]] + delta_II[i]) / (1 + exp(mu[s[i]] + delta_II[i])); // model for linear predictor 
    rhat[i] = p_II[i] * n[i];      // expected value of the numerators
    
     dev[i]   =  2 *(               // Deviance contribution 
			        to_real(r[i]) * log( to_real(r[i]) /rhat[i]) 
		 	        + ( to_real_II(n[i]) - to_real(r[i])) * log( (to_real_II(n[i]) - to_real(r[i])) /(to_real_II(n[i]) - rhat[i])   ) ); 
  }
  
  totalresdev   = sum(dev);        // Total Residual Deviance 
  
   //Provide estimates of treatment effects T[k] on the natural (probability) scale  
  // Given a Mean Effect, meanA, for 'standard' treatment 1, with precision (1/variance) precA 
  T[1:NT] = exp(A + d_II[1:NT]) ./(1 + exp(A + d_II[1:NT])); 

   
}




