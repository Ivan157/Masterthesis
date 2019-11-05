// Program 3(a):  Binomial likelihood, cloglog link, Random Effects (Diabetes example) 


// Binomial likelihood, cloglog link 
// Random effects model for multi-arm trials 




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
  vector[NO] mu; 
  vector[NO] delta; 
  real<lower=0.000001> s_d; 
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
  real p[NO];
  real<lower=0.00000001> p_cloglog[NO];
  vector[NO] delta_II=delta;

  d_II[1] = 0.0;                         // treatment effect is zero for reference treatment 

  // LOOP THROUGH STUDIES 
  for (i in 1:NO) {         
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
 
 
  for (i in 1:NO){ 
    p[i] =    mu[i] + delta[i] + log(time[i]) ;
    p_cloglog[i] = inv_cloglog(p[i]);
  }
}  




model {
  A ~ normal(-4.2, 0.9491579); 
  delta[1:NO] ~ normal(md[1:NO], variance_d[1:NO]);                // trial-specific LOR distributions 
  s_d ~ uniform(0, 1);  
  // target += uniform_lpdf(s_d | 0,1);
  mu[s[1:NO]] ~ normal(0.000001, 1);
  // target += normal_lpdf(mu[s[1:NO]]| 0,1);    
  //  mu[1:NO] ~ normal(0, 1);
  
  d[2:NT] ~ normal(0, 1);
  // target += normal_lpdf(d[t[1:NO]]| 0,1);                        // vague priors for treatment effects 
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



