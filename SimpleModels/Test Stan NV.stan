// Test Normalverteilung


data {
}

parameters {
  real a; 
  real b; 
  real c; 
  real d; 
  // real e; 
  real f; 
  real g; 
  real h;
  real i; 
}



model {
  	a ~ normal(0, inv_sqrt(0.0001)); // Referenz
		
		// variiert wird die Prior - von vage bis informativ
		b ~ normal(0, inv_sqrt(0.001)); 
		c ~ normal(0, inv_sqrt(0.01));
		d ~ normal(0, inv_sqrt(0.1));
	  //	e ~ normal(0, 0); // Exception: normal_lpdf: Scale parameter is 0, but must be > 0!  (in 'model257057d87032_Test_Stan_NV' at line 27)
		f ~ normal(0, inv_sqrt(1));
		
		g ~ normal(0, inv_sqrt(0.00001));
		h ~ normal(0, inv_sqrt(0.000001));
		i ~ normal(0, 1e6);
}

