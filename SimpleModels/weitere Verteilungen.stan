// Test Normalverteilung


data {
  // int l; // posterior für Binomialverteilung darf keine reelle Zahl (real) und somit auch kein Parameter sein
  // int m; 
  // int n; 
  int o; // posterior für Poissonverteilung darf keine reelle Zahl (real) und somit auch kein Parameter sein
  int p; 
  int q; 
}

parameters {
  real a; 
  real b; 
  real c; 
  real<lower=0> d; 
  real<lower=0> e; 
  real<lower=0> f; 
  real<lower=0> g; 
  real<lower=0> h;
  real<lower=0> i; 
  real<lower=0> j; 
  real<lower=0> k; 


}

model{
	// Uniform
	a ~ uniform(0,5);  // Referenz
	// variierte Prior
	b ~ uniform(0,1);  
	c ~ uniform(0,10);  

  // Gamma
	d ~ gamma(.001,.001); // Referenz; shape parameter (r)
	// variierte Prior
	e ~ gamma(.01,.001); 
	f ~ gamma(.1,.001); 
	g ~ gamma(.001,.01); 
  h ~ gamma(.001,.1); 
	i ~ gamma(.01,.01);
	j ~ gamma(.1,.1);
	k ~ gamma(1,1);

	// Binomial
  //	l ~ binomial(1,1); // die wie bei Winbugs verwendeten Werte für theta können nicht verwendet werden => Log probability evaluates to log(0), i.e. negative infinity.
  //	l ~ binomial(1,0); 
  //	m ~ binomial(2,1); 
  //	m ~ binomial(2,0); 
  //	n ~ binomial(3,1); 
	
	// Poisson
	o ~ poisson(1); 
	p ~ poisson(2); 
	q ~ poisson(3); 
}

