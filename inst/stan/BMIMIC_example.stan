data{
  //Data size
  int<lower=1> nitemWorked;
  int<lower=1> nstud;
  int<lower=1> nitem;
  
  // DIF Index
  int<lower=0, upper=1> unidif_idx [nitem];
  int<lower=0> unidifeffect_idx [nitem];
  int<lower=0, upper=1> nondif_idx [nitem];
  int<lower=0> nondifeffect_idx [nitem];

  // number of DIF
  int<lower=0> unidif_n;
  int<lower=0> nondif_n;
  
  // indices
  int<lower=1,upper=nstud> studIdx[nitemWorked];
  int<lower=1,upper=nitem> itemIdx[nitemWorked];

  // index for first item
  int<lower=0, upper=1> firstitem[nitem];

  // data
  int<lower=0,upper=1> group[nstud];
  int<lower=0,upper=1> response[nitemWorked];
}

parameters{
  // IRT model parameters
  vector[nstud] eta;               // individual latent traits
 
  real<lower=0> lambda_free[nitem]; // item slope
  real tau_free[nitem];             // item intercept 

  // Impact parameters
  real impact0;       // intercept in latent mean
  real impact;        // difference in latent mean
  real impactv0;      // intercept in latent variance
  real impactv;       // difference in latent variance
 
  // DIF effect parameters
  vector[unidif_n] difeffect ;
  vector[nondif_n] nondifeffect; 
  
}

transformed parameters {
  real<lower=0> lambda[nitem];
  real tau[nitem];
 
  // Factor loading constraints
  for(jj in 1:nitem) {
    if(firstitem[jj] == 1) { // first loading per factor constrained to 1.
       lambda[jj] = 1;
	   tau[jj] = 0;
    } else {
       lambda[jj] = lambda_free[jj];
	     tau[jj] = tau_free[jj];
    }      
  };
}

model{
  real linPred[nitemWorked];
  vector[nstud] muEta;
  vector[nstud] sigEta;

  for(i in 1:nstud){
    muEta[i] = impact0 + group[i]*impact; 
	sigEta[i] = sqrt(exp(impactv0 + group[i]*impactv));
  };

   eta ~ normal(muEta, sigEta);

//MIMIC
  for(j in 1:nitemWorked) {
    linPred[j] = tau[itemIdx[j]] + lambda[itemIdx[j]] * eta[studIdx[j]];
	  if( unidif_idx[itemIdx[j]] == 1 && nondif_idx[itemIdx[j]] == 0 ) {
	    // Uniform DIF
	    linPred[j] = linPred[j] + group[studIdx[j]]*difeffect[itemIdx[j]-unidifeffect_idx[itemIdx[j]]];
		
	  } else if(unidif_idx[itemIdx[j]] == 0 && nondif_idx[itemIdx[j]] == 1 ) {
	    // non Uniform DIF
	    linPred[j] = linPred[j] + eta[studIdx[j]]*group[studIdx[j]]*nondifeffect[itemIdx[j]-nondifeffect_idx[itemIdx[j]]];
	  
	  } else if(unidif_idx[itemIdx[j]] == 1 && nondif_idx[itemIdx[j]] == 1 ) {
	    // uni + non
        linPred[j] = linPred[j] + group[studIdx[j]]*difeffect[itemIdx[j]-unidifeffect_idx[itemIdx[j]]] + eta[studIdx[j]]*group[studIdx[j]]*nondifeffect[itemIdx[j]-nondifeffect_idx[itemIdx[j]]];	
		
	  }
	  
	response[j] ~ bernoulli_logit(linPred[j]);
  }
  
//priors
  tau_free ~ normal(0, 1);
  for(i in 1:nitem) {
    lambda_free[i] ~ lognormal(0, 1);
  };

  impact0 ~ normal(0, 1);
  impact ~ normal(0, 1.44);
  impactv0 ~ normal(0, 0.1);
  impactv ~ normal(0, 1.44);
  
  difeffect ~ normal(0, 1);
  nondifeffect  ~ normal(0, 1);
}
// last line blank
