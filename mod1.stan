// generated with brms 2.14.4
functions {
}
data {
  int<lower=1> N;  // total number of observations
  int<lower=1> N_x;  // number of observations
  vector[N_x] Y_x;  // response variable
  int<lower=0> Nmi_x;  // number of missings
  int<lower=1> Jmi_x[Nmi_x];  // positions of missings
  int<lower=1> N_y;  // number of observations
  vector[N_y] Y_y;  // response variable
  int<lower=0> Nmi_y;  // number of missings
  int<lower=1> Jmi_y[Nmi_y];  // positions of missings
  int<lower=1> Ksp_y;  // number of special effects terms
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector[Nmi_x] Ymi_x;  // estimated missings
  real Intercept_x;  // temporary intercept for centered predictors
  real<lower=0> sigma_x;  // residual SD
  vector[Nmi_y] Ymi_y;  // estimated missings
  real Intercept_y;  // temporary intercept for centered predictors
  vector<lower=0>[Ksp_y] bsp_y;  // special effects coefficients
  real<lower=0> sigma_y;  // residual SD
}
transformed parameters {
}
model {
  // likelihood including all constants
  if (!prior_only) {
    // vector combining observed and missing responses
    vector[N_x] Yl_x = Y_x;
    // vector combining observed and missing responses
    vector[N_y] Yl_y = Y_y;
    // initialize linear predictor term
    vector[N_x] mu_x = Intercept_x + rep_vector(0.0, N_x);
    // initialize linear predictor term
    vector[N_y] mu_y = Intercept_y + rep_vector(0.0, N_y);
    Yl_x[Jmi_x] = Ymi_x;
    Yl_y[Jmi_y] = Ymi_y;
    for (n in 1:N_y) {
      // add more terms to the linear predictor
      mu_y[n] += (bsp_y[1]) * Yl_x[n];
    }
    target += normal_lpdf(Yl_x | mu_x, sigma_x);
    target += normal_lpdf(Yl_y | mu_y, sigma_y);
  }
  // priors including all constants
  target += student_t_lpdf(Intercept_x | 4,0,1);
  target += normal_lpdf(sigma_x | 0.1,1.0e-5)
    - 1 * normal_lccdf(0 | 0.1,1.0e-5);
  target += student_t_lpdf(Intercept_y | 4,0,1);
  target += lognormal_lpdf(bsp_y | 0.0,1.0)
    - 1 * lognormal_lccdf(0 | 0.0,1.0);
  target += gamma_lpdf(sigma_y |  4, 40 );
}
generated quantities {
  // actual population-level intercept
  real b_x_Intercept = Intercept_x;
  // actual population-level intercept
  real b_y_Intercept = Intercept_y;
}

