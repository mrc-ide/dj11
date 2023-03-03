#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include <vector>
#include <iostream>
using namespace cpp11;
namespace writable = cpp11::writable;

double phi_to_theta(double phi, int transformation_type, double theta_min, double theta_max) {

  double theta;

  if(transformation_type == 0){
    theta = phi;
  }
  if(transformation_type == 1){
    theta = theta_max - exp(phi);
  }
  if(transformation_type == 2){
    theta = exp(phi) + theta_min;
  }
  if(transformation_type == 3){
    theta = (theta_max * exp(phi) + theta_min) / (1 + exp(phi));
  }

  return(theta);
}

double theta_to_phi(double theta, int transformation_type, double theta_min, double theta_max) {

  double phi;

  if(transformation_type == 0){
    phi = theta;
  }
  if(transformation_type == 1){
    phi = log(theta_max - theta);
  }
  if(transformation_type == 2){
    phi = log(theta - theta_min);
  }
  if(transformation_type == 3){
    phi = log(theta - theta_min) - log(theta_max - theta);
  }

  return(phi);
}

double get_adjustment(double theta, double theta_prop, int transformation_type, double theta_min, double theta_max) {

  double adjustment;

  if(transformation_type == 0){
    adjustment = 0.0;
  }
  if(transformation_type == 1){
    adjustment = log(theta_max - theta_prop) - log(theta_max - theta);
  }
  if(transformation_type == 2){
    adjustment = log(theta_prop - theta_min) - log(theta - theta_min);
  }
  if(transformation_type == 3){
    adjustment = log(theta_max - theta_prop) + log(theta_prop - theta_min) - log(theta_max - theta) - log(theta - theta_min);

  }
  return(adjustment);
}


[[cpp11::register]]
list mcmc(
    doubles theta_init,
    integers transform_type,
    doubles theta_min,
    doubles theta_max,
    doubles data,
    int iterations,
    function ll,
    function lp) {

  int n_par = theta_init.size();
  double ll_prop;
  double lp_prop;
  double mh;
  bool mh_accept;
  double adjustment;

  // Initialise log_likelihood vector
  std::vector<double> theta(n_par);
  for(int p = 0; p < n_par; ++p){
    theta[p] = theta_init[p];
    //std::cout << theta[p];
  }
  // Initialise vector for proposal theta
  std::vector<double> theta_prop(n_par);
  // Initialise value for transformed theta: phi
  std::vector<double> phi(n_par);
  for(int p = 0; p < n_par; ++p){
    phi[p] = theta_to_phi(theta[p], transform_type[p], theta_min[p], theta_max[p]);
  }

  // Initialise vector for proposal phi
  std::vector<double> phi_prop(n_par);

  // Initialise log_likelihood vector
  std::vector<double> log_likelihood(iterations);
  log_likelihood[0] = ll(theta, data);
  // Initialise log_prio vector
  std::vector<double> log_prior(iterations);
  log_prior[0] = lp(theta);

  // Initialise output matrix
  writable::doubles_matrix<> out(iterations, n_par);
  for(int p = 0; p < n_par; ++p){
    out(0, p) = theta[p];
  }

  // // Initialise proposal sd
  std::vector<double> proposal_sd(n_par);
  for(int p = 0; p < n_par; ++p){
    proposal_sd[p] = 0.1;
  }
  // Initialise acceptance count vector
  std::vector<double> acceptance(n_par);
  for(int p = 0; p < n_par; ++p){
    acceptance[p] = 0;
  }

  // Initialise vector for proposal theta
  std::vector<double> current(n_par);
  for(int p = 0; p < n_par; ++p){
    current[p] = theta[p];
  }

  for(int i = 1; i < iterations; ++i){
    theta_prop = theta;
    phi_prop = phi;
    for(int p = 0; p < n_par; ++p){
      // Propose new value
      phi_prop[p] = Rf_rnorm(phi[p], proposal_sd[p]);
      theta_prop[p] = phi_to_theta(phi_prop[p], transform_type[p], theta_min[p], theta_max[p]);
      ll_prop = ll(theta_prop, data);
      lp_prop = lp(theta_prop);
      // calculate Metropolis-Hastings ratio
      adjustment = get_adjustment(theta[p], theta_prop[p], transform_type[p], theta_min[p], theta_max[p]);
      mh = (ll_prop - log_likelihood[i - 1]) + (lp_prop - log_prior[i - 1]) + adjustment;
      // accept or reject move
      mh_accept = log(Rf_runif(0, 1)) < mh;
      if(mh_accept){
        proposal_sd[p] = exp(log(proposal_sd[p]) + (1 - 0.24) / sqrt(i));
        acceptance[p] = acceptance[p] + 1;
      } else {
        theta_prop[p] = theta[p];
        phi_prop[p] = phi[p];
        proposal_sd[p] = exp(log(proposal_sd[p]) - 0.24 / sqrt(i));
      }
    }
    theta = theta_prop;
    phi = phi_prop;
    // Record log likelihood
    log_likelihood[i] = ll(theta, data);
    log_prior[i] = lp(theta);
    // Record parameters
    for(int p = 0; p < n_par; ++p){
      out(i, p) = theta[p];
    }

  }

  return writable::list({
    "log_likelihood"_nm = log_likelihood,
      "log_prior"_nm = log_prior,
      "out"_nm = out,
      "proposal_sd"_nm = proposal_sd,
      "acceptance"_nm = acceptance
  });
}
