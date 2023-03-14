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

double sum(std::vector<double> x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}


[[cpp11::register]]
list mcmc(
    doubles theta_init,
    strings theta_names,
    integers transform_type,
    doubles theta_min,
    doubles theta_max,
    integers blocks,
    int n_unique_blocks,
    list data,
    int burnin,
    int samples,
    function ll_f,
    function lp_f,
    double target_acceptance) {


  int iterations = burnin + samples;
  int n_par = theta_init.size();

  // Initialisise variables, ////////////////////////////////////////////////////
  double mh;
  bool mh_accept;
  double adjustment;
  int block;

  // Initialise vector for theta
  writable::doubles theta(n_par);
  for(int p = 0; p < n_par; ++p){
    theta[p] = theta_init[p];
  }
  theta.names() = theta_names;
  // Initialise vector for proposal theta
  writable::doubles theta_prop(n_par);
  theta_prop.names() = theta_names;

  // Initialise value for transformed theta: phi
  std::vector<double> phi(n_par);
  for(int p = 0; p < n_par; ++p){
    phi[p] = theta_to_phi(theta[p], transform_type[p], theta_min[p], theta_max[p]);
  }
  // Initialise vector for proposal phi
  std::vector<double> phi_prop(n_par);

  // Initialise vector to store blocked log likelihood
  std::vector<double> ll(n_unique_blocks);
  for(int b = 0; b < n_unique_blocks; ++b){
    ll[b] = ll_f(theta, data, b);
  }
  // Initialise vector to store proposal blocked log likelihood
  std::vector<double> ll_prop(n_unique_blocks);
  // Initialise log prior
  double lp = lp_f(theta);
  // Initialise proposal log prior
  double lp_prop;
  //////////////////////////////////////////////////////////////////////////////

  // Outputs ///////////////////////////////////////////////////////////////////
  // Initialise log_likelihood output vector
  std::vector<double> log_likelihood(iterations);
  log_likelihood[0] = sum(ll);
  // Initialise log_prior output vector
  std::vector<double> log_prior(iterations);
  log_prior[0] = lp;

  // Initialise output matrix
  writable::doubles_matrix<> out(iterations, n_par);
  for(int p = 0; p < n_par; ++p){
    out(0, p) =  static_cast<double>(theta[p]);
  }
  //////////////////////////////////////////////////////////////////////////////

  // Tuning ////////////////////////////////////////////////////////////////////
  // Initialise proposal sd
  std::vector<double> proposal_sd(n_par);
  for(int p = 0; p < n_par; ++p){
    proposal_sd[p] = 0.1;
  }

  // Initialise acceptance count vector
  std::vector<double> acceptance(n_par);
  for(int p = 0; p < n_par; ++p){
    acceptance[p] = 0;
  }
  //////////////////////////////////////////////////////////////////////////////


  // Run ///////////////////////////////////////////////////////////////////////
  for(int i = 1; i < iterations; ++i){
    for(int p = 0; p < n_par; ++p){
      theta_prop[p] = static_cast<double>(theta[p]);
    }
    // TODO: Why is phi_prop not copied as with theta?
    phi_prop = phi;
    ll_prop = ll;
    lp_prop = lp;
    for(int p = 0; p < n_par; ++p){
      block = blocks[p] - 1;
      // Propose new value
      phi_prop[p] = Rf_rnorm(phi[p], proposal_sd[p]);
      theta_prop[p] = phi_to_theta(phi_prop[p], transform_type[p], theta_min[p], theta_max[p]);
      ll_prop[block] = ll_f(theta_prop, data, block);
      lp_prop = lp_f(theta_prop);

      // get parameter transformation adjustment
      adjustment = get_adjustment(theta[p], theta_prop[p], transform_type[p], theta_min[p], theta_max[p]);
      // calculate Metropolis-Hastings ratio
      mh = (sum(ll_prop) - sum(ll)) + (lp_prop - lp) + adjustment;
      // accept or reject move
      mh_accept = log(Rf_runif(0, 1)) < mh;
      if(mh_accept){
        // Robbins monroe step
        if(i <= burnin){
          proposal_sd[p] = exp(log(proposal_sd[p]) + (1 - target_acceptance) / sqrt(i));
        }
        acceptance[p] = acceptance[p] + 1;
      } else {
        theta_prop[p] =  static_cast<double>(theta[p]);
        phi_prop[p] = phi[p];
        ll_prop[block] = ll[block];
        lp_prop = lp;
        // Robbins monroe step
        if(i <= burnin){
          proposal_sd[p] = exp(log(proposal_sd[p]) - target_acceptance / sqrt(i));
        }
      }
    }
    for(int p = 0; p < n_par; ++p){
      theta[p] = static_cast<double>(theta_prop[p]);
    }
    phi = phi_prop;
    ll = ll_prop;
    lp = lp_prop;
    // Record log likelihood
    log_likelihood[i] = sum(ll);
    log_prior[i] = lp;
    // Record parameters
    for(int p = 0; p < n_par; ++p){
      out(i, p) =  static_cast<double>(theta[p]);
    }
  }

  // Return outputs in a list
  return writable::list({
    "log_likelihood"_nm = log_likelihood,
      "log_prior"_nm = log_prior,
      "out"_nm = out,
      "proposal_sd"_nm = proposal_sd,
      "acceptance"_nm = acceptance
  });
}
