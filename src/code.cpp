#include "cpp11.hpp"
#include "cpp11/matrix.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include "utils.h"
#include "transform.h"
#include <vector>
using namespace cpp11;
namespace writable = cpp11::writable;


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
    double target_acceptance,
    writable::list misc,
    int n_rungs,
    doubles beta_init) {


  int iterations = burnin + samples;
  int n_par = theta_init.size();

  // Initialisise variables, ////////////////////////////////////////////////////
  double mh;
  bool mh_accept;
  double adjustment;
  int block;
  misc.push_back({"block"_nm = 0});

  // Initialise vector for theta
  std::vector<double> theta(n_par);
  for(int p = 0; p < n_par; ++p){
    theta[p] = theta_init[p];
  }
  // Initialise vector for proposal theta
  // Proposal theta is always the theta given to the likelihood function
  // and therefore must be a named vectord, which is why it is writebale::doubles
  writable::doubles theta_prop(n_par);
  for(int p = 0; p < n_par; ++p){
    theta_prop[p] = theta_init[p];
  }
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
    misc["block"] = as_sexp(b);
    ll[b] = ll_f(theta_prop, data, misc);
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
    out(0, p) =  theta[p];
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

  // PT ////////////////////////////////////////////////////////////////////////
  int rung;
  double rung_beta;
  double rung_proposal_sd;

  // Index of rungs, 0 = cold rung, n_rungs - 1 = hot rung
  std::vector<int> rung_index(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    rung_index[r] = r;
  }
  // Betas for each rung
  std:vector<double> beta(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    beta[r] = beta_init[r];
  }
  //////////////////////////////////////////////////////////////////////////////


  // Run ///////////////////////////////////////////////////////////////////////
  for(int i = 1; i < iterations; ++i){
    for(int r = 0; r < n_rungs; ++r){
      rung_beta = beta[rung_index[r]];


      for(int p = 0; p < n_par; ++p){
        theta_prop[p][r] = theta[p][r];
        phi_prop[p][r] = phi[p][r];
      }
      for(int b = 0; b < n_unique_blocks; ++b){
        ll_prop[b][r] = ll[b][r];
      }
      lp_prop[r] = lp[r];

      for(int p = 0; p < n_par; ++p){

        block = blocks[p] - 1;
        misc["block"] = as_sexp(block);
        // Propose new value
        phi_prop[p][r] = Rf_rnorm(phi[p][r], proposal_sd[p, rung_index[r]]);
        theta_prop[p][r] = phi_to_theta(phi_prop[p][r], transform_type[p], theta_min[p], theta_max[p]);
        ll_prop[block][r] = ll_f(theta_prop[,r], data, misc);
        lp_prop[r] = lp_f(theta_prop[,r]);

        // get parameter transformation adjustment
        adjustment = get_adjustment(theta[p][r], theta_prop[p][r], transform_type[p], theta_min[p], theta_max[p]);
        // calculate Metropolis-Hastings ratio
        mh = (sum(ll_prop[,r]) - sum(ll[,r])) + (lp_prop[r] - lp[r]) + adjustment;
        // accept or reject move
        mh_accept = log(Rf_runif(0, 1)) < mh;
        if(mh_accept){
          // Robbins monroe step
          if(i <= burnin){
            proposal_sd[p, rung_index[r]] = exp(log(proposal_sd[p, rung_index[r]]) + (1 - target_acceptance) / sqrt(i));
          }
          acceptance[p][r] = acceptance[p][r] + 1;
        } else {
          theta_prop[p][r] =  theta[p][r];
          phi_prop[p][r] = phi[p][r];
          ll_prop[block][r] = ll[block][r];
          lp_prop[r] = lp[r];
          // Robbins monroe step
          if(i <= burnin){
            proposal_sd[p, rung_index[r]] = exp(log(proposal_sd[p, rung_index[r]]) - target_acceptance / sqrt(i));
          }
        }
      }

      for(int p = 0; p < n_par; ++p){
        theta[p][r] = theta_prop[p][r];
        phi[p][r] = phi_prop[p][r];
      }
      for(int b = 0; b < n_unique_blocks; ++b){
        ll[b][r] = ll_prop[b][r];
      }
      lp[r] = lp_prop[r];

      // Only store values for the cold chain
      if(rung_index[r] == 0){
        // Record log likelihood
        log_likelihood[i] = sum(ll[,r]);
        log_prior[i] = lp[r];
        // Record parameters
        for(int p = 0; p < n_par; ++p){
          out(i, p) =  theta[p][r];
        }
      }
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
