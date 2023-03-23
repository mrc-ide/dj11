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

  // Initialise matrix for theta
  double theta[n_par][n_rungs];
  for(int i = 0; i < n_par; ++i){
    for(int j = 0; j < n_rungs; ++j){
      theta[i][j] = theta_init[i];
    }
  }
  // Initialise vector for proposal theta
  // Proposal theta is always the theta given to the likelihood function
  // and therefore must be a named vectord, which is why it is writeable::doubles
  writable::doubles theta_prop(n_par);
  for(int p = 0; p < n_par; ++p){
    theta_prop[p] = theta_init[p];
  }
  theta_prop.names() = theta_names;

  // Initialise value for transformed theta: phi
  std::vector<double> phi(n_par);
  for(int p = 0; p < n_par; ++p){
    phi[p] = theta_to_phi(theta[p][0], transform_type[p], theta_min[p], theta_max[p]);
  }
  // Initialise vector for proposal phi
  std::vector<double> phi_prop(n_par);

  // Initialise vector to store blocked log likelihood
  double block_ll[n_unique_blocks][n_rungs];
  for(int i = 0; i < n_unique_blocks; ++i){
    misc["block"] = as_sexp(i);
    for(int j = 0; j < n_rungs; ++j){
      block_ll[i][j] = ll_f(theta_prop, data, misc);
    }
  }

  // Initialise vector to store proposal blocked log likelihood
  std::vector<double> block_ll_prop(n_unique_blocks);
  // Initialise vector to store rung log likelihood (summed over blocks)
  std::vector<double> ll(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    double sum_ll = 0;
    for(int b = 0; b < n_unique_blocks; ++b){
      sum_ll += block_ll[b][r];
    }
    ll[r] = sum_ll;
  }

  // Initialise log prior vector
  std::vector<double> lp(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    lp[r] = lp_f(theta_prop);
  }
  // Initialise proposal log prior
  double lp_prop;
  //////////////////////////////////////////////////////////////////////////////

  // Outputs ///////////////////////////////////////////////////////////////////
  // Initialise log_likelihood output vector
  std::vector<double> out_log_likelihood(iterations);
  out_log_likelihood[0] = ll[0];
  // Initialise log_prior output vector
  std::vector<double> out_log_prior(iterations);
  out_log_prior[0] = lp[0];

  // Initialise output matrix
  writable::doubles_matrix<> out_theta(iterations, n_par);
  for(int p = 0; p < n_par; ++p){
    out_theta(0, p) =  theta[p][0];
  }
  //////////////////////////////////////////////////////////////////////////////

  // Tuning ////////////////////////////////////////////////////////////////////
  // Initialise matrix for proposal sd
  double proposal_sd[n_par][n_rungs];
  for(int i = 0; i < n_par; ++i){
    for(int j = 0; j < n_rungs; ++j){
      proposal_sd[i][j] = 0.1;
    }
  }

  // Initialise acceptance count matrix
  int acceptance[n_par][n_rungs];
  for(int i = 0; i < n_par; ++i){
    for(int j = 0; j < n_rungs; ++j){
      acceptance[i][j] = 0;
    }
  }

  // Initialise swap acceptance count vector
  std::vector<int> swap_acceptance(n_rungs - 1);
  for(int i = 0; i < (n_rungs - 1); ++i){
    swap_acceptance[i] = 0;
  }
  //////////////////////////////////////////////////////////////////////////////

  // PT ////////////////////////////////////////////////////////////////////////
  int rung;
  double rung_beta;
  double rung_proposal_sd;
  int index;

  // Index of rungs, 0 = cold rung, n_rungs - 1 = hot rung
  std::vector<int> rung_index(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    rung_index[r] = r;
  }
  // Betas for each rung
  std::vector<double> beta(n_rungs);
  for(int r = 0; r < n_rungs; ++r){
    beta[r] = beta_init[r];
  }
  //////////////////////////////////////////////////////////////////////////////


  // Run ///////////////////////////////////////////////////////////////////////
  for(int i = 1; i < iterations; ++i){
    for(int r = 0; r < n_rungs; ++r){
      rung_beta = beta[r];
      index = rung_index[r];

      // Copy rung theta and phi
      for(int p = 0; p < n_par; ++p){
        theta_prop[p] = theta[p][index];
        phi_prop[p] = phi[p];
      }
      for(int b = 0; b < n_unique_blocks; ++b){
        // TODO I think block_ll also needs to be indexed like theta
        block_ll_prop[b] = block_ll[b][index];
      }
      lp_prop = lp[r];

      for(int p = 0; p < n_par; ++p){
        // Set block for parameter
        block = blocks[p] - 1;
        misc["block"] = as_sexp(block);
        // Propose new value
        phi_prop[p] = Rf_rnorm(phi[p], proposal_sd[p][r]);
        theta_prop[p] = phi_to_theta(phi_prop[p], transform_type[p], theta_min[p], theta_max[p]);
        block_ll_prop[block] = ll_f(theta_prop, data, misc);
        lp_prop = lp_f(theta_prop);

        // get parameter transformation adjustment
        adjustment = get_adjustment(theta[p][index], theta_prop[p], transform_type[p], theta_min[p], theta_max[p]);
        // calculate Metropolis-Hastings ratio
        mh = rung_beta * (sum(block_ll_prop) - ll[r]) + (lp_prop - lp[r]) + adjustment;
        // accept or reject move
        mh_accept = log(Rf_runif(0, 1)) < mh;
        if(mh_accept){
          // Robbins monroe step
          if(i <= burnin){
            proposal_sd[p][r] = exp(log(proposal_sd[p][r]) + (1 - target_acceptance) / sqrt(i));
          }
          acceptance[p][r] = acceptance[p][r] + 1;
        } else {
          theta_prop[p] =  theta[p][index];
          phi_prop[p] = phi[p];
          block_ll_prop[block] = block_ll[block][index];
          lp_prop = lp[r];
          // Robbins monroe step
          if(i <= burnin){
            proposal_sd[p][r] = exp(log(proposal_sd[p][r]) - target_acceptance / sqrt(i));
          }
        }
      }

      for(int p = 0; p < n_par; ++p){
        theta[p][index] = theta_prop[p];
        phi[p] = phi_prop[p];
      }
      for(int b = 0; b < n_unique_blocks; ++b){
        block_ll[b][index] = block_ll_prop[b];
      }
      lp[r] = lp_prop;
      ll[r] = sum(block_ll_prop);

      // Only store values for the cold chain
      if(r == 0){
        // Record log likelihood
        out_log_likelihood[i] = ll[r];
        out_log_prior[i] = lp[r];
        // Record parameters
        for(int p = 0; p < n_par; ++p){
          out_theta(i, p) =  theta[p][index];
        }
      }
    }
  }

  // Return outputs in a list
  return writable::list({
    "log_likelihood"_nm = out_log_likelihood,
      "log_prior"_nm = out_log_prior,
      "out"_nm = out_theta,
      "proposal_sd"_nm = proposal_sd[0][0],
                                       "acceptance"_nm = acceptance[0][0]
  });
}
