#include "cpp11.hpp"
#include "cpp11/doubles.hpp"
#include "Rmath.h"
#include <vector>
using namespace cpp11;
namespace writable = cpp11::writable;

[[cpp11::register]]
double cpp11_ll(doubles params, list data, list misc){
  // extract parameters
  double mu = params["mu"];
  double sigma = params["sigma"];

  // unpack data
  doubles x = data["x"];

  // sum log-likelihood over all data
  double ret = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    ret += Rf_dnorm4(x[i], mu, sigma, 1);
  }

  // return as SEXP
  return(ret);
}

