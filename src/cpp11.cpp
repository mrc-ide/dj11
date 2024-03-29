// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// code.cpp
list mcmc(doubles theta_init, strings theta_names, integers transform_type, doubles theta_min, doubles theta_max, list blocks_list, int n_unique_blocks, list data, int burnin, int samples, function ll_f, function lp_f, double target_acceptance, writable::list misc, int n_rungs, doubles beta_init, bool swap);
extern "C" SEXP _dj11_mcmc(SEXP theta_init, SEXP theta_names, SEXP transform_type, SEXP theta_min, SEXP theta_max, SEXP blocks_list, SEXP n_unique_blocks, SEXP data, SEXP burnin, SEXP samples, SEXP ll_f, SEXP lp_f, SEXP target_acceptance, SEXP misc, SEXP n_rungs, SEXP beta_init, SEXP swap) {
  BEGIN_CPP11
    return cpp11::as_sexp(mcmc(cpp11::as_cpp<cpp11::decay_t<doubles>>(theta_init), cpp11::as_cpp<cpp11::decay_t<strings>>(theta_names), cpp11::as_cpp<cpp11::decay_t<integers>>(transform_type), cpp11::as_cpp<cpp11::decay_t<doubles>>(theta_min), cpp11::as_cpp<cpp11::decay_t<doubles>>(theta_max), cpp11::as_cpp<cpp11::decay_t<list>>(blocks_list), cpp11::as_cpp<cpp11::decay_t<int>>(n_unique_blocks), cpp11::as_cpp<cpp11::decay_t<list>>(data), cpp11::as_cpp<cpp11::decay_t<int>>(burnin), cpp11::as_cpp<cpp11::decay_t<int>>(samples), cpp11::as_cpp<cpp11::decay_t<function>>(ll_f), cpp11::as_cpp<cpp11::decay_t<function>>(lp_f), cpp11::as_cpp<cpp11::decay_t<double>>(target_acceptance), cpp11::as_cpp<cpp11::decay_t<writable::list>>(misc), cpp11::as_cpp<cpp11::decay_t<int>>(n_rungs), cpp11::as_cpp<cpp11::decay_t<doubles>>(beta_init), cpp11::as_cpp<cpp11::decay_t<bool>>(swap)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_dj11_mcmc", (DL_FUNC) &_dj11_mcmc, 17},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_dj11(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
