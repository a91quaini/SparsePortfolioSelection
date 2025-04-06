// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// compute_sr
double compute_sr(const arma::vec& weights, const arma::vec& mu, const arma::mat& sigma, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_compute_sr(SEXP weightsSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sr(weights, mu, sigma, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// compute_mve_sr
double compute_mve_sr(const arma::vec& mu, const arma::mat& sigma, const arma::uvec& selection, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_compute_mve_sr(SEXP muSEXP, SEXP sigmaSEXP, SEXP selectionSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_mve_sr(mu, sigma, selection, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// compute_mve_weights
arma::vec compute_mve_weights(const arma::vec& mu, const arma::mat& second_moment, const arma::uvec& selection, const double gamma, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_compute_mve_weights(SEXP muSEXP, SEXP second_momentSEXP, SEXP selectionSEXP, SEXP gammaSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type second_moment(second_momentSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type selection(selectionSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_mve_weights(mu, second_moment, selection, gamma, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// compute_sparse_mve_sr
Rcpp::List compute_sparse_mve_sr(const arma::vec& mu, const arma::mat& sigma, unsigned int max_card, const double greedy_perc, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_compute_sparse_mve_sr(SEXP muSEXP, SEXP sigmaSEXP, SEXP max_cardSEXP, SEXP greedy_percSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_card(max_cardSEXP);
    Rcpp::traits::input_parameter< const double >::type greedy_perc(greedy_percSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sparse_mve_sr(mu, sigma, max_card, greedy_perc, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// compute_sr_sparsity_loss
Rcpp::List compute_sr_sparsity_loss(const double mve_sr, const arma::vec& mu, const arma::mat& sigma, const arma::vec& mu_sample, const arma::mat& sigma_sample, unsigned int max_card, const double greedy_perc, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_compute_sr_sparsity_loss(SEXP mve_srSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP mu_sampleSEXP, SEXP sigma_sampleSEXP, SEXP max_cardSEXP, SEXP greedy_percSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type mve_sr(mve_srSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu_sample(mu_sampleSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma_sample(sigma_sampleSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type max_card(max_cardSEXP);
    Rcpp::traits::input_parameter< const double >::type greedy_perc(greedy_percSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_sr_sparsity_loss(mve_sr, mu, sigma, mu_sample, sigma_sample, max_card, greedy_perc, do_checks));
    return rcpp_result_gen;
END_RCPP
}
// simulate_sr_loss
Rcpp::List simulate_sr_loss(const double mve_sr, const arma::vec& mu, const arma::mat& sigma, const unsigned int n_obs, const unsigned int max_card, const double greedy_perc, const bool do_checks);
RcppExport SEXP _SparsePortfolioSelection_simulate_sr_loss(SEXP mve_srSEXP, SEXP muSEXP, SEXP sigmaSEXP, SEXP n_obsSEXP, SEXP max_cardSEXP, SEXP greedy_percSEXP, SEXP do_checksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type mve_sr(mve_srSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_obs(n_obsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type max_card(max_cardSEXP);
    Rcpp::traits::input_parameter< const double >::type greedy_perc(greedy_percSEXP);
    Rcpp::traits::input_parameter< const bool >::type do_checks(do_checksSEXP);
    rcpp_result_gen = Rcpp::wrap(simulate_sr_loss(mve_sr, mu, sigma, n_obs, max_card, greedy_perc, do_checks));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SparsePortfolioSelection_compute_sr", (DL_FUNC) &_SparsePortfolioSelection_compute_sr, 4},
    {"_SparsePortfolioSelection_compute_mve_sr", (DL_FUNC) &_SparsePortfolioSelection_compute_mve_sr, 4},
    {"_SparsePortfolioSelection_compute_mve_weights", (DL_FUNC) &_SparsePortfolioSelection_compute_mve_weights, 5},
    {"_SparsePortfolioSelection_compute_sparse_mve_sr", (DL_FUNC) &_SparsePortfolioSelection_compute_sparse_mve_sr, 5},
    {"_SparsePortfolioSelection_compute_sr_sparsity_loss", (DL_FUNC) &_SparsePortfolioSelection_compute_sr_sparsity_loss, 8},
    {"_SparsePortfolioSelection_simulate_sr_loss", (DL_FUNC) &_SparsePortfolioSelection_simulate_sr_loss, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SparsePortfolioSelection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
