// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pairs_freq
Eigen::MatrixXd pairs_freq(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC);
RcppExport SEXP _plCFA_pairs_freq(SEXP YSEXP, SEXP C_VECSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    rcpp_result_gen = Rcpp::wrap(pairs_freq(Y, C_VEC));
    return rcpp_result_gen;
END_RCPP
}
// get_Lam
Eigen::MatrixXd get_Lam(Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> c_vec, const Eigen::VectorXd& theta);
RcppExport SEXP _plCFA_get_Lam(SEXP ASEXP, SEXP c_vecSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type c_vec(c_vecSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_Lam(A, c_vec, theta));
    return rcpp_result_gen;
END_RCPP
}
// get_Sigma_u
Eigen::MatrixXd get_Sigma_u(Eigen::Map<Eigen::MatrixXd> A, const Eigen::VectorXd& theta);
RcppExport SEXP _plCFA_get_Sigma_u(SEXP ASEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(get_Sigma_u(A, theta));
    return rcpp_result_gen;
END_RCPP
}
// matrixprod
Eigen::MatrixXd matrixprod(Eigen::MatrixXd A, Eigen::MatrixXd B);
RcppExport SEXP _plCFA_matrixprod(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixprod(A, B));
    return rcpp_result_gen;
END_RCPP
}
// multiThread_completePairwise
Rcpp::List multiThread_completePairwise(Eigen::Map<Eigen::MatrixXd> Y, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::VectorXd> TAU, Eigen::Map<Eigen::VectorXd> LAMBDA, Eigen::Map<Eigen::VectorXd> TRANSFORMED_RHOS, Eigen::Map<Eigen::MatrixXd> FREQ, int CORRFLAG, int GRFLAG, int SILENTFLAG);
RcppExport SEXP _plCFA_multiThread_completePairwise(SEXP YSEXP, SEXP C_VECSEXP, SEXP ASEXP, SEXP TAUSEXP, SEXP LAMBDASEXP, SEXP TRANSFORMED_RHOSSEXP, SEXP FREQSEXP, SEXP CORRFLAGSEXP, SEXP GRFLAGSEXP, SEXP SILENTFLAGSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type Y(YSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TAU(TAUSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type LAMBDA(LAMBDASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TRANSFORMED_RHOS(TRANSFORMED_RHOSSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type FREQ(FREQSEXP);
    Rcpp::traits::input_parameter< int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< int >::type GRFLAG(GRFLAGSEXP);
    Rcpp::traits::input_parameter< int >::type SILENTFLAG(SILENTFLAGSEXP);
    rcpp_result_gen = Rcpp::wrap(multiThread_completePairwise(Y, C_VEC, A, TAU, LAMBDA, TRANSFORMED_RHOS, FREQ, CORRFLAG, GRFLAG, SILENTFLAG));
    return rcpp_result_gen;
END_RCPP
}
// plCFA
Rcpp::List plCFA(Eigen::Map<Eigen::MatrixXd> DATA, Eigen::Map<Eigen::VectorXd> C_VEC, Eigen::Map<Eigen::MatrixXd> CONSTRMAT, Eigen::Map<Eigen::VectorXd> TAU, Eigen::Map<Eigen::VectorXd> LAMBDA, Eigen::Map<Eigen::VectorXd> TRANSFORMED_RHOS, const unsigned int CORRFLAG, const unsigned int METHODFLAG, const unsigned int PAIRS_PER_ITERATION, const double PROB, double ETA, const unsigned int BURN, const unsigned int MAXT, const double TOLGRAD, const double TOLPAR, const double TOLOBJ, const unsigned int TOLCOUNT, const unsigned int SEED, const unsigned int SILENTFLAG, const bool CHECKCONVERGENCE, const double TOL, const int TOLCOUNTER);
RcppExport SEXP _plCFA_plCFA(SEXP DATASEXP, SEXP C_VECSEXP, SEXP CONSTRMATSEXP, SEXP TAUSEXP, SEXP LAMBDASEXP, SEXP TRANSFORMED_RHOSSEXP, SEXP CORRFLAGSEXP, SEXP METHODFLAGSEXP, SEXP PAIRS_PER_ITERATIONSEXP, SEXP PROBSEXP, SEXP ETASEXP, SEXP BURNSEXP, SEXP MAXTSEXP, SEXP TOLGRADSEXP, SEXP TOLPARSEXP, SEXP TOLOBJSEXP, SEXP TOLCOUNTSEXP, SEXP SEEDSEXP, SEXP SILENTFLAGSEXP, SEXP CHECKCONVERGENCESEXP, SEXP TOLSEXP, SEXP TOLCOUNTERSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type DATA(DATASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type C_VEC(C_VECSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type CONSTRMAT(CONSTRMATSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TAU(TAUSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type LAMBDA(LAMBDASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::VectorXd> >::type TRANSFORMED_RHOS(TRANSFORMED_RHOSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type CORRFLAG(CORRFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type METHODFLAG(METHODFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type PAIRS_PER_ITERATION(PAIRS_PER_ITERATIONSEXP);
    Rcpp::traits::input_parameter< const double >::type PROB(PROBSEXP);
    Rcpp::traits::input_parameter< double >::type ETA(ETASEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type BURN(BURNSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type MAXT(MAXTSEXP);
    Rcpp::traits::input_parameter< const double >::type TOLGRAD(TOLGRADSEXP);
    Rcpp::traits::input_parameter< const double >::type TOLPAR(TOLPARSEXP);
    Rcpp::traits::input_parameter< const double >::type TOLOBJ(TOLOBJSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type TOLCOUNT(TOLCOUNTSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SEED(SEEDSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SILENTFLAG(SILENTFLAGSEXP);
    Rcpp::traits::input_parameter< const bool >::type CHECKCONVERGENCE(CHECKCONVERGENCESEXP);
    Rcpp::traits::input_parameter< const double >::type TOL(TOLSEXP);
    Rcpp::traits::input_parameter< const int >::type TOLCOUNTER(TOLCOUNTERSEXP);
    rcpp_result_gen = Rcpp::wrap(plCFA(DATA, C_VEC, CONSTRMAT, TAU, LAMBDA, TRANSFORMED_RHOS, CORRFLAG, METHODFLAG, PAIRS_PER_ITERATION, PROB, ETA, BURN, MAXT, TOLGRAD, TOLPAR, TOLOBJ, TOLCOUNT, SEED, SILENTFLAG, CHECKCONVERGENCE, TOL, TOLCOUNTER));
    return rcpp_result_gen;
END_RCPP
}
// sampling_step
std::vector<int> sampling_step(std::vector<int>& FULL_POOL, const unsigned int METHODFLAG, const double PROB, const int PAIRS_PER_ITERATION, const unsigned int N_ITEMS, const unsigned int SEED, const int SILENTFLAG, const unsigned int ITER);
RcppExport SEXP _plCFA_sampling_step(SEXP FULL_POOLSEXP, SEXP METHODFLAGSEXP, SEXP PROBSEXP, SEXP PAIRS_PER_ITERATIONSEXP, SEXP N_ITEMSSEXP, SEXP SEEDSEXP, SEXP SILENTFLAGSEXP, SEXP ITERSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int>& >::type FULL_POOL(FULL_POOLSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type METHODFLAG(METHODFLAGSEXP);
    Rcpp::traits::input_parameter< const double >::type PROB(PROBSEXP);
    Rcpp::traits::input_parameter< const int >::type PAIRS_PER_ITERATION(PAIRS_PER_ITERATIONSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type N_ITEMS(N_ITEMSSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type SEED(SEEDSEXP);
    Rcpp::traits::input_parameter< const int >::type SILENTFLAG(SILENTFLAGSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type ITER(ITERSEXP);
    rcpp_result_gen = Rcpp::wrap(sampling_step(FULL_POOL, METHODFLAG, PROB, PAIRS_PER_ITERATION, N_ITEMS, SEED, SILENTFLAG, ITER));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_plCFA_pairs_freq", (DL_FUNC) &_plCFA_pairs_freq, 2},
    {"_plCFA_get_Lam", (DL_FUNC) &_plCFA_get_Lam, 3},
    {"_plCFA_get_Sigma_u", (DL_FUNC) &_plCFA_get_Sigma_u, 2},
    {"_plCFA_matrixprod", (DL_FUNC) &_plCFA_matrixprod, 2},
    {"_plCFA_multiThread_completePairwise", (DL_FUNC) &_plCFA_multiThread_completePairwise, 10},
    {"_plCFA_plCFA", (DL_FUNC) &_plCFA_plCFA, 22},
    {"_plCFA_sampling_step", (DL_FUNC) &_plCFA_sampling_step, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_plCFA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
