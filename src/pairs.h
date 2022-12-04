#ifndef pairs_H
#define pairs_H


#include "bivariateNormal.h"
#include "genericUtils.h"
#include "biProbs.h"

//Individual pair of item contribution to
// 1. Log-likelihood
// 2. Gradient
// 3. Diagonal of Hessian
void pair_contribution (
        // Parameters
        Eigen::Map<Eigen::MatrixXd> A,
        Eigen::Map<Eigen::VectorXd> C_VEC,
        const Eigen::VectorXd &THETA,
        const int CORRFLAG,

        // Input:
        const unsigned int k,
        const unsigned int l,
        const Eigen::MatrixXd &PAIRS_TABLE,

        // Options:
        const unsigned int SILENTFLAG,
        const unsigned int DFLAG,
        const unsigned int GRADFLAG,

        // Outputs:
        double &ll,
        Eigen::VectorXd &gradient,
        Eigen::VectorXd &gradient2

){
    unsigned int p = A.rows();
    unsigned int q = A.cols();
    unsigned int d = THETA.size();
    unsigned int c = C_VEC.sum();
    unsigned int nthr = c-p;
    unsigned int ncorr = q*(q-1)/2;
    unsigned int nload = d-nthr-ncorr;

    double tmp_ll = 0;
    Eigen::VectorXd tmp_gradient(d); tmp_gradient.setZero();
    Eigen::VectorXd tmp_gradient2(d); tmp_gradient2.setZero();

    // rearrange parameters
    Eigen::MatrixXd Lam            = get_Lam(A, C_VEC, THETA);
    Eigen::MatrixXd Sigma_u        = get_Sigma_u(A, THETA);
    Eigen::VectorXd tau            = get_tau(THETA, C_VEC);
    Eigen::VectorXd transformed_rhos=THETA.segment(nthr+nload, d);

    // Identifies quantities related to pair (k,l)
    unsigned int ck = C_VEC(k);
    unsigned int cl = C_VEC(l);
    Eigen::VectorXd lambdak = Lam.row(k);
    Eigen::VectorXd lambdal = Lam.row(l);
    double rho_kl = lambdak.transpose() * Sigma_u * lambdal;
    Eigen::MatrixXd pairs_tab = PAIRS_TABLE;
    pairs_tab.conservativeResize(PAIRS_TABLE.rows() + 1, Eigen::NoChange_t() );
    // identify column index in freq table
    // i1: starting index item k
    unsigned int i1 = 0;
    if(k > 1){
        for(unsigned int u = 1; u < k; u++){
            unsigned int cu = C_VEC(u);
            //if(SILENTFLAG == 0)Rcpp::Rcout << "u: " << u << ", cu: "<< cu << "\n";
            i1 += cu * C_VEC.segment(0,u).sum();
        }
    }

    // i2 starting index from i1 for item l
    unsigned int i2 = 0;
    if(l > 0){
        i2 = C_VEC.segment(0,l).sum() * C_VEC(k);
    }

    if(DFLAG!=1){
        ////////////////////////////
        /* LIKELIHOOD COMPUTATION */
        ////////////////////////////
        for(unsigned int sk = 0; sk < ck; sk ++){

            // i3: starting index from i2 for cat sk
            unsigned int i3 = sk * cl;

            for(unsigned int sl = 0; sl < cl; sl ++){

                // final column index for pairs_tab. Print to check
                unsigned int r = i1 + i2 + i3 + sl;

                // read frequency
                unsigned int n_sksl = PAIRS_TABLE(4, r);

                // identify thresholds
                Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);

                // compute pi
                double pi_sksl = compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
                if(SILENTFLAG == 0)Rcpp::Rcout << "("<<k<<","<<l<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
                pairs_tab(5,r) = pi_sksl;

                // update ll
                tmp_ll += n_sksl * log(pi_sksl+1e-8);
            }
        }

        //////////////////////////
        /* GRADIENT COMPUTATION */
        /////////////////////////

        if(GRADFLAG == 1){
            unsigned int iter = 0;

            /////////////////////////////////////////////////////
            // (k,l)-pair likelihood derivative wrt thresholds //
            /////////////////////////////////////////////////////
            if(SILENTFLAG == 0)Rcpp::Rcout << "- Gradient wrt thresholds: \n";

            // loop: terate over elements of thresholds vector
            for(unsigned int s = 0; s < tau.size(); s++){
                double grs = 0; // temporary location for gradient related to s-th element of tau
                if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ gradient("<< s<< ")\n";

                // List three cases: 1. threshold related to item k, 2. threshold related to item l, 3. threshold non relevant to items couple (k,l)
                if(s >= (C_VEC.segment(0, k).sum()) - (k) & s < C_VEC.segment(0, k + 1).sum() - (k + 1)){
                    // [CASE 1]: threshold related to item k

                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item k:\n";
                    unsigned int sk = s - (C_VEC.segment(0, k).sum()) + (k);

                    // i3: starting index from i2 for cat sk and sk+1
                    unsigned int i3 = sk * cl;
                    unsigned int i3suc = (sk+1) * cl;
                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ sk: " << sk << ". Summing over categories item l: ";

                    // iterate over categories of item l
                    for(unsigned int sl = 0; sl < cl; sl ++){
                        if(SILENTFLAG == 0)Rcpp::Rcout << " ... cat" << sl ;

                        // identify pairs_tab column for (sk,sl) and (sk+1, sl)
                        unsigned int r = i1 + i2 + i3 + sl;
                        unsigned int rsuc = i1 + i2 + i3suc + sl;

                        // read frequences
                        unsigned int n_sksl = pairs_tab(4, r);
                        unsigned int n_sksucsl = pairs_tab(4, rsuc);

                        // read probabilities
                        double pi_sksl = pairs_tab(5, r);
                        double pi_sksucsl = pairs_tab(5, rsuc);

                        // identify tau_sk, tau_sl, tau_sl-1
                        Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
                        double t_sk = pi_thresholds(0); double t_sl = pi_thresholds(1); double t_sk_prev = pi_thresholds(2); double t_sl_prev = pi_thresholds(3);

                        // compute gradient
                        double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_sksucsl/(pi_sksucsl+1e-8)));
                        double tmp2 = R::dnorm(t_sk, 0, 1, 0);
                        double tmp3 = R::pnorm((t_sl-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
                        double tmp4 = R::pnorm((t_sl_prev-rho_kl*t_sk)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
                        grs += tmp1 * tmp2 * (tmp3 - tmp4);
                    }
                    if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

                }else if(s >= (C_VEC.segment(0, l).sum())-(l) & s<C_VEC.segment(0, l + 1).sum()-(l + 1)){
                    // [CASE 2]: threshold related to item l

                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_ tau item l\n";
                    unsigned int sl = s - (C_VEC.segment(0, l).sum()) + (l);

                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  sl: " << sl << ". Summing over categories item k: ";

                    // iterate over categories item k
                    for(unsigned int sk = 0; sk < ck; sk ++){

                        // i3: starting index from i2 for cat sk
                        unsigned int i3 = sk * cl;

                        // identify pairs_tab column for (sk,sl) and (sk, sl + 1)
                        unsigned int r = i1 + i2 + i3 + sl;
                        unsigned int rsuc = i1 + i2 + i3 + sl + 1;

                        // read frequences
                        unsigned int n_sksl = pairs_tab(4, r);
                        unsigned int n_skslsuc = pairs_tab(4, rsuc);

                        // read probabilities
                        double pi_sksl = pairs_tab(5, r);
                        double pi_skslsuc = pairs_tab(5, rsuc);

                        // identify tau_sk, tau_sl, tau_sl-1
                        Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
                        double t_sk = pi_thresholds(0); double t_sl = pi_thresholds(1); double t_sk_prev = pi_thresholds(2); double t_sl_prev = pi_thresholds(3);


                        if(SILENTFLAG == 0)Rcpp::Rcout<<"\n  |    |   |_ sk:"<< sk << ", r: "<< r<<", n_sksl:"
                                                      << n_sksl<< ", n_sksl+1:" << n_skslsuc << ", pi_sksl:"
                                                      << pi_sksl << ", pi_sksl+1:"<< pi_skslsuc << ", t_sk:"
                                                      << t_sk<< ", t_sl:" << t_sl << "t_sk-1:"<< t_sk_prev;
                        // compute gradient
                        double tmp1 = ((n_sksl/(pi_sksl+1e-8))-(n_skslsuc/(pi_skslsuc+1e-8)));
                        double tmp2 = R::dnorm(t_sl, 0, 1, 0);
                        double tmp3 = R::pnorm((t_sk-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
                        double tmp4 = R::pnorm((t_sk_prev-rho_kl*t_sl)/(pow(1-pow(rho_kl,2), .5)), 0, 1, 1, 0);
                        if(SILENTFLAG == 0)Rcpp::Rcout<<" => out" << sk << ":" << tmp1 * tmp2 * (tmp3 - tmp4);
                        grs += tmp1 * tmp2 * (tmp3 - tmp4);
                    }
                    if(SILENTFLAG == 0)Rcpp::Rcout << "\n";

                }else{
                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  tau of other item\n";
                }
                if(SILENTFLAG == 0)Rcpp::Rcout << "  |    |_  Thresholds:: " << iter<< "/"<< tau.size()
                                               <<". Tot:: "<< iter << "/"<< d <<". Val ="<< grs <<"\n";
                tmp_gradient(iter) += grs;
                iter ++;
            }
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ Done. \n";

            ///////////////////////////////////////////////////////////
            // (k,l)-pair likelihood derivative wrt URV correlation: //
            // intermediate step for derivatives wrt                 //
            // loadings and factor correlations                      //
            ///////////////////////////////////////////////////////////
            if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Intermediate derivative for loadings and correlation: \n";
            double tmp_kl = 0; // temporary location of the gradient

            // double loop: iterate over each combination of categories of items k and l
            for(unsigned int sk = 0; sk < ck; sk ++){
                for(unsigned int sl = 0; sl < cl; sl ++){
                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ sk: "<< sk << ", sl: " << sl << ": \n";

                    // identify pairs_tab column for (sk,sl)
                    unsigned int i3 = sk * cl;
                    unsigned int r = i1 + i2 + i3 + sl;

                    // read freq
                    unsigned int n_sksl = pairs_tab(4, r);

                    // read prob
                    double pi_sksl = pairs_tab(5, r);

                    // identify tau_sk, tau_sl, tau_sk-1, tau_sl-1
                    Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);
                    double t_sk = pi_thresholds(0); double t_sl = pi_thresholds(1); double t_sk_prev = pi_thresholds(2); double t_sl_prev = pi_thresholds(3);


                    // phi(t_sk, t_sl; rho_kl)
                    double d1 = pbv_rcpp_dbvnorm0( t_sk, t_sl, rho_kl, 0);

                    // phi(t_sk, t_sl-1; rho_kl)
                    double d2 = pbv_rcpp_dbvnorm0( t_sk, t_sl_prev, rho_kl, 0);

                    // phi(t_sk-1, t_sl; rho_kl)
                    double d3 = pbv_rcpp_dbvnorm0( t_sk_prev, t_sl, rho_kl, 0);

                    // phi(t_sk-1, t_sl-1; rho_kl)
                    double d4 = pbv_rcpp_dbvnorm0( t_sk_prev, t_sl_prev, rho_kl, 0);

                    tmp_kl += (n_sksl/(pi_sksl+1e-8)) * ( d1 - d2 - d3 + d4);
                }
            }
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ tmp_kl:" << tmp_kl << "\n";

            ///////////////////////////////////////////////////
            // (k,l)-pair likelihood derivative wrt loadings //
            ///////////////////////////////////////////////////
            if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Gradient wrt loadings: \n";

            // double loop: iterate over elements of loadings matrix
            for(unsigned int j = 0; j < p; j++){
                for(unsigned int v = 0; v < q; v++){
                    if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ visiting lambda_"<< j << v <<":\n";

                    // elicit three cases: 1. free loading item k, 2. free loading l, 3. other
                    if(j == k){
                        if(A(j,v)!=0 ){
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item k, free loading:\n";
                            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
                            double d_rho_kl = ev.transpose() * Sigma_u * lambdal;
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
                            tmp_gradient(iter) += tmp_kl * d_rho_kl;
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << "\n";

                            iter ++;
                        }
                    }else if (j == l){
                        if(A(j,v)!=0 ){
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ item l, free loading:\n";
                            Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
                            double d_rho_kl = lambdak.transpose() * Sigma_u * ev;
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ d_rho_kl:" << d_rho_kl << "\n";
                            if(SILENTFLAG == 0)Rcpp::Rcout << "  |   |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << "\n";
                            tmp_gradient(iter) += tmp_kl * d_rho_kl;

                            iter ++;
                        }
                    }else if(A(j,v)!=0){

                        if(SILENTFLAG == 0)Rcpp::Rcout << "  |  |_ Loadings:: " << iter - tau.size() << "/"<< nload <<". Tot:: "<< iter << "/"<< d << " [not included]\n";
                        iter ++;
                    }
                }
            }
            if(SILENTFLAG == 0)Rcpp::Rcout << "  |_ Done. \n";

            //////////////////////////////////////////////////////////////
            // (k,l)-pair likelihood derivative wrt latent correlations //
            //////////////////////////////////////////////////////////////
            if(SILENTFLAG == 0)Rcpp::Rcout << "\n- Gradient wrt correlations: \n";

            if(CORRFLAG == 1){
                // double loop: iterate over each non-redundant latent correlation
                for(unsigned int v = 1; v < q; v++){
                    for(unsigned int  t = 0; t < v; t++){
                        Eigen::VectorXd ev(q); ev.fill(0.0); ev(v) = 1;
                        Eigen::VectorXd et(q); et.fill(0.0); et(t) = 1;
                        double trho = transformed_rhos(iter - tau.size() - nload);
                        double drho = 2*exp(2*trho) * pow((exp(2*trho) + 1),-1) * ( 1 - ( exp(2*trho) - 1) * pow((exp(2*trho) + 1),-1) );

                        // impose symmetric structure
                        Eigen::MatrixXd Jvt = ev * et.transpose();
                        Eigen::MatrixXd Jtv = et * ev.transpose();
                        Eigen::MatrixXd Svt = Jvt + Jtv - Jvt*Jvt;

                        double d_rho_kl = lambdak.transpose() * (Svt * drho) * lambdal;

                        //if(SILENTFLAG == 0)
                        tmp_gradient(iter) += tmp_kl * d_rho_kl;
                        iter ++;
                    }
                }
            }

            if(SILENTFLAG == 0)Rcpp::Rcout << "\n=====> gradient r-th pair:\n" << gradient << "\n";
        }
    }else{
        for(unsigned int sk = 0; sk < ck; sk ++){

            // i3: starting index from i2 for cat sk
            unsigned int i3 = sk * cl;

            for(unsigned int sl = 0; sl < cl; sl ++){

                // final column index for pairs_tab. Print to check
                unsigned int r = i1 + i2 + i3 + sl;

                // read frequency
                unsigned int n_sksl = pairs_tab(4, r);

                // identify thresholds
                Eigen::VectorXd pi_thresholds = extract_thresholds(tau, C_VEC, k, l, sk, sl);

                // compute pi
                double pi_sksl = compute_pi(C_VEC, pi_thresholds, rho_kl, k, l, sk, sl);
                //if(SILENTFLAG == 0)Rcpp::Rcout << "("<<k<<","<<l<<","<<sk<<","<<sl<<"), rho_kl:"<<rho_kl<<", t_sk:"<< pi_thresholds(0)<<", t_sl:"<< pi_thresholds(1)<<", t_sk-1:"<< pi_thresholds(2)<<", t_sl-1:"<< pi_thresholds(3)<<", pi: "<< pi_sksl<< "\n";
                pairs_tab(5,r) = pi_sksl;
                Eigen::VectorXd pi_grad = compute_pi_grad(A, C_VEC, pi_thresholds, Sigma_u, Lam, THETA, rho_kl, k, l, sk, sl, CORRFLAG);
                Eigen::VectorXd pi_grad2 = compute_pi_grad2(A, C_VEC, pi_thresholds, Sigma_u, Lam, THETA, rho_kl, k, l, sk, sl, CORRFLAG);
                Eigen::VectorXd pi_grad_sq = pi_grad.array().square();

                tmp_gradient += (n_sksl/(pi_sksl+1e-8))*pi_grad;
                tmp_gradient2+= (n_sksl/(pi_sksl+1e-8))*pi_grad2-(n_sksl/(pow(pi_sksl+1e-8,2)))*pi_grad_sq;

                // update ll
                tmp_ll += n_sksl * log(pi_sksl+1e-8);
            }
        }
    }

    ll += tmp_ll;
    gradient += tmp_gradient;
    gradient2 += tmp_gradient2;

}

// RcppParallel Worker to compute pair_contribution
// on multiple pairs in parallel
struct SubsetWorker : public RcppParallel::Worker{
    // Declaration parameters:
    //// Global:
    const Eigen::Map<Eigen::MatrixXd> constrmat;
    const Eigen::Map<Eigen::VectorXd> c_vec;
    const Eigen::MatrixXd &pairs_table;
    const Eigen::MatrixXd &items_pairs;
    const unsigned int corrFLAG;
    const unsigned int silentFLAG;
    const unsigned int DFLAG;
    const unsigned int gradFLAG;

    //// Iteration:
    const Eigen::VectorXd &theta;
    const std::vector<int> &index_vector;

    // Output quantities:
    double subset_ll;
    Eigen::VectorXd subset_gradient, subset_gradient2;

    // Constructor 1:
    SubsetWorker(
        const Eigen::Map<Eigen::MatrixXd> CONSTRMAT,
        const Eigen::Map<Eigen::VectorXd> C_VEC,
        const Eigen::MatrixXd &PAIRS_TABLE,
        const Eigen::MatrixXd &ITEMS_PAIRS,
        const unsigned int CORRFLAG,
        const unsigned int SILENTFLAG,
        const unsigned int DFLAG,
        const unsigned int GRADFLAG,
        const Eigen::VectorXd &THETA,
        const std::vector<int> &INDEX_VECTOR
    ):
        constrmat(CONSTRMAT), c_vec(C_VEC), pairs_table(PAIRS_TABLE), items_pairs(ITEMS_PAIRS),
        corrFLAG(CORRFLAG), silentFLAG(SILENTFLAG), DFLAG(DFLAG), gradFLAG(GRADFLAG),
        theta(THETA), index_vector(INDEX_VECTOR), subset_ll(0.0),
        subset_gradient(Eigen::VectorXd::Zero(THETA.size())),
        subset_gradient2(Eigen::VectorXd::Zero(THETA.size())){}

    // Constructor 2:
    SubsetWorker(const SubsetWorker &OBJ, RcppParallel::Split):
        constrmat(OBJ.constrmat), c_vec(OBJ.c_vec), pairs_table(OBJ.pairs_table), items_pairs(OBJ.items_pairs),
        corrFLAG(OBJ.corrFLAG), silentFLAG(OBJ.silentFLAG), DFLAG(OBJ.DFLAG), gradFLAG(OBJ.gradFLAG),
        theta(OBJ.theta), index_vector(OBJ.index_vector), subset_ll(0.0),
        subset_gradient(Eigen::VectorXd::Zero(theta.size())),
        subset_gradient2(Eigen::VectorXd::Zero(theta.size())){}

    // MEMBER FUNCTIONS:
    //// Overload operator () for main computation:
    void operator()(std::size_t BEGIN, std::size_t END);
    void join(const SubsetWorker &RHS);

};

void SubsetWorker::operator()(std::size_t BEGIN, std::size_t END){
    // Run along the pairs identified by indexes in index_vector
    // computing their contribution to nll, gradient and eventually Hessian.
    unsigned int d = theta.size();
    //subset_gradient.resize(d); subset_gradient.setZero();
    for (unsigned int h = BEGIN; h < END; h++){

        // identify corresponding column in pairs_table
        unsigned int col = index_vector[h];

        // identify the pair
        unsigned int k = items_pairs(0, col);
        unsigned int l = items_pairs(1, col);

        // initialize empty pair-output
        double pair_ll = 0;
        Eigen::VectorXd pair_gradient(d); pair_gradient.setZero();
        Eigen::VectorXd pair_gradient2(d); pair_gradient2.setZero();

        // computation of log-likrlihood, gradient and Hessian diagonal
        pair_contribution(constrmat, c_vec, theta, corrFLAG, k, l, pairs_table, silentFLAG, DFLAG, gradFLAG, pair_ll, pair_gradient, pair_gradient2);

        // update
        {
            subset_ll += pair_ll;
            subset_gradient += pair_gradient;
            subset_gradient2 += pair_gradient2;
        }
    }
}

void SubsetWorker::join(const SubsetWorker &RHS){
    subset_ll += RHS.subset_ll;
    for(unsigned i = 0; i < subset_gradient.size(); i++){
        subset_gradient(i) += RHS.subset_gradient(i);

    }
    for(unsigned i = 0; i < subset_gradient.size(); i++){
        subset_gradient2(i) += RHS.subset_gradient2(i);

    }
}


#endif
