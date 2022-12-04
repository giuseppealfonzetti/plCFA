#ifndef optimisationUtils_H
#define optimisationUtils_H

#include "genericUtils.h"

/* SAMPLING PAIRS */
//' @export
// [[Rcpp::export]]
std::vector<int> sampling_step(
        std::vector<int> &FULL_POOL,                // Address of pairs indices
        const unsigned int METHODFLAG,                    // 0 hypergeometric, 1 Bernoulli
        const double PROB,
        const int PAIRS_PER_ITERATION,
        const unsigned int N_ITEMS,
        const unsigned int SEED,
        const int SILENTFLAG,
        const unsigned int ITER
){
    // Set-up the randomizer
    std::mt19937 randomizer(SEED + ITER);

    // empty vector to collect chosen indices
    std::vector<int> chosen_index;

    if(METHODFLAG == 0){
        std::shuffle(FULL_POOL.begin(), FULL_POOL.end(), randomizer);

        for(unsigned int draw = 0; draw < PAIRS_PER_ITERATION; draw++){
            if(SILENTFLAG == 0)Rcpp::Rcout << "Draw "<< draw + 1 <<":\n |_ Shuffling the pairs...\n";

            // Collect the index drawn
            unsigned int pair_index = FULL_POOL[draw];
            chosen_index.push_back(pair_index);
            if(SILENTFLAG == 0)Rcpp::Rcout << " |_ Position drawn:" << pair_index << "\n";

        }

    }else if(METHODFLAG == 1){
        for(unsigned int pair = 0; pair < FULL_POOL.size(); pair++){
            if(R::runif(0,1) < PROB ) chosen_index.push_back(pair);
;
        }

    }

    return chosen_index;
}

void stabilise_loadings(
        Eigen::Map<Eigen::MatrixXd> A,
        Eigen::Map<Eigen::VectorXd> C_VEC,
        Eigen::VectorXd &THETA,          // Address used to update projected theta
        bool &CHECKEVENT                 // Turn on the indicator whether there was a projection or not
){
    //dimensions
    const unsigned int p = A.rows(); // number of items
    const unsigned int q = A.cols(); // number of latents
    const unsigned int d = THETA.size(); // number of parameters
    const unsigned int c = C_VEC.sum();

    // Extract parameters
    const Eigen::MatrixXd Sigma_u = get_Sigma_u(A, THETA);
    Eigen::MatrixXd Lam = get_Lam(A, C_VEC, THETA);
    const double eps = 1e-8;


    // Normalisation if needed
    for(unsigned int k = 0; k < p; k++){
        const Eigen::VectorXd lambdak = Lam.row(k);
        const double vk = lambdak.transpose()*Sigma_u*lambdak;
        if(vk>=1){
            Lam.row(k) = lambdak/(vk+eps);
            CHECKEVENT = true;
        }
    }


    if(CHECKEVENT){
        Eigen::VectorXd norm_lambda(d-c+p-q*(q-1)/2);
        unsigned int iter = 0;
        for(unsigned int h = 0; h < q; h++){
            for(unsigned int j = 0; j < p; j++){
                if(A(j, h) != 0.0)
                {
                    norm_lambda(iter) = Lam(j, h);
                    iter ++;
                };
            };
        }

        THETA.segment(c-p, d-c+p-q*(q-1)/2) = norm_lambda;
    }
}

#endif
