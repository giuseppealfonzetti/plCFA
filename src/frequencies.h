#ifndef frequencies_H
#define frequencies_H

/* COMPUTE FREQUENCIES */
// It returns a 5-rows matrix with each combination of items and categories as columns.
// Row0: item k, Row1: item l, Row2; category item k, Row3: category item l, Row4: freq
// It is computed just once, before the optimization of the complete pairwise
//' @export
// [[Rcpp::export]]
Eigen::MatrixXd pairs_freq(
        Eigen::Map<Eigen::MatrixXd> Y,
        Eigen::Map<Eigen::VectorXd> C_VEC
){

    const unsigned int n = Y.rows(); // number of units
    const unsigned int p = Y.cols(); // number of items

    const unsigned int R = p*(p-1)/2; // number of pairs
    const unsigned int c = C_VEC.sum();
    unsigned int iter; // helper iterator

    Eigen::MatrixXd freq(5,1);

    // Find how many possible pairs and setup freq matrix
    iter = 0;
    for(int k = 1; k < p; k++){
        int ck = C_VEC(k);
        for(int l = 0; l < k; l ++){
            int cl = C_VEC(l);
            for(int sk = 0; sk < ck; sk++){
                for(int sl = 0; sl < cl; sl ++){
                    freq.conservativeResize(5, iter + 1);
                    Eigen::VectorXd freq_coord(5);
                    freq_coord << k, l, sk, sl, 0;
                    freq.col(iter) = freq_coord;
                    iter++;
                }
            }
        }
    }

    auto task = [&freq, &n, &Y](unsigned int r){

        // read data
        int k  = freq(0,r);
        int l  = freq(1,r);
        int sk = freq(2,r);
        int sl = freq(3,r);
        Eigen::MatrixXd obs_resp(n,2);
        obs_resp.col(0) = Y.col(k); obs_resp.col(1) = Y.col(l);
        //if(silentFLAG == 0)Rcpp::Rcout << "k: " << k << " , l: "<< l << " ,sk: " << sk << " , sl: " << sl << "\n";
        //if(silentFLAG == 0)Rcpp::Rcout << "obs responses : \n" << obs_resp << "\n";
        // compute frequency
        int n_sksl = 0;
        for(int i = 0; i < n; i++){
            if(sk == obs_resp(i,0) & sl == obs_resp(i,1)) {
                //if(silentFLAG == 0)Rcpp::Rcout << "Match on unit " << i <<"!\n";
                n_sksl++;
            }
        }

        // update
        freq(4, r) = n_sksl;
    };

    // Parallel loop to compute frequencies
    for(int r = 0; r < freq.cols(); r++){
        task(r);
    }


    return freq;
}

/* FREQUENCIES */
// 1. Compute frequencies for each pair-pattern and store it in "freq". Each pair is a identified by a column index.
// 2. Build a pool of indexes related to each item and store it in "items_pools".
// 3. Build up a matrix with all possible pairs o items in "items_pairs".
// Note: the function is not exported to R session but only used internally since it's only
//       used for the stochastic optimization.
void pairs_freq_cpp(
        // input
        const Eigen::MatrixXd Y,
        const Eigen::VectorXd C_VEC,

        // output
        Eigen::MatrixXd &FREQ,
        Eigen::MatrixXd &ITEMS_PAIRS,
        std::vector< std::vector<int> > &ITEMS_POOLS
){

    int n = Y.rows(); // number of units
    int p = Y.cols(); // number of items

    int R = p*(p-1)/2; // number of pairs
    int c = C_VEC.sum();
    int iter; // helper iterator
    int iter1;

    // Serial loop:
    // Build up items_pairs and items_pools, while setting up freq
    iter = 0;
    iter1 = 0;

    for(int k = 1; k < p; k++){
        int ck = C_VEC(k);
        for(int l = 0; l < k; l ++){
            int cl = C_VEC(l);
            for(int sk = 0; sk < ck; sk++){
                for(int sl = 0; sl < cl; sl ++){
                    FREQ.conservativeResize(5, iter + 1);
                    Eigen::VectorXd freq_coord(5);
                    freq_coord << k, l, sk, sl, 0;
                    FREQ.col(iter) = freq_coord;

                    iter++;
                }
            }

            ITEMS_PAIRS(0, iter1) = k;
            ITEMS_PAIRS(1, iter1) = l;
            ITEMS_POOLS[k].push_back(iter1);
            ITEMS_POOLS[l].push_back(iter1);
            iter1++;
        }
    }

    // lambda function to count frequencies
    auto task = [&FREQ, &n, &Y](unsigned int r){

        // read data
        int k  = FREQ(0,r);
        int l  = FREQ(1,r);
        int sk = FREQ(2,r);
        int sl = FREQ(3,r);
        Eigen::MatrixXd obs_resp(n,2);
        obs_resp.col(0) = Y.col(k); obs_resp.col(1) = Y.col(l);
        // compute frequency
        int n_sksl = 0;
        for(int i = 0; i < n; i++){
            if(sk == obs_resp(i,0) & sl == obs_resp(i,1)) {
                //if(silentFLAG == 0)Rcpp::Rcout << "Match on unit " << i <<"!\n";
                n_sksl++;
            }
        }

        // update
        FREQ(4, r) = n_sksl;
    };

    // Loop to compute frequencies
    for(int r = 0; r < FREQ.cols(); r++){
        task(r);
    }
}

#endif
