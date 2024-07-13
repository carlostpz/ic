// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include <Rcpp.h>
using namespace Rcpp;

//////////////////////////////////////////////
//              Main function
//////////////////////////////////////////////


///////////////////////////////
//  indexes for the documents
///////////////////////////////

// [[Rcpp::export]]
arma::mat test_sum_over_rows( arma::mat x){
  return sum(x, 1);
}


// [[Rcpp::export]]
arma::mat mini_index( int minibatch_size, int n_minibatch){
  
  arma::mat result( n_minibatch, 2);
  
  for ( int line=1; line < n_minibatch + 1.0; line++ ){
    result( line-1, 0 ) = minibatch_size * ( line - 1.0 );
    result( line-1, 1 ) = minibatch_size * line - 1.0;
  }
  
  return result;
}

// [[Rcpp::export]]
arma::mat select_doc( arma::mat data ){
  
  int D = data( data.n_rows - 1, 0);
  arma::mat doc_indexes( D, 2);
  doc_indexes.fill(0);
  
  int line_1 = 0;
  int line_2 = 0;
  int current_line = 0;
  int current_doc = data( current_line, 0 );
  
  for (int d = 1; d < D; d++){
    
    // finding the beginning of document d
    while ( current_doc < d ){
      current_line = current_line + 1;
      current_doc = data( current_line, 0 );
    }
    
    line_1 = current_line;
    
    // finding the end of document d
    while ( current_doc < d+1 ){
      current_line = current_line + 1;
      current_doc = data( current_line, 0 );
    }
    line_2 = current_line - 1;
  
    // beginning and end of document d
    doc_indexes(d-1, 0) = line_1 ;
    doc_indexes(d-1, 1) = line_2 ;
    
  }
  
  // the last document is a bit different
  doc_indexes(D-1, 0) = line_2 + 1;
  doc_indexes(D-1, 1) = data.n_rows - 1;
  
  return doc_indexes;
    
}

// [[Rcpp::export]]
arma::mat digamma_mat(arma::mat x){
  
  int c = x.n_cols;
  int r = x.n_rows;
  
  arma::mat res( r, c);
  
  for ( int i=0; i<r; i++){
    for ( int j=0; j<c; j++){
      res(i,j) = R::digamma( x(i,j) ); 
    }
  }
  
  return res;
  
}

// [[Rcpp::export]]
Rcpp::List lda_vi_cpp ( arma::mat original_data, 
                    int n_topics, 
                    int n_local_iter, 
                    int minibatch_size,
                    double step_t,
                    double tau,
                    double kappa,
                    arma::mat gamma_init,
                    arma::mat lambda_init,
                    arma::mat alpha,
                    arma::mat eta,
                    int when_calc_lik,
                    int n_loops_through_data){
  
  // initializing lambda
  arma::mat lambda = lambda_init;
  
  // initializing gamma
  arma::mat gamma = gamma_init;
  
  // documents index
  arma::colvec index = original_data.col(0);
  
  // excluding column corresponding to the documents from the dataset
  arma::mat data = original_data.cols(1, 2);
  data.col(0) = data.col(0) - 1; // accounting for the fact that c++ starts to count from 0
  
  // Number of documents
  int D = max( original_data.col(0) );
  
  // Number of unique words in the entire corpse
  int N = max( data.col(0) ) + 1;     // be careful with calulations. N and D are INTEGERS!
                                      // the data may contain words that never appear in the docs!
                                      // +1 because data.col(0) now starts from 0 instead of 1
  
  // Number of minibatches that fit in the data
  int n_minibatch = D/minibatch_size;
  
  // defining psi_lambda = digamma(lambda)
  arma::mat psi_lambda( n_topics, N);
  
  // defining psi( sum_n( lambda_kn ) )
  arma::colvec psi_sum_lambda( n_topics );
  
  // pre-cashing psi( sum_n( lambda_kn ) )
  psi_sum_lambda = digamma_mat( sum( lambda, 1 ) );
  
  // matrix with minibatch indexes
  arma::mat mini_id_matrix( n_minibatch, n_loops_through_data);
  mini_id_matrix = mini_index( minibatch_size, n_minibatch );
  
  // matrix with document indexes
  arma::mat doc_id_matrix( D, 2);
  doc_id_matrix = select_doc( original_data );
  
  // Creating list for phi and initializing dimensions
  Rcpp::List phi_list(D);
  for(int d=0; d<D; d++){
    
    // indices corresponding to document d
    arma::rowvec doc_d_index = doc_id_matrix.row(d);
    // how many UNIQUE words there are in document d
    int Nd = doc_d_index(1) - doc_d_index(0) + 1;
    // creating object phi
    arma::mat phi( Nd, n_topics); phi.fill(1.0/n_topics);
    // storing phi in the list
    phi_list(d) = phi;
    
  }
  
  
  // vector that will store the log likelihood
  arma::rowvec log_lik_vec( D * n_loops_through_data / when_calc_lik );
  
  // auxiliar variable that says which entry should be updated in log_lik_vec
  int entry_to_update = -1;
  
  arma::mat document_d;
  double rho_t;
  
  for ( int iter_data = 0; iter_data < n_loops_through_data;  iter_data++ ){
    
    // Loop through different minibatches     (n_minibatch)
    for ( int mini = 0; mini < n_minibatch; mini++){
      
      //if (mini == 261){
      //  cout << " sum(lambda_261) inicio " << sum(lambda, 1) << endl;
      //}
      
      // pre-cashing psi(lambda) and psi( sum_over_n( lambda_k ) ) to be used in phi step
      psi_lambda = digamma_mat( lambda );
      
      // pre-cashing psi( sum_n( lambda_kn ) )
      psi_sum_lambda = digamma_mat( sum( lambda, 1 ) ); // sum over each row. Yes, row is 1 and column is 0.
      
      arma::rowvec mini_id = mini_id_matrix.row( mini );
      step_t = step_t + 1;
      
      // cashing E( log(beta_kn) )
      arma::mat exp_E_log_beta(n_topics, N);
      for ( int k=0; k<n_topics; k++){
        exp_E_log_beta.row( k ) = exp( psi_lambda.row(k) - psi_sum_lambda(k) );
      }
  
      /////////////////////////////////////////////////
      //  Loop trhough documents within the minibatch
      /////////////////////////////////////////////////
      for( int d = mini_id_matrix(mini, 0) ; d < mini_id_matrix(mini, 1) + 1 ; d++ ){
        
        cout << "d= " << d << endl;
        
        //////////////////////
        //  update for phi_d
        //////////////////////
        
        // indices corresponding to document d
        arma::rowvec doc_d_index = doc_id_matrix.row(d);
        
        // document d
        document_d = data.rows( doc_d_index(0), doc_d_index(1) );
        
        // unique words and freq in document d
        arma::colvec uniq_words_doc_d = document_d.col(0); 
        arma::colvec freq_words_doc_d = document_d.col(1);
        
        // how many UNIQUE words there are in document d
        int Nd = document_d.n_rows;
  
        // creating object phi
        arma::mat phi( Nd, n_topics);
        
        // Notice: phi is indexed by phi_nk, where n goes from 1 to Nd. 
        // Opposing the remaining variables, n here reffers to the n-th UNIQUE word *in document d*
        
        
        for (int local_iter=0; local_iter < n_local_iter; local_iter++){
         
          //cashing digamma( gamma[d, k] ) for every topic k    (vector)
          arma::rowvec psi_gamma_d = digamma_mat( gamma.row(d) );
          
          // cashing digamma( sum_over_k( gamma[d, k] ) )    (scalar)
          double psi_sum_gamma_d = R::digamma( sum(gamma.row(d) ) );
          
          // cashing exp( E_q[ log(theta_d) ] ) (vecor of dimension n_topics)
          arma::rowvec exp_E_log_theta_d = exp( psi_gamma_d - psi_sum_gamma_d );

          //////////////////////////
          // main update for phi
          //////////////////////////
          for( int n=0; n<Nd; n++){
            // cout << "updating phi, n= "<< n << " k=0 " << endl;
            for ( int k=0; k<n_topics; k++){
          
              phi(n, k) = exp_E_log_theta_d(k) * exp_E_log_beta( k, uniq_words_doc_d(n) );
              
            }
            //normalize phi
            phi.row(n) = phi.row(n)/ sum( phi.row(n) );
          }
          // end of update for phi_d
  
          /////////////////////
          // update for gamma
          /////////////////////
          for (int k = 0; k < n_topics; k++){
            gamma(d,k) = alpha(d, k) + dot( phi.col(k), freq_words_doc_d);
          }
          
        }
        // end of loop for local variables in document d
        
        //saving phi_d to be accessed latter when updating phi
        phi_list(d) = phi;
        phi.fill(0);
        
      }
      // end of loop through documents within the minibatch
      
      // for each minibatch, update lambda once
      
      arma::mat lambda_hat( n_topics, N); lambda_hat.fill(0);
      
      ///////////////////////////////////
      //       update for lambda
      ///////////////////////////////////
      
      for( int d = mini_id_matrix(mini, 0) ; d < mini_id_matrix(mini, 1) + 1 ; d++  ){
        
        arma::mat phi = phi_list(d);
        
        // indices corresponding to document d
        arma::rowvec doc_d_index = doc_id_matrix.row(d);
        
        // restricting ourselves to document d
        document_d = data.rows( doc_d_index(0), doc_d_index(1) );
        
        // unique words and freq in document d
        arma::colvec uniq_words_doc_d = document_d.col(0); 
        arma::colvec freq_words_doc_d = document_d.col(1);
        
        // how many UNIQUE words there are in document d
        int Nd = document_d.n_rows;
        
        for (int n=0; n<Nd; n++){
          int word = uniq_words_doc_d(n);
          for ( int k=0; k<n_topics; k++){
            lambda_hat( k, word ) = lambda_hat( k, word ) + freq_words_doc_d(n) * phi(n, k);
          }
        }
  
      }
      // we already went through all documents in the minibatch to accumulate the sum for updating phi

      // finishing update for phi
      lambda_hat = ( ( D + 0.0 ) / ( minibatch_size + 0.0 ) ) * lambda_hat;
      lambda_hat = lambda_hat + eta;
      
      // Robins-Monroe stepzise rule
      rho_t = pow( 1/( step_t + tau ), kappa );
      for ( int n=0; n < N; n++ ){
        for ( int k=0; k<n_topics; k++){
          lambda(k, n) = ( 1- rho_t) * lambda(k, n) + rho_t * lambda_hat(k, n);
        }
      }

      ///////////////////////////////
      //       Log likelihood
      ///////////////////////////////
      
      // calculating the estimated log likelihood (only once at each "when_calc_lik" documents)
      if ( ( (mini + 1) * minibatch_size) % ( when_calc_lik ) == 0 ){
        
        // which entry of lok_lik should I update
        entry_to_update = entry_to_update + 1;
        
        // estimating beta
        arma::mat beta_hat( n_topics, N);
        for (int k=0; k < n_topics; k++ ){
          beta_hat.row(k) = lambda.row(k)/sum( lambda.row(k) );
        }
        
        // initializing log likelihood
        double log_lik = 0;
        
        for (int d=0; d<D; d++){
          
          //if( d == 452){
          //  continue;
          //}
          
          // indices corresponding to document d
          arma::rowvec doc_d_index = doc_id_matrix.row(d);
          
          // restricting ourselves to document d
          document_d = data.rows( doc_d_index(0), doc_d_index(1) );
          
          // unique words and freq in document d
          arma::colvec uniq_words_doc_d = document_d.col(0); 
          arma::colvec freq_words_doc_d = document_d.col(1);
          
          
          // extracting the right phi
          arma::mat phi = phi_list(d);
          
          // how many UNIQUE words there are in document d
          int Nd = document_d.n_rows;
          
          for ( int n=0; n < Nd; n++){
            
            int z_hat_dn = index_max( phi.row(n) );
            log_lik = log_lik + log( beta_hat( z_hat_dn, uniq_words_doc_d(n) ) ) * freq_words_doc_d(n) ;
          }
          
        }
        // end of loop through documents to calculate the log likelihood

        log_lik_vec( entry_to_update ) = log_lik; 
        
      }
      // end of one log likelihood calculation
      
    }
    // end of loop through minibatches 
    
  }
  //end of loop trhough the entire dataset
  
  Rcpp::List ans = Rcpp::List::create( Rcpp::Named("rho") = rho_t,
                                       Rcpp::Named("phi") = phi_list,
                                       Rcpp::Named("lambda") = lambda,
                                       Rcpp::Named("gamma") = gamma,
                                       Rcpp::Named("doc_id_matrix") = doc_id_matrix,
                                       Rcpp::Named("log_likelihood") = log_lik_vec);
                                       
  return ans;
  
}
