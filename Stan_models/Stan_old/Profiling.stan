functions {

   matrix matern_ns_corr_mat(matrix Locations, matrix[] Aniso, vector dets, int  nu)
  { 
    int N = cols(Locations);
    vector[N] sqrsqrdets = sqrt(sqrt(dets));
    real norm_const = 2^(1-nu)/exp(lgamma(nu));
    matrix[N,N] C;
    for ( i in 1:N){
      C[i,i] = 1;
      for ( j in (i+1):N)
      { 
        matrix[2,2] metric = 0.5 * (Aniso[i] + Aniso[j]);
        real normdet = sqrt(1/(metric[1,1]*metric[2,2] - square(metric[1,2])));
        real Q = sqrt( (col(Locations,i)  - col(Locations,j))' * (metric\(col(Locations,i)  - col(Locations,j))));
        C[i,j] = sqrsqrdets[i]*sqrsqrdets[j] * normdet * norm_const * Q^nu * modified_bessel_second_kind(nu,Q);
        C[j,i] = C[i,j];
      }
    }
        return (C);
  }
  
    vector matern_ns_corr_vec(vector ref_loc, matrix Locations, matrix ref_Aniso, matrix[] Aniso, real ref_det, vector dets, int nu )
  {
    int N = cols(Locations);
    vector[N] sqrsqrdets = sqrt(sqrt(dets));
    real sqrsqrdet = sqrt(sqrt(ref_det));
    real norm_const = 2^(1-nu)/exp(lgamma(nu));
    vector[N] v;
    for ( j in 1:N)
      { 
        matrix[2,2] metric = 0.5 * (ref_Aniso + Aniso[j]);
        real normdet = sqrt(1/(metric[1,1]*metric[2,2] - square(metric[1,2])));
        real Q = sqrt( (ref_loc  - col(Locations,j))' * (metric\(ref_loc  - col(Locations,j))));
        v[j] = sqrsqrdet*sqrsqrdets[j] * normdet * norm_const * Q^nu * modified_bessel_second_kind(nu,Q);
      }
        return (v);
  }
  
  
  
  real nngp_lpdf(vector Y, vector X_beta, real sigmasq, real tausq, vector theta, vector lambda1, vector lambda2,
                matrix Locations, int[,] NN_ind,int N, int M, int nu){

          vector[N] V;
          vector[N] YXb = Y - X_beta;
          vector[N] U = YXb;
          real kappa_p_1 = tausq / sigmasq + 1;

          matrix[2,2] Anis[N];
          vector[N] dets;
          
          dets = lambda1 .* lambda2;
          for (i in 1:N){
            matrix[2,2] rot = [[cos(theta[i]), -sin(theta[i]) ] , [sin(theta[i]), cos(theta[i])]];
            matrix[2,2] eig = [ [lambda1[i],0] , [0,lambda2[i]]];
            Anis[i] = rot*eig*rot';
          }



          for (i in 2:N) {
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNdistM;
              matrix[ i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1): M]
              iNNCholL;
              vector[ i < (M + 1) ? (i - 1) : M] iNNcorr;
              vector[ i < (M + 1) ? (i - 1) : M] v;
              row_vector[i < (M + 1) ? (i - 1) : M] v2;
              matrix[rows(iNNdistM)+1,cols(iNNdistM)+1] t_iNNdistM;
              int dim = (i < (M + 1))? (i - 1) : M;
              int idxs[dim];
              idxs= NN_ind[i-1,1:dim];
              if(dim == 1){iNNdistM[1, 1] = kappa_p_1;}
              else{
                 
                  iNNdistM = matern_ns_corr_mat(Locations[1:2,idxs],Anis[idxs],dets[idxs],nu);
                  for(jj in 1:dim){
                      iNNdistM[jj, jj] = kappa_p_1;
                  }
              }

            iNNCholL = cholesky_decompose(iNNdistM);
            iNNcorr = matern_ns_corr_vec(Locations[1:2,i],Locations[1:2,idxs],Anis[i],Anis[idxs],dets[i],dets[idxs],nu);  

             v = mdivide_left_tri_low(iNNCholL, iNNcorr);

             V[i] = kappa_p_1 - dot_self(v);

             v2 = mdivide_right_tri_low(v', iNNCholL);

             for (j in 1:dim){
                  U[i] = U[i] - v2[j] * YXb[NN_ind[(i - 1), j]];
              }
          }
          V[1] = kappa_p_1;
          return - 0.5 * ( 1 / sigmasq * dot_product(U, (U ./ V)) +
                          sum(log(V)) + N * log(sigmasq));
      }

  
}



      data {
      int<lower=1> N;
      int<lower=1> M;
      int<lower=0> P;
      int<lower=0> P_lambdas;
      int<lower=0> P_theta;
      int<lower=1> nu;
      vector[N] Y;
      matrix[N, P + 1] X;
      matrix[N, P_lambdas + 1] X_lambdas;
      matrix[N, P_theta + 1] X_theta;
      int NN_ind[N - 1, M];
      matrix[N,2] Locations;
//      vector[P + 1] uB;
      matrix[P + 1, P + 1] VB_beta;
      matrix[P_lambdas + 1, P_lambdas + 1] VB_lambdas;
      matrix[P_theta + 1, P_theta + 1] VB_theta;
//      real ss;
//      real st;
  }

  transformed data {
      cholesky_factor_cov[P + 1] L_VB_beta;
      cholesky_factor_cov[P_lambdas + 1] L_VB_lambdas;
      cholesky_factor_cov[P_theta + 1] L_VB_theta;

      matrix[2,N] t_Locations;
      L_VB_beta = cholesky_decompose(VB_beta);
      L_VB_lambdas = cholesky_decompose(VB_lambdas);
      L_VB_theta = cholesky_decompose(VB_theta);
      t_Locations = Locations';
  }

  parameters{
      vector[P + 1] alfa_beta;
      vector[P_lambdas + 1] alfa_beta_lambda1;
      vector[P_lambdas + 1] alfa_beta_lambda2;
      vector[P_theta + 1] alfa_beta_theta;
  //    real<lower = 0> sigma;
  //    real<lower = 0> tau;
  }

  transformed parameters {
  //    real sigmasq = square(sigma);
  //    real tausq = square(tau);
      vector[P + 1] beta = L_VB_beta*alfa_beta;
      vector[P_lambdas + 1] beta_lambda1 = 1 + L_VB_lambdas*alfa_beta_lambda1;
      vector[P_lambdas + 1] beta_lambda2 = 1 + L_VB_lambdas*alfa_beta_lambda2; 
      vector[P_theta + 1] beta_theta = 1 + L_VB_theta*alfa_beta_theta; 
  }

  model{
      alfa_beta ~ std_normal();
      alfa_beta_lambda1 ~ std_normal();
      alfa_beta_lambda2 ~ std_normal();
      alfa_beta_theta ~ std_normal();
  //    sigma ~ normal(0, ss);
  //    tau ~ normal(0, st);
      Y ~ nngp(X * beta, 1, 0,pi()/2 * inv_logit(X_theta*beta_theta), exp(X_lambdas*beta_lambda1) , exp(X_lambdas*beta_lambda2),t_Locations, NN_ind, N, M, nu);
  }
  