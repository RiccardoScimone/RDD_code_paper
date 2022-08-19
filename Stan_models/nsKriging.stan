functions {
 real besselK(real nu, real d);

   matrix matern_ns_corr_mat(matrix Locations, matrix[] Aniso, vector dets, real  nu, vector sigma)
  { 
    int N = cols(Locations);
    vector[N] sqrsqrdets = sqrt(sqrt(dets));
    real norm_const = 2^(1-nu)/exp(lgamma(nu));
    matrix[N,N] C;
    for ( i in 1:N){
      C[i,i] = (1+1e-8)*square(sigma[i]);
      for ( j in (i+1):N)
      { 
        matrix[2,2] metric = 0.5 * (Aniso[i] + Aniso[j]);
        real normdet = sqrt(1/(metric[1,1]*metric[2,2] - square(metric[1,2])));
        real Q = sqrt( (col(Locations,i)  - col(Locations,j))' * (metric\(col(Locations,i)  - col(Locations,j))));
        if(Q > 0)
          C[i,j] = sigma[i]*sigma[j]*sqrsqrdets[i]*sqrsqrdets[j] * normdet * norm_const * Q^nu * besselK(nu,Q);
        else
          C[i,j] = C[i,i];
        C[j,i] = C[i,j];
      }
    }
        return (C);
  }
  
    vector matern_ns_corr_vec(vector ref_loc, matrix Locations, matrix ref_Aniso, matrix[] Aniso, real ref_det, vector dets, int nu, real ref_sigma, vector sigma )
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
        v[j] = ref_sigma*sigma[j]*sqrsqrdet*sqrsqrdets[j] * normdet * norm_const * Q^nu * besselK(nu,Q);
      }
        return (v);
  }
  matrix simple_kriging(matrix Sigma, int N_obs, int N_pred)
{ 
  matrix[N_obs,N_obs] L = cholesky_decompose(Sigma[1:N_obs,1:N_obs]);
  matrix[N_pred,N_obs] Weights = mdivide_right_tri_low(mdivide_left_tri_low(L,Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)])',L);
//  return(append_col(Weights,Sigma[(N_obs+1):(N_pred + N_obs),(N_obs+1):(N_pred + N_obs)]-Weights*Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)] ));
  return(Weights);
}

matrix sum_row_mat_vec(matrix A, vector v)
{ 
  matrix[rows(A), cols(A)] Temp;
  for(i in 1:rows(A))
  {
    Temp[i,] = A[i,] + v[i];
  }
  return(Temp);
}

vector colsum(matrix A)
{
  int n = rows(A);
  vector[n] ret;
  for (i in 1:n)
  {
    ret[i] = sum(A[i,]);
  }
  return(ret);
}

matrix ordinary_kriging(matrix Sigma,  int N_obs, int N_pred)
{ 
  matrix[N_obs,N_obs] Precision = mdivide_left_spd(Sigma[1:N_obs,1:N_obs],diag_matrix(rep_vector(1,N_obs)));
  matrix[N_pred,N_obs] Simple_Weights = Sigma[(N_obs+1):(N_pred + N_obs),1:N_obs]*Precision;
  real m_hat_variance = 1/sum(Precision);
//  return(append_col( (sum_row_mat_vec(Sigma[(N_obs+1):(N_pred + N_obs),1:N_obs], -m_hat_variance*colsum(Simple_Weights)) + m_hat_variance)*Precision,
//                      Sigma[(N_obs+1):(N_pred + N_obs),(N_obs+1):(N_pred + N_obs)] -Simple_Weights*Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)] +
//                      m_hat_variance * tcrossprod( to_matrix(1 - colsum(Simple_Weights)) ) ));
  return ((sum_row_mat_vec(Sigma[(N_obs+1):(N_pred + N_obs),1:N_obs], -m_hat_variance*colsum(Simple_Weights)) + m_hat_variance)*Precision);
}

matrix[] Compute_Aniso (vector lambda_1, vector lambda_2, vector theta)
{ 
  int N = num_elements(lambda_1);
  matrix[2,2] Anis[N];
  for (i in 1:N){
            matrix[2,2] rot = [[cos(theta[i]), -sin(theta[i]) ] , [sin(theta[i]), cos(theta[i])]];
            matrix[2,2] eig = [ [lambda_1[i],0] , [0,lambda_2[i]]];
            Anis[i] = rot*eig*rot';
          }
  return(Anis);
          
}

vector Compute_Determinants(vector lambda_1, vector lambda_2)
{
  return(lambda_1 .* lambda_2);
}

}

data {
  int N_prediction;
  int N_obs;
  matrix[2,N_prediction] x_prediction;
  matrix[2,N_obs] x_obs;
  vector[N_obs + N_prediction] lambda_1;
  vector[N_obs + N_prediction] lambda_2;
  vector[N_obs + N_prediction] theta;
  vector[N_obs + N_prediction] sigma;  
  real<lower=0> nu;
}

transformed data{
  matrix[2,2] Aniso[N_prediction + N_obs] = Compute_Aniso(lambda_1,lambda_2,theta);
  matrix[2,N_prediction + N_obs] x = append_col(x_obs, x_prediction);
  vector[N_obs + N_prediction] dets = Compute_Determinants(lambda_1,lambda_2);
  matrix[N_prediction+N_obs,N_prediction + N_obs] cov = matern_ns_corr_mat(x,Aniso,dets,nu,sigma); 
  cov[1:N_obs, 1:N_obs] = add_diag(cov[1:N_obs, 1:N_obs],0.05);
}

generated quantities{
  matrix[N_prediction + N_obs,N_prediction + N_obs] C = cov; 
}









