matrix simple_kriging(matrix Sigma, int N_obs, int N_pred)
{ 
  matrix[N_obs,N_obs] L = cholesky_decompose(Sigma[1:N_obs,1:N_obs]);
  matrix[N_pred,N_obs] Weights = mdivide_right_tri_low(mdivide_left_tri_low(L,Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)])',L);
  return(append_col(Weights,Sigma[(N_obs+1):(N_pred + N_obs),(N_obs+1):(N_pred + N_obs)]-Weights*Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)] ));
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
  return(append_col( (sum_row_mat_vec(Sigma[(N_obs+1):(N_pred + N_obs),1:N_obs], -m_hat_variance*colsum(Simple_Weights)) + m_hat_variance)*Precision,
                      Sigma[(N_obs+1):(N_pred + N_obs),(N_obs+1):(N_pred + N_obs)] -Simple_Weights*Sigma[1:N_obs,(N_obs+1):(N_pred + N_obs)] +
                      m_hat_variance * tcrossprod( to_matrix(1 - colsum(Simple_Weights)) ) ));
}

