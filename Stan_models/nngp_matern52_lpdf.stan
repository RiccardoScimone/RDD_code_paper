real nngp_matern52_lpdf(vector Y, vector estimated_means, vector[] Locations, int[,] NN_ind, int N, int M, real sigma, real rho)
{
          vector[N] V;
          vector[N] Y_centered = Y - estimated_means;
          vector[N] U = Y_centered;
          for (i in 2:N) {
              int dim = (i < (M + 1))? (i - 1) : M;
              matrix[dim,dim] iNNMcovariance;
              matrix[dim,dim]   iNNMcovarianceCholL;
              vector[dim] iNNMcov_vector;
              vector[dim] v;
              row_vector[dim] v2;                            
              int indexes[dim] = NN_ind[i - 1, 1:dim];
              int dummy[1];
              matrix[dim+1,dim+1] tcov;
              dummy[1] = i;
              tcov = gp_matern52_cov(Locations[append_array(dummy,indexes)],sigma, rho);
              iNNMcovariance = tcov[2:,2:];
              iNNMcovarianceCholL = cholesky_decompose(iNNMcovariance);
              iNNMcov_vector = tcov[2:,1];
              v = mdivide_left_tri_low(iNNMcovarianceCholL, iNNMcov_vector);
              V[i] = square(sigma) - dot_self(v);
              v2 = mdivide_right_tri_low(v', iNNMcovarianceCholL);
              U[i] = U[i] - sum(v2' .* Y_centered[indexes]);
          }
          V[1] = square(sigma);
          return - 0.5 * ( dot_product(U, (U ./ V)) + sum(log(V)));
}