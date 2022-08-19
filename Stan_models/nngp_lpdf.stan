real nngp_lpdf(vector Y, vector estimated_means, vector sigma, matrix Locations, matrix[] Aniso, vector dets, int[,] NN_ind, 
int N, int M, int nu)
{

          vector[N] V;
          vector[N] Y_centered = Y - estimated_means;
          vector[N] U = Y_centered;
          // real kappa_p_1 = tausq / sigmasq + 1;
         //  int dim;

          for (i in 2:N) {
              int dim = (i < (M + 1))? (i - 1) : M;
              matrix[dim,dim] iNNMcovariance;
              matrix[dim,dim]   iNNMcovarianceCholL;
              vector[dim] iNNMcov_vector;
              vector[dim] v;
              row_vector[dim] v2;                            
              int indexes[dim] = NN_ind[i - 1, 1:dim];

              if(dim == 1) iNNMcovariance[1, 1] = square(sigma[i]);
              
              else  iNNMcovariance = matern_ns_corr_mat(Locations[,indexes],Aniso[indexes],dets[indexes], nu, sigma[indexes]);
              
              iNNMcovarianceCholL = cholesky_decompose(iNNMcovariance);
              iNNMcov_vector = matern_ns_corr_vec(Locations[,i],Locations[,indexes],Aniso[i],Aniso[indexes],dets[i],dets[indexes],nu,sigma[i],sigma[indexes]);
              v = mdivide_left_tri_low(iNNMcovarianceCholL, iNNMcov_vector);
              V[i] = square(sigma[i]) - dot_self(v);
              v2 = mdivide_right_tri_low(v', iNNMcovarianceCholL);
              U[i] = U[i] - sum(v2' .* Y_centered[indexes]);
          }
          V[1] = square(sigma[1]);
          return - 0.5 * ( dot_product(U, (U ./ V)) + sum(log(V)));
}