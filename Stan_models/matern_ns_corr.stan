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
        C[i,j] = sigma[i]*sigma[j]*sqrsqrdets[i]*sqrsqrdets[j] * normdet * norm_const * Q^nu * besselK(nu,Q);
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
  