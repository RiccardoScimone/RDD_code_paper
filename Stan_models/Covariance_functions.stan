real besselK(real nu, real d);

real stationary_matern_1d(real x_1, real x_2, real nu, real rho)
{ 
  real d = (x_1 - x_2)^2;
  if(d>0)
  return ( 2^(1-nu/2) / tgamma(nu) * (sqrt(nu*d)/rho)^nu*besselK(nu,sqrt2()*sqrt(nu*d)/rho) );
  return 1;
}

real stationary_matern_vec(vector x_1, vector x_2, real nu, real rho)
{ 
  real d = (x_1 - x_2)'*(x_1 - x_2);
  if(d>0)
  return ( 2^(1-nu/2) / tgamma(nu) * (sqrt(nu*d)/rho)^nu*besselK(nu,sqrt2()*sqrt(nu*d)/rho) );
  return 1;
}

matrix cov_stationary_matern_multi(matrix Locations,real nu, real rho, real sigma)
  { 
    int N = cols(Locations);
    matrix[N,N] C;
    real sigmasq = square(sigma);
    for ( i in 1:N){
      C[i,i] = (1+1e-8)*sigmasq;
      for ( j in (i+1):N)
      { 
        C[i,j] = sigmasq*stationary_matern_vec(Locations[,i],Locations[,j],nu,rho);
        C[j,i] = C[i,j];
      }
    }
        return (C);
  }
  
matrix cov_stationary_matern_uni(real[] Locations,real nu, real rho, real sigma)
  { 
    int N = num_elements(Locations);
    matrix[N,N] C;
    real sigmasq = square(sigma);
    for ( i in 1:N){
      C[i,i] = (1+1e-8)*sigmasq;
      for ( j in (i+1):N)
      { 
        C[i,j] = sigmasq*stationary_matern_1d(Locations[i],Locations[j],nu,rho);
        C[j,i] = C[i,j];
      }
    }
        return (C);
  }
  
real brownian_1d(real x_1, real x_2)
{ 
  return ( fmin(x_1,x_2) );
}

matrix cov_brownian_uni(real[] Locations)
  { 
    int N = num_elements(Locations);
    matrix[N,N] C;
    for ( i in 1:N){
      for ( j in i:N)
      { 
        C[i,j] = brownian_1d(Locations[i],Locations[j]);
        C[j,i] = C[i,j];
      }
      C[i,i] = (1 + 1e-8)*C[i,i];
    }
        return (C);
  }
  
  