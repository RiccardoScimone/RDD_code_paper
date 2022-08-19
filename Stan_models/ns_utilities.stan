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
