% Calculate transfer matrix from overlap matrix S
function X = transf_matrix(S)
  [U,s] = eig(S);
  X = U*pinv(sqrt(s))*conj(U');
  
  
end
