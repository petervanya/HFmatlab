% Perform MC integration and return the matrix
% !!! CURRENTLY WORKS ONLY FOR MINIMAL BASIS H2 !!!
% arguments:
% -- nuclei coords Ri (or later xyz)
% -- atomic orbital function phi
% -- constant C
function S = overlap_matrix(R1,R2,phi,C);
  N=2;
  S=zeros(N);
  % make this function arguments
  nMC=1e5;
  sigma=5;
  bound=sigma/C;
  
  % MC integration of the lower half of S
  for i=2:N
    for j=1:i-1
	    for k=1:nMC
		    p = rand(1,3)*2*bound - bound;
		    S(i,j) = S(i,j) + phi(p,R1,C)*phi(p,R2,C);
      end
    end
  end
  
  S = S * (2*bound)^3/nMC;    % MC normalisation
  S = S + conj(S)';           % get the upper half of S
  S = S + diag(ones(N,1));    % get the diagonal (trivially ones)
end
