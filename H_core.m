% Hamiltonian matrix for core electrons
% H_core = T + V
function Hc = H_core(R1,R2,phi,C)
  N=2;
  Z=1;
  T=zeros(N);
  V=zeros(N);
  h=0.01;
  
  nMC=1e5;
  sigma=5;
  bound=sigma/C;
  
  % Laplacian
  del_phi=@(r,h,R,C) ...
      (phi(r+h*[1 0 0],R,C) + phi(r-h*[1 0 0],R,C) + ...
       phi(r+h*[0 1 0],R,C) + phi(r-h*[0 1 0],R,C) + ...
       phi(r+h*[0 0 1],R,C) + phi(r-h*[0 0 1],R,C) + ...
      -6*phi(r,R,C))/h^2;
  
  % --- UNFINISHED
  % MC integration of lower triang. part of T and V
  for ii=1:N
    for jj=1:ii
	    for k=1:nMC
		    %p = rand(1,3)*2*bound - bound;
		    %T(ii,jj) -= phi(p,R1,C)*del_phi(p,h,R2,C)/2;
		    %V(ii,jj) -= phi(p,R1,C)*(Z/norm(p-R1)+Z/norm(p-R2))*phi(p,R2,C);
      end
    end
  end
  
  % MC int'n of upper triang. part of T for i=1:N
  for ii=1:N-1
	  for jj=ii+1:N
	    for k=1:nMC
		    %p = rand(1,3)*2*bound - bound;
		    % provisional, now only T(2,2) is computed
		    %T(ii,jj) -= phi(p,R2,C)*del_phi(p,h,R1,C)/2;
      end
    end
  end
  
  % symmetrise V
  for ii=1:N-1
    for jj=ii+1:N
      %V(ii,jj) = conj(V(jj,ii));
    end
  end
  % --- END OF UNFINISHED
  
  % minimal basis hydrogen atom, requires only 5 computations
  for k=1:nMC
		p = rand(1,3)*2*bound - bound;
		T(1,1) = T(1,1) - phi(p,R1,C)*del_phi(p,h,R1,C)/2;
		%T(2,2) = T(2,2) - phi(p,R2,C)*del_phi(p,h,R2,C)/2; % redundant
		T(1,2) = T(1,2) - phi(p,R2,C)*del_phi(p,h,R1,C)/2;		
  	V(1,1) = V(1,1) + phi(p,R1,C)^2*(-Z/norm(p-R1)-Z/norm(p-R2));
  	%V(2,2) = V(2,2) + phi(p,R2,C)^2*(-Z/norm(p-R1)-Z/norm(p-R2)); % redundant
  	V(1,2) = V(1,2) + phi(p,R2,C)*phi(p,R1,C)*(-Z/norm(p-R1)-Z/norm(p-R2));
  end
  T(2,1) = T(1,2);
  T(2,2) = T(1,1);
  V(2,1) = V(1,2);
  V(2,2) = V(1,1);
  
  T = T * (2*bound)^3/nMC;
  V = V * (2*bound)^3/nMC;
    
  Hc = T + V;
end


% ----- functions ----- DOES NOT WORK
%function y=d2f(f,x,h)
%  y = (f(x+h)-2*f(x)+f(x-h))/h**2;
%end
