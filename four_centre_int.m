% Four-centre integrals without symmetries
% DO LATER THE SYMMETRIES
function fci = four_centre_int(R1,R2,phi,C)
  N=2;
  fci=zeros(N,N,N,N);
  
  n=1e5;
  sigma=5;
  bound=sigma/C;
  
  % MC integration of T and V
  for ii=1:N
    for jj=1:N
      for kk=1:N
        for ll=1:N
        % NOT IMPLEMENTED YET
					%for k=1:N
						%p1 = rand(1,3)*2*bound - bound*ones(1,3);
						%p2 = rand(1,3)*2*bound - bound*ones(1,3);
						% general form
						% have to include the symmetries! 
						% (mu nu|lambda sigma) = (nu mu | lambda sigma) = (lambda sigma | mu nu)
						%fci(ii,jj,kk,ll) += phi(p1,R1,C)*phi(p1,R2,C)*phi(p2,R1,C)*phi(p2,R2,C)/norm(p1-p2);
				  %end
				  %fci(ii,jj) *= (2*bound)^6/n;
				end
			end
    end
  end
  
  % minimal basis hydrogen atom, requires only 5 computations
  % (11|11), (22|22), (11|22), (11|12), (12|22)
  % (mu nu | la si) = psi_mu(1) psi_nu(1) psi_la(2) psi_si(2)
	for k=1:n
	  % ALL WRONG!
		p1 = rand(1,3)*2*bound - bound;
		p2 = rand(1,3)*2*bound - bound;
		fci(1,1,1,1) = fci(1,1,1,1) + phi(p1,R1,C)^2*phi(p2,R1,C)^2/norm(p1-p2);
		fci(1,1,2,2) = fci(1,1,2,2) + phi(p1,R1,C)^2*phi(p2,R2,C)^2/norm(p1-p2);
		fci(1,2,2,1) = fci(1,2,2,1) + phi(p1,R1,C)*phi(p1,R2,C)*phi(p2,R2,C)*phi(p2,R1,C)/norm(p1-p2);
		%fci(1,1,1,2) += phi(p1,R1,C)^2*phi(p2,R1,C)*phi(p2,R2,C)/norm(p1-p2);
		fci(1,2,2,2) = fci(1,2,2,2) + phi(p1,R1,C)*phi(p1,R2,C)*phi(p2,R2,C)^2/norm(p1-p2);
  end
  
  fci(2,2,2,2)=fci(1,1,1,1);
  fci(2,2,1,1)=fci(1,1,2,2);
  fci(1,2,1,2)=fci(1,2,2,1);
  fci(2,1,1,2)=fci(1,2,2,1);
  fci(2,1,2,1)=fci(1,2,2,1);
  
  fci(1,1,1,2)=fci(1,2,2,2);
  fci(2,1,1,1)=fci(1,1,1,2);
  fci(1,2,1,1)=fci(1,1,1,2);
  fci(1,1,2,1)=fci(1,1,1,2);
  fci(2,2,2,1)=fci(1,2,2,2);
  fci(2,2,1,2)=fci(1,2,2,2);
  fci(2,1,2,2)=fci(1,2,2,2);
  
  fci = fci * (2*bound)^6/n;  % MC normalisation
    
end



