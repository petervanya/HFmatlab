% construct matrix G
% in the future this will reguire a far mor efficient routine
function G = get_G(P,fci)
  N=2;
  G=zeros(2);
  for mu=1:N
    for nu=1:N
      for la=1:N
        for si=1:N
          G(mu,nu) = G(mu,nu) + P(la,si)*(fci(mu,nu,la,si)-fci(mu,la,si,nu)/2);
        end
      end
    end
  end
end
