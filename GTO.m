% Gaussian-type orbitals
% r,r0 are row vectors
function phi=GTO(r,r0,A)
  phi = (2*A/pi)^0.75*exp(-A*(r-r0)*(r-r0)');  
end
