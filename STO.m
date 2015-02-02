% Slater-type orbitals
% r,r0 are row vectors
function phi=STO(r,r0,C)
  phi = sqrt(C^3/pi)*exp(-C*norm(r-r0));  
end
