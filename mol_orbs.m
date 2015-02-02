% function to produce molecular orbitals
% from atomic orbitals and nuclear coords
% !!! CURRENTLY WORKS ONLY FOR MINIMAL BASIS H2 !!!
% input:
% -- position r
% -- nuclear coords matrix xyz
% -- function to produce atomic orbitals

% THIS IS A USELESS FUNCTIONS MOS ARE DETERMINED IN HFA
function psi = mol_orbs(r,xyz,AO)
  % N=size(xyz,1); % number of atoms, NOT USED NOW
  phi1=AO(r,xyz(1,:),C);
  phi2=AO(r,xyz(2,:),C);
  
  psi=zeros(2,1);
  psi(1) = phi1 + phi2;
  psi(2) = phi1 - phi2;
  % normalise by MC integration  
  n=10e4;
  sigma=5;
  bound=sigma/C;
  A=zeros(2,1);
  for k=1:n
    p = rand(3,1)*2*bound-bound*ones(3,1);
    A(1) += orb(p,R1,C)*orb(p,R2,C);
  end
  
  psi./=A;  
end

function psi=normaliseMOs(
