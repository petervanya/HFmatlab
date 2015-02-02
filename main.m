% HF approximation to H2 molecule
% 2 nuclei, 2 electrons
% 19/11/14
% ========
% PROCEDURE
% 1. load file with coordinates of nuclei, create the basis set
% 2. calculate S_ij, H^core_ij and (mu nu | lambda sigma)
% 3. diagonalise S and get X
% 4. guess density matrix P, start with zero
% 5. calculate G from 
% 6. create Fock operator F = H_core + G
% 7. calc F'=X^T*F*X
% 8. diagonalise F', get C' and e
% 9. calc C = XC'
% 10. get a new P, check if converged, repeat
% 
% constants
zeta = 1.24; % Slater orbital exponent
Rnucl = 1.4;

% 1. step
% =======
% read in r1...rN from file (this is provisional)
R1 = [0 0 0];
R2 = [0 0 Rnucl];
xyz = [R1; R2];

% form the basis
%basis = [1  1; 1 -1];
phi = @(r,R,zeta) STO(r,R,zeta);

t0=tic;
% Get core Hamiltonian matrix
tic
Hc=H_core(R1,R2,phi,zeta)
%Hc = [-1.10016  -0.94869; -0.94869  -1.10016];
toc

% get four-centre integrals
tic
mnls=four_centre_int(R1,R2,phi,zeta)
%mnls=zeros(2,2,2,2);
%mnls(:,:,1,1) = [11.731   10.794; 10.794   10.794];
%mnls(:,:,2,1) = [10.79384    0.00000; 0.00000   10.79384];
%mnls(:,:,1,2) = [10.79384    0.00000; 0.00000   10.79384];
%mnls(:,:,2,2) = [10.794   10.794; 10.794   10.794];
toc

% get overlap and transformation matrix
tic
S = overlap_matrix(R1,R2,phi,zeta)
%S = [1.00000  0.67368; 0.67368  1.00000];
X = transf_matrix(S);
toc

P = zeros(2);
G = zeros(2);
G = get_G(P,mnls);
F = Hc + G;
Fp = conj(X')*F*X;
[Cp,E] = eig(Fp);
C = X*Cp;
% new density matrix
Pn = density_mat(C);

fprintf('=== Results:\n');
E
C
P
Pn
fprintf('Total time: %.2f s.\n',toc(t0));




