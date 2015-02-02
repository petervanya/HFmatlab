function rho = density_mat(C)
  rho = zeros(2);
  rho(1,1) = C(1,1)^2;
  rho(2,2) = C(2,1)^2;
  rho(1,2) = C(1,1)*C(2,1);
  rho(2,1) = rho(1,2);
  rho = rho * 2;
end
