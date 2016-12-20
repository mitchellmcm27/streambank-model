function [Q,U,Cf,tau] = model_average_velocity(B,H,S)
% Sefick et al JAWRA 2015

Abkf = B.*H;
Rh = Abkf./(B + H + H);

Q = Abkf .* Rh.^0.6906 .* S.^0.1216;
U = Q./Abkf;

g = 9.81; % gravity m/s2
rho = 999.9; % density of water at 20 C kg/m3
tau = g*rho*Rh.*S;
shear_velocity = sqrt(tau./rho);

Cf = (shear_velocity./U).^2;

end

