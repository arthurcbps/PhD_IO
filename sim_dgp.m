% Simulate DGP

% productivity terms
xi = normrnd(0, 1, 1000, 10);
omega_0 = normrnd(0, 1, 1000, 1);
omega = zeros(1000, 10);
omega(:, 1) = 0.9*omega_0 + xi(:, 1);

for i = 2:10
    omega(:, i) = 0.9*omega(:, i-1) + xi(:, i);
end

epsilon = normrnd(0, 3, 1000, 10);
u = normrnd(0, 0.2, 1000, 10);

% capital and labor are determined one period in advance

L = zeros(1000, 10);
L(:, 1) = (0.108)^3*(0.73)^0.8*exp(4.75 + 3.6*omega_0);
L(:, 2:10) = (0.108)^3*(0.73)^0.8*exp(4.75 + 3.6*omega(:, 1:9));

K = 0.73*L;

% materials are determined on the same period, and can be conditioned on
% epsilon

M = (0.4)^(5/3)*(exp(0.1 + 0.3*epsilon + omega).*L.^(0.3).*K.^(0.2)).^(4/3);

% Realized revenue depend on the actual realizations of u

R = (exp(u + omega + 0.3*epsilon).*L.^(0.3).*K.^(0.2).*M.^(0.5)).^(0.8); 

