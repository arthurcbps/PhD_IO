% Simulate DGP


T = 10;
I = 1000;

rng(2)
% simulate productivity 50 periods past just to get stationarity
omega_past =  zeros(I, 50);

for t = 2:50
  omega_past(:, t) = 0.9*omega_past(:, t-1) + normrnd(0, 1, I, 1);
end


% productivity terms
xi = normrnd(0, 1, 1000, 10);

omega = zeros(I, T);
omega(:, 1) = 0.9*omega_past(:, 50) + xi(:, 1);

for t = 2:T
  omega(:, t) = 0.9*omega(:, t-1) + xi(:, t);
end

epsilon = normrnd(0, 3, I, T);
u = normrnd(0, 0.2, I, T);

%# capital and labor are determined one period in advance

L = zeros(I, T);
L(:, 1) =  (0.108)^3*(0.73)^(0.8)*exp(3.69 + 3.6*omega_past(:,50));
L(:, 2:10) = (0.108)^3*(0.73)^(0.8)*exp(3.69 + 3.6*omega(:,1:9));

K = 0.73*L;

% materials are determined on the same period, and can be conditioned on epsilon

M = (0.4)^(5/3)*exp(8/75)*(exp(0.3*epsilon + omega).*L.^(0.3).*K.^(0.2)).^(4/3);

% Realized revenue depends on the actual realizations of u

R = (exp(u + omega + 0.3*epsilon).*L.^(0.3).*K.^(0.2).*M.^(0.5)).^(0.8); 


% can back out quantities using elasticity of demand:

eta = -5;
theta = 1/eta + 1;

Y = R.^(1/theta);


y = log(Y);
k = log(K);
l = log(L);
m = log(M);

inputX = cat(3,k,l,m);
inputY = y;
ndegree = 3;

[inputN,inpuT,nvar] = size(inputX);


    
    Phi = xx_flexible_poly(inputY,inputX,ndegree);



% Reshaping

stackR = reshape(R, I*T, 1, []);
stackK = reshape(K, I*T, 1, []);
stackL = reshape(L, I*T, 1, []);
stackM = reshape(M, I*T, 1, []);


lagR = 


