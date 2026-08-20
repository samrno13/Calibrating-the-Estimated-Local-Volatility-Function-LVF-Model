%%% Samuel Reano
%%% Calibrating the Estimated Local Volatility Function (LVF) Model from
%%% Market Prices of the European Call Option Project

%%% This is for the Calibration problem part of the project.

function [F,J] = myfun(x)

T = 0.25; % time expiry in 3 months
K = [0.80 0.85 0.90 0.95 1.00 1.05 1.10]; % a vector of strike prices
V0 = [0.3570 0.2792 0.2146 0.1747 0.1425 0.1206 0.0676]; % a vector of initial market European call option prices

%Vector of Market Priced European Call Options
V_MC = Eur_Call_LVF_MC(1,K,0.25,0.03,x,10000,100);
F = V_MC - V0; % vector of model option value errors from the market prices (objective function values at x)

%Applying the finite difference approximation
dx = 0.1; %delta x can be any small number like -0.1, 0.01, -0.01, etc.
x1 = [x(1)+dx,x(2),x(3)];
x2 = [x(1),x(2)+dx,x(3)];
x3 = [x(1),x(2),x(3)+dx];

%Column vectors for the Jacobian matrix using the Monte Carlo Method
V_MCx1 = Eur_Call_LVF_MC(1,K,0.25,0.03,x1,10000,100);
Fx1 = (V_MCx1 - V0)/dx; % output is the first row vector
fdx1 = transpose(Fx1); % output is the first column vector 

V_MCx2 = Eur_Call_LVF_MC(1,K,0.25,0.03,x2,10000,100);
Fx2 = (V_MCx2 - V0)/dx; % output is the second row vector
fdx2 = transpose(Fx2); % output is the second column vector

V_MCx3 = Eur_Call_LVF_MC(1,K,0.25,0.03,x3,10000,100);
Fx3 = (V_MCx3 - V0)/dx; % output is the third row vector  
fdx3 = transpose(Fx3); % output is the third column vector

if nargout > 1 % Two output arguments
    J = [fdx1,fdx2,fdx3]; % Jacobian of the function evaluated at x
end

