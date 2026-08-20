%%% Samuel Reano
%%% Calibrating the Estimated Local Volatility Function (LVF) Model from
%%% Market Prices of the European Call Option Project


%%% The following chunk of code solves the stochastic differential equation 
%%% version of the LVF model by the Monte Carlo method. This chunk of code 
%%% is a function file in MATLAB, which means we can use the created 
%%% function in another MATLAB script file.



function V = Eur_Call_LVF_MC(S0,K,T,r,x,M,N)
%
% Price the European call option of the LVF model by the Monte Carlo method
%
% Input
%   S0 - initial stock price
%   K - strike price
%   T - maturity
%   r - risk free interest rate
%   x - vector parameters for the LVF sigma, [x(1),x(2),x(3)]
%   M - number of simulated paths
%   N - number of time steps, i.e., deltaT = T/N
%
% Output
%   V - European call option price at t = 0 and S0

deltaT = T/N;

%Monte Carlo method
for i=1:M
    A(1)=0;
    for j=1:N
        S(1)=S0;
        sigma = max(0.0,x(1)+x(2)*S(j)+x(3)*(S(j))^2);
        S(j+1) = S(j)*exp(r-((sigma^2)/2)*deltaT+sigma*sqrt(deltaT)*randn);
        A(j+1) = (j/(j+1))*A(j)+(1/(j+1)*S(j+1));
    end
payoff = max(A(j)-K,0);
end

%European call option price at t=0
V = exp(-r*T)*(1/M)*sum(payoff);

%% An explanation 

%%% A(j) is the Monte Carlo simulation of the stock price. A(j) is an 
%%% approximation of the stock price when implementing the Monte Carlo 
%%% simulation. Line 37 shows the implementation of Ito's Lemma. Line 38 
%%% shows path dependent option pricing.