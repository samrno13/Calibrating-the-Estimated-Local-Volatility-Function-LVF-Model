%%% Samuel Reano (500973322)
%%% Calibrating the Estimated Local Volatility Function (LVF) Model from
%%% Market Prices of the European Call Option Project

%%% The following chunk of code solves the partial differential equation 
%%% version of the LVF model by the explicit finite difference method. 

function V0 = Eur_Call_LVF_FD(S0,K,T,r,x,Smax,M,N)
%
% Price the European call option of the LVF model by the explicit finite
% difference method.
%
% Input
%   S0 - initial stock price
%   K - strike price
%   T - maturity
%   r - risk free interest rate
%   x - vector parameters for the LVF sigma, [x(1),x(2),x(3)]
%   Smax - upper bound of the stock price
%   M - number of stock price difference, i.e., dS = Smax/M
%   N - number of time steps, i.e., dtau = T/N
%
% Output
%   V0 - European call option price at t = 0 and S0

%Construct a grid of S and tau
sigma = max(0, x(1) + x(2)*S0 + x(3)*(S0^2));
dS = Smax/M;
dt = T/N;

%Initialize V, alpha and beta
V = zeros(N,M);
alpha = zeros(1,M);
beta = zeros(1,M);

%Subinterval ends of S
S = 0:dS:Smax;

%Initial condition of the PDE
for i = 1:M
    V(1,i) = max(S(i)-K, 0);
end

%Boundary condition on Smax
V(1,M) = Smax;

%Explicit method
for n=1:N-1
    V(n+1,1) = V(n,1)*((1-r)*dt); % boundary condition on S=0
    V(n+1,M) = V(n,M); % boundary condition on Smax

    for i = 2:M-1

        %Central differencing
        alpha_central = ((sigma^2)*(S(i)^2))/(((S(i)-S(i-1))*(S(i+1)-S(i-1)))) - (r*S(i))/(S(i+1)-S(i-1));
        beta_central = ((sigma^2)*(S(i)^2))/(((S(i+1)-S(i))*(S(i+1)-S(i-1)))) + (r*S(i))/(S(i+1)-S(i-1));
        if (alpha_central >= 0) && (beta_central >= 0)
            alpha(i) = alpha_central;
            beta(i) = beta_central;
        else

            %Forward differencing
            alpha_forward = ((sigma^2)*(S(i)^2))/(((S(i)-S(i-1))*(S(i+1)-S(i-1))));
            beta_forward = ((sigma^2)*(S(i)^2))/(((S(i+1)-S(i))*(S(i+1)-S(i-1)))) + (r*S(i))/(S(i+1)-S(i-1));
            alpha(i) = alpha_forward;
            beta(i) = beta_forward;
        end

        %Update V(n+1) by V(n)
        V(n+1,i) = V(n,i)*(1-(alpha(i)+beta(i)+r)*dt) + alpha(i)*dt*V(n,i-1) + beta(i)*dt*V(n,i+1);
    end
end

% Evaluate the current option price V0
% Find the smallest interval including S0
indx1 = max(find(S<=S0));
indx2 = min(find(S>=S0));

if indx1 == indx2
    V0 = V(N,indx1);
else
    %Estimate V0 by the linear interpolation
    w = (S0-S(indx1))/(S(indx2)-S(indx1));
    V0 = V(N,indx1)*w + (1-w)*V(N,indx2);
end


%% An explanation

%%% Lines 47 through 72 is the implementation of the explicit finite 
%%% difference method. Upstream weighting is shown on lines 54 through 67. 
%%% The MATLAB function, find, searches indices and values of nonzero 
%%% elements (hence the name of the variables, indx1 and indx2). In 
%%% line 76, it searches the elements in the S vector from line 37 and 
%%% compares which element(s) is(are) less than or equal to the initial 
%%% stock price, S0, and finds the highest number among the comparison 
%%% (same idea for line 77). Implementing the initial and boundary 
%%% conditions are shown on lines 39 through 50. In line 83, w is the 
%%% slope of the linear interpolation equation ( where  is the slope). 