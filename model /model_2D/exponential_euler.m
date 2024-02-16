function [x,t]=exponential_euler ( num_comp, x0,y0, T, dt , weight, alpha, omega,sigma)
%*****************************************************************************
% Applies the order 1 Exponential Euler method to the Ornstein-Uhlenbeck SDE.
% The stochastic differential equation (SDE) is:
%
%      dx(t) = -alpha * x(t)  dt -omega * y(t)  dt + sigma dW,
%      dy(t) = omega *x(t) dt - alpha * y(t) dt

%      x(0) = x0
%      y(0) = y0
%
%The discretized Brownian path uses a constant stepsize:    dW=sqrt(dt)*rand_number
%
%The method has the form:
%
%      s(j) = s(j-1) + f(s(j-1)) * dt + g(s(j-1)) * dW(j-1) with s,f,g,W vectors

%  For x equation:
%    f(s(j-1))= - alpha * x - omega * y
%    g(s(j-1))= sigma
%      
%  For y equation:
%    f(s(j-1))= omega * x - alpha * y
%    g(s(j-1))= 0
%     
% Each component of the superposition of SDE is individually integarted and
% summed according to the weights only at the end.


%  INPUT:
%   - num_comp, number of compnents of the fit; 
%   - x0, x initial point;
%   - y0, y initial point;
%   - T, time interval;
%   - dt, integration steps;
%   - weight, alpha, omega, sigma,  SDE coefficients.
%
% OUTPUT:
%   - x, wind signal values;
%   - t, time values of the signal.
  
%******************************************************************************
%Number of steps
n_steps = T / dt;

%Initial points
X = zeros (num_comp, n_steps + 1 );
Y = zeros (num_comp, n_steps + 1 );

X(:,1) = x0;
Y(:,1) = y0;

%Integration cycle
    for j = 1:n_steps
       X(:,j+1) = X(:,j) -  alpha .* X(:,j) * dt -  omega .* Y(:,j) * dt + (exp(-alpha*dt) + exp(-omega*dt)) .* sigma * sqrt ( dt ) .* randn ( num_comp,1 );
       Y(:,j+1) = Y(:,j) +  omega .* X(:,j) * dt -  alpha .* Y(:,j) * dt ;
    end

%Final weighted sum
x=sum(sqrt(weight) .* X,1);
y=sum(sqrt(weight) .* Y,1);
t = linspace ( 0, T, n_steps + 1 );

end
