function [x,t]=taylor1_5( num_comp, x0,y0, T, dt , weight, alpha, omega,sigma)
%*****************************************************************************
% Applies the 1.5 order Taylor method to the Ornstein-Uhlenbeck SDE.
% The stochastic differential equation (SDE) is:
%
%      dx(t) = -alpha * x(t) * dt -omega * y(t) * dt + sigma * dW
%      dy(t) = omega * x(t) * dt - alpha * y(t) * dt
%
%      x(0) = x0
%      y(0) = y0
%
% The sysytem can be written
%
%          dS(t)= A * s(t) * dt + B * dW
% with 
%   S=[x,y]
%   A=[-alpha, -omega; omega, -alpha]
%   B=[sigma; 0]

%The discretized Brownian path uses a constant stepsize:    dW=sqrt(dt)*rand_number
%
%The method has the form:
%
%      s(j) = s(j-1) + dt * A * s(j-1) + B * dW + 1/2 * dt^2 * A^2 + B * dW + A * B * dt/2 * (dW + sqrt(dt/6) * R_2)

%   with dW=sqrt(dt)*R1
%        R=random number

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
       R1=randn ( num_comp,1 );
       R2=randn ( num_comp,1 );
       X(:,j+1) = X(:,j) -  alpha .* X(:,j) * dt -  omega .* Y(:,j) * dt + sigma * sqrt ( dt ) .* R1 + 0.5 * dt^2 * alpha.^2 .* X(:,j) + 0.5 * dt^2 * omega.^2 .* Y(:,j) + sigma * sqrt ( dt ) .* R1 - alpha .*sigma * dt/2 .* (sqrt(dt) * R1 + sqrt(dt/6) * R2)- omega .* sigma * dt/2 .* (sqrt(dt)* R1 + sqrt(dt/6) * R2);
       Y(:,j+1) = Y(:,j) +  omega .* X(:,j) * dt -  alpha .* Y(:,j) * dt + 0.5 * dt^2 *omega.^2 .* X(:,j) + 0.5 * dt^2 * alpha.^2 .* Y(:,j);
    end

%Final weighted sum
x=sum(sqrt(weight) .* X,1);
y=sum(sqrt(weight) .* Y,1);
t = linspace ( 0, T, n_steps + 1 );

end
