function [x,t]=heun( num_comp, sigma, x0,y0, T, dt, weight, alpha, omega)
% Applies the Heun method to the Ornstein-Uhlenbeck SDE.
% The stochastic differential equation (SDE) is:
%
%      dx(t) = -alpha * x(t) * dt -omega * y(t) * dt + sigma * dW
%      dy(t) = omega * x(t) * dt - alpha * y(t) * dt

%      x(0) = x0
%      y(0) = y0
%
%The discretized Brownian path uses a constant stepsize:    dW=sqrt(dt)*rand_number
%
%The method has the form:
%
%      s(j)_p = s(j-1) + f(s(j-1)) * dt + g(s(j-1)) * dW(j-1) 
%      s(j) = s(j-1) + dt/2*  [f(s(j-1)) * dt + g(s(j-1)) * dW(j-1) + f(s(j-1)_p) * dt + g(s(j-1)_p) * dW(j-1) ]
%

%  For x equation:
%    f()= - alpha * x - omega * y
%    g()= sigma
%      
%  For y equation:
%    f()= omega * x - alpha * y
%    g()= 0
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
num_steps=T/dt;

%Initial points
x = zeros(num_comp,num_steps + 1);
y = zeros(num_comp,num_steps + 1);
t = linspace(0, T, num_steps + 1);

x(:,1) = x0;
y(:,1) = y0;

%Integration cycle
for i = 1:num_steps
    % Initial slopes
    k1x = - alpha .* x(:,i) * dt - omega .* y(:,i) * dt+ sigma * sqrt(dt) .* randn(num_comp,1);
    k1y = - alpha .* y(:,i) * dt + omega .* x(:,i) * dt;

    %Predictor step
    x_temp = x(:,i) + k1x;
    y_temp = y(:,i) + k1y;

    % Final slopes
    k2x = - alpha .* x_temp * dt - omega .* y_temp * dt+ sigma * sqrt(dt) .* randn(num_comp,1);
    k2y = - alpha .* y_temp * dt + omega .* x_temp * dt;

    %Corrector step
    x(:,i + 1) = x(:,i) + 0.5 * (k1x + k2x);
    y(:,i + 1) = y(:,i) + 0.5 * (k1y + k2y);
end

%Final weighted sum
x=sum(sqrt(weight) .* x,1);
y=sum(sqrt(weight) .* y,1);

end