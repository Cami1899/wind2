function [x,t]=OU1D_euler( alpha, mu, sigma, x0, tmax, n_steps, weight )

%*****************************************************************************
%
%% ornstein_uhlenbeck_EULER applies the Euler method to the Ornstein-Uhlenbeck SDE.
%
%  Discussion:
%
%    The stochastic differential equation (SDE) is:
%
%      dx(t) = alpha * ( mu - x(t) ) dt + sigma dW,   
%      x(0) = x0.
%
%    The discretized Brownian path uses a constant stepsize.
%
%    For an SDE of the form:
%
%      dx = f(x(t)) dt + g(x(t)) dW(t),
%
%    the Euler method has the form:
%
%      x(j) = x(j-1) + f(x(j-1)) * dt + g(x(j-1)) * dW(j-1)
%
%    Note that if SIGMA is zero, the problem becomes deterministic,
%    with solution:
%
%      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
%
% 
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Desmond Higham,
%    An Algorithmic Introduction to Numerical Simulation of
%    Stochastic Differential Equations,
%    SIAM Review,
%    Volume 43, Number 3, September 2001, pages 525-546
%
%  Input:
%
%    real THETA, MU, SIGMA, the value of problem parameters.
%
%    real X0, the initial condition.  When studying many
%    realizations of this problem, it is usual for X0 to be chosen
%    from a normal distribution.
%
%    real TMAX, the final time.
%
%    integer N, the number of time steps.
%
  
%******************************************************************************
%  Set the discrete time stepsize.
%
  dt = tmax / n_steps;
%
%  Compute the Brownian increments.
% 
seed = 42; 
rng(seed);
  dw = sqrt ( dt ) * randn ( length(alpha), n_steps );
%
%  Carry out the Euler approximate integration process.
%
 
  x = zeros ( length(alpha), n_steps + 1 );

  x(:,1) = x0;

for i = 1:length(alpha)
    for j = 1:n_steps
        x(i, j+1) = x(i, j) + dt * alpha(i) * (mu(i) - x(i, j)) + sigma(i) * dw(i, j);
    end
end
x = sum(x .* sqrt(weight'), 1);
 t = linspace ( 0, tmax, n_steps + 1 );



  
end