function [X_f,X,t]=sovrapposition_OU_2D_euler ( num_comp, mu, sigma, x0,y0, tmax, n_steps, weight, alpha, omega)

%*****************************************************************************
%
%% ornstein_uhlenbeck_EULER applies the Euler method to the Ornstein-Uhlenbeck SDE.
%
%  Discussion:
%
%    The stochastic differential equation (SDE) is:
%
%      dx(t) = -alpha * x(t)  dt -omega * y(t)  dt + sigma dW,
%      dy(t) = omega *x(t) dt - alpha * y(t) dt

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

  dW= sqrt ( dt ) * randn (num_comp, n_steps );
%
%  Carry out the Euler approximate integration process.
%
 
  X = zeros (num_comp, n_steps + 1 );
  Y = zeros ( num_comp, n_steps + 1 );

  X(:,1) = x0;
  Y(:,1) = y0;
for i = 1:num_comp
    for j = 1:n_steps
       X(i,j+1) = -alpha(i) * X(i,j) * dt - omega(i) * Y(i,j) * dt + sigma(i) * dW(i,j)+X(i,j);
       Y(i,j+1) = omega(i) * X(i,j) * dt - alpha(i) * Y(i,j) * dt + Y(i,j);
    end
end
X_f = sum(X .* sqrt(weight), 1);
 t = linspace ( 0, tmax, n_steps + 1 );

%
%  Plot the approximate solution.
%


  return
end
