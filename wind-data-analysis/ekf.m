function [x]=ekf(n,f,h,x0,P0,Q,R,J_h,J_f,z)


% 0)Initial values:

        % x0 - initial state
        % P0 - initial error covariance
[a,b]=size(z);
x=zeros(n,b);
x(:,1)=x0;
P=P0;

for k=1:b


% 1) Predict state and error covariance
x_p=f(x(:,k));
P_p=J_f*P*J_f'+Q;

% 2) Kalman gain

K_k=P_p*J_h'*(J_h*P_p*J_h'+R)^(-1);

% 3) Estimate state

x(:,k+1)=x_p+K_k*(z(:,k)-h(x_p));

% 4) Error covariance

P=P_p-K_k*J_h*P_p';

end



