function [estimate_next, covariance_sqrt] = cdekf_update_phase(R,P,C,estimate,measurement)

K = P * C' / (C * P * C' + R);
estimate_next = estimate + K * (measurement - C * estimate);
covariance_sqrt = (eye(size(P)) - K * C) * P;
end