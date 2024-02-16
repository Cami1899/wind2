function [estimate_next, covariance_sqrt] = cdekf_update_phase(t_val,J_h_num,P_p,x_p,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,t_G,t_B,t_wind,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance)
[h_x_p,y,R]=observation_function(x_p,t_val,t_G,t_B,t_wind,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM)
[TT, a, pp, rhoG_val] = atmosisa(double(x_p(3)));
[TT, a, pp, rhoB_val] = atmosisa(double(x_p(6)));
K = P_p * J_h_num(x_p,rhoG_val,rhoB_val)' / (J_h_num(x_p,rhoG_val,rhoB_val) * P_p * J_h_num(x_p,rhoG_val,rhoB_val)' + R);
estimate_next = x_p + K * (y - h_x_p );
covariance_sqrt = (eye(size(P_p)) - K * J_h_num(x_p,rhoG_val,rhoB_val)) * P_p;
end