function [estimate_next, covariance_sqrt] = cdekf_update_phase(t_val,P_p,x_p,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,t_G,t_B,t_wind,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num)
[h_x_p,y,R,J_h]=observation_function(x_p,t_val,t_G,t_B,t_wind,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num);
[~, ~, ~, rhoG_val] = atmosisa(double(x_p(3)));
[~, ~, ~, rhoB_val] = atmosisa(double(x_p(6)));
K = P_p * J_h(x_p,rhoG_val,rhoB_val)' / (J_h(x_p,rhoG_val,rhoB_val) * P_p * J_h(x_p,rhoG_val,rhoB_val)' + R);
estimate_next = x_p + K * (y - h_x_p );
covariance_sqrt = P_p - K * J_h(x_p,rhoG_val,rhoB_val) * P_p;
end