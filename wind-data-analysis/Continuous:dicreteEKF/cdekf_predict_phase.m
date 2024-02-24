function [x_p, Pp,P_diag,t] = cdekf_predict_phase(jacobian_func,dt_between_measurements,start_time,P_update,x_0,Q,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,h)

    state_count = length(x_0);
    [~, ~, ~, rhoG_val] = atmosisa(double(x_0(3)));
    [~, ~, ~, rhoB_val] = atmosisa(double(x_0(6)));

	finish_time = start_time + dt_between_measurements;

   
    [t,x_p]=ode45(@(t,x) model_function(x,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma),[start_time:h:finish_time],x_0);
    x_p=x_p';

%     x_p_end=x_p(:,end);

%     jacobian_func(x_0,rhoG_val,rhoB_val)

dot_SE_covariance= @(x,Pp_old,rhoG,rhoB) jacobian_func(x,rhoG,rhoB) * reshape(Pp_old,[state_count,state_count]) +  reshape(Pp_old,[state_count,state_count]) * jacobian_func(x,rhoG,rhoB)'+ Q;
%     dot_SE_covariance= @(x,Pp_old,rhoG,rhoB) jacobian_func(x,rhoG,rhoB) * reshape(Pp_old,[state_count,state_count]) * jacobian_func(x,rhoG,rhoB)'+ Q;
    [t,Pp_vector]=ode113 (@(t,Pp_vector) reshape(dot_SE_covariance(x_0,Pp_vector,rhoG_val,rhoB_val),[state_count^2,1]),[start_time:h:finish_time],reshape(P_update,[state_count^2,1]));
    Pp=reshape(Pp_vector(end,:),size(P_update));
	P_diag=Pp_vector(:,1:state_count+1:end)';
end
