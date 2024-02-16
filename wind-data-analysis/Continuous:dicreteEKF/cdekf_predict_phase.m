function [estimate, SE_covariance,P_diag,t] = cdekf_predict_phase(jacobian_func,dt_between_measurements,start_time,P_update,x_0,Q,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,h)

    state_count = length(x_0);
    [TT, a, pp, rhoG_val] = atmosisa(double(x_0(3)));
    [TT, a, pp, rhoB_val] = atmosisa(double(x_0(6)));

	finish_time = start_time + dt_between_measurements;

    
    [t,estimate]=ode45(@(t,x) model_function(t,x,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma),[start_time:h:finish_time],x_0);
    estimate=estimate';

    x_p=estimate(:,end);

    
	dot_SE_covariance= @(x,Pp_old,rhoG,rhoB) jacobian_func(x,rhoG,rhoB) * reshape(Pp_old,[state_count,state_count]) + jacobian_func(x,rhoG,rhoB) * reshape(P_update,[state_count,state_count])'+ Q;
    [t,Pp_vector]=ode45 (@(t,Pp_vector) reshape(dot_SE_covariance(x_p,Pp_vector,rhoG_val,rhoB_val),[state_count^2,1]),[start_time:h:finish_time],reshape(P_update,[state_count^2,1]));
    SE_covariance=reshape(Pp_vector(end,:),size(P_update));
	P_diag=Pp_vector(:,1:state_count+1:end)';
end
