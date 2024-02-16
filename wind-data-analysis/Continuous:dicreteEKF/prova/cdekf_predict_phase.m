function [estimate, SE_covariance,P_diag,t] = cdekf_predict_phase(jacobian_func,dt_between_measurements,start_time,P_update,x_0,Q,h)

    state_count = length(x_0);
    

	finish_time = start_time + dt_between_measurements;

    
    [t,estimate]=ode45(@(t, x) [x(2); x(3); 0],[start_time:h:finish_time],x_0);
    estimate=estimate';
%     plot(t_x,estimate(6,:))

    
    x_p=estimate(:,end);

    
	dot_SE_covariance= @(Pp_old) jacobian_func * reshape(Pp_old,[state_count,state_count]) + jacobian_func * reshape(P_update,[state_count,state_count])'+ Q;
    [t,Pp_vector]=ode45 (@(t,Pp_vector) reshape(dot_SE_covariance(Pp_vector),[state_count^2,1]),[start_time:h:finish_time],reshape(P_update,[state_count^2,1]));
    SE_covariance=reshape(Pp_vector(end,:),size(P_update));

    P_diag=Pp_vector(:,1:state_count+1:end)'; 
	
end
