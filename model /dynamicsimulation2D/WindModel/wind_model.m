
function [wind_speed_int,t_wind,v_wind]=wind_model(new_magnitude,new_time)
%Create the training and testing data sample
percent=70;
number_of_data=round(percent/100*length(new_magnitude));
index_training=sort(randperm(length(new_magnitude),number_of_data));
training_ws=new_magnitude(index_training);
training_time=new_time(index_training);

testing_wd=setdiff(new_magnitude,training_ws);
testing_time=setdiff(new_time,training_time);

mean_value=mean(training_ws);
%*************************************************
% STEP 1: Creata a superposition of SDE with the desired ACF of data
%1.1 Calculate the ACF of data
maxLag=100;
acf=autocorr(training_ws,'NumLag', maxLag);
lag=(0:maxLag)';


%1.2 Fit an exponential function to data ACF
[fitresult, gof] = fit_acf_2exp(lag, acf); %fit the acf with a duble sum of exponential functions


%1.3 Create and integrate the superposition of SDE 
eq = formula(fitresult); 
values = coeffvalues(fitresult);
weight=values(1:2:end);
alpha=abs(values(2:2:end));
%mu=mean*ones(length(alpha));
mu=zeros(length(alpha));
sigma=sqrt(2*alpha);
v0=mean_value;
% tmax=10000;
n_steps=10000;
tmax=100;
[v_acf,t_wind]=sovrapposition_ornstein_uhlenbeck1D_euler ( alpha, mu, sigma, v0, tmax, n_steps, weight ); %Euler method

%% 
% ****************************************************

%STEP 2: impose the desired PDF

%2.1 Calculate the inverse CDF of data
[ecdf_value,x_ecdf]=ecdf(training_ws); %Use the Kaplan-Meier estimator to estiamte ecdf
inv_ECDF=@(z) interp1(ecdf_value,x_ecdf,z,'linear','extrap'); 

%2.2 SDE gaussian CDF
gaussianCDF=normcdf(v_acf);

%2.3 Memoryless transformation
v_wind=inv_ECDF(gaussianCDF);
wind_speed_int=@(t) interp1(t_wind,v_wind,t,'linear','extrap');

figure(4)
plot(t_wind,v_wind, 'r-', 'LineWidth', 3 )
xlabel ( 't [s]', 'FontSize', 16 )
ylabel ( 'v [m/s]', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
title ( 'SDE withe the desired PDF and ACF', 'FontSize', 16 )
grid ( 'on' );

%**************************************************

%TEST the results


% figure;
% subplot(2, 1, 1);
% histogram(v_pdf, 'Normalization', 'pdf', 'EdgeColor', 'none');
% xlabel('v [m/s]', 'FontSize', 16)
% ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
% title('PDF of generated wind speed', 'FontSize', 16)
% grid('on');
% subplot(2, 1, 2);
% histogram(training_ws, 'Normalization', 'pdf', 'EdgeColor', 'none');
% xlabel('v [m/s]', 'FontSize', 16)
% ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
% title('PDF of wind speed data', 'FontSize', 16)
% grid('on');
% 
% acf_final=autocorr(v_pdf,'NumLag', maxLag);
% 
% figure
% subplot(2,1,1)
% stem(acf_final, 'LineWidth', 1.5);
% title('Auto-correlation function of generated wind speed','FontSize', 16);
% xlabel('Lag [s]','FontSize', 16);
% ylabel('ACF','FontSize', 16);
% grid on;
% subplot(2, 1, 2);
% stem(acf, 'LineWidth', 1.5);
% title('Auto-correlation function of wind speed data','FontSize', 16);
% xlabel('Lag [s]','FontSize', 16);
% ylabel('ACF','FontSize', 16);
% grid on;

end





