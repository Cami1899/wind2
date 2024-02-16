clc
clear
close all
%Import data from flight test
load('wind_value_restricted.mat');
load('time_restricted.mat');
load('angles.mat')

wind_speed=new_magnitude;
time=new_time;
angle=new_angle;

%Create the training and testing data sample
percent=70;
number_of_data=round(percent/100*length(angle));
index_training=sort(randperm(length(angle),number_of_data));
training_ws=wind_speed(index_training);
training_time=time(index_training);
training_angle=angle(index_training);

testing_wd=setdiff(wind_speed,training_ws);
testing_time=setdiff(time,training_time);
testing_angle=setdiff(angle,training_angle);

mean_value_angle=mean(training_angle);
mean_value_ws=mean(training_ws);
%*************************************************
% STEP 1: Creata a superposition of SDE with the desired ACF of data
%1.1 Calculate the ACF of data
maxLag=100;
acf_angle=autocorr(training_angle,'NumLag', maxLag);
acf_ws=autocorr(training_ws,'NumLag', maxLag);
lag=(0:maxLag)';
% figure;
% stem(acf, 'LineWidth', 1.5);
% title('Auto-correlation function','FontSize', 16);
% xlabel('Lag [s]');
% ylabel('Auto-correlation of angles data');
% grid on;

%1.2 Fit an exponential function to data ACF
fitresult_ws = fit(lag,acf_ws,'exp2');
fitresult_angle = fit(lag,acf_angle,'exp2'); %fit the acf with a duble sum of exponential functions



values_angle = coeffvalues(fitresult_angle);
weight_angle=values_angle(1:2:end);
alpha_angle=abs(values_angle(2:2:end));
%mu=mean*ones(length(alpha));
mu_angle=zeros(length(alpha_angle));
sigma_angle=sqrt(2*alpha_angle);
a0=mean_value_angle;
tmax=10000;
n_steps=10000;
[angle_acf,t_angle]=OU1D_euler( alpha_angle, mu_angle, sigma_angle, a0, tmax, n_steps, weight_angle );

values_ws = coeffvalues(fitresult_ws);
weight_ws=values_ws(1:2:end);
alpha_ws=abs(values_ws(2:2:end));
%mu=mean*ones(length(alpha));
mu_ws=zeros(length(alpha_ws));
sigma_ws=sqrt(2*alpha_ws);
v0=mean_value_ws;
[v_acf,t_ws]=OU1D_euler( alpha_ws, mu_ws, sigma_ws, v0, tmax, n_steps, weight_ws );

%TEST the results
%Test the gaussian PDF of SDE

figure;
histogram(angle_acf, 'Normalization', 'pdf');
hold on;
pd = fitdist(angle_acf', 'Normal');
x_values = linspace(min(angle_acf), max(angle_acf), 1000);
y_values = pdf(pd, x_values);
plot(x_values, y_values, 'LineWidth', 2);
title('Angle');
legend('Istogramma', 'PDF Gaussiana');

figure;
histogram(v_acf, 'Normalization', 'pdf');
hold on;
pd = fitdist(v_acf', 'Normal');
x_values = linspace(min(v_acf), max(v_acf), 1000);
y_values = pdf(pd, x_values);
plot(x_values, y_values, 'LineWidth', 2);
title('Wind speed');
legend('Istogramma', 'PDF Gaussiana');

%Test the ACF
angle_acf_step1=autocorr(angle_acf,'NumLag', maxLag);

figure
subplot(2,1,1)
stem(angle_acf_step1, 'LineWidth', 1.5);
title('Angle','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;
subplot(2, 1, 2);
stem(acf_angle, 'LineWidth', 1.5);
title('Angle','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;

v_acf_step1=autocorr(angle_acf,'NumLag', maxLag);

figure
subplot(2,1,1)
stem(v_acf_step1, 'LineWidth', 1.5);
title('Wind speed','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;
subplot(2, 1, 2);
stem(acf_ws, 'LineWidth', 1.5);
title('Wind speed','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;


% ****************************************************

%STEP 2: impose the desired PDF

%2.1 Calculate the CDF of data using interpolation of ECDF

[ecdf_angle,x_ecdf_angle]=ecdf(training_angle); %Use the Kaplan-Meier estimator to estiamte ecdf
inv_ECDF_angle=@(z) interp1(ecdf_angle,x_ecdf_angle,z,'linear','extrap'); 
[ecdf_ws,x_ecdf_ws]=ecdf(training_ws); %Use the Kaplan-Meier estimator to estiamte ecdf
inv_ECDF_ws=@(z) interp1(ecdf_ws,x_ecdf_ws,z,'linear','extrap'); 

%SDE gaussian CDF
gaussianCDF_angle=normcdf(angle_acf);
angle_pdf=inv_ECDF_angle(gaussianCDF_angle); %Final generated directions

gaussianCDF_ws=normcdf(v_acf);
ws_pdf=inv_ECDF_ws(gaussianCDF_ws);%Final generated speeds

figure
polarscatter(angle_pdf,ws_pdf,[],'filled')
 
%**************************************************
%Wind component along axes
v_x=ws_pdf.*cos(angle_pdf);
v_y=ws_pdf.*sin(angle_pdf);

%**************************************************

%TEST the results


figure;
subplot(2, 1, 1);
histogram(angle_pdf, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('angle [rad]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of generated wind direction', 'FontSize', 16)
grid('on');
subplot(2, 1, 2);
histogram(training_angle, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('angle [rad]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of wind direction data', 'FontSize', 16)
grid('on');

figure;
subplot(2, 1, 1);
histogram(ws_pdf, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('speed [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of generated wind speed', 'FontSize', 16)
grid('on');
subplot(2, 1, 2);
histogram(training_ws, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('speed [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of wind speed data', 'FontSize', 16)
grid('on');

acf_final_angle=autocorr(angle_pdf,'NumLag', maxLag);

figure
subplot(2,1,1)
stem(acf_final_angle, 'LineWidth', 1.5);
title('Auto-correlation function of generated wind directions','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;
subplot(2, 1, 2);
stem(acf_angle, 'LineWidth', 1.5);
title('Auto-correlation function of wind direction data','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;

acf_final_ws=autocorr(ws_pdf,'NumLag', maxLag);

figure
subplot(2,1,1)
stem(acf_final_ws, 'LineWidth', 1.5);
title('Auto-correlation function of generated wind speed','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;
subplot(2, 1, 2);
stem(acf_ws, 'LineWidth', 1.5);
title('Auto-correlation function of wind speed data','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;


