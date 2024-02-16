clc
clear
close all
%Import data from flight test
load('wind_value_restricted.mat');
load('time_restricted.mat');

% Define your data
wind_speed=new_magnitude;
time=new_time;


%Create the training and testing data sample
percent=70;
number_of_data=round(percent/100*length(wind_speed));
index_training=sort(randperm(length(wind_speed),number_of_data));
training_ws=wind_speed(index_training);
training_time=time(index_training);


testing_wd=setdiff(wind_speed,training_ws);
testing_time=setdiff(time,training_time);


mean_value_ws=mean(training_ws);

%Calculate the ACF of data
maxLag=100;

acf_ws=autocorr(training_ws,'NumLag', maxLag);
t=(0:maxLag)';

% Define the sum of product of exponential and cosine functions
f = @(p, t) p(1) * exp(-p(2) * t) .* cos( p(3) * t)+p(4) * exp(-p(5) * t) .* cos( p(6) * t);

% Set initial guess for parameters
p0 = [0.099266     0.11401     0.99977     0.97131   0.0052164 -1.3472e-07 ];
% Use curve fitting to find the optimal parameters
options = optimoptions('lsqnonlin', 'Display', 'none');
[p, fval] = lsqnonlin(@(p) (acf_ws - f(p, t)).^2, p0, [],[], options);

% Print the optimal parameters
disp(['The optimal parameters are: ', num2str(p)]);

plot(t,f(p,t),t,acf_ws)