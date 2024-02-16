%Weilbull distribution 2-parameter
clc
% new_wind_magnitude=
figure()
% plot(new_time, new_magnitude)

wind_velocity_rescaled = new_magnitude;
pd_w= fitdist(wind_velocity_rescaled, 'Weibull'); %create probability distribution fitting the Weibull distribution to data 
%methods(pd) %show supporteded function
x_w = 0:0.1:max(wind_velocity_rescaled); % defining value steps 
pdf_w = pdf(pd_w, x_w); %probability density function

figure()
histogram(new_magnitude,'BinMethod','scott','Normalization','pdf','FaceColor','c')
hold on
plot(x_w, pdf_w,'r', LineWidth=2); 
hold off
xlabel('wind speed [m/s]');
ylabel('probability density');


params = paramci(pd_w); % Calculate parameters 
k = params(1); % shape parameter 
c = params(2); % scale speed [m/s]
P=cdf('Weibull',x_w,c,k); %cumulative probability
figure()
plot(x_w,P);
hold on
xlabel('wind velocity [m/s]');
ylabel('cumulative probability P(<V9)');
hold off
%% 

%%kernel

[F,xi]=ksdensity(new_magnitude);
figure()
histogram(new_magnitude,'Normalization','pdf','BinMethod','scott','FaceColor','c')
hold on
plot(xi, F,'r', LineWidth=2); 
hold off
xlabel('wind speed [m/s]');
ylabel('probability density');

pd_k=fitdist(wind_velocity_rescaled,'Kernel', 'Kernel','normal','Width',0.04);
x_k=min(wind_velocity_rescaled):0.01:max(wind_velocity_rescaled);
pdf_k=pdf(pd_k,x_k);

figure()
histogram(new_magnitude,'Normalization','pdf','BinMethod','scott','FaceColor','c')
hold on
plot(x_k, pdf_k,'r', LineWidth=2); 
hold off
xlabel('Wind speed [m/s]');
ylabel('Probability density');
legend('Data', 'Kernel model')
title('Probability density function of wind speed and Kernel model')

% num_samples = 1000; % Numero di campioni casuali desiderati
% wind_random = random('Weibull', c, k, [1, num_samples]);
% 
% figure()
% plot([1:num_samples],wind_random)
% % wind_random=random('Weibull',c,k)
% figure()
% histfit(x) %Create a histogram
% qqplot(x,pd) %create a quantile-quantile plot of the quantiles of the sample data x versus the theoretical quantile values of the fitted distribution.