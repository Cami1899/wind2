clc
clear

%Import data from flight test
load('wind_value_restricted.mat');
load('time_restricted.mat');

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

mean_value=mean(training_ws);
%*************************************************
% STEP 1: Creata a superposition of SDE with the desired ACF of data
%1.1 Calculate the ACF of data
maxLag=100;
acf=autocorr(training_ws,'NumLag', maxLag);
lag=(0:maxLag)';
% figure;
% stem(acf, 'LineWidth', 1.5);
% title('Auto-correlation function','FontSize', 16);
% xlabel('Lag [s]');
% ylabel('Auto-correlation of wind speed data');
% grid on;

%1.2 Fit an exponential function to data ACF
[fitresult, gof] = fit_acf_2exp(lag, acf); %fit the acf with a duble sum of exponential functions

% %fit the ACF of data to a weighted sum of a dumped sinusoidal function
% (this is for 2D-O.-U. process)

% numTerms = 1;  
% 
% params0 = ones(1, 3 * numTerms);
% options = optimoptions(@lsqcurvefit, 'MaxIterations', 10000, 'MaxFunctionEvaluations', 10000);
% params = lsqcurvefit(@weightedExpCosSumModel, params0, lag, acf, [], [], options);
% 
% 
% yFit = weightedExpCosSumModel(lag, params);
% 

% figure;
% scatter(lag, acf, 'DisplayName', 'Data');
% hold on;
% plot(xData, yFit, 'r', 'DisplayName', 'Fit');
% legend();
% xlabel('X');
% ylabel('Y');
% title('The ACF of Data Set  and the fitted sum.');



eq = formula(fitresult); 
values = coeffvalues(fitresult);
weight=values(1:2:end);
alpha=abs(values(2:2:end));
%mu=mean*ones(length(alpha));
mu=zeros(length(alpha));
sigma=sqrt(2*alpha);
v0=mean_value;
tmax=10000;
n_steps=100000;
[v_acf,t]=sovrapposition_ornstein_uhlenbeck1D_euler( alpha, mu, sigma, v0, tmax, n_steps, weight );


%TEST the results
%Test the gaussian PDF of SDE

figure;
histogram(v_acf, 'Normalization', 'pdf');
hold on;
pd = fitdist(v_acf', 'Normal');
x_values = linspace(min(v_acf), max(v_acf), 1000);
y_values = pdf(pd, x_values);
plot(x_values, y_values, 'LineWidth', 2);
title('Confronto Istogramma con PDF Gaussiana');
legend('Istogramma', 'PDF Gaussiana');

% h=lillietest(v_acf); %lillie test method
% if h == 1
%     disp('The distribution of SDE is not Gaussian');
%     return;
% end

acf_step1=autocorr(v_acf,'NumLag', maxLag);

% figure
% subplot(2,1,1)
% stem(acf_step1, 'LineWidth', 1.5);
% title('Auto-correlation function','FontSize', 16);
% xlabel('Lag [s]');
% ylabel('Auto-correlation of generated wind speed at the end of dtep1');
% grid on;
% subplot(2, 1, 2);
% stem(acf, 'LineWidth', 1.5);
% title('Auto-correlation function','FontSize', 16);
% xlabel('Lag [s]');
% ylabel('Auto-correlation of wind speed data');
% grid on;

%% 

% ****************************************************

%STEP 2: impose the desired PDF

%2.1 Calculate the CDF of data using interpolation of ECDF
% n_ecdf=length(training_ws);
% v_sort = sort(training_ws);
% [Fi,xi] = ecdf(v_sort);


% figure
% stairs(xi,Fi,'r');
% xlim([2.000 2.001]); xlabel('x'); ylabel('F(x)'); %zoom to see the steps
% 
% xj = xi(2:end);
% Fj = (Fi(1:end-1)+Fi(2:end))/2;

% figure
% stairs(xi,Fi,'r');
% hold on
% plot(xj,Fj,'b.', xj,Fj,'b-');
% xlim([2.000 2.001]); xlabel('x'); ylabel('F(x)'); %zoom to see the steps
% hold off
% legend({'ECDF' 'Breakpoints' 'Piecewise Linear Estimate'},'location','NW');

% figure
% stairs(xi,Fi,'r');
% hold on
% plot(xj,Fj,'b.', xj,Fj,'b-');
%  xlabel('x'); ylabel('F(x)'); 
% hold off
% legend({'ECDF' 'Breakpoints' 'Piecewise Linear Estimate'},'location','NW');

%Since the smallest data point corresponds to a height of 1/2n, and the largest to 1-1/2n, the first and last linear segments must be extended beyond the data, to make the function reach 0 and 1
% xj = [xj(1)-Fj(1)*(xj(2)-xj(1))/((Fj(2)-Fj(1)));
%       xj;
%       xj(n_ecdf)+(1-Fj(n_ecdf))*((xj(n_ecdf)-xj(n_ecdf-1))/(Fj(n_ecdf)-Fj(n_ecdf-1)))];
% Fj = [0; Fj; 1];
% % hold on
% plot(xj,Fj,'b-','HandleVisibility','off');
% hold off

% cdf = @(y) interp1(xj,Fj,y,'linear','extrap');
% y = linspace(min(xj),max(xj),n_ecdf/1000);
% u =sort(min(training_ws)+(max(training_ws)-min(training_ws)) *rand(1,4000));
% figure
% stairs(xj,Fj,'r');
% hold on
% plot(u,cdf(u),'k-');
% xlabel('x'); ylabel('F(x)');
% legend({'ECDF'  'CDF'},'location','NW');
% title('ECDFvsCDF','FontSize', 16)
%(https://www.mathworks.com/help/stats/nonparametric-estimates-of-cumulative-distribution-functions-and-their-inverses.html?searchHighlight=cumulative%20distribution%20function%20&s_tid=srchtitle_support_results_1_cumulative%20distribution%20function%20)

%calculte the inverse CDF
% figure
% stairs(Fi,[xi(2:end); xi(end)],'r');
% hold on
% plot(Fj,xj,'b-');
% hold off
% ylabel('x'); xlabel('F(x)');
% legend({'ECDF' 'Piecewise Linear Estimate'},'location','NW');
% title('Inverse CDF','FontSize', 16)
% g = rand(1, n_ecdf);
% Finv = @(l) interp1(Fj,xj,g,'linear','extrap');

%alternative method
[ecdf,x_ecdf]=ecdf(training_ws); %Use the Kaplan-Meier estimator to estiamte ecdf
inv_ECDF=@(z) interp1(ecdf,x_ecdf,z,'linear','extrap'); 

%SDE gaussian CDF
gaussianCDF=normcdf(v_acf);
v_pdf=inv_ECDF(gaussianCDF);



% wind_speed=matlabFunction(v_final,'File', 'wind_model');

figure
plot(t,v_pdf, 'r-', 'LineWidth', 3 )
xlabel ( 't [s]', 'FontSize', 16 )
ylabel ( 'v [m/s]', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
title ( 'SDE withe the desired PDF and ACF', 'FontSize', 16 )
grid ( 'on' );

%**************************************************

%TEST the results


figure;
subplot(2, 1, 1);
histogram(v_pdf, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('v [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of generated wind speed', 'FontSize', 16)
grid('on');
subplot(2, 1, 2);
histogram(training_ws, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('v [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of wind speed data', 'FontSize', 16)
grid('on');

acf_final=autocorr(v_pdf,'NumLag', maxLag);

figure
subplot(2,1,1)
stem(acf_final, 'LineWidth', 1.5);
title('Auto-correlation function of generated wind speed','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;
subplot(2, 1, 2);
stem(acf, 'LineWidth', 1.5);
title('Auto-correlation function of wind speed data','FontSize', 16);
xlabel('Lag [s]');
ylabel('ACF');
grid on;





