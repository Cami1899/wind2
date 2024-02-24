clc
clear
close 
%------------------------------------------------------------------------
%                     WIND SPEED MODEL WITH ARBITRARY PDF and ACF

%The model generates a wind signal that has the same auto-correaltion function (ACF) and probability density function (PDF) as real wind data. 
%The model is based on a weighted superposition of 2-dimensional Ornstein-Uhlenbeck processes, in the form:

%       dx(t) = weight * sum( -alpha * x(t)  dt -omega * y(t)  dt + sigma dW)
%       dy(t) = weight * sum( omega *x(t) dt - alpha * y(t) dt)

%whose ACF is in the form

%       R = sum( w*exp( -alfa * t ) * cos( omega * t ) )

%Weight, alpha, omega and sigma are the model coefficents.

%Only the component in x is evaluated as the wind signal

%The model is composed by 2 steps:
%1)The coefficents of the model are obtained fitting the ACF of real data with a superposition of the model ACF.
%  How the fit is realized is explained in 'function_fitACF'
%  The SDE whith the desired coefficient is integrated to obtain the wind signal.
%  4 different integration methods are used: Euler-Maruyama, Exponential Euler,order-1.5 Taylor and Heun method

%2)A memoryless transformation to impose the real data PDF to the generatedsignal is applied:
%  the Gaussian CDF ɸ  is applied to the inverse of CDF of real data F^-1
 
%------------------------------------------------------------------------

% %Import data from flight test
% 
% load('wind_matrix_finale1.mat')
% load('wind_matrix_finale2.mat')
% 
%         % %Flight1 
%         
%         wind_speed1=wind_matrix_final_1(:,2);
%         time1=wind_matrix_final_1(:,1);
%         
%         % %Flight2
%         wind_speed2=wind_matrix_final_2(:,2); 
%         time2=wind_matrix_final_2(:,1);
% 
% %Fix the seed of random generator to reproduce same results 
% seed_value = 42;
% rng(seed_value);
% 
% %Choose between wind data of Flight 1 and Flight 2 and the period of time
% 
% choice_flight=input('Which flight test? (Digit 1 for flight 1 or 2 for flight 2) ')
% if choice_flight==1
% 
%     % %Flight1
%     %index (1:17961)=10min
%         wind_speed=wind_speed1(1:17961);
%         time=time1(1:17961);
% else
% 
%     % %Flight2
%         %index (1:17533)=0-10 min  
%         %index (1:3593)=0-2 min
%         %index (89806:107767)=50-60 min
%         wind_speed=wind_speed2(1:17533);
%         time=time2(1:17533);
% end
% 
% %Resampling of values with uniformly time intervals
% choice_resample=input('Which uniformly time interval? (In seconds) ')
% time_interval=choice_resample;
% Fs=1/time_interval;
% [wind_speed,time]=resample(wind_speed,time,Fs);

%------------------------------------------------------------------------





% Noise parameters
signal_length = 6000; % Signal length (e.g., 10 minutes at 10 Hz)
base_mean = 5; % Base mean of the noise
base_std_dev = 1; % Base standard deviation of the noise

% Define the time intervals and corresponding mean and variance values
time_intervals = [1, 4000, 6000]; % Start and end of each interval
mean_intervals = [5, 9, 7.55]; % Mean values for each interval
variance_intervals = [1, 2, 1]; % Variance values for each interval

% Linearly interpolate the mean and variance values between intervals
interpolated_mean = interp1(time_intervals, mean_intervals, 1:signal_length, 'linear', 'extrap');
interpolated_variance = interp1(time_intervals, variance_intervals, 1:signal_length, 'linear', 'extrap');

% Generate noise with varying mean and variance
gust_signal_generated= interpolated_variance .* randn(1, signal_length) + interpolated_mean;


% Define sinusoidal parameters
frequency = 0.5; % Frequency of the sinusoid (e.g., 0.1 Hz)
amplitude = 0.5; % Amplitude of the sinusoid

% Generate sinusoidal component
t = 0:1/10:600; % Time vector at 10 Hz (6000 samples for 10 minutes)
sinusoidal_component = amplitude * sin(2 * pi * frequency * t);

% Add sinusoidal component to the gust signal
gust_signal_with_sine = gust_signal_generated + sinusoidal_component(1:6000);

% Define second sinusoidal parameters
frequency2 =0.092; % Frequency of the second sinusoid (e.g., 0.05 Hz)
amplitude2 = 1; % Amplitude of the second sinusoid
phase_shift = pi/4; % Phase shift of the second sinusoid (e.g., pi/2)

% Generate second sinusoidal component
sinusoidal_component2 = amplitude2 * cos(2 * pi * frequency2 * t + phase_shift);

% Add second sinusoidal component to the gust signal
gust_signal = gust_signal_with_sine + sinusoidal_component2(1:6000);

% Plot the resulting signal
figure;
plot(gust_signal);
xlabel('Time [s]');
ylabel('Wind Speed (m/s)');
title('Wind gust Signal with Two Sinusoidal Components');

wind_speed=gust_signal;


%*************************************************

% STEP 1: Creata a superposition of SDE with the desired ACF of data

%1.1 Calculate the ACF of data

% choice_nLag=input('How many lags? ')
% maxLag=choice_nLag;
maxLag=100;
[acf,lag]=autocorr(wind_speed,maxLag); 
% lag=(0:maxLag)';


        % %ACF plots
        % %1)Discrete sequence data plot 

        figure;
        stem(acf, 'LineWidth', 1.5);
        title('Auto-correlation function','FontSize', 16);
        xlabel('Lag [s]');
        ylabel('Auto-correlation of wind speed data');
        grid on;

        % %2)Continue sequence data plot

        % figure;
        % plot(lag,acf, 'LineWidth', 1.5);
        % title('Auto-correlation function','FontSize', 16);
        % xlabel('Lag [s]');
        % ylabel('Auto-correlation of wind speed data');
        % grid on;
 
 
  
  
 

%1.2 Fit a weighted sum of exponential and/or damped sinusoidal functions to data ACF

%1.2. a)Evaluate the fit

PopSz = 500;
[coef,num_comp,fit_function] = function_fitACF(acf,lag,PopSz) %coef=(w,alfa,omega)

    save('f2_coef_2min_lag100_res1',"coef")
    save('f2_fit_2min_lag100_res1',"fit_function")
    save('f2_ncomp_2min_lag100_res1',"num_comp")

%1.2. b)Load a fit

%Flight 1 

%%This data are for wind_speed 10, resample 0.1, lag 100;
% load('f1_coef_10min_lag100_res01.mat')
% load('f1_fit_10min_lag100_res01.mat')
% load('f1_ncomp_10min_lag100_res01.mat')

%%This data are for wind_speed 10, resample 0.5, lag 100;
% load('f1_coef_10min_lag100_res05.mat')
% load('f1_fit_10min_lag100_res05.mat')
% load('f1_ncomp_10min_lag100_res05.mat')

% %%This data are for wind_speed 10, resample 1, lag 100;
% load('f1_coef_10min_lag100_res1.mat')
% load('f1_fit_10min_lag100_res1.mat')
% load('f1_ncomp_10min_lag100_res1.mat')

%Flight 2 

%%This data are for wind_speed(1:3593)=2min, resample 0.01, lag 100;
% load('coef_2min_lag100_res01.mat')
% load('fit_2min_lag100_res01.mat')
% load('ncomp_2min_lag100_res01.mat')

%%This data are for wind_speed(1:3593)=2min, resample 0.5, lag 100;
% load('coef_2min_lag100_res5.mat')
% load('fit_2min_lag100_res5.mat')
% load('ncomp_2min_lag100_res5.mat')

%%This data are for wind_speed 50-60 min, resample 1, lag 100;
% load('f2_coef_50min_lag100_res1.mat')
% load('f2_fit_50min_lag100_res1.mat')
% load('f2_ncomp_50min_lag100_res1.mat')


%%This data are for wind_speed 10min, resample 1, lag 100;
% load('f2_coef_10min_lag100_res1.mat')
% load('f2_fit_10min_lag100_res1.mat')
% load('f2_ncomp_10min_lag100_res1.mat')


    %Plot  fit results
    Cfit = fit_function(coef,lag);
    
    figure
    hd = plot(lag, acf, 'p', 'DisplayName','Data');
    for k1 = 1:size(acf,2)
        CV(k1,:) = hd(k1).Color;
        hd(k1).MarkerFaceColor = CV(k1,:);
    end
    hold on
    hlp = plot(lag, Cfit, 'DisplayName',['Sum of ' num2str(num_comp) ' weighted dumped exponential ']);
    title('Fit of ACF with weighted sum of exponential and damped sinusoidal functions ')
    for k1 = 1:size(acf,2)
        hlp(k1).Color = CV(k1,:);
    end
    hold off
    grid
    legend('Location','best')

%Estrapolate model coefficients
weight=coef(1:3:end);
alpha=coef(2:3:end);
omega=coef(3:3:end);
mean_value=mean(wind_speed);
%mu=mean*ones(length(num_comp));
mu=zeros(num_comp);
sigma=sqrt(2*abs(alpha));
%% 

%1.3 Integrate the SDE 

x0=mean_value;  %mean_value; %x initial point
y0=0;   %mean_value; %y initial poin
T= 700; %time interval
dt=0.0001; %integration steps

%Integretion methods

[v_acf_exp,t_exp]=exponential_euler ( num_comp, x0,y0, T, dt , weight, alpha, omega,sigma); %Exponential Euler
[v_acf_eul,t_eul]=euler_maruyama ( num_comp, x0,y0, T, dt , weight, alpha, omega,sigma); %Euler-Maruyama
[v_acf_tay,t_tay]=taylor1_5( num_comp, x0,y0, T, dt , weight, alpha, omega,sigma); %order-1.5 Taylor
[v_acf_heun,t_heun]=heun( num_comp, sigma, x0,y0, T, dt, weight, alpha, omega); %Heun 

    %Wind signal plots
    figure;
    subplot(2, 2, 1);
    plot(t_exp, v_acf_exp);
    title('Exponetial Euler');
    
    subplot(2, 2, 2);
    plot(t_eul, v_acf_eul);
    title('Euler-Maruyama');
    
    subplot(2, 2, 3);
    plot(t_tay, v_acf_tay);
    title('Taylor 1.5');
    
    subplot(2, 2, 4);
    plot(t_heun, v_acf_heun);
    title('Heun');


%1.4 ACF comparison
%To compare the ACF of real data with the one of the generated data it's necessary to compare values at same time istants


v_acf_r_exp=v_acf_exp(1:1/Fs/dt:end);
t_r_exp=t_exp(1:1/Fs/dt:end);
[acf_exp,lag_exp]=autocorr(v_acf_r_exp,'NumLag', maxLag);

v_acf_r_eul=v_acf_eul(1:1/Fs/dt:end);
t_r_eul=t_eul(1:1/Fs/dt:end);
[acf_eul,lag_eul]=autocorr(v_acf_r_eul,'NumLag', maxLag);


v_acf_r_tay=v_acf_tay(1:1/Fs/dt:end);
t_r_tay=t_tay(1:1/Fs/dt:end);
[acf_tay,lag_tay]=autocorr(v_acf_r_tay,'NumLag', maxLag);


v_acf_r_heun=v_acf_heun(1:1/Fs/dt:end);
t_r_heun=t_heun(1:1/Fs/dt:end);
[acf_heun,lag_heun]=autocorr(v_acf_r_heun,'NumLag', maxLag);



    %ACF of generated signal
    figure 
    plot(lag_exp,acf_exp,'-',LineWidth=1)
    hold on 
    plot(lag_eul,acf_eul,'-',LineWidth=1)
    plot(lag_tay,acf_tay,'-',LineWidth=1)
    plot(lag_heun,acf_heun,'-',LineWidth=1)
    plot(lag,acf,'-',LineWidth=2)
    legend('Exponential euler','Euler-Maruyama','Taylor 1.5','Heun','Real data')
    xlabel('Lag [s]')
    ylabel('ACF')
    title('Integration methods','FontSize',20)
    hold off

% figure 
% plot(lag_exp,mean((acf_exp-acf).^2),'-',LineWidth=1)
% hold on 
% plot(lag_eul,mean((acf_eul-acf).^2),'-',LineWidth=1)
% plot(lag_tay,mean((acf_tay-acf).^2),'-',LineWidth=1)
% plot(lag_heun,mean((acf_heun-acf).^2),'-',LineWidth=1)
% legend('Exponential euler error','Euler-Maruyama error','Taylor 1.5 error','Heun')
% xlabel('Lag [s]')
% ylabel('Mean square error')
% title('Mean square error','FontSize',20)
% hold off
% 
% figure 
% plot(lag_exp,mean((acf_exp-acf)),'-',LineWidth=1)
% hold on 
% plot(lag_eul,mean((acf_eul-acf)),'-',LineWidth=1)
% plot(lag_tay,mean((acf_tay-acf)),'-',LineWidth=1)
% plot(lag_heun,mean((acf_heun-acf)),'-',LineWidth=1)
% legend('Exponential euler error','Euler-Maruyama error','Taylor 1.5 error','Heun')
% xlabel('Lag [s]')
% ylabel('Mean absolute error')
% title('Mean absolute error','FontSize',20)
% hold off


%TEST the results

% Test the gaussian PDF of SDE

figure;
histogram(v_acf_heun, 'Normalization', 'pdf');
hold on;
pd = fitdist(v_acf_heun', 'Normal');
x_values = linspace(min(v_acf_heun), max(v_acf_heun), 100);
y_values = pdf(pd, x_values);
plot(x_values, y_values, 'LineWidth', 2);
title('Confronto Istogramma con PDF Gaussiana');
legend('Istogramma', 'PDF Gaussiana');





[acf_step1,lag_step1]=autocorr(v_acf_r_heun,'NumLag', maxLag);

    %Compare the ACF of real data with the ACF of generated data
    figure
    subplot(2,1,1)
    stem(acf_step1, 'LineWidth', 1.5);
    title('Auto-correlation function','FontSize', 16);
    xlabel('Lag [s]');
    ylabel('Auto-correlation of generated wind speed at the end of dtep1');
    grid on;
    subplot(2, 1, 2);
    stem( acf, 'LineWidth', 1.5);
    title('Auto-correlation function','FontSize', 16);
    xlabel('Lag [s]');
    ylabel('Auto-correlation of wind speed data');
    grid on;


% ****************************************************

%STEP 2: impose the desired PDF

%2.1 Calculate the inverse CDF of real data 

[ecdf,x_ecdf]=ecdf(wind_speed); %Use the Kaplan-Meier estimator to estiamte ecdf
inv_ECDF=@(z) interp1(ecdf,x_ecdf,z,'linear','extrap'); 

%2.2 Calculate the gaussian CDF of generated data 
gaussianCDF=normcdf(v_acf_heun);

%2.3 Memoryless transformation 
v_pdf=inv_ECDF(gaussianCDF);

    %Plot final wind generated signal
    figure
    plot(t_heun,v_pdf, 'r-', 'LineWidth', 3 )
    xlabel ( 't [s]', 'FontSize', 16 )
    ylabel ( 'v [m/s]', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
    title ( 'SDE withe the desired PDF and ACF', 'FontSize', 16 )
    grid ( 'on' );

%**************************************************

%TEST the results

%Comparison plot of PDF
figure;
subplot(2, 1, 1);
histogram(v_pdf, 100,'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('v [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of generated wind speed', 'FontSize', 16)
xlim([0,0.55])
ylim([0,15])
grid('on');
subplot(2, 1, 2);
histogram(wind_speed, 100,'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('v [m/s]', 'FontSize', 16)
ylabel('PDF', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right')
title('PDF of wind speed data', 'FontSize', 16)
xlim([0,0.55])
ylim([0,15])
grid('on');

%Comparison plot of ACF
v_pdf_r=v_acf_heun(1:1/Fs/dt:end);
acf_final=autocorr(v_pdf_r,'NumLag', maxLag);

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





