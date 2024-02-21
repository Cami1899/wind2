%VENUSIAN GUSTS
%***********************************************************************************
%Having no available signals of horizontal wind gusts in the Venus atmosphere, 
% a speculative simulation is carried out by going to impose a PDF on the model 
% based on the measurements of horizontal winds detected by Venera 8. 
% Instead, the gust dynamics are maintained using an ACF function from a gust signal on Earth 


%***********************************************************************************


PopSz = 700;
[coef,num_comp,fit_function] = function_fitACF(acf,lag,PopSz);


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

%1.3 Integrate the SDE 
x0=mean_value;  %mean_value; %x initial point
y0=0;   %mean_value; %y initial poin
T= 700; %time interval
dt=0.0001; %integration steps

[v_acf_heun,t_heun]=heun( num_comp, sigma, x0,y0, T, dt, weight, alpha, omega); %Heun 

v_acf_r_heun=v_acf_heun(1:1/Fs/dt:end);
t_r_heun=t_heun(1:1/Fs/dt:end);
[acf_heun,lag_heun]=autocorr(v_acf_r_heun,'NumLag', maxLag);
 
figure 
plot(lag_heun,acf_heun,'-',LineWidth=1)
plot(lag,acf,'-',LineWidth=2)
legend('Heun','Real data')
xlabel('Lag [s]')
ylabel('ACF')
hold off

%STEP 2: impose the desired PDF
%Wind pdf of the real gust it's rescaled with the value of horizontal wind of Venera 8 measurments.

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
