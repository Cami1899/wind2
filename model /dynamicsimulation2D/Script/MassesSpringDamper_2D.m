clear all
close all
clc


%% MSD Earth

load('wind_value_restricted.mat');
load('time_restricted.mat');
tmax=100;
n=3;
l_t = 10;         %[m] Tether lenght
l_s = l_t/n;
z_eq = 56000;     %[m] Design equilibrium altitude
h=56100; %[m] Initial altitude at the begin of the simulation
dt = 0.01; % [s] Sampling time for the output results
% Initial states [x_b;Vx_b,z_b;Vz_b;(x_t;Vx_t;z_t;Vz_t;)*n;x_G;Vx_G;z_G;Vz_G]

offset = h - z_eq;
% s_0 = [0;0;z_eq+offset;0;...
%        l_s*sin(20*pi/180);0;z_eq+offset-l_s*cos(20*pi/180);0;...
%        l_s*(sin(20*pi/180)+sin(10*pi/180)); 0 ; z_eq+offset-l_s*(cos(20*pi/180)+cos(10*pi/180));0];
%wind model 
[wind_speed_int,t_wind,v_wind]=wind_model(new_magnitude,new_time);


s_0 = zeros((n-1)*4+2*4,1);
s_0(1:4,1) = [ 0 ; 0 ; z_eq+offset  ; 0 ];
for i = 1:(n-1)
    s_0(4+4*(i-1)+1,1) = 0; 
    s_0(4+4*(i-1)+2,1) = 0; 
    s_0(4+4*(i-1)+3,1) = z_eq-i*l_s+offset;
    s_0(4+4*(i-1)+4,1) = 0;
end
s_0((end-3):end,1) = [ 0 ; 0 ; z_eq-l_t+offset  ; 0 ];

atm_data = readmatrix('vira-venus-atmosphere-45.csv');
[t,s] = ode23s(@(t,s)odeMSD_2D(t,s,atm_data,z_eq,n,wind_speed_int),[0:dt:100], s_0);

for i=1:length(t)
    [dsdt(i,:),wind_value(i)]=odeMSD_2D(t(i),s(i,:),atm_data,z_eq,n,wind_speed_int);
end
%% Plots

figure(1)
plot(t,s(:,3),t,s(:,7:4:(end-1-4)),t,s(:,end-1))
xlabel('t [s]','Interpreter','latex')
ylabel('z [m]','Interpreter','latex')
legend('Balloon','$T_1$','$T_2$','Gondola','Interpreter','latex')
grid on

figure(2)
plot(t,s(:,1),t,s(:,5:4:(end-3-4)),t,s(:,end-3))
% xlim([0 200])
xlabel('t [s]','Interpreter','latex',FontSize=15)
ylabel('x [m]','Interpreter','latex',FontSize=15)
legend('Balloon','$T_1$','$T_2$','Gondola','Interpreter','latex',FontSize=10)
title('Motiono along X',FontSize=16)
grid on

for i = 0:(n-1)
    theta(:,i+1) = atan((s(:,4*i+5)-s(:,4*i+1))./abs(s(:,4*i+3)-s(:,4*i+7)));
end

figure(3)
subplot(2,1,1);
plot(t, (180/pi)*theta(:,1),'b','LineWidth',0.1);
% xlim([0 20]);
grid on
ylabel('$\theta_1$ [deg]','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');

subplot(2,1,2);
plot(t, (180/pi)*theta(:,2),'r','LineWidth',0.1);
% xlim([0 20]);
grid on
ylabel('$\theta_2$ [deg]','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex');

figure(4)
plot(t_wind,v_wind, 'r-' )
xlabel ( 't [s]', 'Interpreter','latex' )
ylabel ( 'v [m/s]', 'Interpreter','latex', 'Rotation', 0, 'HorizontalAlignment', 'right' )
title ( 'wind value generated', 'FontSize', 16 )
grid ( 'on' );

figure(5)
plot(t,wind_value)
xlabel ( 't [s]', 'Interpreter','latex' ,'FontSize', 14)
ylabel ( 'v [m/s]', 'Interpreter','latex', 'Rotation', 0, 'HorizontalAlignment', 'right', 'FontSize', 14 )
title ( 'Wind value evaluated in t of integration', 'FontSize', 16 )
grid ( 'on' );

figure
plot(t_wind,wind_speed_int(t_wind))
xlabel ( 't [s]', 'Interpreter','latex' )
ylabel ( 'v [m/s]', 'Interpreter','latex', 'Rotation', 0, 'HorizontalAlignment', 'right' )
title ( 'wind value evaluated in t of integration', 'FontSize', 16 )
grid ( 'on' );

figure(6)
plot(t,s(:,2),t,s(:,6:4:end-2-4),t,s(:,end-2))
xlabel('t [s]','Interpreter','latex',FontSize=15)
ylabel('$\dot{x}$ [m/s]', 'Interpreter', 'latex',FontSize=15)
legend('Balloon','$T_1$','$T_2$','Gondola','Interpreter','latex',FontSize=12)
title ( 'Velocity along x', 'FontSize', 16 )
grid on

figure(7)
plot(t,dsdt(:,2),t,dsdt(:,6:4:end-2-4),t,dsdt(:,end-2))
xlabel('t [s]','Interpreter','latex',FontSize=15)
ylabel('$\ddot{x} [m/s^2]$ ', 'Interpreter', 'latex',FontSize=15)
legend('Balloon','$T_1$','$T_2$','Gondola','Interpreter','latex',FontSize=12)
title ( 'Acceleration along x', 'FontSize', 16 )
grid on


figure(7)
% subplot(3,1,3);
% plot(t, (180/pi)*theta(:,2)-(180/pi)*theta(:,1),'c');
% xlim([0 20]);
% grid on
% ylabel('$\theta_1-\theta_2$ [deg]','Interpreter','latex')
% xlabel('Time (s)','Interpreter','latex');

%% Animation

figure(20)
for i = 1:length(t)

    for j = 1:n-2
        plot(s(i,+1),s(i,+3),'or','Markersize',30)
        hold on
        plot(s(i,4*j+1),s(i,4*j+3),'ok','Markersize',5)
        plot(s(i,4*n+1),s(i,4*n+3),'squareb','Markersize',10)
        plot([s(i,4*(j-1)+1) s(i,4*(j-1+1)+1)],[s(i,4*(j-1)+3) s(i,4*(j-1+1)+3)],'-k','Linewidth',1.5)
        plot([s(i,4*j+1) s(i,4*(j+1)+1)],[s(i,4*j+3) s(i,4*(j+1)+3)],'-k','Linewidth',1.5)
        plot([s(i,4*(j+1)+1) s(i,4*(j+1+1)+1)],[s(i,4*(j+1)+3) s(i,4*(j+1+1)+3)],'-k','Linewidth',1.5)
        xlim([-5+s(i,1)  5+s(i,1)])
        ylim([-15+s(i,3) 5+s(i,3)])
        xlabel('x [m]','Interpreter','latex')
        ylabel('z [m]','Interpreter','latex')
    end
    plot(s(i,4*(j+1)+1),s(i,4*(j+1)+3),'ok','Markersize',5)
    grid on
    legend('Balloon','Tether','Gondola','Interpreter','latex')
    hold off

    pause(dt)
end
 

%% Frequency analysis

N = length(t_wind);
Fs = 1/dt;
frequency_query = (0:N-1)*Fs/N;
X = stft(v_wind);
X= X(1:N/2+1);
magnitude = abs(fftshift(X));
Pxx = abs(X).^2 / (Fs * N);  
frequencies = linspace(0, Fs/2, N/2 + 1);  
% frequency = (-Fs/2):(Fs/length(fftshift(X))):(Fs/2 - Fs/length(fftshift(X)));
figure(8)
plot(frequencies,magnitude,'Linewidth',1.2)
grid on
xlabel('Frequency [Hz]', 'FontSize',14)
ylabel('Magnitude', 'FontSize',14)
title('Fourier transform', 'FontSize',16)


figure(9);
plot(frequencies, 10*log10(Pxx(1:N/2 + 1)), 'LineWidth', 1.2); 
grid on;
xlabel('Frequency [Hz]', 'FontSize',14);
ylabel('Power [dB]', 'FontSize',14);
title('Power Spectral Density (PSD)', 'FontSize',16);

% Calcolo dello spettrogramma
window_size = 1;  % window dimension in sec
overlap = 0.5;   %  50% sovrapposition
spectrogram(v_wind, window_size*Fs, round(overlap*window_size*Fs), [], Fs, 'yaxis');
title('Spectrogram ', 'FontSize', 15)


% %%spectrum

frequencyLimits = [0 1]*pi; % Normalized frequency (rad/sample)

timeLimits = [1 10001]; 
v_wind_ROI = v_wind(:);
v_wind_ROI = v_wind_ROI(timeLimits(1):timeLimits(2));

[Pv_wind_ROI, Fv_wind_ROI] = pspectrum(v_wind_ROI, ...
    'FrequencyLimits',frequencyLimits);
figure(10)
pspectrum(v_wind_ROI, ...
    'FrequencyLimits',frequencyLimits);

%spectrogram
overlapPercent = 50;

[P,F,T] = pspectrum(v_wind_ROI, ...
    'spectrogram', ...
    'FrequencyLimits',frequencyLimits, ...
    'OverlapPercent',overlapPercent);
figure(11)
pspectrum(v_wind_ROI, ...
    'spectrogram', ...
    'FrequencyLimits',frequencyLimits, ...
    'OverlapPercent',overlapPercent);