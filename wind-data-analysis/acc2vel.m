clc
clear all


%GPS didn't work at high altitudes, bcm velocities are not measured and gondola speeds are wrongly measured.
%Linear velocity are calculated integrating linear acceleration

% load("allPinsGondola1.mat")
% load("allPinsGondola2.mat")
load("diffWind1.mat")
load("diffWind2.mat")
load("gondolaAllPimu1.mat")
load("gondolaAllPimu2.mat")
load('gondolaAlt1.mat')
% load("bcmAllPins1.mat")
% load("bcmAllPins2.mat")
load("bcmAllPimu1.mat")
load("bcmAllPimu2.mat")

gondolaAllPimu1=gondolaPimu1;
gondolaAllPimu2=allPimu;

%allPins contains velocity along x,y,z in body frame  in columnn 11,12,13, the 1st coloumn is CPU time, 3rd column is GPS time
%diffWind contains GPS time, differential wind 3-1 ((positive datapoint means windspeed at 3 is higher),diff wind 4-2 (positive datapoint means windspeed at 4 is higher)
%allPimu contains angular acceleration in columns 4,5,6 and linear
%acceleration in columns 7,8,9


%Resempling of data with constant frequency 
f_nnuni=1./diff(CPU2GPS('bcm1',bcmAllPimu1(:,1)));
plot(f_nnuni,f_nnuni)
Fs=mean(f_nnuni);
[acc_res_x,t_res]=resample(bcmAllPimu1(:,7),CPU2GPS('bcm1',bcmAllPimu1(:,1)),Fs);
[acc_res_y,t_res]=resample(bcmAllPimu1(:,8),CPU2GPS('bcm1',bcmAllPimu1(:,1)),Fs);
%spectrogram
figure
spectrogram(acc_res_x,[],[],[],Fs,'yaxis')

%Maximum frequency 
L=length(t_res);                      
T = 1/Fs;                               
Y = fft(acc_res_x);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure
semilogx(f,P1)
hold on
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('FFT analysis of x accelerations')
hold off
% [max_freq,index_max]=max(f);


%Maximum frequency 
                              
Y = fft(acc_res_y);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure
semilogx(f,P1)
hold on
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('FFT analysis of y accelerations')
hold off
% [max_freq,index_max]=max(f);
% 
% %%Design High Pass Filter
% fc = 0.6;  % Cut off Frequency
% order = 2; % 6th Order Filter
% %%Filter  Acceleration Signals
% [b1 a1] = butter(order,fc,'low');
% accxf=filter(b1,a1,mag);
% plot(time,accxf,'r'); hold on
% plot(time,accxf)
% xlabel('Time (sec)')
% ylabel('Acceleration')
% %%First Integration (Acceleration - Veloicty)
% velocity=cumtrapz(time,accxf);
% figure (3)
% plot(time,velocity)
% xlabel('Time (sec)')
% ylabel('Velocity')


