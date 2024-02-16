%Process Wind Data from Flights

clc
clear all
load('aprsAlt2.mat')
load('flight2Wind.mat')
load('gondolaAlt1.mat')
load('pressures1.mat')
load('pressures2.mat')
load("gondolaAllPimu1.mat")
load('gondolaAllPimu2.mat')
gondolaPimu2=allPimu;

allTime1 = pressures1(:, 1);
allTime2 = pressures2(:, 1);


%pressures1 & pressures2 are matrices of time, and dynamic pressures
%measures at each of the 4 pitot tubes: SW1`, SW2, SW3, LW1

%Important: pressures1 & pressures2 are internal gondola pressure minus
%dynamic pressure. Atmospheric effects caused observed internal pressure to exceed
%observed dynamic pressure 

diff13Pressure1 = pressures1(:, 2)-pressures1(:, 4);
diff24Pressure1 = pressures1(:, 3)-pressures1(:, 5); 

diff13Pressure2 = pressures2(:, 2)-pressures2(:, 4);
diff24Pressure2 = pressures2(:, 3)-pressures2(:, 5);


%Atmosisa atmospheric model to calculate density as a function of altitude
%(atmospheric model has less variation at very high altitudes, near apex of flight)

[T1, A1, P1, Rho1] = atmosisa(gondolaAlt1(:, 2));

%Use APRS Altitude for density calculation in flight 2: fewer datapoints but more
%consistent than INS altitude
[T2, A2, P2, Rho2] = atmosisa(aprsAlt(:, 2));

%sync12 is called to compute operations using points samples at different
%times.

syncedRho1 = sync12(gondolaAlt1(:, 1), Rho1, pressures1(:, 1));
syncedRho2 = sync12(aprsAlt(:, 1), Rho2, allTime2);


%Approximate differential windspeed
%Important: A positive datapoint in diff13Wind1 or diff13Wind2 means that windspeed at transducer 3 is
%higher than windspeed at transducer 1. Likewise, a positive datapoint in diff24Wind1 or diff24Wind2 means that windspeed at transducer 4 is
%higher than windspeed at transducer 2.

diff13Wind1 = pressure2Wind(diff13Pressure1, syncedRho1);
    
diff24Wind1 = pressure2Wind(diff24Pressure1, syncedRho1);

diff13Wind2 = pressure2Wind(diff13Pressure2, syncedRho2);
    
diff24Wind2 = pressure2Wind(diff24Pressure2, syncedRho2);


windMagnitude1 = (diff13Wind1.^2 + diff24Wind1.^2).^0.5;
windMagnitude2 = (diff13Wind2.^2 + diff24Wind2.^2).^0.5;

% Plot (transducer 3 windspeed - transducer 1)
% plot(allTime1-allTime1(1), diff13Wind1)
% title("(Windspeed 3 - Windspeed 1) vs Time Since Startup, Flight 1", 'FontSize', 20)
% ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% set(gcf,'Position',[0 0 1500 500])
teta1= atan2(diff24Wind1,diff13Wind1);
teta2= atan2(diff24Wind2,diff13Wind2);

% figure()
% plot(allTime2-allTime2(1),windMagnitude2)
% hold on
% title("(Wind magnitude) vs Time Since Startup, Flight 2", 'FontSize')
% ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% %set(gcf,'Position',[0 0 1500 500])
% hold off

% figure()
% plot(allTime1-allTime1(1),windMagnitude1)
% hold on
% title("(Wind magnitude) vs Time Since Startup, Flight 1", 'FontSize')
% ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% %set(gcf,'Position',[0 0 1500 500])
% hold off
% 
% figure()
% plot(allTime2-allTime2(1),teta2)
% hold on
% title("(Wind direction: wind vector's angle in gongola x-y plane) vs Time Since Startup, Flight 2", 'FontSize')
% ylabel("Angle in x-y plane(rad)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% %set(gcf,'Position',[0 0 1500 500])
% hold off
% 
%Extrappolation of wind data from launch to burst
    % For Flight 1, launch happens at GPS Time = 398613.3
    % For Flight 1, burst happens at GPS Time = 402713.3
    % For Flight 2, launch happens at GPS Time = 399283.6
    % For Flight 2, burst happens at GPS Time =  403638.1
%Launch 1
absoluteDifferencesL1 = abs(pressures1(:,1) - 398613.3); 
[minDifferenceL1, linearIndexL1] = min(absoluteDifferencesL1(:));
%Burst 1
absoluteDifferencesB1 = abs(pressures1(:,1) - 402713.3);
[minDifferenceB1, linearIndexB1] = min(absoluteDifferencesB1(:))
%Launch 2
absoluteDifferencesL2 = abs(pressures2(:,1) - 399283.6); 
[minDifferenceL2, linearIndexL2] = min(absoluteDifferencesL2(:));
%Burst 1
absoluteDifferencesB2 = abs(pressures2(:,1) - 403638.1);
[minDifferenceB2, linearIndexB2] = min(absoluteDifferencesB2(:));
% index1 = find(pressures1(:,1) == 398613.3)
ascentPeriod1 = linearIndexL1:linearIndexB1;
ascentPeriod2 = linearIndexL2:linearIndexB2;

allTimeplot1=allTime1-allTime1(1);
new_magnitude1=windMagnitude1(linearIndexL1:linearIndexB1);
wind_matrix=[pressures1(ascentPeriod1,1),diff13Wind1(ascentPeriod1),diff24Wind1(ascentPeriod1)];
%new_time=1:size(new_magnitude)';
new_time1=allTime1(linearIndexL1:linearIndexB1);
new_angle1=teta1(linearIndexL1:linearIndexB1);
v_x1=diff13Wind1(linearIndexL1:linearIndexB1);
v_y1=diff24Wind1(linearIndexL1:linearIndexB1);
frequency1=1./diff(new_time1);
new_data1=diff(new_magnitude1);

wind_matrix_final_1=[new_time1,new_magnitude1]; %wind vector from launch to burst of flight 1: 1st column is GPS time from the first value after the GPS istant of launch 398613.3 to the last value before burst istanat at 402713.3, second column values are the wind magnitude
wind_matrix_components_1=[new_time1,v_x1,v_y1];

new_magnitude2=windMagnitude2(linearIndexL2:linearIndexB2);
%new_time=1:size(new_magnitude)';
new_time2=allTime2(linearIndexL2:linearIndexB2);
new_angle2=teta2(linearIndexL2:linearIndexB2);
v_x2=diff13Wind2(linearIndexL2:linearIndexB2);
v_y2=diff24Wind2(linearIndexL2:linearIndexB2);
frequency2=1./diff(new_time2);
new_data2=diff(new_magnitude2);

wind_matrix_final_2=[new_time2,new_magnitude2]; %wind vector from launch to burst of flight 2: 1st column is GPS time from the first value after the GPS istant of launch 399283.6 to the last value before burst istanat at 403638.1
wind_matrix_components_2=[new_time2,v_x2,v_y2];
%from 5 min 20 sec to 6 min 20 sec since launch
w1=398613.3+5*60+20;
w2=398613.3+6*60+20;
[val,ind1]=min(abs(pressures1(:,1) - w1));
[val,ind2]=min(abs(pressures1(:,1) - w2));
w1_time_mag_dir=[allTime1(ind1:ind2),teta1(ind1:ind2),windMagnitude1(ind1:ind2)];
w1_restricted=[allTime1(ind1:ind2),diff13Wind1(ind1:ind2),diff24Wind1(ind1:ind2)];
%%


% Plot wind magnitude 

figure()
plot(new_time1,new_magnitude1)
hold on
title("Wind speed magnitude, Flight 1", 'FontSize', 20)
ylabel(" Wind speed (m/s)", 'FontSize', 15)
xlabel("Time Since Launch (s)", 'FontSize', 15)
%set(gcf,'Position',[0 0 1500 500])
hold off

figure()
hist(new_magnitude1)
title('Histogram of wind speed magnitude', 'FontSize', 20)
xlabel('Wind velocity [m/s]', 'FontSize', 15)
ylabel('Amount','FontSize', 15)



figure()
plot(new_time2,new_magnitude2)
hold on
title("Wind speed magnitude, Flight 2", 'FontSize', 20)
ylabel(" Wind speed (m/s)", 'FontSize', 15)
xlabel("Time Since Launch (s)", 'FontSize', 15)
%set(gcf,'Position',[0 0 1500 500])
hold off

figure()
hist(new_magnitude2)
title('Histogram of wind speed magnitude', 'FontSize', 20)
xlabel('Wind velocity [m/s]', 'FontSize', 15)
ylabel('Amount','FontSize', 15)
% 
% figure
% hist(new_angle2)
% title('distribution velocity angle ')
% xlabel('angle [rad]')
% ylabel('Amount')
% 
figure
polarhistogram(new_angle1,100)
title('Wind angle variation, Flight 1', 'FontSize', 20)

figure
polarhistogram(new_angle2,100)
title('Wind angle variation, Flight 2', 'FontSize', 20)

% % 
% varFs=1./diff(new_time2);
% Fs_mean=mean(varFs);
% d=diff(new_time2);
% 
% 
% figure
% plot(varFs,varFs,'o',LineWidth=5)
% hold on
% plot(Fs_mean,Fs_mean,'*',LineWidth=5)
% xlabel('Sample frequencies [Hz]')
% ylabel('Sample frequencies [Hz]')
% hold off
% 
% figure
% hist(varFs)
% xlabel('Sample frequencies [Hz]')
% 
% Frequency extrapolation
% ft=nufft(new_magnitude2, new_time2);
% [ft,freq]=stft(new_magnitude2);
% n = length(new_time2);
% f = (0:n-1)/n;
% 
% 
% figure
% semilogx(freq,ft, LineWidth=2)
% title('Wind speed frequencies', 'FontSize', 20)
% xlabel('f [HZ]', 'FontSize', 15)
% ylabel('Amplitude','FontSize', 15)
% 
% 
% % Fs=Fs_mean;
% L=length(new_magnitude2);
% n=floor(L/2);
% % frequencies = (-N/2:N/2-1) * f_s / N; %includes both negative and positive frequencies
% frequencies = (0:n-1)*(Fs/L); %only positive frequencies
% power_spectrum = abs(ft).^2/L;
% amplitude=abs(ft(1:n) / L);
% semilogx(frequencies,amplitude );
% xlabel('log10(f) [HZ]')
% ylabel('Amplitude')
% 
% % 
% 
% Polar plot variation of wind vector with quiver
% figure()
% h=polarplot(new_angle2(1),new_magnitude2(1),'or', LineWidth=500);
% title("(Wind direction) vs Time Since Startup, Flight 2")
% thetalim('auto')
%     for i=2:size(new_angle2)
%     set(h, 'ThetaData', new_angle2(i), 'RData', new_magnitude2(i));
%     pause(0.001);
%     drawnow;
%     end
%% 

%Plot variation of wind vector with quiver
figure()
num_punti = length(new_magnitude2);
h = quiver(0, 0, new_magnitude2(1) * cos(new_angle2(1)), new_magnitude2(1) * sin(new_angle2(1)), 'r');
title('Variation of wind vector over time');
xlim([-1, 1]);
ylim([-1, 1]);
xlabel('asse 3-1 [m/s]')
ylabel('asse 2-4 [m/s]')
for i = 2:num_punti
      set(h, 'UData', v_x2(i), 'VData', v_y2(i));
    pause(0.001);
    drawnow;
end



close all




dpMagnitude1 = hypot(diff13Pressure1, diff24Pressure1);
dpMagnitude2 = hypot(diff13Pressure2, diff24Pressure2);


plot(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2)
absoluteCombinedRoll = (diff13Wind1.^2 + diff24Wind1.^2).^0.5;

syncWindMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncWindMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), windMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));

syncPMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), dpMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncPMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), dpMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));





% % Calculate pressure from voltage
% for i = 1:length(flight2Wind)
%     for j = 2:5
%         if flight2Wind(i, j) < 2.5
%             newPressure = -133*abs((flight2Wind(i, j)/5.0 - 0.5))*(flight2Wind(i, j)/2.0 - 1.25)^2;
%         else
%             newPressure = 133*abs((flight2Wind(i, j)/5.0 - 0.5))*(flight2Wind(i, j)/2.0 - 1.25)^2;
%         end
%         pressures_total(i, j) = newPressure;
%     end
% end




