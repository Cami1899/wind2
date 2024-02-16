%Process Wind Data from Flights

clc
clear 
load('aprsAlt2.mat')
load('flight2Wind.mat')
load('gondolaAlt1.mat')
load('pressures1.mat')
load('pressures2.mat')

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

% % Plot (transducer 3 windspeed - transducer 1)
% plot(allTime1-allTime1(1), diff13Wind1)
% title("(Windspeed 3 - Windspeed 1) vs Time Since Startup, Flight 1", 'FontSize', 20)
% ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% set(gcf,'Position',[0 0 1500 500])

teta= atan2(diff24Wind2,diff13Wind2);

% figure()
% plot(allTime2-allTime2(1),windMagnitude2)
% hold on
% title("(Wind magnitude) vs Time Since Startup, Flight 2", 'FontSize')
% ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% %set(gcf,'Position',[0 0 1500 500])
% hold off

figure()
plot(allTime1-allTime1(1),windMagnitude1)
hold on
title("(Wind magnitude) vs Time Since Startup, Flight 1", 'FontSize')
ylabel("Approximate Windspeed (m/s)", 'FontSize', 15)
xlabel("Time Since Startup (s)", 'FontSize', 15)
%set(gcf,'Position',[0 0 1500 500])
hold off

% figure()
% plot(allTime2-allTime2(1),teta)
% hold on
% title("(Wind direction: wind vector's angle in gongola x-y plane) vs Time Since Startup, Flight 2", 'FontSize')
% ylabel("Angle in x-y plane(rad)", 'FontSize', 15)
% xlabel("Time Since Startup (s)", 'FontSize', 15)
% %set(gcf,'Position',[0 0 1500 500])
% hold off

%Extrappolation of wind data from launch to land

allTimeplot=allTime2-allTime2(1);
index1=66125;
index2=181929;

new_magnitude=windMagnitude2(index1:index2);
%new_time=1:size(new_magnitude)';
new_time=allTime2(index1:index2)-allTime2(index1);
new_angle=teta(index1:index2);
v_x=diff13Wind2(index1:index2);
v_y=diff24Wind2(index1:index2);
frequency=1./diff(new_time);
new_data=diff(new_magnitude);



%Plot wind magnitude 

figure()
plot(new_time,new_magnitude)
hold on
title("Wind magnitude vs Time Since Launch, Flight 2", 'FontSize', 20)
ylabel(" Wind speed (m/s)", 'FontSize', 15)
xlabel("Time Since Launch (s)", 'FontSize', 15)
%set(gcf,'Position',[0 0 1500 500])
hold off

figure()
hist(new_magnitude)
title('Histogram of wind speed magnitude', 'FontSize', 20)
xlabel('Wind velocity [m/s]', 'FontSize', 15)
ylabel('Amount','FontSize', 15)

figure
hist(new_angle)
title('distribution velocity angle ')
xlabel('angle [rad]')
ylabel('Amount')

figure
polarhistogram(new_angle,100)
title('Wind angle variation', 'FontSize', 20)

polarscatter(new_angle,new_magnitude,[],'filled')
%% 
varFs=1./diff(new_time);
Fs_mean=mean(varFs);
d=diff(new_time);


figure
plot(varFs,varFs,'o',LineWidth=5)
hold on
plot(Fs_mean,Fs_mean,'*',LineWidth=5)
xlabel('Sample frequencies [Hz]')
ylabel('Sample frequencies [Hz]')
hold off

% figure
% hist(varFs)
% xlabel('Sample frequencies [Hz]')

%Frequency extrapolation
% ft=nufft(new_magnitude, new_time);
[ft,freq]=stft(new_magnitude);
n = length(new_time);
f = (0:n-1)/n;


figure
semilogx(freq,ft, LineWidth=2)
title('Wind speed frequencies', 'FontSize', 20)
xlabel('f [HZ]', 'FontSize', 15)
ylabel('Amplitude','FontSize', 15)


% % Fs=Fs_mean;
% L=length(new_magnitude);
% n=floor(L/2);
% % frequencies = (-N/2:N/2-1) * f_s / N; %includes both negative and positive frequencies
% frequencies = (0:n-1)*(Fs/L); %only positive frequencies
% power_spectrum = abs(ft).^2/L;
% amplitude=abs(ft(1:n) / L);
% semilogx(frequencies,amplitude );
% xlabel('log10(f) [HZ]')
% ylabel('Amplitude')

%% 

%Polar plot variation of wind vector with quiver
figure(10)
h=polarplot(new_angle(1),new_magnitude(1),'or', LineWidth=500);
title("(Wind direction) vs Time Since Startup, Flight 2")
thetalim('auto')
    for i=2:size(new_angle)
    set(h, 'ThetaData', new_angle(i), 'RData', new_magnitude(i));
    pause(0.001);
    drawnow;
    end
%% 

%Plot variation of wind vector with quiver
figure(20)
num_punti = length(new_magnitude);
h = quiver(0, 0, new_magnitude(1) * cos(new_angle(1)), new_magnitude(1) * sin(new_angle(1)), 'r');
title('Variation of wind vector over time');
xlim([-1, 1]);
ylim([-1, 1]);
xlabel('asse 3-1 [m/s]')
ylabel('asse 2-4 [m/s]')
for i = 2:num_punti
      set(h, 'UData', v_x(i), 'VData', v_y(i));
    pause(0.001);
    drawnow;
end







%{
dpMagnitude1 = hypot(diff13Pressure1, diff24Pressure1);
dpMagnitude2 = hypot(diff13Pressure2, diff24Pressure2);


plot(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2)
absoluteCombinedRoll = (diff13Wind1.^2 + diff24Wind1.^2).^0.5;

syncWindMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), windMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncWindMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), windMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));

syncPMag2 = sync12(CPU2GPS('gondola2', pressures2(:, 1)), dpMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));
syncPMag1 = sync12(CPU2GPS('gondola1', pressures1(:, 1)), dpMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));
%}




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




