
%****************************************************************
%% Absolute wind velocity
%Discussion
%The script is used to compute the absolute velocity of wind. 
% Initially the script align the wind value, provided in the refernce
% system composed by the orthogonal axis of Pitot tubes, to the INS
% reference system.

%Pitot tube RS: axis 24, axis 13, origin in CoM of the system
%INS RS: axis x, axis y, origin in the INS

%axis 24 will be align to axis x
%axis 13 will be align to axis y


%****************************************************************+
clc
clear all
load("allPinsGondola1.mat")
load("allPinsGondola2.mat")
load("diffWind1.mat")
load("diffWind2.mat")
load("gondolaAllPimu1.mat")
load("gondolaAllPimu2.mat")
load('gondolaAlt1.mat')
load("bcmAllPins1.mat")
load("bcmAllPins2.mat")
load("bcmAllPimu1.mat")
load("bcmAllPimu2.mat")

gondolaAllPimu1=gondolaPimu1;
gondolaAllPimu2=allPimu;


%allPins contains velocity along x,y,z in body frame  in columnn 11,12,13, the 1st coloumn is CPU time, 3rd column is GPS time
%diffWind contains GPS time, differential wind 3-1 ((positive datapoint means windspeed at 3 is higher),diff wind 4-2 (positive datapoint means windspeed at 4 is higher)
%allPimu contains angular acceleration in column 4,5,6

%%*************************************************************************

%GPS didn't work at high altitudes, bcm velocities are not measured and gondola speeds are wrongly measured.
%Integration of linear acceleration



vel_x1BCM=cumtrapz(CPU2GPS('bcm1',bcmAllPimu1(:,1)),bcmAllPimu1(:,7));
vel_y1BCM=cumtrapz(bcmAllPimu1(:,1),CPU2GPS('bcm1',bcmAllPimu1(:,8)));
figure
plot(CPU2GPS('bcm1',bcmAllPimu1(:,1)),bcmAllPimu1(:,7))
hold on 
title("Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
ylabel("bcm a_x  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

figure
plot(CPU2GPS('bcm1',bcmAllPimu1(:,1)),vel_x1BCM)
hold on 
title("Acceleration vs Time Since Startup, Flight 1", 'FontSize', 20)
ylabel("bcm Vx  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
%****************************************************************+

%Syncronization of GPS time of gondola velocity to GPS time of wind


vel_x1G=sync12(allPinsgondola1(:,3), allPinsgondola1(:,11),diffWind1(:,1));
vel_y1G=sync12(allPinsgondola1(:,3), allPinsgondola1(:,12),diffWind1(:,1));
vel_z1G=sync12(allPinsgondola1(:,3), allPinsgondola1(:,13),diffWind1(:,1));

vel_INS1=[vel_x1G'; vel_y1G'; vel_z1G'];

vel_x2G=sync12(allPinsgondola2(:,3), allPinsgondola2(:,11),diffWind2(:,1));
vel_y2G=sync12(allPinsgondola2(:,3), allPinsgondola2(:,12),diffWind2(:,1));
vel_z2G=sync12(allPinsgondola2(:,3), allPinsgondola2(:,13),diffWind2(:,1));

vel_INS2=[vel_x2G'; vel_y2G'; vel_z2G'];


omega_x1G=sync12(CPU2GPS('gondola1', gondolaAllPimu1(:,1)), gondolaAllPimu1(:,4),diffWind1(:,1));
omega_y1G=sync12(CPU2GPS('gondola1',gondolaAllPimu1(:,1)), gondolaAllPimu1(:,5),diffWind1(:,1));
omega_z1G=sync12(CPU2GPS('gondola1', gondolaAllPimu1(:,1)), gondolaAllPimu1(:,6),diffWind1(:,1));

omega_INS1=[omega_x1G'; omega_y1G'; omega_z1G'];

omega_x2G=sync12(CPU2GPS('gondola2', gondolaAllPimu2(:,1)), gondolaAllPimu2(:,4),diffWind2(:,1));
omega_y2G=sync12(CPU2GPS('gondola2', gondolaAllPimu2(:,1)), gondolaAllPimu2(:,5),diffWind2(:,1));
omega_z2G=sync12(CPU2GPS('gondola2', gondolaAllPimu2(:,1)), gondolaAllPimu2(:,6),diffWind2(:,1));

omega_INS2=[omega_x2G'; omega_y2G'; omega_z2G'];

%Geometric features of Gondola
r=[0.80;0;0];
alpha=50*(pi/180);

for i=1:length(omega_INS1)
%Align the origin of INS rf with the origin of Pitot tube system which is also the system CoM 
v_INS_CoM1(:,i)=vel_INS1(:,i)-cross(omega_INS1(:,i),r);

%Transform wind vector from pitot tubes' rs to INS rs with a rotational matrix
v_wind_INS1(:,i)=[-sin(alpha), -cos(alpha);-cos(alpha), sin(alpha)]*[diffWind1(i,2); diffWind1(i,3)];

%Compute the absolute
v_wind_abs1(:,i)=v_wind_INS1(:,i)-v_INS_CoM1(1:2,i);

absolute_wind_mag1(i)=hypot(v_wind_abs1(1,i),v_wind_abs1(2,i));

teta1(i)= atan2(v_wind_abs1(2,i),v_wind_abs1(1,i));
end



for i=1:length(omega_INS2)
v_INS_CoM2(:,i)=vel_INS2(:,i)-cross(omega_INS2(:,i),r);
v_wind_INS2(:,i)=[-sin(alpha), -cos(alpha);-cos(alpha), sin(alpha)]*[diffWind2(i,2); diffWind2(i,3)];
v_wind_abs2(:,i)=v_wind_INS2(:,i)-v_INS_CoM2(1:2,i);
absolute_wind_mag2(i)=hypot(v_wind_abs2(1,i),v_wind_abs2(2,i));
teta2(i)= atan2(v_wind_abs2(2,i),v_wind_abs2(1,i));
end

figure
plot(diffWind1(:,1),absolute_wind_mag1)
hold on 
title("Absolute Gust Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
ylabel("Gust Speed  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

figure
plot(diffWind2(:,1),absolute_wind_mag2)
hold on
title("Absolute Gust Speed  vs Time Since Startup, Flight 2", 'FontSize', 20)
ylabel("Gust Speed  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
%% Comparison between velocities of gondola and velocities of bcm to understand if we can apply the same wind data to balloon
figure
plot(allPinsgondola1(:,3), allPinsgondola1(:,11), bcmAllPins1(:,3),bcmAllPins1(:,11))
hold on 
title("Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
ylabel("Speed  (m/s)", 'FontSize', 15)
legend('gondola Vx','bcm Vx')
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

figure
plot(allPinsgondola1(:,3), allPinsgondola1(:,12), bcmAllPins1(:,3),bcmAllPins1(:,12))
hold on 
title("Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
ylabel("Speed  (m/s)", 'FontSize', 15)
legend('gondola Vy','bcm Vy')
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

figure

title("BCM Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
subplot(3,1,1)
plot( bcmAllPins1(:,3),bcmAllPins1(:,11))
hold on 
ylabel("Vx   (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
subplot(3,1,2)
plot( bcmAllPins1(:,3),bcmAllPins1(:,12))
hold on 
ylabel("Vy  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
subplot(3,1,3)
plot(bcmAllPins1(:,3),bcmAllPins1(:,13))
hold on 
ylabel("Vz  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

figure

title("Gondola Speed  vs Time Since Startup, Flight 1", 'FontSize', 20)
subplot(3,1,1)
plot(allPinsgondola1(:,3), allPinsgondola1(:,11))
hold on 
ylabel("Vx   (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
subplot(3,1,2)
plot(allPinsgondola1(:,3), allPinsgondola1(:,12))
hold on 
ylabel("Vy  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off
subplot(3,1,3)
plot(allPinsgondola1(:,3), allPinsgondola1(:,13))
hold on 
ylabel("Vz  (m/s)", 'FontSize', 15)
xlabel("GPS Time Since Startup (s)", 'FontSize', 15)
hold off

%% 
var(gondolaAllPimu1(1:10000,7))
var(gondolaAllPimu1(10000:20000,7))
var(gondolaAllPimu1(300000:330000,7))