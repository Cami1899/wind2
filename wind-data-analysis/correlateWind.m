clc
load("gondolaAllPimu1.mat")
load("gondolaAllPimu2.mat")
%This script correlates rolling averages of XY accelerations with XY Calculated Dynamic Pressure 

windMagnitude1 = hypot(diff13Wind1, diff24Wind1);
syncWindMag1 = sync12( pressures1(:, 1), windMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));

windMagnitude2 = hypot(diff13Wind2, diff24Wind2);
syncWindMag2 = sync12(pressures2(:, 1), windMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));

%Absolute XY Acceleration magnitude
xyGondAcc1 = hypot(gondolaPimu1(:, 7), gondolaPimu1(:, 8));
xyGondAcc2 = hypot(gondolaPimu2(:, 7), gondolaPimu2(:, 8));

ascentPeriod1 = 97400:261600;
ascentPeriod2 = 97070:271250;


ascentAMag1 = xyGondAcc1(ascentPeriod1); %Absolute acceleration magnitude, flight 1
ascentAMag2 = xyGondAcc2(ascentPeriod2);
avgAMag1 = zeros(size(ascentAMag1)); %Rolling Average Acceleration Magnitude, Flight 1
avgAMag2 = zeros(size(ascentAMag2));


ascentPMag1 = syncWindMag1(ascentPeriod1); %Absolute Dynamic Pressure Magnitude, Flight 1
ascentPMag2 = syncWindMag2(ascentPeriod2);
avgPMag1 = zeros(size(ascentPMag1)); %Rolling Average Dynamic Pressure Magnitude, Flight 1
avgPMag2 = zeros(size(ascentPMag2));

windowSize = 100; %100 datapoints either side corresponds to 5 seconds (41.05 sec)



% cross correlation using moving average 
avgAMag1 = movmean(ascentAMag1, windowSize, 'Endpoints', 'discard');
avgAMag2 = movmean(ascentAMag2, windowSize, 'Endpoints', 'discard');

avgPMag1 = movmean(ascentPMag1, windowSize, 'Endpoints', 'discard');
avgPMag2 = movmean(ascentPMag2, windowSize, 'Endpoints', 'discard');

[avg_corr1, avg_lags1]=xcorr(avgPMag1,avgAMag1,'normalized');


[avg_corr2, avg_lags2]=xcorr(avgPMag2,avgAMag2);

figure
stem(avg_lags1, avg_corr1)

%cross correlation
[corr1, lags1]=xcorr(ascentPMag1,ascentAMag1,'normalized');
figure
stem(lags1, corr1)

[max_corr1,index_max_corr1]=max(corr1);
lags1_max=lags1(index_max_corr1);
al_ascentPMag1 = circshift(ascentPMag1, lags1_max);

[r_pearson, p_pearson] = corr(al_ascentPMag1,ascentAMag1, 'Type', 'Pearson');
%% 

%Create rolling averages, using window size either side of selected point
for i = 1:length(ascentAMag1)
    if i > windowSize && i < length(ascentAMag1) - windowSize
        avgAMag1(i) = mean(ascentAMag1(i-windowSize:i+windowSize));
    elseif i <= windowSize
        avgAMag1(i) = mean(ascentAMag1(1:i));
        
    else
        avgAMag1(i) = mean(ascentAMag1(i:end));
    end
    q = i;
end

for i = 1:length(ascentAMag2)
    if i > windowSize && i < length(ascentAMag2) - windowSize
        avgAMag2(i) = mean(ascentAMag2(i-windowSize:i+windowSize));
    elseif i <= windowSize
        avgAMag2(i) = mean(ascentAMag2(1:i));
    else
        avgAMag2(i) = mean(ascentAMag2(i:end));
    end
    q = i;
end



for i = 1:length(ascentPMag1)
    if i > windowSize && i < length(ascentPMag1) - windowSize
        avgPMag1(i) = mean(ascentPMag1(i-windowSize:i+windowSize));
    elseif i <= windowSize
        avgPMag1(i) = mean(ascentPMag1(1:i));
        
    else
        avgPMag1(i) = mean(ascentPMag1(i:end));
    end
    q = i;
end

for i = 1:length(ascentPMag2)
    if i > windowSize && i < length(ascentPMag2) - windowSize
        avgPMag2(i) = mean(ascentPMag2(i-windowSize:i+windowSize));
    elseif i <= windowSize
        avgPMag2(i) = mean(ascentPMag2(1:i));
    else
        avgPMag2(i) = mean(ascentPMag2(i:end));
    end
    q = i;
end

figure
plot(syncWindMag1(ascentPeriod1), xyGondAcc1(ascentPeriod1), '.r', 'MarkerSize', 0.01)
makeLabels("Calculated Windspeed 1 vs Combined XY Acceleration, Flight 1 Ascent", "Calculated Windspeed (m/s)", "Combined XY Acceleration (m/s^2)")

figure
plot(syncWindMag2(ascentPeriod2), xyGondAcc2(ascentPeriod2), '.r', 'MarkerSize', 0.01)
makeLabels("Calculated Windspeed 2 vs Combined XY Acceleration, Flight 2 Ascent", "Calculated Windspeed (m/s)", "Combined XY Acceleration (m/s^2)")

figure
plot(syncPMag1(ascentPeriod1), xyGondAcc1(ascentPeriod1), '.r', 'MarkerSize', 0.2)
makeLabels("Calculated Dynamic Pressure vs Combined XY Acceleration, Flight 1 Ascent", "Differential Pressure (Pa)", "Combined XY Acceleration (m/s^2)")
%plot(syncPMag2(ascentPeriod2), xyGondAcc2(ascentPeriod2), '.r', 'MarkerSize', 0.2)
%makeLabels("Calculated Dynamic Pressure vs Combined XY Acceleration, Flight 2 Ascent", "Differential Pressure (Pa)", "Combined XY Acceleration (m/s^2)")

figure
plot(avgPMag1, xyGondAcc1(ascentPeriod1), '.r', 'MarkerSize', 0.2)
plot(avgPMag2, xyGondAcc2(ascentPeriod2), '.r', 'MarkerSize', 0.2)
makeLabels("Rolling Average Dynamic Pressure vs Combined XY Acceleration, Flight 1 Ascent", "5s Rolling Average Dynamic Pressure (Pa)", "Combined XY Acceleration (m/s^2)")

figure
plot(gondolaTime1(ascentPeriod1)-gondolaTime1(ascentPeriod1(1)), avgPMag1)
makeLabels("Rolling Average Dynamic Pressure vs Time Since Launch, Flight 1 Ascent", "Time Since Launch(s)", "5s Rolling Average Dynamic Pressure (Pa)")
plot(gondolaTime2(ascentPeriod2)-gondolaTime2(ascentPeriod2(1)), avgPMag2)
makeLabels("Rolling Average Dynamic Pressure vs Time Since Launch, Flight 2 Ascent", "Time Since Launch(s)", "5s Rolling Average Dynamic Pressure (Pa)")
