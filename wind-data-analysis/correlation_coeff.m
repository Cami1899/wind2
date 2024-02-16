clc
% this script need to run processWindReworked before

windMagnitude1 = hypot(diff13Wind1, diff24Wind1);
syncWindMag1 = sync12(pressures1(:, 1), windMagnitude1, CPU2GPS('gondola1', gondolaPimu1(:, 1)));

windMagnitude2 = hypot(diff13Wind2, diff24Wind2);
syncWindMag2 = sync12( pressures2(:, 1), windMagnitude2, CPU2GPS('gondola2', gondolaPimu2(:, 1)));

%Absolute XY linear Acceleration magnitude
xyGondAcc1 = hypot(gondolaPimu1(:, 7), gondolaPimu1(:, 8));
xyGondAcc2 = hypot(gondolaPimu2(:, 7), gondolaPimu2(:, 8));


%only the ascent flight
ascentPeriod1 = 97400:261600;
ascentPeriod2 = 97070:271250;

ascentAMag1 = xyGondAcc1(ascentPeriod1); 
ascentAMag2 = xyGondAcc2(ascentPeriod2);
ascentWindMag1 = syncWindMag1(ascentPeriod1);
ascentWindMag2 = syncWindMag2(ascentPeriod2);

%Normalized vector
vec1_normalized = (ascentAMag2 - mean(ascentWindMag2)) / std(ascentAMag2);
vec2_normalized = (ascentWindMag2 - mean(ascentWindMag2)) / std(ascentWindMag2);

% Pearson correlation coefficient
[r_pearson1, p_pearson1] = corr(ascentWindMag1,ascentAMag1, 'Type', 'Pearson');

[r_pearson2, p_pearson2] = corr(ascentWindMag2,ascentAMag2, 'Type', 'Pearson');

% [r_pearson, p_pearson] = corr(vettore1_normalized,vettore2_normalized, 'Type', 'Pearson');

% Spearman correlation coefficient
[r_spearman1, p_spearman1] = corr(ascentWindMag1,ascentAMag1, 'Type', 'Spearman');
[r_spearman2, p_spearman2] = corr(ascentWindMag2,ascentAMag2, 'Type', 'Spearman');

%Cross correlation: in order to find the delay between the gust application
%and the consequential acceleration 
[corr1, lags1]=xcorr(ascentWindMag1,ascentAMag1,'normalized');
figure
stem(lags1, corr1)

[max_corr1,index_max_corr1]=max(corr1);
lags1_max=lags1(index_max_corr1);
al_ascentAMag1 = circshift(ascentAMag1, lags1_max);

[r_pearson_shift, p_pearson_shift] = corr(ascentWindMag1(1:end-2),al_ascentAMag1(3:end), 'Type', 'Pearson')
%pearson coefficient doesn't change moving the values 






