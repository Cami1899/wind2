%Fourier transform script - currently in development, not yet complete
%September 2023 - Milo Eirew
%Stack does not use uniform sampling time, so use nufft rather than fft

%Vector indices corresponding to ascent: between launch and burst
ascentPeriod1 = 97400:261350;
ascentPeriod2 = 97070:271250;

%Gondola 1 Angular Rates
gondolaTime1 = CPU2GPS('gondola1', gondolaPimu1(:, 1));
gondolaRoll1 = gondolaPimu1(:, 4);
gondolaPitch1 = gondolaPimu1(:, 5);
gondolaYaw1 = gondolaPimu1(:, 6);

%BCM 1 Angular Rates (time-synced & frame-aligned)
bcmRoll1 = bcmSyncARates1(:, 2);
bcmPitch1 = bcmSyncARates1(:, 3);
bcmYaw1 =  bcmSyncARates1(:, 4);

%Gondola 2 Angular Rates (time-synced & frame-aligned)
gondolaTime2 = CPU2GPS('gondola2', gondolaPimu2(:, 1));
gondolaRoll2 = gondolaPimu2(:, 4);
gondolaPitch2 = gondolaPimu2(:, 5);
gondolaYaw2 = gondolaPimu2(:, 6);

%BCM 2 Angular Rates (time-synced & frame-aligned)
bcmRoll2 = bcmSyncARates2(:, 2);
bcmPitch2 = bcmSyncARates2(:, 3);
bcmYaw2 =  bcmSyncARates2(:, 4);

%Differential Angular Rates
dRoll1 = bcmRoll1 - gondolaRoll1;
dPitch1 = bcmPitch1 - gondolaPitch1;
dYaw1 = bcmYaw1 - gondolaYaw1;

dRoll2 = bcmRoll2 - gondolaRoll2;
dPitch2 = bcmPitch2 - gondolaPitch2;
dYaw2 = bcmYaw2 - gondolaYaw2;

len1 = length(gondolaTime1(ascentPeriod1));
gondolaFreq1 = 800*(0:len1-1)/len1;
%gondolaFreq1 = 0:0.01:1000;

gondolaFreq2 = 0:0.01:200;

%Fourier Transforms - incomplete so far
%Gondola Roll Rate 1, Pitch Rate 1, Yaw Rate 1
GR1 = nufft(gondolaRoll1(ascentPeriod1), gondolaTime1(ascentPeriod1), gondolaFreq1);
GP1 = nufft(gondolaPitch1(ascentPeriod1), gondolaTime1(ascentPeriod1), gondolaFreq1);
GY1 = nufft(gondolaYaw1(ascentPeriod1), gondolaTime1(ascentPeriod1), gondolaFreq1);

GR2 = nufft(gondolaRoll2(ascentPeriod2), gondolaTime1(ascentPeriod2), gondolaFreq2);
GP2 = nufft(gondolaPitch2(ascentPeriod2), gondolaTime1(ascentPeriod2), gondolaFreq2);
GY2 = nufft(gondolaYaw2(ascentPeriod2), gondolaTime1(ascentPeriod2), gondolaFreq2);

plot(gondolaFreq1, GR1)