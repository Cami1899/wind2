%Script to get differential absolute roll, pitch, & yaw for the two flights

bcmTime1 = CPU2GPS('bcm1', bcmPimu1(:, 1));
bcmTime2 = CPU2GPS('bcm2', bcmPimu2(:, 1));

gondolaTime1 = CPU2GPS('gondola1', gondolaPimu1(:, 1));
gondolaTime2 = CPU2GPS('gondola2', gondolaPimu2(:, 1));

dataMatrix = bcmPins1;

allQuat = zeros(length(dataMatrix), 4);
allAngles = zeros(length(dataMatrix), 3);

newW = dataMatrix(:,7);
newX = dataMatrix(:,8);
newY = dataMatrix(:,9);
newZ = dataMatrix(:,10);

newQuat = quaternion(newW, newX, newY, newZ);
%allQuat(i, :) = newQuat;

[xAbs, yAbs, zAbs] = quat2angle(quatinv(newQuat), 'XYZ');
allAngles(:, :) = [xAbs yAbs zAbs];
allAngles = [dataMatrix(:, 3) allAngles];

for i=1:length(bcmAngles2)
    if bcmAngles2(i, 2) < 0
        bcmAngles2(i, 2) = bcmAngles1(i, 2) + pi;
    else
        bcmAngles2(i, 2) = bcmAngles1(i, 2) - pi;
    end
end