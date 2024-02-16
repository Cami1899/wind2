%This script extracts absolute angles from the quaternions contained in the
%PINS2 Data from the Lucerne Valley flights on July 20th & 27th, 2023. The
%script also acts as a visualization tool to animate the motion of the
%payloads.
%Milo Eirew September 2023
dataMatrix = bcmPins1;
dataMatrix2 = gondolaPins1;

windData = [CPU2GPS('gondola1', pressures1(:, 1)) diff13Wind1, diff24Wind1];
syncWindData = [dataMatrix(:, 3), sync12(windData(:, 1), windData(:, 2), dataMatrix(:, 3)), sync12(windData(:, 1), windData(:, 3), dataMatrix(:, 3))];


allQuat = zeros(length(dataMatrix), 4);

b = 0;

[x,y,z] = cylinder([0,1,1,0],10);

for i = 1:length(z)
    z(:, i) = [0 0 0.2 0.2]';
end

%Extract quaternion data from matrices
newW = dataMatrix(:,7);
newX = dataMatrix(:,8);
newY = dataMatrix(:,9);
newZ = dataMatrix(:,10);

newW2 = dataMatrix2(:,7);
newX2 = dataMatrix2(:,8);
newY2 = dataMatrix2(:,9);
newZ2 = dataMatrix2(:,10);

newQuat = quaternion(newW, newX, newY, newZ);
newQuat2 = quaternion(newW2, newX2, newY2, newZ2);
%allQuat(i, :) = newQuat;


[xAbs, yAbs, zAbs] = quat2angle(quatinv(newQuat), 'XYZ');
allAngles = [xAbs yAbs zAbs];
allAngles = [dataMatrix(:, 3) allAngles];

[xAbs2, yAbs2, zAbs2] = quat2angle(quatinv(newQuat2), 'XYZ');
allAngles2 = [xAbs2 yAbs2 zAbs2];
allAngles2 = [dataMatrix2(:, 3) allAngles2];

allAnglesDeg = [allAngles(:, 1), rad2deg(allAngles(:, 2:4))];
allAnglesDeg2 = [allAngles2(:, 1), rad2deg(allAngles2(:, 2:4))];

windAxis = 2;

startPlotTime = 3000;
timeDivider = 10;
viewPoint = [-2 -2 0.5];
unchangedMesh = mesh(x,y,z);



for i = 1:length(dataMatrix)
    timeElapsed = syncWindData(i, 1)-syncWindData(1,1);
    if mod(i, timeDivider) == 0 && timeElapsed >= startPlotTime 
        clf
        subplot(1, 3, 1)
        %{
        rP = rPMatrix(:, :, i);
        plot3([rP(:, 1)]', [rP(:, 2)]', [rP(:, 3)]');
        axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
        title("Gondola, Flight 1: "+(dataMatrix(i, 3)-dataMatrix(1, 3)) + " Seconds Since Startup")
        grid
        view(viewPoint)
        %}
        
        
        originalMesh = mesh(x,y,z);
        

        rotate(originalMesh, [1 0 0], allAnglesDeg(i, 2))
        rotate(originalMesh, [0 1 0], allAnglesDeg(i, 3))
        rotate(originalMesh, [0 0 1], allAnglesDeg(i, 4))

        axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
        title("Gondola, Flight 1: "+(dataMatrix(i, 3)-dataMatrix(1, 3)) + " Seconds Since Startup")
        grid
        view(viewPoint)
        %{
        subplot(1, 3, 2)
        originalMesh2 = mesh(x,y,z);
        rotate(originalMesh2, [1 0 0], allAnglesDeg2(i, 2))
        rotate(originalMesh2, [0 1 0], allAnglesDeg2(i, 3) + 180)
        rotate(originalMesh2, [0 0 1], allAnglesDeg2(i, 4) + 75)

        axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
        title("BCM, Flight 1")
        grid
        view(viewPoint)
        %}
        subplot(1, 3, 3)
           
        invisPlot=compass([0 windAxis 0 -windAxis], [windAxis 0 -windAxis 0]);
        set(invisPlot,'Visible','off');
        hold on;
        compass(syncWindData(i, 2), syncWindData(i, 3));
        title("Time Elapsed = " + timeElapsed)
        view(-allAnglesDeg(i, 4), 90)
        axis vis3d
        drawnow

    end
end

%{
for i = 1:length(dataMatrix)
    elapsedTime
    if(mod(i, 20) == 0 && dataMatrix(i, 3)-dataMatrix(1, 3) > 5000)
        rP = rPMatrix(:, :, i);
        plot3([rP(:, 1)]', [rP(:, 2)]', [rP(:, 3)]');
        axis([-1.5 1.5 -1.5 1.5 -1.5 1.5]);
        title("Gondola, Flight 1: "+(dataMatrix(i, 3)-dataMatrix(1, 3)) + " Seconds Since Startup")
        grid
        drawnow
    end
end
%}


%end

%{
for i = 1:length(allAngles)
    if(allAngles(i, 2)) < 0
        allAngles(i, 2) = allAngles(i, 2) + pi;
    else
        allAngles(i, 2) = allAngles(i, 2) - pi;
    end
end
%}

plot(allAngles(:, 1)-allAngles(1, 1), allAngles(:, 3), '.r', 'MarkerSize', 0.1, 'Color', 'red')
title("Gondola Abolute Roll (Body Frame to NED) Since Startup, Flight 1", 'FontSize', 20)
xlabel("Time Since Startup (s)", 'FontSize', 15)
ylabel("BCM Abolute Roll (rad)", 'FontSize', 20)
set(gcf,'Position',[0 0 1500 500])
ylim([-3.5 3.5])


%{
rPMatrix = zeros(42, 3, length(dataMatrix));
for i = 1:length(dataMatrix)
%if b > 2500
    rP = rotatepoint(newQuat(i), [X Y Z]);
    rPMatrix(:,:, i) = rP;
    a = i
end
%}


%{
plot(allAngles(:, 1)-allAngles(1, 1), allAngles(:, 2)-allAngles(1, 2), '-o', 'MarkerSize', 1)
title("BCM Abolute Roll From Original Position vs Time Since Startup, Flight 2", 'FontSize', 20)
xlabel("Time Since Startup (s)", 'FontSize', 15)
ylabel("BCM Abolute Roll (rad)", 'FontSize', 20)
set(gcf,'Position',[0 0 1500 500])
%}
