
clc
clear 
close 
%This script fit the ACF of data with a weighted sum of exponential and/or damped sinusoidal functions 
% using a non linera regression model. The max number of
%components of the sum evaluate for the fit is 10, using a while cycle with
%a stop method based on the Root Mean Squared Error and R-squared. In order to find the best initial
%coefficients for the fit a gentic algorithm is used.
tic
%Import data from flight test
load('windvalues1.mat');
load("windvalues2.mat")
load("time1.mat");
load("time2.mat");
wind_speed=new_magnitude2(1:3593); 
% time=new_time2(1:17533); %index 17533=10min  %index 3593=2min
%Compute the ACF of data
maxLag=1000;
acf=autocorr(wind_speed,'NumLag', maxLag);
lag=(0:maxLag)'; 
allfunction={@(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5)*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x)+coef(16).*exp(-coef(17).*x).*cos(coef(18).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x)+coef(16).*exp(-coef(17).*x).*cos(coef(18).*x)+coef(19).*exp(-coef(20).*x).*cos(coef(21).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x)+coef(16).*exp(-coef(17).*x).*cos(coef(18).*x)+coef(19).*exp(-coef(20).*x).*cos(coef(21).*x)+coef(22).*exp(-coef(23).*x).*cos(coef(24).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x)+coef(16).*exp(-coef(17).*x).*cos(coef(18).*x)+coef(19).*exp(-coef(20).*x).*cos(coef(21).*x)+coef(22).*exp(-coef(23).*x).*cos(coef(24).*x)+coef(25).*exp(-coef(26).*x).*cos(coef(27).*x),
    @(coef,x)coef(1).*exp(-coef(2).*x).*cos(coef(3).*x)+coef(4).*exp(-coef(5).*x).*cos(coef(6).*x)+coef(7).*exp(-coef(8).*x).*cos(coef(9).*x)+coef(10).*exp(-coef(11).*x).*cos(coef(12).*x)+coef(13).*exp(-coef(14).*x).*cos(coef(15).*x)+coef(16).*exp(-coef(17).*x).*cos(coef(18).*x)+coef(19).*exp(-coef(20).*x).*cos(coef(21).*x)+coef(22).*exp(-coef(23).*x).*cos(coef(24).*x)+coef(25).*exp(-coef(26).*x).*cos(coef(27).*x)+coef(28).*exp(-coef(29).*x).*cos(coef(30).*x)}; %all the weighted sum are defined
Rsq=inf; 
R_squared=0;
i=1;
while Rsq>0.002||R_squared<0.99
%Genetic algorith evaluate the min of the difference between real acf
%data and the fit function. Implementing more simulations of the genetic 
% algorith (kga) can improove the best extimation of the minimum.  
% To modify the velocity of ga chanche the population size
normfunction=@(coef) norm(acf - allfunction{i}(coef,lag));
Parms = 3*i;
optsAns = optimoptions(@ga, 'PopulationSize', PopSz, 'InitialPopulationMatrix', randi(1E+4, PopSz, Parms)*1E-3, 'MaxGenerations', 5E3, 'FunctionTolerance', 1E-10); %define the option of genetic algorithm 
for kga = 1:1
    [theta,fval,exitflag,output,population,scores] = ga(normfunction, Parms, [],[],[],[],zeros(Parms,1),Inf(Parms,1),[],[],optsAns);
    GAdata(kga,1:i*3+1) = [fval theta(:).'];
end
[MaxFtns,idx] = min(GAdata(:,1)); 
coef0 = GAdata(idx, 2:end) ;   %initial coefficients from ga
%Non linear regression fit
functions=@(coef,x) allfunction{i}(coef,x);
ACFmodel = fitnlm(lag, acf, functions, coef0) 
coef(1:3*i,i) = ACFmodel.Coefficients.Estimate; %the coefficient of the fit
Rsq=ACFmodel.RMSE;
R_squared = ACFmodel.Rsquared.Ordinary;
Rsq_vec(i)=ACFmodel.RMSE;
R_squared_vec(i) = ACFmodel.Rsquared.Ordinary;

%Evaluate the fit results
Cfit1 = allfunction{i}(coef(1:3*i,i), lag);
figure(i)
hd = plot(lag, acf, 'p', 'DisplayName','Data');
for k1 = 1:size(acf,2)
    CV(k1,:) = hd(k1).Color;
    hd(k1).MarkerFaceColor = CV(k1,:);
end
hold on
hlp = plot(lag, Cfit1, 'DisplayName',['Sum of ' num2str(i-1) ' weighted dumped exponential '],'Color','r','LineWidth',2);
xlabel('lag')
ylabel('ACF')
hold off
grid
legend('Location','best')
i=i+1
t=clock;
if i==11
    [val,index]=min(Rsq);
    coef_f=coef(:,index)
    break;    
end 
end

if i~=11
    coef_f=coef(:,i-1);
end



