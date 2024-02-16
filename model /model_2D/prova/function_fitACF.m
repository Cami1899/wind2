function [coef_f,num_comp,fit_function] = function_fitACF(acf,lag,PopSz)
%------------------------------------------------------------------------
%This function fit the ACF of data with a weighted sum of exponential and/or damped sinusoidal 
%functions using a non linera regression model: 
% 
% R = sum( w*exp( -alfa * t ) * cos( omega * t ) )

% The max number of components of the sum evaluated for the fit is 10, using a while cycle.
% For each iteration the best initial set of coefficients for the fit is identified using a gentic algorithm GA.
% To identify the best fit a stop method based on the Root Mean Squared Error and R-squared is implemented. 
%Finally, only the coefficients that constitute stable SDE equations are selected for integration 


%INPUT
% - acf, Auto-correlation values of real wind data;
% - lag, lag in which the ACF is evaluated;
% - PopSz, initial population size that the GA evaluates as initial coefficients for the fit.
%OUTPUT:
% - coef_f=(weight_1,alpha_1,omega_1,...,weight_n,alpha_n,omega_n),coefficients of the fit with 1,...,n number of components of the fit;
% - num_comp, number of components of the fit;
% - fit_function, final best function for the fit.
%------------------------------------------------------------------------

%1.Define all the fit superposition functions from 1 to 10 components
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

%2.Identify the best initial point in terms of coefficents for the fit using a Genetic algorithm GA

Rsq=inf; 
R_squared=0;
i=1;
while i<=10 && (Rsq>0.002||R_squared<0.99)
%Genetic algorithm evaluate the min of the difference between real acf
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



%3.Fit the ACF using a Non Linear Regression method
functions=@(coef,x) allfunction{i}(coef,x);
ACFmodel = fitnlm(lag, acf, functions, coef0); 
coef(1:3*i,i) = ACFmodel.Coefficients.Estimate; 
Rsq=ACFmodel.RMSE; %Root Mean Squared Error
R_squared = ACFmodel.Rsquared.Ordinary; %R-squared
Rsq_vec(i)=ACFmodel.RMSE;
R_squared_vec(i) = ACFmodel.Rsquared.Ordinary;

i=i+1;
 

if i==11 %If none of the fits meet the stopping criterion the one with the least error is chosen
    [val,index]=min(Rsq_vec);
    coef_f=coef(1:3*index,index)
    num_comp=index;
    fit_function=allfunction{num_comp};
    break;    
end 
end

if i~=11
    coef_f=coef(:,i-1);
    num_comp=i-1;
    fit_function=allfunction{num_comp};
end

weight=coef_f(1:3:end);
alpha=coef_f(2:3:end);
omega=coef_f(3:3:end);


%4.Filter the coefficients for integration stability

j = 1;
while j <= length(alpha)
eigenvalues=eig([-alpha(j) -omega(j); omega(j) -alpha(j)]);
if any(real(eigenvalues) > 0)
    weight(j)=[];
    alpha(j)=[];
    omega(j)=[];
    num_comp=num_comp-1;
    fit_function=allfunction{num_comp};
    coef_f(3*j-2:3*j)=[];
    else
        j = j + 1;

end
end

end