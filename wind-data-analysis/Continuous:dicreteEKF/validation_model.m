clc
clear
m_G=2.418+2.2355;    %gondola  mass (plate + avionic)[kg]
l_t=sqrt(0.066^2+2.55^2);    %tethers lenghts
k=83e4* 0.00635^2*pi/(4*l_t);      %elastic coefficent tethers %83e9
c=0.1;      %damping coefficent tethers
radius_G=0.204; %Radius of attachment point from center of gondola plate [m]
radius_BCM=0.138; %Radius of attachment point from center of BCM plate [m]
C_D_G=1.3;  %gondola drag coefficient  (considering as a rectangle with d/h=3)
S_G=[2*radius_G *0.134;2*radius_G *0.134;radius_G^2*pi] ;    %gondola normal to flow surface [m^2] (outer diameter of the plate * height of the plate)
m_BCM=0.880+1.822;    %BCM mass (plate + avionic) [kg]
C_D_B=0.47; %0.38 %1.1;  %BCM drag coefficient (considering as a rectangle with d/h=2.8)
    %BCM normal to flow surface   [m^2] (outer diameter of the plate * height of the plate) 
lambda=0; %0.15;    %linear coffcient of evolution of absolute vertical wind on gondola 
sigma=0;
rho_gas=0.166; %Helium density [kg/m^3]
[~, ~, ~, rho_atm_balance] = atmosisa(double(10000)); %balance point at 32 km
VolB=abs((m_G+m_BCM)/(rho_atm_balance-rho_gas)); %45.0045;  %balloon as a sphere with a radius of 12 m  
r_balloon=(3*VolB/(4*pi))^(1/3);
S_BCM=[0.375*0.134+r_balloon^2*pi;0.375*0.134+r_balloon^2*pi;r_balloon^2*pi;]; %BCM normal to flow surface   [m^2] (outer diameter of the plate * height of the plate+ sphere surface normal to flux) 
g=[0; 0; -9.81]; %gravity vector [m/s^2]
Id=eye(4);
v_ref=1.5;
z_ref=700;



% s(:,1)=[0;0;0;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2);0;0;0;0.0000001;0.0000001;0.0000001;0;0;0.000000001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
% l_t_ref=norm((MC(1,1:3)'+s(4:6,1)-(MC(1,4:6)'+s(1:3,1)) - (MC(1,1:3)'+s(4:6,1))-(MC(1,4:6)'+s(1:3,1))));
dt=0.001;


% syms z v_ref z_ref gamma
% g_sym=v_ref*(z/z_ref).^gamma;
% g_dot_sym=diff(g_sym,z)
% g_dot=matlabFunction(g_dot_sym)


% [t,s]=ode45(@(t,s) model_function3D(t,s,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,l_t_ref,MC,I_G,I_B,Id),[1:dt:100],s(:,1));

%model 1D traslation dynamic
 
s(:,1)=[0;0;0;0.000000000;0.000000000;l_t;0;0;0;0;0.000000;0.000000;0;0;0];
[t,s]=ode45(@(t,s) model_function(s,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,v_ref,z_ref),[1:dt:2000],s(:,1));
s=s';
 


a_G=(s(7:9,2:end)-s(7:9,1:end-1))./dt;
a_BCM=(s(10:12,2:end)-s(10:12,1:end-1))./dt;


% [TT, a, pp, rhoG_val] = atmosisa(double(s(3,i)));
% [TT, a, pp, rhoB_val] = atmosisa(double(s(6,i)));
% if (norm(s(4:6,i)-s(1:3,i))-l_t <= 0)   
%                 F_el = [0;0;0];
%                 F_d = [0;0;0];
% else  
% 
% F_el = k * (s(4:6,i)-s(1:3,i)); 
% F_d= c * s(10:12,i)-s(7:9,i); 
% end
% 
% F_t = (F_el + F_d);
% Fae_b = -sign(s(10:12,i) - s(13:15,i)) * 0.5 * rhoB_val * C_D_B .* S_BCM .* (s(10:12,i) - s(13:15,i)).^2;
% Fae_G = -sign(s(7:9,i) - s(13:15,i)) * 0.5 * rhoG_val * C_D_G .* S_G .* (s(7:9,i) - s(13:15,i)).^2;
% Fa= -(rhoB_val-rho_gas)*VolB*g;
% F_ext_BCM = Fa -  F_t + Fae_b + m_BCM*g;
% F_ext_G = F_t + Fae_G + m_G*g;
% w_G=s(7:9,i) - s(13:15,i);
% a_G(:,i)=F_ext_G /m_G;
% a_BCM(:,i)=F_ext_BCM/m_BCM;
% end

% figure
% plot(t,a_G(3,:))
figure
plot(t,s(6,:));
legend('z_bcm')
figure
plot(t,s(12,:));
legend('v_z_bcm')
figure
plot(t(1:end-1),a_BCM(3,:))
legend('acc_z_bcm')
figure
plot(t,s(15,:));
legend('wind_z')
%% 

sigma=[0.1;0.1;0.1];
measure_G=a_G+sigma.*rand(3,length(t));
measure_BCM=a_BCM+sigma.*rand(3,length(t));
measure_wind=w_G+0.3*rand(3,length(t));

% figure
% plot(t(1:1000),w_G(1,1:1000),'k')
% hold on 
% plot(t(1:1000),measure_wind(1,1:1000),'--')

% figure
% plot(t(1:1000),a_G(3,1:1000),'k')
% hold on 
% plot(t(1:1000),measure_G(3,1:1000),'--')

% x_0=[0;0;0;0.000000000;0.000000000;l_t;0;0;0;0;0.000000;0.000000;0;0;0];
% for i=1:length(t)
% measures(:,i)=model_function(x_0,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,v_ref,z_ref)
% x_0=measures(:,i)
% end 


figure
plot(t,s(1:2,:));
legend('x_gondola','y_gondola')
figure
plot(t,s(4:5,:));
legend('x_bcm','y_bcm')

figure
plot(t,s(3,:));
legend('z_gondola')
figure
plot(t,s(6,:));
legend('z_bcm')

figure
plot(t,s(9,:));
legend('v_z_gondola')
figure
plot(t,s(12,:));
legend('v_z_bcm')

figure
plot(t,s(15,:));
legend('wind_z')


% g=@(x)1.5*(x/700).^0.3
% plot(g(linspace(0,32000)),linspace(0,32000))

% figure
% plot(t,s(13,:))
% legend('wind_x')
% figure
% plot(t,s(14,:))
% legend('wind_y')
% figure
% plot(t,s(15,:))
% legend('wind_z')

