clc
clear
m_G=2.418+2.2355;    %gondola  mass (plate + avionic)[kg]
l_t=sqrt(0.066^2+2.55^2);    %tethers lenghts
k=83e2* 0.00635/l_t;      %elastic coefficent tethers %83e9
c=0.1;      %damping coefficent tethers
C_D_G=1.3;  %gondola drag coefficient  (considering as a rectangle with d/h=3)
S_G=0.5 *0.134 ;    %gondola normal to flow surface [m^2] (outer diameter of the plate * height of the plate)
m_BCM=0.880+1.822;    %BCM mass (plate + avionic) [kg]
C_D_B=1.1;  %BCM drag coefficient (considering as a rectangle with d/h=2.8)
S_BCM=0.375*0.134 ;    %BCM normal to flow surface   [m^2] (outer diameter of the plate * height of the plate) 
lambda=0.15;    %linear coffcient of evolution of absolute vertical wind on gondola 
sigma=0.1;
I_G=  [0.075830966	  0.000025218	-0.001299338;
	   0.000025218	  0.150474758    0.000286682;
	  -0.001299338    0.000286682	 0.076836793];  %gondola inertia matrix [kg*m^2]
I_B=   [0.012183291   0.000044346   -0.000054387;
        0.000044346   0.020870987    0.000048206;
       -0.000054387   0.000048206    0.012163675];  %BCM inertia matrix [kg*m^2]
radius_G=0.204; %Radius of attachment point from center of gondola plate [m]
radius_BCM=0.138; %Radius of attachment point from center of BCM plate [m]
alfa= 2*asin(radius_BCM*sin(60*pi/180)/l_t);         %angle betwen the v shape of tether
delta= acos(0.066/l_t);         %angle betwen  tether and gondola plane
rho_gas=0.166; %Helium density [kg/m^3]
VolB=45.0057; %4/3*pi*r; %balloon as a sphere with a radius of 12 m 
g=[0; 0; -9.81]; %gravity vector [m/s^2]
Id=eye(4);

MC=[radius_BCM 0 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0 ;
    radius_BCM 0 0 radius_G*cos(300*pi/180) radius_G*sin(300*pi/180) 0;
    radius_BCM*cos(120*pi/180) radius_BCM*sin(120*pi/180) 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0; 
    radius_BCM*cos(120*pi/180) radius_BCM*sin(120*pi/180) 0 radius_G*cos(180*pi/180) radius_G*sin(180*pi/180) 0;
    radius_BCM*cos(240*pi/180) radius_BCM*sin(240*pi/180) 0 radius_G*cos(180*pi/180) radius_G*sin(180*pi/180) 0;
    radius_BCM*cos(240*pi/180) radius_BCM*sin(240*pi/180) 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0 ];

% s(:,1)=[0;0;0;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2);0;0;0;0.0000001;0.0000001;0.0000001;0;0;0.000000001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
% l_t_ref=norm((MC(1,1:3)'+s(4:6,1)-(MC(1,4:6)'+s(1:3,1)) - (MC(1,1:3)'+s(4:6,1))-(MC(1,4:6)'+s(1:3,1))));
dt=0.001;



% [t,s]=ode45(@(t,s) model_function3D(t,s,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,l_t_ref,MC,I_G,I_B,Id),[1:dt:100],s(:,1));

%model 1D traslation dynamic
s(:,1)=[0;0;0;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2);0;0;0;0.0000001;0.0000001;0.0000001;0;0;0.1];
[t,s]=ode45(@(t,s) model_function(t,s,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma),[1:dt:100],s(:,1));
s=s';

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
% figure
% plot(t(1:i),F(7:9,:))
% legend('F_ae_G_x','F_ae_G_y','F_ae_G_z')
% figure
% plot(t(1:i),F(10:12,:))
% legend('F_ae_B_x','F_ae_B_y','F_ae_B_z')
% 

% figure
% plot(t,s(13,:))
% legend('wind_x')
% figure
% plot(t,s(14,:))
% legend('wind_y')
% figure
% plot(t,s(15,:))
% legend('wind_z')


