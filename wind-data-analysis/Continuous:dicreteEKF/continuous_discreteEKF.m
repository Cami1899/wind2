%EXTENDED KALMAN FILTER APPLIED TO AN AEROBOT, COMBOSED BY A BALLOON AND A
%GONDOLA
%*********************************************************************************

% MEASURES/OBSERVATIONS

%  -Gondola linear accleration                             ac_mG
%  -BCM linear acceleration                                ac_mB
%  -Realtive wind velocity on Gondola only along x e y     w_mG



% STATES

%  -Gondola position                                  r_G
%  -BCM position                                      r_B
%  -Gondola linear velocity                           v_G
%  -BCM linear velocity                               v_B
%  -Absolute wind velocity on Gondola along all axis  W_G

%*********************************************************************************
clc
clear all

%MEASURES
%Load data and create the measures vectors

load("gondolaAllPimu1.mat")
gondolaAllPimu1=gondolaPimu1;
load("gondolaAllPimu2.mat")
gondolaAllPimu2=allPimu;
load("bcmAllPimu1.mat")
load("bcmAllPimu2.mat")
load("wind_matrix1.mat")


%allPimu contains angular acceleration in columns 4,5,6 and linear acceleration in columns 7,8,9
%FLIGHT 1

[v,i_l]=min(abs(CPU2GPS('bcm1',bcmAllPimu1(:,1))-398613.3)); %For Flight 1, launch happens at GPS Time = 398613.3
[v,i_b]=min(abs(CPU2GPS('bcm1',bcmAllPimu1(:,1))-402713.3)); %For Flight 1, burst happens at GPS Time = 402713.3

[v,i_lG]=min(abs(CPU2GPS('gondola1',gondolaAllPimu1(:,1))-398613.3)); %For Flight 1, launch happens at GPS Time = 398613.3
[v,i_bG]=min(abs(CPU2GPS('gondola1',gondolaAllPimu1(:,1))-402713.3)); %For Flight 1, burst happens at GPS Time = 402713.3

t_B=CPU2GPS('bcm1',bcmAllPimu1(i_l:i_b,1));

ac_mB=[bcmAllPimu1(i_l:i_b,7)'; bcmAllPimu1(i_l:i_b,8)' ; bcmAllPimu1(i_l:i_b,9)' ];
ar_mB=[bcmAllPimu1(i_l:i_b,4)'; bcmAllPimu1(i_l:i_b,5)' ; bcmAllPimu1(i_l:i_b,6)' ];


ac_mG=[gondolaAllPimu1(i_lG:i_bG,7)';gondolaAllPimu1(i_lG:i_bG,8)' ; gondolaAllPimu1(i_lG:i_bG,9)' ];
ar_mG=[gondolaAllPimu1(i_lG:i_bG,4)'; gondolaAllPimu1(i_lG:i_bG,5)' ; gondolaAllPimu1(i_lG:i_bG,6)' ];
t_G=CPU2GPS('gondola1',gondolaAllPimu1(i_lG:i_bG,1));


windG=[wind_matrix(:,2)';wind_matrix(:,3)'];
t_wind=wind_matrix(:,1);

% time_vectors = {t_B, t_G, t_wind};
% time_vectors = {t_G};
time_vectors = {t_B, t_G};
times= sort(cat(1, time_vectors{:}));
t=unique(times-times(1));
t_B=t_B-times(1);
t_G=t_G-times(1);
t_wind=t_wind-times(1);
start_time=t(1);

%*************************************************************************
load("allPinsGondola1.mat")
load('bcmAllPins1.mat')


% [v,index_lB]=min(abs(CPU2GPS('bcm1',bcmAllPins1(:,1))-398613.3)); %For Flight 1, launch happens at GPS Time = 398613.3
% [v,index_bB]=min(abs(CPU2GPS('bcm1',bcmAllPins1(:,1))-402713.3)); %For Flight 1, burst happens at GPS Time = 402713.3
% 
% [v,index_lG]=min(abs(CPU2GPS('gondola1',allPinsgondola1(:,1))-398613.3)); %For Flight 1, launch happens at GPS Time = 398613.3
% [v,index_bG]=min(abs(CPU2GPS('gondola1',allPinsgondola1(:,1))-402713.3)); %For Flight 1, burst happens at GPS Time = 402713.3
% 
% t_velB=CPU2GPS('bcm1',bcmAllPins1(index_lB:index_bB,1));
% t_velB=t_velB-t_velB(1);
% t_velG=CPU2GPS('gondola1',allPinsgondola1(index_lG:index_bG,1));
% t_velG=t_velG-t_velG(1);
% vel_G_measure=[allPinsgondola1(index_lG:index_bG,11)';allPinsgondola1(index_lG:index_bG,12)';allPinsgondola1(index_lG:index_bG,13)'];
% vel_B_measure=[bcmAllPins1(index_lB:index_bB,11)';bcmAllPins1(index_lB:index_bB,12)';bcmAllPins1(index_lB:index_bB,13)'];


% figure
% plot(t_velB,vel_B_measure(3,:))
%*************************************************************************

%SYSTEM DATA

m_G=2.418+2.2355;    %gondola  mass (plate + avionic)[kg]
l_t=sqrt(0.066^2+2.55^2);    %tethers lenghts
k=883e6* 0.00635^2*pi/l_t;      %elastic coefficent tethers %83e9
c=0.1;      %damping coefficent tethers
radius_G=0.204; %Radius of attachment point from center of gondola plate [m]
radius_BCM=0.138; %Radius of attachment point from center of BCM plate [m]
C_D_G=1.3;  %gondola drag coefficient  (considering as a rectangle with d/h=3)
S_G=2*radius_G *0.134;   %gondola normal to flow surface [m^2] (outer diameter of the plate * height of the plate)
% C_D_G=1.3;  %gondola drag coefficient  (considering as a rectangle with d/h=3)
% S_G=[2*radius_G *0.134;2*radius_G *0.134;radius_G^2*pi] ;    %gondola normal to flow surface [m^2] (outer diameter of the plate * height of the plate)
m_BCM=0.880+1.822;    %BCM mass (plate + avionic) [kg]
C_D_B= 0.38 %0.47;  %BCM drag coefficient (considering as a rectangle with d/h=2.8)
lambda=0 %0.15;  
sigma=0 %0.1; %linear coffcient of evolution of absolute wind on gondola 
I_G=  [0.075830966	 0.000025218	    -0.001299338;
	    0.000025218	   0.150474758   0.000286682;
	   -0.001299338   0.000286682	  0.076836793];  %gondola inertia matrix [kg*m^2]
I_B=   [0.012183291 0.000044346 -0.000054387;
        0.000044346 0.020870987 0.000048206;
        -0.000054387 0.000048206 0.012163675];  %BCM inertia matrix [kg*m^2]
alfa= 2*asin(radius_BCM*sin(60*pi/180)/l_t);         %angle betwen the v shape of tether
delta= acos(0.066/l_t);         %angle betwen  tether and gondola plane
rho_gas=0.166; %Helium density [kg/m^3]
[~, ~, ~, rho_atm_balance] = atmosisa(double(32000)); %balance point at 32 km
VolB=abs((m_G+m_BCM)/(rho_atm_balance-rho_gas)); %45.0045;  %balloon as a sphere with a radius of 12 m 
r_balloon=(3*VolB/(4*pi))^(1/3);
S_BCM=0.375*0.134
% S_BCM=[0.375*0.134+r_balloon^2*pi;0.375*0.134+r_balloon^2*pi;r_balloon^2*pi;]; %BCM normal to flow surface   [m^2] (outer diameter of the plate * height of the plate+ sphere surface normal to flux) 
g=[0; 0; -9.81]; %gravity vector [m/s^2]
Id=eye(4);

%*************************************************************************

state_count=3*5;
measurement_count= 100 %length(t);
sensor_count=3*2+2;

% SYMBOLIC FUNCTIONS (symbolic is used just to generate the Jacobian matrix)

%MODEL FUNCTION


syms r_G1 r_G2 r_G3 r_B1 r_B2 r_B3 v_G1 v_G2 v_G3 v_B1 v_B2 v_B3 W_G1 W_G2 W_G3 real
syms rhoG rhoB  real


rG = [r_G1; r_G2; r_G3];
vG = [v_G1; v_G2; v_G3];
rB = [r_B1; r_B2; r_B3];
vB = [v_B1; v_B2; v_B3];
WG = [W_G1; W_G2; W_G3];

%p=0.0001



F_el = k * (rB-rG);                 %elastic force
F_d=c * (vB-vG);                          %dumping force
F_t = 6*(F_el + F_d);                                   %tether force
F_a_b = -sign(vB - WG) * 0.5 * rhoB * C_D_B .* S_BCM .* (vB - WG).^2;  %aerodinamic force on gondola
F_a_G =-sign(vG - WG) * 0.5 * rhoG * C_D_G .* S_G .* (vG - WG).^2;    %aerodinamic force on bcm
Fa=-(rhoB-rho_gas)*VolB*g;                                 %buoyancy force
F_ext_BCM = Fa -  F_t + F_a_b + m_BCM*g;          %externak forces on bcm
F_ext_G = F_t + F_a_G + m_G*g;                    %externak forces on gondola

f_sym=[vG;
       vB;
       F_ext_G /m_G;
       F_ext_BCM /m_BCM;
       sigma  ;
       lambda ;
       lambda ];

J_f=jacobian(f_sym,[rG;rB;vG;vB;WG]);
J_f_num= matlabFunction(J_f, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3],  rhoG, rhoB});


%MODEL FUNCTION




h_sym=[F_ext_G /m_G;
       F_ext_BCM /m_BCM;
       v_G1 - W_G1;
       v_G2 - W_G2 ];



%Not all measures are available at each time istant, 7 possibilities:

%1)all observations
h_sym1=h_sym;
J_h1=jacobian(h_sym1,[rG;rB;vG;vB;WG]);
J_h1_num= matlabFunction(J_h1, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%2)only gondola IMU observations 
h_sym2=F_ext_G /m_G;

J_h2=jacobian(h_sym2,[rG;rB;vG;vB;WG]);
J_h2_num= matlabFunction(J_h2, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%3)only bcm IMU observations
h_sym3=F_ext_BCM /m_BCM;

J_h3=jacobian(h_sym3,[rG;rB;vG;vB;WG]);
J_h3_num= matlabFunction(J_h3, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%4)only wind observations
h_sym4=[v_G1 - W_G1;
       v_G2 - W_G2 ];
J_h4=jacobian(h_sym4,[rG;rB;vG;vB;WG]);
J_h4_num= matlabFunction(J_h4, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%5)only wind and gondola IMU
h_sym5=[F_ext_G /m_G;
       v_G1 - W_G1;
       v_G2 - W_G2 ];
J_h5=jacobian(h_sym5,[rG;rB;vG;vB;WG]);
J_h5_num= matlabFunction(J_h5, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%6)only wind and bcm IMU
h_sym6=[F_ext_BCM /m_BCM;
       v_G1 - W_G1;
       v_G2 - W_G2 ];
J_h6=jacobian(h_sym6,[rG;rB;vG;vB;WG]);
J_h6_num= matlabFunction(J_h6, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});

%7)only gondola and bcm IMU
h_sym7=[F_ext_G /m_G;
       F_ext_BCM /m_BCM];
J_h7=jacobian(h_sym7,[rG;rB;vG;vB;WG]);
J_h7_num= matlabFunction(J_h7, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3], rhoG, rhoB});



 





%PROCESS NOISE COVARIANCE MATRIX
 q=[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,8,8,8];
%  q=[0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,];
% q=[5,5,5,5,5,5,5,5,5,5,5,5,8,8,8];
 %q=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
%  q=zeros(state_count,1);
 Q=diag(q);








%MEASUREMENT NOISE COVARIANCE MATRIX

angular_velocity_noise_density = 5e-3;  % 5 mdps/√Hz
acceleration_noise_density = 60e-6;  % 60 μg/√Hz

% Sampling rate
angular_velocity_sampling_frequency = 8e3;  % 8 KHz
acceleration_sampling_frequency = 4e3;  % +4 KHz

% Variance
ar_variance = (angular_velocity_noise_density) * sqrt(angular_velocity_sampling_frequency);
ac_variance = (acceleration_noise_density) * sqrt(acceleration_sampling_frequency)*9.81;


err_wind=[0.01,0.01];



%*************************************************************************

%EKF

%State and P initialization

vel_0G=(ac_mG(:,2)-ac_mG(:,1))*(t_G(2)-t_G(1)); 
vel_0B=(ac_mB(:,2)-ac_mB(:,1))*(t_B(2)-t_B(1)); 
vel_1G=(ac_mG(:,3)-ac_mG(:,2))*(t_G(3)-t_G(2)); 
vel_1B=(ac_mB(:,3)-ac_mB(:,2))*(t_B(3)-t_B(2));
x_0G=(vel_1G-vel_0G)*(t_G(3)-t_G(2));
x_0B=(vel_1B-vel_0B)*(t_B(3)-t_B(2));
z_0G=5*1.60;
% x_0=[0;0;z_0G;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2)+z_0G;vel_0G;vel_0B;vel_0G(1,1)-windG(1);vel_0G(1,1)-windG(2);0.1];
% x_0=[x_0G;x_0B;vel_0G;vel_0B;vel_0G(1,1)-windG(1);vel_0G(1,1)-windG(2);0.1];
x_0=[0;0;z_0G;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2)+z_0G;3;3;3;3;3;3;vel_0G(1,1)-windG(1);vel_0G(1,1)-windG(2);0.1];
wind_std_xy=std(windG,0,2);
acc_std_G=std(ac_mG,0,2);
acc_std_B=std(ac_mB,0,2);

p=[0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,0.000001,];
% p=[acc_std_G'*4,acc_std_B'*4,acc_std_G'*2,acc_std_B'*2,acc_std_G(1:2)'-wind_std_xy',0.1];
% p=[acc_std_G'*4,acc_std_B'*4,acc_std_G',acc_std_B',2,2,2];
P=diag(p);
diag_P(:,1) = p;

%Observability study 
% 
% [TT, a, pp, rhoG_val] = atmosisa(double(x_0(3)));
% [TT, a, pp, rhoB_val] = atmosisa(double(x_0(6)));
% eig_values=eig(J_f_num(x_0,rhoG_val,rhoB_val))
% O=obsv(J_f_num(x_0,rhoG_val,rhoB_val),J_h1_num(x_0,rhoG_val,rhoB_val));
% rank(O)
% null_basis=null(O)
% [U, S, V] = svd(O);
% 
% % Find the unobservable states
% unobservable_states = [];
% for i = 1:size(V, 2)
%     if S(i, i) < 1e-5 
%         unobservable_states = [unobservable_states, V(:, i)];
%     end
% end
% 
% unobservable_states_indices = find(diag(S) < 1e-10);
% disp("Unobservable states indices: ");
% disp(unobservable_states_indices);

[~, ~, ~, rhoG_val0] = atmosisa(double(x_0(3)));
[~, ~, ~, rhoB_val0] = atmosisa(double(x_0(6)));
J_0=J_f_num(x_0,rhoG_val0,rhoB_val0);

P_dot_0=J_0*P+P*J_0'+Q;
% P_P=integral(@(t)P_dot_0,t(1),0.5,'ArrayValued', true)



fun = @(t, P) reshape((J_0*reshape(P, size(J_0)) + reshape(P, size(J_0))*J_0' + Q)', [], 1);


[t, P_P] = ode45(fun, [0:0.000001:0.5], P);
P_P = reshape(P_P(end, :), size(P));
rap=max(P_P,[],'all')/min(P_P,[],'all')
sqrt(P_P(9,9))


%% 


	

x_old = x_0;
P_old = P;


   
	
	x_p(:,1) = x_0; 

	dt_int=0.0001;

%eliminate the interval between t element that is less than h (integration step ode)
   index_t = length(t);
   new_t = t(1); 
for i = 1:(index_t - 1)
    if t(i + 1) - t(i) > 2 * dt_int
        new_t = [new_t, t(i + 1)];  
    end
end 

t = new_t;  
l=0;
 

for j=1:measurement_count

        start_time=t(j);

        dt_between_measurements= t(j+1)-t(j);
        
        [x_p,P_p,P_diag,t_p] = cdekf_predict_phase(J_f_num,dt_between_measurements,start_time,P_old,x_old,Q,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma,dt_int);
        
        t_new(l+1:l+length(t_p))=t_p;
       
        [x_u,P_u] =  cdekf_update_phase(t(j+1),P_p,x_p(:,end),C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,t_G,t_B,t_wind,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num);
        x_p(:,end) = x_u;
        state(:,l+1:l+length(t_p))=x_p;
        P_diag(:,end)=diag(P_u);
        P_diag_matrix(:,l+1:l+length(t_p))=P_diag;

		x_old = x_u;
		P_old = P_u;

        diag_P(:,j+1) = diag(P_u);
        
         l=length(t_new);
end

% measurement_count=449;
std=sqrt(P_diag_matrix);



index_G=find(t_G<t(measurement_count));
index_B=find(t_B<t(measurement_count));
figure
plot(t(1:measurement_count),state(3,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(3,1:measurement_count)+3*std(3,1:measurement_count))
plot(t(1:measurement_count),state(3,1:measurement_count)-3*std(3,1:measurement_count))
legend('z-gondola')
hold off
figure
plot(t(1:measurement_count),state(6,1:measurement_count))
hold on
plot(t_B(1:index_B(end)),zeros(1,length(index_B(end))),'*')
plot(t(1:measurement_count),state(6,1:measurement_count)+3*std(6,1:measurement_count))
plot(t(1:measurement_count),state(6,1:measurement_count)-3*std(6,1:measurement_count))
legend('z-bcm')
hold off

figure
plot(t(1:measurement_count),state(2,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(2,1:measurement_count)+3*std(2,1:measurement_count))
plot(t(1:measurement_count),state(2,1:measurement_count)-3*std(2,1:measurement_count))
legend('y-gondola')
hold off
figure
plot(t(1:measurement_count),state(5,1:measurement_count))
hold on
plot(t_B(1:index_B(end)),zeros(1,length(index_B(end))),'*')
plot(t(1:measurement_count),state(5,1:measurement_count)+3*std(5,1:measurement_count))
plot(t(1:measurement_count),state(5,1:measurement_count)-3*std(5,1:measurement_count))
legend('y-bcm')
hold off

figure
plot(t(1:measurement_count),state(1,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(1,1:measurement_count)+3*std(1,1:measurement_count))
plot(t(1:measurement_count),state(1,1:measurement_count)-3*std(1,1:measurement_count))
legend('x-gondola')
hold off
figure
plot(t(1:measurement_count),state(4,1:measurement_count))
hold on
plot(t_B(1:index_B(end)),zeros(1,length(index_B(end))),'*')
plot(t(1:measurement_count),state(4,1:measurement_count)+3*std(4,1:measurement_count))
plot(t(1:measurement_count),state(4,1:measurement_count)-3*std(4,1:measurement_count))
legend('x-bcm')
hold off

figure
plot(t(1:measurement_count),state(9,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(9,1:measurement_count)+3*std(9,1:measurement_count))
plot(t(1:measurement_count),state(9,1:measurement_count)-3*std(9,1:measurement_count))
legend('v_z-gondola')
hold off

figure
plot(t(1:measurement_count),state(8,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(8,1:measurement_count)+3*std(8,1:measurement_count))
plot(t(1:measurement_count),state(8,1:measurement_count)-3*std(8,1:measurement_count))
legend('v_y-gondola')
hold off

figure
plot(t(1:measurement_count),state(7,1:measurement_count))
hold on
plot(t_G(1:index_G(end)),zeros(1,length(index_G(end))),'*')
plot(t(1:measurement_count),state(7,1:measurement_count)+3*std(7,1:measurement_count))
plot(t(1:measurement_count),state(7,1:measurement_count)-3*std(7,1:measurement_count))
legend('v_x-gondola')
hold off

index_wind=find(t_wind<t(measurement_count) );

figure
plot(t(1:measurement_count),state(15,1:measurement_count))
hold on
plot(t(1:measurement_count),state(15,1:measurement_count)+3*std(15,1:measurement_count))
plot(t(1:measurement_count),state(15,1:measurement_count)-3*std(15,1:measurement_count))
legend('wind_z')
hold off

figure
plot(t(1:measurement_count),state(8,1:measurement_count)-state(13,1:measurement_count))
hold on
plot(t_wind(1:index_wind(end)),windG(1,1:index_wind(end)),'*')
% plot(t(1:measurement_count),state(13,1:measurement_count)+3*std(13,1:measurement_count))
% plot(t(1:measurement_count),state(13,1:measurement_count)-3*std(13,1:measurement_count))
legend('wind_x')
hold off
figure
plot(t(1:measurement_count),state(7,1:measurement_count)-state(14,1:measurement_count))
hold on
plot(t_wind(1:index_wind(end)),windG(2,1:index_wind(end)),'*')
% plot(t(1:measurement_count),state(14,1:measurement_count)+3*std(14,1:measurement_count))
% plot(t(1:measurement_count),state(14,1:measurement_count)-3*std(14,1:measurement_count))
legend('wind_y')
hold off




