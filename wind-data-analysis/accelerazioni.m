%EXTENDED KALMAN FILTER APPLIED TO AN AEROBOT, COMBOSED BY A BALLOON AND A
%GONDOLA
%*********************************************************************************

% MEASURES/OBSERVATIONS

%  -Gondola linear accleration                             ac_mG
%  -BCM linear acceleration                                ac_mB
%  -Gondola angular rates                                  ar_mG
%  -BCM angular rates                                      ar_mB
%  -Realtive wind velocity on Gondola only along x e y     w_mG

%Sampling time isn't costant, data need to be resampled

% STATES

%  -Gondola position                                  r_G
%  -BCM position                                      r_B
%  -Gondola linear velocity                           v_G
%  -BCM linear velocity                               v_B
%  -Absolute wind velocity on Gondola along all axis  W_G
%  -Gondola angular rates                             ar_G
%  -BCM angular rates                                 ar_B
%  -Gondola angular acceleration                      aa_G
%  -BCM angular acceleration                          aa_B
%  -Gondola quaternions                               qG
%  -BCM quaternions                                   qB



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

% i_b=90000;
% i_bG=90000;

t_B=CPU2GPS('bcm1',bcmAllPimu1(i_l:i_b,1));

ac_mB=[bcmAllPimu1(i_l:i_b,7)'; bcmAllPimu1(i_l:i_b,8)' ; bcmAllPimu1(i_l:i_b,9)' ];
ar_mB=[bcmAllPimu1(i_l:i_b,4)'; bcmAllPimu1(i_l:i_b,5)' ; bcmAllPimu1(i_l:i_b,6)' ];


ac_mG=[gondolaAllPimu1(i_lG:i_bG,7)';gondolaAllPimu1(i_lG:i_bG,8)' ; gondolaAllPimu1(i_lG:i_bG,9)' ];
ar_mG=[gondolaAllPimu1(i_lG:i_bG,4)'; gondolaAllPimu1(i_lG:i_bG,5)' ; gondolaAllPimu1(i_lG:i_bG,6)' ];
t_G=CPU2GPS('gondola1',gondolaAllPimu1(i_lG:i_bG,1));


windG=[wind_matrix(:,2)';wind_matrix(:,3)'];
t_wind=wind_matrix(:,1);

time_vectors = {t_B, t_G, t_wind};
times= sort(cat(1, time_vectors{:}));
t=times-times(1);
t_B=t_B-times(1);
t_G=t_G-times(1);
t_wind=t_wind-times(1);



%*************************************************************************

%SYSTEM DATA

m_G=2.418+2.2355;    %gondola  mass (plate + avionic)[kg]
l_t=sqrt(0.066^2+2.55^2);    %tethers lenghts
k=83e2* 0.00635/l_t;      %elastic coefficent tethers %83e9
c=0.1;      %damping coefficent tethers
C_D_G=1.3;  %gondola drag coefficient  (considering as a rectangle with d/h=3)
S_G=0.5 *0.134 ;    %gondola normal to flow surface [m^2] (outer diameter of the plate * height of the plate)
m_BCM=0.880+1.822;    %BCM mass (plate + avionic) [kg]
C_D_B=1.1;  %BCM drag coefficient (considering as a rectangle with d/h=2.8)
S_BCM=0.375*0.134 ;    %BCM normal to flow surface   [m^2] (outer diameter of the plate * height of the plate) 
lambda=1.00001; 
sigma=1.00003;%linear coffcient of evolution of absolute wind on gondola 
I_G=  [0.075830966	 0.000025218	    -0.001299338;
	    0.000025218	   0.150474758   0.000286682;
	   -0.001299338   0.000286682	  0.076836793];  %gondola inertia matrix [kg*m^2]
I_B=   [0.012183291 0.000044346 -0.000054387;
        0.000044346 0.020870987 0.000048206;
        -0.000054387 0.000048206 0.012163675];  %BCM inertia matrix [kg*m^2]
radius_G=0.204; %Radius of attachment point from center of gondola plate [m]
radius_BCM=0.138; %Radius of attachment point from center of BCM plate [m]
alfa= 2*asin(radius_BCM*sin(60*pi/180)/l_t);         %angle betwen the v shape of tether
delta= acos(0.066/l_t);         %angle betwen  tether and gondola plane
rho_gas=0.166; %Helium density [kg/m^3]
VolB=45.0057;  %balloon as a sphere with a radius of 12 m 
g=[0; 0; -9.81]; %gravity vector [m/s^2]
Id=eye(4);

%*************************************************************************

% SYMBOLIC FUNCTIONS (symbolic is used just to generate the Jacobian matrix)

syms dt real positive
syms r_G1 r_G2 r_G3 r_B1 r_B2 r_B3 v_G1 v_G2 v_G3 v_B1 v_B2 v_B3 W_G1 W_G2 W_G3 ar_G1 ar_G2 ar_G3  ar_B1 ar_B2 ar_B3 aa_G1 aa_G2 aa_G3 aa_B1 aa_B2 aa_B3 q0G q1G q2G q3G q0B q1B q2B q3B real
syms rhoG rhoB  real


rG = [r_G1; r_G2; r_G3];
vG = [v_G1; v_G2; v_G3];
WG = [W_G1; W_G2; W_G3];
arG = [ar_G1; ar_G2; ar_G3];
aaG = [aa_G1; aa_G2; aa_G3];
qG = [q0G; q1G; q2G; q3G];

rB = [r_B1; r_B2; r_B3];
vB = [v_B1; v_B2; v_B3];
arB = [ar_B1; ar_B2; ar_B3];
aaB = [aa_B1; aa_B2; aa_B3];
qB = [q0B; q1B; q2B; q3B];


OmegaG=[0 -ar_G1 -ar_G2 -ar_G3; ar_G1 0 ar_G3 -ar_G2; ar_G2 -ar_G3 0 ar_G1; ar_G3 ar_G2 -ar_G1 0];
OmegaB=[0 -ar_B1 -ar_B2 -ar_B3; ar_B1 0 ar_B3 -ar_B2; ar_B2 -ar_B3 0 ar_B1; ar_B3 ar_B2 -ar_B1 0];
DCM_B_I = [ 1-2*(q2B^2+q3B^2), 2*(q1B*q2B+q3B*q0B), 2*(q1B*q3B-q2B*q0B);
        2*(q2B*q1B-q3B*q0B), 1-2*(q1B^2+q3B^2), 2*(q2B*q3B+q1B*q0B);
        2*(q3B*q1B+q2B*q0B), 2*(q3B*q2B-q1B*q0B), 1-2*(q1B^2+q2B^2)];
DCM_I_B = DCM_B_I';
DCM_G_I = [ 1-2*(q2G^2+q3G^2), 2*(q1G*q2G+q3G*q0G), 2*(q1G*q3G-q2G*q0G);
        2*(q2G*q1G-q3G*q0G), 1-2*(q1G^2+q3G^2), 2*(q2G*q3G+q1G*q0G);
        2*(q3G*q1G+q2G*q0G), 2*(q3G*q2G-q1G*q0G), 1-2*(q1G^2+q2G^2)];
MC=[radius_BCM 0 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0;
    radius_BCM 0 0 radius_G*cos(300*pi/180) radius_G*sin(300*pi/180) 0;
    radius_BCM*cos(120*pi/180) radius_BCM*sin(120*pi/180) 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0; 
    radius_BCM*cos(120*pi/180) radius_BCM*sin(120*pi/180) 0 radius_G*cos(180*pi/180) radius_G*sin(180*pi/180) 0;
    radius_BCM*cos(240*pi/180) radius_BCM*sin(240*pi/180) 0 radius_G*cos(180*pi/180) radius_G*sin(180*pi/180) 0;
    radius_BCM*cos(240*pi/180) radius_BCM*sin(240*pi/180) 0 radius_G*cos(60*pi/180) radius_G*sin(60*pi/180) 0 ];

%MODEL FUNCTION
tau=(rB-rG)/norm(rB-rG);

F_el = k * (norm(rB-rG)-l_t) * tau; 
F_d=c * norm(vB-vG) .*tau; 
F_t = F_el + F_d;
F_a_b = 0.5 * rhoB * C_D_B * S_BCM * (vB - WG).^2;
F_a_G = 0.5 * rhoG * C_D_G * S_G * (vB - WG).^2;
Fa=-rho_gas*VolB*g;
F_ext_BCM = Fa - 1* F_t + F_a_b + m_BCM*g;
F_ext_G = 1*F_t + F_a_G + m_G*g;

T_G=cross(MC(1,1:3)',DCM_B_I*F_t)+cross(MC(2,1:3)',DCM_B_I*F_t)+cross(MC(3,1:3)',DCM_B_I*F_t)+cross(MC(4,1:3)',DCM_B_I*F_t)+cross(MC(5,1:3)',DCM_B_I*F_t)+cross(MC(6,1:3)',DCM_B_I*F_t);
T_B=cross(MC(1,4:6)',-DCM_B_I*F_t)+cross(MC(2,4:6)',-DCM_B_I*F_t)+cross(MC(3,4:6)',-DCM_B_I*F_t)+cross(MC(4,4:6)',-DCM_B_I*F_t)+cross(MC(5,4:6)',-DCM_B_I*F_t)+cross(MC(6,4:6)',-DCM_B_I*F_t);


f_sym=[ rG + dt * vG;
    rB + dt * vB;
    vG + F_ext_G * dt/m_G;
    vB + F_ext_BCM * dt/m_BCM;
    sigma  * W_G1;
    lambda * W_G2;
    lambda * W_G3;
    arG + dt * aaG;
    arB + dt * aaB;
    I_G \ (T_G - cross(arG, I_G*arG));
    I_B \ (T_B - cross(arB, I_B*arB));
    (Id + dt/2 * OmegaG) * qG;
    (Id + dt/2 * OmegaB) * qB];

J_f=jacobian(f_sym,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_f_num= matlabFunction(J_f, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

% m_added_balloon = 0.5*rho_gas*VolB;
% m_b_final = m_BCM + m_added_balloon;




%PROCESS NOISE COVARIANCE MATRIX

q=[0.001,0.001,0.001,0.001,0.001,0.001,0.01,0.01,0.01,0.05,0.05,0.05,0.5,0.5,0.5,0.001,0.001,0.001,0.001,0.001,0.001,0.01,0.01,0.01,0.01,0.01,0.01,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001];
Q=diag(q);

%OBSERVATION FUNCTION

h_sym=[vG/dt;
    vB/dt;
    arG;
    arB;
    vG(1:2)-WG(1:2)];

%Not all measures are available at each time istant, 7 possibilities:

%1)all observations
h_sym1=h_sym;
J_h1=jacobian(h_sym1,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h1_num= matlabFunction(J_h1, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%2)only gondola IMU observations 
h_sym2=[vG/dt;
arG];

J_h2=jacobian(h_sym2,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h2_num= matlabFunction(J_h2, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%3)only bcm IMU observations
h_sym3=[vB/dt;
arB];

J_h3=jacobian(h_sym3,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h3_num= matlabFunction(J_h3, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%4)only wind observations
h_sym4=vG(1:2)-WG(1:2);
J_h4=jacobian(h_sym4,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h4_num= matlabFunction(J_h4, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%5)only wind and gondola IMU
h_sym5=[vG/dt;
arG;
vG(1:2)-WG(1:2)];
J_h5=jacobian(h_sym5,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h5_num= matlabFunction(J_h5, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%6)only wind and bcm IMU
h_sym6=[ vB/dt;
arB;
vG(1:2)-WG(1:2)];
J_h6=jacobian(h_sym6,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h6_num= matlabFunction(J_h6, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});

%7)only gondola and bcm IMU
h_sym7=[vG/dt;
vB/dt;
arG;
arB];
J_h7=jacobian(h_sym7,[rG;rB;vG;vB;WG;arG;arB;aaG;aaB;qG;qB]);
J_h7_num= matlabFunction(J_h7, 'Vars', {[r_G1; r_G2; r_G3; r_B1; r_B2; r_B3; v_G1; v_G2; v_G3; v_B1; v_B2; v_B3; W_G1; W_G2; W_G3; ar_G1; ar_G2; ar_G3; ar_B1; ar_B2; ar_B3; aa_G1; aa_G2; aa_G3; aa_B1; aa_B2; aa_B3; q0G; q1G; q2G; q3G; q0B; q1B; q2B; q3B], dt, rhoG, rhoB});



%MEASUREMENT NOISE COVARIANCE MATRIX

angular_velocity_noise_density = 5e-3;  % 5 mdps/√Hz
acceleration_noise_density = 60e-6;  % 60 μg/√Hz

% Sampling rate
angular_velocity_sampling_frequency = 8e3;  % 8 KHz
acceleration_sampling_frequency = 4e3;  % +4 KHz

% Variance
ar_variance = (angular_velocity_noise_density^2) * (1/angular_velocity_sampling_frequency);
ac_variance = (acceleration_noise_density^2) * (1/acceleration_sampling_frequency);

err_wind=[0.0001,0.0001];




% F_elG(:,1)=ones(3,1);
% F_elB(:,1)=ones(3,1);
% F_dampG=ones(3,1);
% F_dampB=ones(3,1);


%NUMERIC FUNCTION 

%Initial state
s(:,1)=[0;0;0;0.0000000001;0.0000000001;sqrt(0.066^2+2.55^2);0;0;0;0.0000001;0.0000001;0.0000001;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
P=diag(ones(1,35));
dt_int=0.001;
t_int=[0:dt_int:t(10000)]';
t_tot=sort([t;t_int]);
%% 

for i=2:10000 %l0ength(t)
    
    
    t_val=t_tot(i);
    dt_val=t_tot(i)-t_tot(i-1);

    % 1) Prediction state
   
    [TT, a, pp, rhoG_val] = atmosisa(double(s(3,i-1)));
    [TT, a, pp, rhoB_val] = atmosisa(double(s(6,i-1)));

   [s_p,F(:,i)]=model_function(s(:,i-1),dt_val,rhoG_val,rhoB_val,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,I_G,I_B,Id,MC,lambda,sigma);
    
   %error covariance
    P_p=J_f_num(s(:,i-1),dt_val,rhoG_val,rhoB_val)*P*J_f_num(s(:,i-1),dt_val,rhoG_val,rhoB_val)'+Q;
    uuuuu=J_f_num(s(:,i-1),dt_val,rhoG_val,rhoB_val);
    
    if any(ismember(t_val,t))==1
    % 2) Estimate state
    [h_xp,y,J_h,R]=observation_function(s_p,t_val,t,dt_val,t_G,t_B,t_wind,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance);

    % Kalman gain
    K_k=P_p*J_h(s_p,dt_val,rhoG_val,rhoB_val)'*(J_h(s_p,dt_val,rhoG_val,rhoB_val)*P_p*J_h(s_p,dt_val,rhoG_val,rhoB_val)'+R)^(-1);
    
    s(:,i)=s_p+K_k*(y-h_xp);
    
    % 4) Error covariance
    
    P=P_p-K_k*J_h(s_p,dt_val,rhoG_val,rhoB_val)*P_p';
    else

    s(:,i)= s_p;  
    P=P_p;
    end
i;
end

figure
plot(t_tot(1:i),s(1:3,:));
legend('x_gondola','y_gondola','z_gondola')
figure
plot(t_tot(1:i),s(4:6,:));
legend('x_bcm','y_bcm','z_bcm')

