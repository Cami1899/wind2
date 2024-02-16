function [s_0] = Initial_Conditions(z_0,z_eq,l_t1,l_t2,lG,n)
% This function allow to determine the Bifilar tethered system states
% Inputs:
% z_0 = Balloon starting altitude [m]
% z_eq = Balloon design altitude [m]
% l_t1 = Length of the Tether unifilar portion [m]
% l_t2 = Length of the Tether bifilar portion [m]
% lG = Cubic gondola side [m]
% n = number of unifilar segments 

l0 = lG/2;
l_s2 = sqrt(l_t2.^2+l0.^2);
offset = z_0 - z_eq;

l_t = l_t1 + l_t2;
l_s1 = (l_t1 - 2*l_s2)/n;     % We want the node of the bifilar tether to coincide with the last tether mass 

% Gondola attitude 
theta_0 = 0; %[rad]
phi_0 = 0; %[rad]
psi_0 = 0; %[rad]
DCM3 = EulAxAngl2DCM([0;0;1],psi_0);
DCM2 = EulAxAngl2DCM([0;-1;0],phi_0);
DCM1 = EulAxAngl2DCM([1;0;0],theta_0);
DCM_0 = DCM1*DCM2*DCM3;
% Gondola attitude quaternion
quat_0 = DCM_quat(DCM_0);

% State vector inizialization
s_0 = zeros(6*(n+2)+1,1); 
s_0(1:6,1) = [ 0 ; 0 ; z_0; 0 ; 0 ; 0];     % Balloon states
for i = 1 : n
    s_0(6+6*(i-1)+1,1) = 0; 
    s_0(6+6*(i-1)+2,1) = 0; 
    s_0(6+6*(i-1)+3,1) = z_0 - (1+2*(i-1))*l_s1/2;
    s_0(6+6*(i-1)+4,1) = 0;
    s_0(6+6*(i-1)+5,1) = 0;
    s_0(6+6*(i-1)+6,1) = 0;
end
% Bifilar tether states
s_0(6+6*n+1,1) = 0; 
s_0(6+6*n+2,1) = 0; 
s_0(6+6*n+3,1) = z_0 - l_t1;
s_0(6+6*n+4,1) = 0;
s_0(6+6*n+5,1) = 0;
s_0(6+6*n+6,1) = 0;


l0_0 = DCM_0*[0;0;l0];
s_0((6*(n+2)+1):(6*(n+2)+6),1) = [ l0_0(1) ; l0_0(2) ; z_0 - l_t - l0_0(3) ; 0.1 ; 0 ; 0];

s_0(6*(n+3)+1:6*(n+3)+4) = quat_0; % Initial base quaternion
s_0(6*(n+3)+5:6*(n+3)+7) = [0;0;0.1]; % Initial angular velocities
s_0(6*(n+3)+8) = psi_0; % Initial torsional spring displacement

end