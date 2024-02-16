function [dsdt,wind_value] = odeMSD_2D(t,s,atm_data,z_eq,n,wind_speed_int)

  % States

  % s(1) = x balloon
  % s(2) = v_x balloon
  % s(3) = z balloon
  % s(4) = v_z balloon

  % s(5:4:4*(n-1)) = x tether n
  % s(6:4:4*(n-1)) = v_x tether n
  % s(7:4:4*(n-1)) = z tether n
  % s(8:4:4*(n-1)) = v_z tether n

  % s(end-4) = x gondola
  % s(end-3) = v_x gondola
  % s(end-2) = z gondola
  % s(end) = v_z gondola

  % Venus charateristics
  g = 8.87;             % [m/s^2] 
  V_v = 0*1.81;         % Venus equatorial velocity [m/s]
  r_v = 6051.8e3;       % Venus radius [m]
  omega_v = V_v/r_v;    % Venus angular velocity [rad/s]
  
  % Kevlar tether
  AS = 47312;           % Axial stiffness [N]
  l_t = 10;             % Tether length [m]
  l_s = l_t/n;          % Portion tether length [m]
  m_t = 0.1*l_t;        % Tether mass [kg]
  m_s = m_t/(n-1);      % Portion tether mass [kg]
  d_t = 1e-2;           % Tether diameter [m]
  A_t = pi*(d_t/2)^2;   % Tether area [m^2]
  E = 83e9; %83e9       % Young modulus [Pa]
  k = A_t*E/l_t;        % Spring stiffness [N/m]
  c = 47.15/l_s;               % Damping [Ns]
  CD_t = 1.15;          % Tether drag coefficient
  S_t = (l_s * pi*d_t)/2;      % Exposed tether surface [m^2]
    
  % System charateristics (Balloon + Gondola)
  m_He = 22.3;                  % [kg] Helium total mass, distributed between the the balloons 
  m_env_ZP = 79.2;              % [kg] Zero Pressure envelope mass 
  m_env_SP = 28;                % [kg] Super Pressure envelope mass
  m_gon = 75;                   % [kg] Gondola mass
  M_tot = (m_env_ZP + m_env_SP) + m_He + m_gon + m_t;   %[kg] System total mass
  m_b = M_tot - m_gon - m_t;    % [kg] Total mass wo Gondola
  [T_eq,p_eq,rho_eq] = density(z_eq,atm_data);
  V_b = M_tot/rho_eq;           % [m^3] Volume of the equivalent Balloon
  d_b = (6*(V_b)/(pi))^(1/3);   % [m] Diameter of the equivalent Balloon
  S_b = pi*(d_b^2/4);           % [m^2] Balloon exposed surface 
  CD_b = 0.5;
  
  l_G = 0.5;    %[m] gondola side
  S_G = l_G^2;  %[m^2] gondola area exposed to the air flux
  CD_G = 1.05;
  
  % Wind model (Impulse along the x axis and z axis at different times)
 
 w_x=wind_speed_int(t);
 
 wind_value=w_x;

 

%       w_x = -3; %Positive along x axis [m/s]
  
  w_z = 0;
  if ((t>300)&&(t<325))
     
      w_z = 0; %Positive along z axis [m/s]
  end
 
  % Forces
  x_rel = s(1)-s(5);
  z_rel = s(3)-s(7);
  theta = atan(x_rel/z_rel);
  x_dot_rel = s(2)-s(6);
  z_dot_rel = s(4)-s(8);

  % Atmospheric data 
  [T_b,p_b,rho_b] = density(s(3),atm_data);
  T_t = zeros(n-1,3);
  for i = 1 : n-1
     [T_t(i) p_t(i) rho_t(i)] = density(s(3+i*4),atm_data);
  end
  [T_G,p_G,rho_G] = density(s(end-1),atm_data);
  

  
  %Aerostatic force
  Fa = rho_b*V_b*g;
    
  %Tether forces
  for i = 0:n-1   
      F_el(i+1) = k*(sqrt((s(4*i+3)-s(4*i+7)).^2+(s(4*i+1)-s(4*i+5)).^2)-l_s);   % [N] Elastic force
      F_d(i+1) = c*(sqrt((s(4*i+4)-s(4*i+8)).^2+(s(4*i+2)-s(4*i+6)).^2));        % [N] Damping force
      theta(i+1) = atan((s(4*i+5)-s(4*i+1))./abs(s(4*i+3)-s(4*i+7)));            % [rad] Angular displacement btw b e G
      F_t(i+1) = F_el(i+1) + F_d(i+1);
      endq
  
  %Aerodynamic forces

  Dx_b = 0.5*rho_b*CD_b*S_b*(s(2)+omega_v*(s(1)+r_v)-w_x).^2*sign(s(2)+omega_v*(s(1)+r_v)-w_x); 
  Dz_b = 0.5*rho_b*CD_b*S_b*(s(4)-w_z).^2*sign(s(4)-w_z);

  % Dx_b = 0.5*rho_b*CD_b*S_b*s(2).^2*sign(s(2));
  % Dz_b = 0.5*rho_b*CD_b*S_b*s(4).^2*sign(s(4));

  for i = 0:n-2
      Dx_t(i+1) = 0.5*rho_t(i+1)*CD_t*S_t*(s(4*i+6)+omega_v*(s(4*i+7)+r_v)-w_x).^2*sign(s(4*i+6)+omega_v*(s(4*i+7)+r_v)-w_x);   % [N] Drag on the tether (Rotation of the planet involved)
      Dz_t(i+1) = 0.5*rho_t(i+1)*CD_t*S_t*(s(4*i+6)-w_z).^2*sign(s(4*i+6)-w_z);   % [N] Drag on the tether along z (no rotation of the planet involved)
  end 
  
  Dx_G = 0.5*rho_G*CD_G*S_b*(s(end-2)+omega_v*(s(end-3)+r_v)-w_x).^2*sign(s(end-2)+omega_v*(s(end-3)+r_v)-w_x);
  Dz_G = 0.5*rho_G*CD_G*S_b*(s(end)-w_z).^2*sign(s(end)-w_z);
  
  % %Wind forces
  % Fw_x_b = 0.5*rho_b*CD_b*S_b*w_x^2*sign(w_x);
  % Fw_z_b = 0.5*rho_b*CD_b*S_b*w_z^2*sign(w_z);
  % 
  % % HERE ADD WIND FORCES ACTING ON THE TETHER MASSES
  % 
  % Fw_x_G = 0.5*rho_G*CD_G*S_G*w_x^2*sign(w_x);
  % Fw_z_G = 0.5*rho_G*CD_G*S_G*w_z^2*sign(w_z);


  % Ode function definition

  dsdt = [s(2);
          (F_t(1)*sin(theta(1)) - Dx_b) / m_b;
          s(4);
          (Fa - m_b*g - F_t(1)*cos(theta(1)) - Dz_b) / m_b];
            
          j = 1;
          for i = 5:4:5+(n-2)*4 %6:n-2+6
                dsdt = [dsdt;
                        s(i+1); 
                        (-(F_t(j))*sin(theta(j)) + (F_t(j+1))*sin(theta(j+1)))/m_s;
                        s(i+3); 
                        ((F_t(j))*cos(theta(j)) - (F_t(j+1))*cos(theta(j+1))-m_s*g)/m_s];
       
                j = j+1;
          end

  dsdt = [dsdt; 
          s(end-2);
          (-(F_t(end))*sin(theta(end))- Dx_G) / m_gon;
          s(end);
          (F_t(end)*cos(theta(end)) - m_gon*g - Dz_G) / m_gon]; 
   
end

