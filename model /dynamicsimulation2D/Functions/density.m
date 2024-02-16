function [T_h,p_h,rho_h] = density(z,atm_data)
    
  z = z/1000; %[km]
  flag = 0;
  
  for i = 1:size(atm_data,1)-1
      
      if (z>=atm_data(i,1) && z<atm_data(i+1,1))
          z_1 = atm_data(i,1); 
          z_2 = atm_data(i+1,1);
          T_1 = atm_data(i,2); 
          T_2 = atm_data(i+1,2);
          p_1 = atm_data(i,3); 
          p_2 = atm_data(i+1,3);
          rho_1 = atm_data(i,4); 
          rho_2 = atm_data(i+1,4);
          
          % Values at altitude z
          T_h = T_1 + ((T_2-T_1)/(z_2-z_1))*(z-z_1);
          p_h = p_1 + ((p_2-p_1)/(z_2-z_1))*(z-z_1);
          rho_h = rho_1 + ((rho_2-rho_1)/(z_2-z_1))*(z-z_1);
          
          flag = 1;
          break;
      end

  end
        
  if flag == 0
          
          % Values when altitude exceeds the maximum value in the table 
          T_h = atm_data(end,2);
          p_h = atm_data(end,3);
          rho_h = atm_data(end,4);
          
  end
          
%   Earth atmposphere

%   T_0 = 273.15 + 25; %[K] ISA
%   m = 5.2561;
%   h = 0.0065; %[K/m]  
%   rho_0 = 1.225;
%   rho_h = rho_0*((T_0 - h*z)/T_0)^(m-1);
end

