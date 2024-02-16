function [x_p]=model_function(t,x,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,lambda,sigma)

[TT, a, pp, rhoG_val] = atmosisa(double(x(3)));
[TT, a, pp, rhoB_val] = atmosisa(double(x(6)));



% tau=(x(4:6)-x(1:3))/norm(x(4:6)-x(1:3));
% 
% if (norm(x(4:6)-x(1:3))-l_t <= 0)   
%                 F_el = [0;0;0];
%                 F_d = [0;0;0];
% else  
% 
% 
% F_el = k * (norm(x(4:6)-x(1:3))-l_t) * tau; 
% F_d=c * norm(x(10:12)-x(7:9)) .* tau; 

F_el = k * (x(4:6)-x(1:3)); 
F_d= c * x(10:12)-x(7:9); 
% F_d=[0;0;0];
% end



F_t = F_el + F_d;
Fae_b = 0.5 * rhoB_val * C_D_B * S_BCM * (x(10:12) - x(13:15)).^2;
Fae_G = 0.5 * rhoG_val * C_D_G * S_G * (x(7:9) - x(13:15)).^2;
%  Fae_b =[0;0;0];
%  Fae_G =[0;0;0];
Fa=-rho_gas*VolB*g;
F_ext_BCM = Fa - 1* F_t + Fae_b + m_BCM*g;
F_ext_G = 1*F_t + Fae_G + m_G*g;



x_p=[ x(7:9);
      x(10:12);
      F_ext_G /m_G;
      F_ext_BCM/m_BCM;
      sigma;
      sigma;
      lambda];



end