function [x_p,F]=model_function(x,dt,rhoG_val,rhoB_val,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,I_G,I_B,Id,MC,lambda,sigma)
DCM_G_I = [ 1-2*(x(30)^2+x(31)^2), 2*(x(29)*x(30)+x(31)*x(28)), 2*(x(29)*x(31)-x(30)*x(28));
        2*(x(30)*x(29)-x(31)*x(28)), 1-2*(x(29)^2+x(31)^2), 2*(x(30)*x(31)+x(29)*x(28));
        2*(x(31)*x(29)+x(30)*x(28)), 2*(x(31)*x(30)-x(29)*x(28)), 1-2*(x(29)^2+x(30)^2)];
DCM_B_I = [ 1-2*(x(34)^2+x(35)^2), 2*(x(33)*x(34)+x(35)*x(32)), 2*(x(33)*x(35)-x(34)*x(32));
        2*(x(34)*x(33)-x(35)*x(32)), 1-2*(x(33)^2+x(35)^2), 2*(x(34)*x(35)+x(33)*x(32));
        2*(x(35)*x(33)+x(34)*x(32)), 2*(x(35)*x(34)-x(33)*x(32)), 1-2*(x(33)^2+x(34)^2)];
DCM_I_B = DCM_B_I';
DCM_I_G=DCM_G_I'; 

tau=(x(4:6)-x(1:3))/norm(x(4:6)-x(1:3));
if (norm(x(4:6)-x(1:3))-l_t <= 0)   
                F_el = [0;0;0];
                F_d = [0;0;0];
else   
for j=1:6
F_el = k * (DCM_I_B*MC(j,1:3)'+x(4:6)-DCM_I_G*MC(j,4:6)'+x(1:3))* ((DCM_I_B*MC(j,1:3)'+x(4:6)-DCM_I_G*MC(j,4:6)'+x(1:3))/norm((DCM_I_B*MC(j,1:3)'+x(4:6)-DCM_I_G*MC(j,4:6)'+x(1:3))))*norm((DCM_I_B*MC(j,1:3)'+x(4:6)-DCM_I_G*MC(j,4:6)'+x(1:3)));
end
% F_el = k * (norm(x(4:6)-x(1:3))-l_t) * tau; 
F_d=c * norm(x(10:12)-x(7:9)) .* tau; 

end



F_t = F_el + F_d;
Fae_b = 0.5 * rhoB_val * C_D_B * S_BCM * (x(10:12) - x(13:15)).^2;
Fae_G = 0.5 * rhoG_val * C_D_G * S_G * (x(7:9) - x(13:15)).^2;

Fa=-rho_gas*VolB*g;
F_ext_BCM = Fa - 1* F_t + Fae_b + m_BCM*g;
F_ext_G = 1*F_t + Fae_G + m_G*g;

F=[F_el;F_d;Fae_G;Fae_b];

OmegaG=[0 -x(16) -x(17) -x(18); x(16) 0 x(18) -x(17); x(17) -x(18) 0 x(16); x(18) x(17) -x(16) 0];
OmegaB=[0 -x(19) -x(20) -x(21); x(19) 0 x(21) -x(20); x(20) -x(21) 0 x(19); x(21) x(20) -x(19) 0];

T_G=cross(MC(1,1:3)',DCM_B_I*F_t)+cross(MC(2,1:3)',DCM_B_I*F_t)+cross(MC(3,1:3)',DCM_B_I*F_t)+cross(MC(4,1:3)',DCM_B_I*F_t)+cross(MC(5,1:3)',DCM_B_I*F_t)+cross(MC(6,1:3)',DCM_B_I*F_t);
T_B=cross(MC(1,4:6)',-DCM_B_I*F_t)+cross(MC(2,4:6)',-DCM_B_I*F_t)+cross(MC(3,4:6)',-DCM_B_I*F_t)+cross(MC(4,4:6)',-DCM_B_I*F_t)+cross(MC(5,4:6)',-DCM_B_I*F_t)+cross(MC(6,4:6)',-DCM_B_I*F_t);

x_p=[x(1:3) + dt * x(7:9);
        x(4:6) + dt * x(10:12);
        x(7:9) + F_ext_G * dt/m_G;
        x(10:12) + F_ext_BCM * dt/m_BCM;
        sigma * x(13:14);
        lambda * x(15);
        x(16:18)+ dt * x(22:24);
        x(19:21)+ dt * x(25:27);
        I_G \ (T_G - cross(x(16:18), I_G*x(16:18)));
        I_B \ (T_B - cross(x(19:21), I_B*x(19:21)));
        (Id + dt/2 * OmegaG) * x(28:31);
        (Id + dt/2 * OmegaB) * x(32:35)];
x(10:12) - x(13:15);
end