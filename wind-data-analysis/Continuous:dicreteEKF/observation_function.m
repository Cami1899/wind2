function [h,y,R,J_h]=observation_function(x,t_val,t_G,t_B,t_wind,ac_mG,ac_mB,windG,ac_variance,err_wind,C_D_B,S_BCM,C_D_G,S_G,rho_gas,VolB,g,k,c,l_t,m_G,m_BCM,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num)



[TT, a, pp, rhoG_val] = atmosisa(double(x(3)));
[TT, a, pp, rhoB_val] = atmosisa(double(x(6)));


if (norm(x(4:6)-x(1:3))-l_t <= 0)   
                F_el = [0;0;0];
                F_d = [0;0;0];
else  

F_el = k * (x(4:6)-x(1:3)); 
F_d= c * x(10:12)-x(7:9); 
end

F_t = (F_el + F_d);
Fae_b = -sign(x(10:12) - x(13:15)) * 0.5 * rhoB_val * C_D_B .* S_BCM .* (x(10:12) - x(13:15)).^2;
Fae_G = -sign(x(7:9) - x(13:15)) * 0.5 * rhoG_val * C_D_G .* S_G .* (x(7:9) - x(13:15)).^2;
Fa= -(rhoB_val-rho_gas)*VolB*g;
F_ext_BCM = Fa -  F_t + Fae_b + m_BCM*g;
F_ext_G = F_t + Fae_G + m_G*g;



if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==1
           
        h=[F_ext_G /m_G;
           F_ext_BCM/m_BCM;
           x(7:8) - x(13:14)];
        j_B=find(ismember(t_B,t_val));
        j_G=find(ismember(t_G,t_val));
        j_W=find(ismember(t_wind,t_val));
        y=[ac_mG(:,j_G);ac_mB(:,j_B);windG(:,j_W)];
        r=[ac_variance*ones(1,6),err_wind];
        R=diag(r);
        J_h=J_h1_num;
%         diffG=y(1:3)-h(1:3);
%         diffB=y(4:6)-h(4:6);
%         diffW=y(7:8)-h(7:8);
%         disp('ALL MEASURES')
end

    if any(ismember(t_val,t_B))==0 && any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==0
        h=F_ext_G /m_G;   
        j_G=find(ismember(t_G,t_val));   
        y=[ac_mG(:,j_G)];
        r=ac_variance*ones(1,3);
        R=diag(r);
        J_h=J_h2_num;
%         diffG=y(1:3)-h(1:3);
%         disp('ACC GONDOLA MEASURES')
    end
    if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==0 && any(ismember(t_val,t_wind))==0
         h=F_ext_BCM/m_BCM;
        j_B=find(ismember(t_B,t_val));
        y=[ac_mB(:,j_B)];
        r=[ac_variance*ones(1,3)];
        R=diag(r);
        J_h=J_h3_num;
%         diffB=y(4:6)-h(4:6);
%         disp('ACC BALLOON MEASURES')
    end
    if any(ismember(t_val,t_B))==0 && any(ismember(t_val,t_G))==0 && any(ismember(t_val,t_wind))==1
        h= x(7:8) - x(13:14);
        j_W=find(ismember(t_wind,t_val));
        y=windG(:,j_W);
        r=[err_wind];
        R=diag(r);
        J_h=J_h4_num;
%          diffW=y(7:8)-h(7:8);
         disp('WIND')
    end
    if  any(ismember(t_val,t_B))==0 && any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==1
        h=[F_ext_G /m_G;
           x(7:8) - x(13:14)];
        j_G=find(ismember(t_G,t_val));
        j_W=find(ismember(t_wind,t_val));
        y=[ac_mG(:,j_G);windG(:,j_W)];
        r=[ac_variance*ones(1,3),err_wind];
        R=diag(r);
        J_h=J_h5_num;
%         diffG=y(1:3)-h(1:3);
%         diffW=y(7:8)-h(7:8);
         disp('ACC GONDOLA MEASURES+WIND');
    end
    if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==0 && any(ismember(t_val,t_wind))==1
       h=[F_ext_BCM/m_BCM;
          x(7:8) - x(13:14)];
         j_B=find(ismember(t_B,t_val));
        j_W=find(ismember(t_wind,t_val));
        y=[ac_mB(:,j_B);windG(:,j_W)];
        r=[ac_variance*ones(1,3),err_wind];
        R=diag(r);
        J_h=J_h6_num;
%         diffB=y(4:6)-h(4:6);
%         diffW=y(7:8)-h(7:8);
         disp('ACC BALLOON MEASURES+WIND');
    end
    if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==0
        h=[F_ext_G /m_G;
           F_ext_BCM/m_BCM];
        j_B=find(ismember(t_B,t_val));
        j_G=find(ismember(t_G,t_val));
        y=[ac_mG(:,j_G);ac_mB(:,j_B)];
        r=[ac_variance*ones(1,6)];
        R=diag(r);
        J_h=J_h7_num;diffG=y(1:3)-h(1:3);
%         diffB=y(4:6)-h(4:6);
     

%          disp('ACC GONDOLA MEASURES+ACC BALLOON MEASURES');


    end




end