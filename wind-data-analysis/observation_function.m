function [h_xp,y,J_h,R]=observation_function(x,t_val,t,dt,t_G,t_B,t_wind,J_h1_num,J_h2_num,J_h3_num,J_h4_num,J_h5_num,J_h6_num,J_h7_num,ac_mG,ac_mB,ar_mG,ar_mB,windG,ac_variance,err_wind,ar_variance)



h_num=[x(7:9)/dt;
x(10:12)/dt;
x(16:18);
x(19:21);
x(7:8)-x(13:14)];


if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==1
       
        h_xp=[x(7:9)/dt;
        x(10:12)/dt;
        x(16:18);
        x(19:21);
        x(7:8)-x(13:14)];
        J_h=J_h1_num;
        j_B=find(ismember(t_val,t_B));
        j_G=find(ismember(t_val,t_G));
        j_W=find(ismember(t_val,t_wind));
        y=[ac_mG(:,j_G);ac_mB(:,j_B);ar_mG(:,j_G);ar_mB(:,j_B);windG(:,j_W)];
        r=[ac_variance*ones(1,6),ar_variance*ones(1,6),err_wind];
        R=diag(r);
end

    if any(ismember(t_val,t_G))==1 
        h_xp=[x(7:9)/dt;
        x(16:18)];
        J_h=J_h2_num;
       
        j_G=find(ismember(t_val,t_G));
      
        y=[ac_mG(:,j_G);ar_mG(:,j_G)];
        r=[ac_variance*ones(1,3),ar_variance*ones(1,3)];
    R=diag(r);
    end
    if any(ismember(t_val,t_B))==1
        h_xp=[x(10:12)/dt;
        x(19:21)];
        J_h=J_h3_num;
        j_B=find(ismember(t_val,t_B));
        y=[ac_mB(:,j_B);ar_mB(:,j_B)];
        r=[ac_variance*ones(1,3),ar_variance*ones(1,3)];
        R=diag(r);
    end
    if any(ismember(t_val,t_wind))==1
        h_xp=[x(7:8)-x(13:14)];
        J_h=J_h4_num;
        j_W=find(ismember(t_val,t_wind));
        y=windG(:,j_W);
        r=[err_wind];
        R=diag(r);
    end
    if  any(ismember(t_val,t_G))==1 && any(ismember(t_val,t_wind))==1
        h_xp=[x(7:9)/dt;
        x(16:18);
        x(7:8)-x(13:14)];
        J_h=J_h5_num;
        j_G=find(ismember(t_val,t_G));
        j_W=find(ismember(t_val,t_wind));
        y=[ac_mG(:,j_G);ar_mG(:,j_G);windG(:,j_W)];
        r=[ac_variance*ones(1,3),ar_variance*ones(1,3),err_wind];
        R=diag(r);
    end
    if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_wind))==1
        h_xp=[x(10:12)/dt;
        x(19:21);
        x(7:8)-x(13:14)];
        J_h=J_h6_num;
        j_B=find(ismember(t_val,t_B));
        j_W=find(ismember(t_val,t_wind));
        y=[ac_mB(:,j_B);ar_mB(:,j_B);windG(:,j_W)];
        r=[ac_variance*ones(1,3),ar_variance*ones(1,3),err_wind];
        R=diag(r);
    end
    if any(ismember(t_val,t_B))==1 && any(ismember(t_val,t_G))==1
        h_xp=[x(7:9)/dt;
        x(10:12)/dt;
        x(16:18);
        x(19:21)];
        J_h=J_h7_num;
        j_B=find(ismember(t_val,t_B));
        j_G=find(ismember(t_val,t_G));
        y=[ac_mG(:,j_G);ac_mB(:,j_B);ar_mG(:,j_G);ar_mB(:,j_B)];
        r=[ac_variance*ones(1,6),ar_variance*ones(1,6)];
        R=diag(r);
    end




end