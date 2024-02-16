clc
clear
% Inizializzazione
x0 = [0; 0; 2]; % Stima iniziale dello stato [posizione; velocità]
P_0 = eye(3); % Covarianza iniziale della stima
dt = 0.5; % Intervallo di tempo tra le misure
t_end = 100; % Tempo finale
t = 0:dt:t_end;
R = 1; 
Q = diag([0.1, 0.1, 0.1]);
% Simulazione
x_measure = zeros(3, length(t));
x_measure(:, 1) = x0;
x_true = x_measure;

for i = 2:length(t)
    [~, x_next] = ode45(@(t, x) [x(2); x(3); 0] , [t(i-1), t(i)], x_true(:, i-1));
    x_next = x_next(end, :)';
    x_true(:, i) = x_next;
    x_measure(:, i) = x_next + mvnrnd(zeros(3,1), Q)';
end

% Simulazione delle misure

for i = 1:length(t)
    % Simula la misura come la posizione, la velocità e l'accelerazione vere con rumore gaussiano aggiunto
    y(i) = x_measure(3, i) +normrnd(0, R)  ;
end

J_f_num=[0 1 0;
         0 0 1;
         0 0 0];
C=[0 0 1];

P_old=P_0;
x_old=x0;
h=0.1;
dim=t_end/h;
l=0;
% state=zeros(3,dim+3);
for i =1:length(y)-1
    start_time=t(i);
   dt_between_measurements=t(i+1)-t(i);
   dim_int=dt_between_measurements/h;

[x_p,P_p,P_diag,t_p] = cdekf_predict_phase(J_f_num,dt_between_measurements,start_time,P_old,x_old,Q,h);


t_new(l+1:l+length(t_p))=t_p;

[x_new, P_new] = cdekf_update_phase(R,P_p,C,x_p(:,end),y(i+1));

x_p(:,end)=x_new;
state(:,l+1:l+length(t_p))=x_p;
P_diag(:,end)=diag(P_new);
P_diag_matrix(:,l+1:l+length(t_p))=P_diag;
x_old=x_new;
P_old=P_new;
l=length(t_new);
end
% state=state(:,dim_int+1:end);
% P_diag_matrix=P_diag_matrix(:,dim_int+1:end);
std=sqrt(P_diag_matrix);
% t_new=t_new(dim_int+1:end);





figure
plot(t_new,state(1,:))
hold on 
title('position')
plot(t_new,state(1,:)+2*std(1,:))
plot(t_new,state(1,:)-2*std(1,:))
plot(t,x_true(1,:))
hold off

figure
plot(t_new,state(2,:))
hold on 
title('velocity')
plot(t_new,state(2,:)+2*std(2,:))
plot(t_new,state(2,:)-2*std(2,:))
plot(t,x_true(2,:))
hold off

figure
plot(t_new,state(3,:))
hold on 
title('acceleration')
plot(t_new,state(3,:)+2*std(3,:))
plot(t_new,state(3,:)-2*std(3,:))
plot(t,x_true(3,:))
hold off


