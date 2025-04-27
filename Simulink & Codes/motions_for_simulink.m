clc
close all
%% parametes 
a = -20;
b = 70;
c = -84;
d = 35;
e = 0;
f = 0;
%% initial and final joint anglesn ( You could change this )
theta_i = [0 ; 0 ; 0 ; 0 ] ; 
theta_f = [0.5 ; pi/3 ; -pi/6 ; -0.05 ];
%%
tf = 1;
Ts_M = 0.00001;
Ts = 0.001;
T = 0 : Ts_M : tf;
t_n_U = 0.8;
t_n_D = 0.2;
%%
p=0;

% 4567 Motion
for j = 0 : Ts_M : tf
    p = p + 1;
    t_n = j/tf;
    
    s(p) = a * t_n^ 7 + b * t_n ^ 6 + c * t_n ^ 5 + d * t_n ^ 4;
    s_prime(p) = 7*a*t_n^6 + 6*b*t_n^5 + 5*c*t_n^4 + 4*d*t_n^3;
    s_second(p) = 42*a*t_n^5 + 30*b*t_n^4 + 20*c*t_n^3 + 12*d*t_n^2;

% % %     
    theta = theta_i + (theta_f - theta_i) * s(p);
    theta = (theta - round( theta /2 / pi ) *2*pi);
    theta_dot = ((theta_f - theta_i) * s_prime(p) / tf);
    theta_ddot = ((theta_f - theta_i) * s_second(p) / tf^2);
%     
    
    b1_num(p) = theta(1);
    theta2_num(p) = theta(2);
    theta3_num(p) = theta(3);
    b4_num(p) = theta(4);
    
    b1_dot_num(p) = theta_dot(1);
    theta2_dot_num(p) = theta_dot(2);
    theta3_dot_num(p) = theta_dot(3);
    b4_dot_num(p) = theta_dot(4);
    
    b1_ddot_num(p) = theta_ddot(1);
    theta2_ddot_num(p) = theta_ddot(2);
    theta3_ddot_num(p) = theta_ddot(3);
    b4_ddot_num(p) = theta_ddot(4);
    
end
%%

b1_timeseries = timeseries(b1_num,T);
theta2_timeseries = timeseries(theta2_num,T);
theta3_timeseries = timeseries(theta3_num,T);
b4_timeseries = timeseries(b4_num,T);

b1_dot_timeseries = timeseries(b1_dot_num,T);
theta2_dot_timeseries = timeseries(theta2_dot_num,T);
theta3_dot_timeseries = timeseries(theta3_dot_num,T);
b4_dot_timeseries = timeseries(b4_dot_num,T);

b1_ddot_timeseries = timeseries(b1_ddot_num,T);
theta2_ddot_timeseries = timeseries(theta2_ddot_num,T);
theta3_ddot_timeseries = timeseries(theta3_ddot_num,T);
b4_ddot_timeseries = timeseries(b4_ddot_num,T);

%%
close all
markers = {'default','-.','--',':'};
figure()
hold on
plot(0 : Ts_M : tf , s ,'LineStyle',markers{1},'LineWidth',2)
plot(0 : Ts_M : tf , s_prime/max(s_prime) ,'LineStyle',markers{3},'LineWidth',2)
plot(0 : Ts_M : tf , s_second/max(s_second) ,'LineStyle',markers{4},'LineWidth',2)
axis( [0 1 -1.5 2])
xlabel( '$t_n$(s)','interpreter','latex')
legend('$s$','$s''/s''_{max}$','$s''''/s''''_{max}$','Interpreter','latex')
