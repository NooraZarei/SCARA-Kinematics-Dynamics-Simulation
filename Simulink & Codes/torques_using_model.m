%%
Ts_M = 0.01;
tf = 1;
bias = [ 0 ; atan(0.2) ; -atan(0.7); 0 ] ;
theta_i_m = theta_i + bias;
theta_f_m = theta_f;
p = 0;
Ts_M = 0.001;
TAU = zeros(4, round(1 + tf / Ts_M));

% 4567 Motion
for j = 0 : Ts_M : tf
    p = p + 1;
    t_n = j/tf;
    
    s = a * t_n^ 7 + b * t_n ^ 6 + c * t_n ^ 5 + d * t_n ^ 4;
    s_prime = 7*a*t_n^6 + 6*b*t_n^5 + 5*c*t_n^4 + 4*d*t_n^3;
    s_second = 42*a*t_n^5 + 30*b*t_n^4 + 20*c*t_n^3 + 12*d*t_n^2;
    
    theta = theta_i_m + (theta_f_m - theta_i_m) * s;
    theta = theta - round( theta /2 / pi ) *2*pi;

    theta_dot = (theta_f_m - theta_i_m) * s_prime / tf;
    theta_ddot = (theta_f_m - theta_i_m) * s_second / tf^2; 
    
    b1 = theta(1);
    theta2 = theta(2);
    theta3 = theta(3);
    b4 = theta(4);

    
    b1_dot = theta_dot(1);
    theta2_dot = theta_dot(2);
    theta3_dot = theta_dot(3);
    b4_dot = theta_dot(4);

    b1_ddot = theta_ddot(1);
    theta2_ddot = theta_ddot(2);
    theta3_ddot = theta_ddot(3);
    b4_ddot = theta_ddot(4);

    % For example, if your function name is "OpenMan_torquesandforce" :
    TAU(:,p) = OpenMan_torquesandforce(b1,b1_dot,b1_ddot,b4_ddot,theta2,theta3,theta2_dot,theta3_dot,theta2_ddot,theta3_ddot);
    % This will calculate the dynamic model that you wrote
    % And put it into TAU Matrix ( Use this matrix to plot the data )
end
