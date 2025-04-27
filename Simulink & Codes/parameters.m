syms b1 theta2 theta3 b4
syms b1_dot theta2_dot theta3_dot b4_dot b1_ddot theta2_ddot theta3_ddot b4_ddot

%% Physical Data
% Link0
theta = [ b1 ; theta2 ; theta3 ; b4 ];
theta_dot = [ b1_dot ; theta2_dot; theta3_dot ; b4_dot];
theta_ddot = [ b1_ddot ; theta2_ddot; theta3_ddot ; b4_ddot];

theta1=0;
theta4=0;
theta1_0=0;
theta4_0=0;
b1_0=0.08;
theta2_0 = 0;
theta3_0 = 0;
b4_0 = 0;

% com, mass, I
com0=[0;0;0];
com1=[0.15;0.08;0.1];
com2=[0.35;0.1;0.25];
com3=[0.45;0.1;0.3];
com4=[0.65;0.02;0.4];
mass_link0=0.8;
mass_link1=0.8;
mass_link2=1.2;
mass_link3=0.85;
mass_link4=0.15;
I0=[0.44,0,0;0,0.003,0;0,0,0.44];
I1=[2000,2,160;2,3000,1;160,1,2000]*1e-06;
I2=[1600,-0.25,-1200;-0.25,30000,0.05;-1200,0.05,30000]*1e-06;
I3=[6000,-0.3,1300;-0.3,15000,-0.3;1300,-0.3,12000]*1e-06;
I4=[3000,0.07,0.04;0.07,3000,0.1;0.04,0.1,20]*1e-06;
%% D-H Form
% T = [Q,P;0,1];
% Com_dh = T*Com_link
% I_dh = Q*I*Q'

% Write your code here
T0=[round(Qx(-pi/2)),[0;0;73.5]/1000;[0,0,0,1]]*[eye(3),[0;0;0.08];[0,0,0,1]];
% link 1
T=T0*[round(Qx(pi/2)),[-0.2;0;0];[0,0,0,1]];
com1_dh=T*[com1;1];
com1_dh=com1_dh(1:3);
Q=T(1:3,1:3);
b=T(1:3,4);
I1_dh=Q*(I1 + mass_link1*((b'*b)*eye(3) - b*b'))*inv(Q);
% link 2
T=T*[eye(3),[0.2;0.08;0.078];[0,0,0,1]]*[eye(3),[-0.073;-0.083;-0.079];[0,0,0,1]];
com2_dh=T*[com2;1];
com2_dh=com2_dh(1:3);
Q=T(1:3,1:3);
I2_dh=Q*(I2 + mass_link2*((b'*b)*eye(3) - b*b'))*inv(Q);
% link 3
T=T*[eye(3),[0.47;0.09;0.2];[0,0,0,1]]*[eye(3),[-0.4;-0.09;-0.2];[0,0,0,1]];
com3_dh=T*[com3;1];
com3_dh=com3_dh(1:3);
b=T(1:3,4);
Q=T(1:3,1:3);
I3_dh=Q*(I3 + mass_link3*((b'*b)*eye(3) - b*b'))*inv(Q);
% link 4
T=T*[eye(3),[0.65;0.07;0.28];[0,0,0,1]]*[eye(3),[-0.65;0;-0.28];[0,0,0,1]];
com4_dh=T*[com4;1];
com4_dh=com4_dh(1:3);
b=T(1:3,4);
Q=T(1:3,1:3);
I4_dh=Q*(I4 + mass_link4*((b'*b)*eye(3) - b*b'))*inv(Q);

%% First dot ( Write W and N vectors)

% Write your code here
Q0=Qx(-pi/2);
Q1=[1,0,0;0,0,-1;0,1,0];
Q2=[cos(theta2),-sin(theta2),0;sin(theta2),cos(theta2),0;0,0,1];
Q3=[cos(theta3),-sin(theta3),0;sin(theta3),cos(theta3),0;0,0,1];
Q4=[1,0,0;0,1,0;0,0,1];

zero=[0;0;0];
z=[0;0;1];
e1=Q0*z;
e2=[0;0;1];
e3=e2;
e4=e2;

av1=[0;0;b1];
av2=[0.3971*cos(theta2);0.3971*sin(theta2);0.121];
av3=[0.1254*cos(theta3);0.1254*sin(theta3);0.08];
av4=[0;0;b4];

r11 = Q0*(com1_dh+[0;0;b1] );

r12 =   Q0*(av1 + Q1*Qz( theta2 -theta2_0 )*(com2_dh) );
r22 =   Q0*Q1*Qz( theta2 - theta2_0)*(com2_dh) ;

r13 =   Q0*(av1 + Q1*av2 + Q1*Q2*Qz( theta3 - theta3_0)*(com3_dh) ) ;
r23 =   Q0*(Q1*av2 + Q1*Q2*Qz( theta3 - theta3_0)*(com3_dh)) ;
r33 =   Q0*(Q1*Q2*Qz( theta3 - theta3_0)*(com3_dh)) ;

r14 =   Q0*(av1 + Q1*av2 + Q1*Q2*av3 + Q1*Q2*Q3*(com4_dh+[0;0;b4])) ;
r24 =   Q0*(Q1*av2 + Q1*Q2*av3 + Q1*Q2*Q3*(com4_dh+[0;0;b4])) ;
r34 =   Q0*(Q1*Q2*av3 + Q1*Q2*Q3*(com4_dh+[0;0;b4]) );
r44 =   Q0*(Q1*Q2*Q3*(com4_dh+[0;0;b4]) ) ;

% Ni in world frame
N1 = [e1 zero zero zero   ];
N2 = [e1  cross(e2,r22)  zero            zero          ];
N3 = [e1  cross(e2,r23)  cross(e3,r33)   zero          ];
N4 = [e1  cross(e2,r24)  cross(e3,r34)   e4 ];
% Wi in i th frame <=== Ii in i th frame

W1 = [zero     zero               zero     zero];
W2 = [zero     z                     zero     zero];
W3 = [zero     Q2'*z             z           zero];
W4 = [zero     Q3'*Q2'*z     Q3'*z   zero];

%% Final calculation ( Energy eq. - Use the paper )

% Write your code here
M1 = mass_link1 * N1'*N1 + W1'*I1_dh*W1;
M2 = mass_link2 * N2'*N2 + W2'*Qz(theta2 - theta2_0)*I2_dh*Qz(theta2 - theta2_0)'*W2;
M3 = mass_link3 * N3'*N3 + W3'*Qz(theta3 - theta3_0)*I3_dh*Qz(theta3 - theta3_0)'*W3;
M4 = mass_link4 * N4'*N4 + W4'*I4_dh*W4;
M = M1 + M2 + M3 + M4  ;
T = 0.5*theta_dot'*M*theta_dot;
%hi all in world frame
h1 =  (Q0*(com1_dh+[0;0;b1]) )'*z;
h2 =   (Q0*(av1 + Q1*Qz(theta2 - theta2_0)*com2_dh ))'*z;
h3 =   (Q0*(av1 + Q1*av2 + Q1*Q2*Qz(theta3 - theta3_0)*com3_dh ))'*z;
h4 =   (Q0*(av1 + Q1*av2 + Q1*Q2*av3 + Q1*Q2*Q3*(com4_dh+[0;0;b4])))'*z;

g = 9.81;
v1 = mass_link1 *g*h1;
v2 = mass_link2 *g*h2;
v3 = mass_link3 *g*h3;
v4 = mass_link4 *g*h4;
V = v1 + v2 + v3 + v4 ;
M_dot = M_dot_generator(M,theta,theta_dot);

%% Euler-Lagrange equations (jacobian etc. )

% Write your code here
n = M_dot*theta_dot ...
    -jacobian(T , theta)' + jacobian(V,theta)';

tav = (M*theta_ddot + n);

tf = 1;
Ts_M = 0.01;
Ts = 0.001;
T = 0 : Ts_M : tf;
p=0;
%% plynomianl 4 5 6 7  
a = -20;
b = 70;
c = -84;
d = 35;
e = 0;
f = 0;
theta_i = [0.08 ; 0 ; 0 ; 0 ] ;
theta_f = [0.5 ; pi/3 ; -pi/6 ; -0.05 ];
%%
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
    
    TAU(:,p)=eval(tav);
end
%%
figure(1)
% plot(out.f1.Time,out.f1.Data);
plot(0:Ts_M:tf,TAU(1,:));
xlabel('t (s)')
ylabel('f1 (n.m)')

figure(2)
% plot(out.tau2.Time,out.tau2.Data);
plot(0:Ts_M:tf,TAU(2,:));
xlabel('t (s)')
ylabel('\tau2 (n.m)')

figure(3)
% plot(out.tau3.Time,out.tau3.Data);
plot(0:Ts_M:tf,TAU(3,:));
xlabel('t (s)')
ylabel('\tau3 (n.m)')

figure(4)
% plot(out.f4.Time,out.f4.Data);
plot(0:Ts_M:tf,TAU(4,:));
xlabel('t (s)')
ylabel('f4 (n.m)')

%% Matlab Function

% For example, this code, convert your ForceandTorques Matrix equations to "OpenMan_torquesandforce" function to use later
% This will generate a function that output your calculations based on yuour symbolic inputs
ForceandTorques = [Fb1;Tt2;Tt3;Fb4];
ForceandTorques = simplify(ForceandTorques);
matlabFunction(ForceandTorques,'File','OpenMan_torquesandforce');