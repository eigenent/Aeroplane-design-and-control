% Dimitrios Nentidis AEM: 6821

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Control %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% define full system

A = [-0.0352,  0.1070,       0, -32.2;
     -0.2140, -0.4400,     305,     0;
  1.1980e-04, -0.0154, -0.4498,     0;
           0,      0 ,       1,    0];
  
B = [      0; 
    -22.1206; 
     -4.6580;
           0]; 

C = eye(4,4);
D = [0;0;0;0];

sys_full =  ss(A,B,C,D);

%% to transfer functions
sys_tf_full  = tf(sys_full);

damp(sys_tf_full)

% do various things such as rlocus, step etc in command window

%% Control q

Gq = sys_tf_full(3);

% this seems to be the optimal Kq
Kq = -0.42;

% Create the closed-loop transfer function
Gq_cl= feedback(Gq, Kq);

damp(Gq_cl)

figure
step(Gq)
hold on
step(Gq_cl)
legend('open loop', 'closed loop')
ylabel('Amplitude q')
hold off
grid on

%% reduced state space

A_red = A(2:3, 2:3);
B_red = B(2:3, :);
C_red = C(:, 2:3);
D_red = 0;  

% define the reduced state space model
sys_red = ss(A_red, B_red, C_red, D_red);
sys_tf_red = tf(sys_red);

figure
step(sys_tf_red(2))
hold on
step(TF(2))
legend('reduced', 'full')
ylabel('Amplitude w')
hold off
grid on

figure
step(sys_tf_red(3))
hold on
step(TF(3))
legend('reduced', 'full')
ylabel('Amplitude q')
hold off
grid on

% Calculate the desired poles
zeta = 0.6;
omega_s = 2;
poles = [-zeta * omega_s + 1i * omega_s * sqrt(1 - zeta^2), real_part - 1i * omega_s * sqrt(1 - zeta^2)];

% Calculate the state feedback matrix K for the reduced system
K_red = place(A_red, B_red, poles);

% Check the poles of the closed-loop system
A_cl_red = A_red - B_red * K_red;
cl_red_poles = eig(A_cl_red);
disp('Reduced close loop poles:');
disp(cl_red_poles);

% Define the new state-space model with feedback
C_red = [1 0; 0 1];
sys_red_cl = ss(A_cl_red, B_red, C_red, D_red);
sys_tf_red_cl = tf(sys_red_cl);

% uncomment for output

% figure
% step(sys_tf_red(2))
% hold on
% step(sys_tf_red_cl(1))
% legend('open loop', 'closed loop')
% ylabel('Amplitude w')
% hold off
% grid on
% 
% figure
% step(sys_tf_red(3))
% hold on
% step(sys_tf_red_cl(2))
% legend('open loop', 'closed loop')
% ylabel('Amplitude q')
% hold off
% grid on


% full model closed loop

K_full = [0, K_red, 0];
A_full_cl = A - B*K_full;
sys_full_cl = ss(A_full_cl , B ,C,D);
sys_tf_full_cl = tf(sys_full_cl);

% uncomment for output

% figure
% step(sys_tf_full_cl(2))
% hold on
% step(sys_tf_full(2))
% legend('open loop', 'closed loop')
% ylabel('Amplitude w')
% hold off
% grid on

figure
step(sys_tf_full_cl(3))
hold on
step(sys_tf_full(3))
legend('open loop', 'closed loop')
ylabel('Amplitude q')
hold off
grid on
