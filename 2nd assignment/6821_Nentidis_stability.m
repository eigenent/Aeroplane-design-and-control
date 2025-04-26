% Dimitrios Nentidis AEM: 6821

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Stability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%               Approach                  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C_m_α and static margin

% NACA 0006  symmerical
x_cg_w_approach = 1.443;
x_ac_w_approach = 1.4;
c_sound = 2.91;
C_l_a = 2*pi;

% for main wing
e_oswald = 0.99;
b = 6.7^2;
S = 18.22;
AR = b/S;

C_L_a_wing = 2.98; % by hand, Λ is taken into account
C_m_a_wing = C_L_a_wing*(x_cg_w_approach - x_ac_w_approach)/c_sound;


% for horizontal tail wing
Tail_arm = 5.3;
St = 4.2;
V_H = Tail_arm*St/(S*c_sound);
h = 1; %%% tail above wings
b_tail = 3.6;
AR_tail = (b_tail^2)/St;
e_tail =  0.99;

C_L_a_tail = C_l_a/(1 + (C_l_a/(pi*e_oswald*AR_tail))); %3D correction

deda = 2*C_L_a_wing/(e_tail*pi*AR);
C_m_a_tail = - h*V_H*C_L_a_tail*(1-deda);

C_m_a_total = C_m_a_tail + C_m_a_wing;
C_L_a_total = C_L_a_tail*(h*(St/S)) + C_L_a_wing;

% neutral point and static margin
x_np = (x_ac_w_approach/c_sound) + h*V_H*(C_L_a_tail/C_L_a_wing)*(1 - deda);
SM = x_np - (x_cg_w_approach/c_sound);


%% C_L_0 and C_m_0

% data
CL_ds = 0.68;
Cm_ds = -1.46;
AoA_approach = 2.3;
V = 87.5; 
mass = 6400;
g = 9.81;


Lift = mass*g;
C_L_approach = Lift/(0.5*1.2*(V^2)*S); %CL_approach = Cl0_wing + AoA*dCl_wing_da
C_L_0_wing = C_L_approach - C_L_a_wing*deg2rad(AoA_approach);
C_m_ac_wing = -1; % by hand
C_m_0_Wing = C_m_ac_wing + C_L_0_wing*((x_cg_w_approach - x_ac_w_approach)/c_sound);
e_0 = 2*C_L_a_tail/(pi*e_tail*AR_tail);

% for tail wing
C_m_0_tail =  h*V_H*C_L_a_tail*(e_0 );

% total
C_m_0_total = C_m_0_tail + C_m_0_Wing;

C_m_O = -C_m_a_total*(deg2rad(AoA_approach)) - Cm_ds*(deg2rad(-7.1));
C_L_O = CL_ds*(deg2rad(-7.1));


%% plotting C_m vs α

A = -5:0.1:5;

figure
plot(A , C_m_O + C_m_a_total*deg2rad(A) + Cm_ds*(deg2rad(-7.1)), 'r', 'LineWidth', 2);
hold on
plot(A, C_m_O + C_m_a_total*deg2rad(A), 'g', 'LineWidth', 2);
yline(0)
xline(0)
legend('δς = -7.1','δς = 0')
ylabel('Cm')
xlabel('α')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%                 Cruise                  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define C_L_α, C_m_α, C_L_δς, C_m_δς

% Define the coefficient matrix
coef = [ -5.524,  34.344, -78.940,  76.698, -21.539;  % C_L_a
         -4.485,  23.523, -41.934,  29.286,  -8.085;  % C_m_a
        0.00331, 0.878,  -3.694,   4.432,  -0.760;    % C_L_ds
          0.086,  -1.912,   6.986,  -7.934,   1.115   % C_m_ds
];

% Define a function to calculate the coefficient for a given Mach number
calculate_coefficient = @(coeffs, M) ...
    coeffs(1) * M^4 + coeffs(2) * M^3 + coeffs(3) * M^2 + coeffs(4) * M + coeffs(5);

% Define the individual coefficient functions
C_L_alpha = @(M) calculate_coefficient(coef(1, :), M);
C_m_alpha = @(M) calculate_coefficient(coef(2, :), M);
C_L_deltas = @(M) calculate_coefficient(coef(3, :), M);
C_m_deltas = @(M) calculate_coefficient(coef(4, :), M);


%% C_L_0 and C_m_0 for α = 3 and δς = -1.5, speed = 961 km/h

% by hand

%% Mach 0.7 - 2 trim points
c_sound = (1.4*287.05*218.808)^(0.5);
density = 0.38;
X = zeros(131,2);
i=1;
for M = 0.7:0.01:2

% Calculate the coefficients for the given Mach number
C_Lalpha_v = C_L_alpha(M);
C_malpha_v = C_m_alpha(M);
C_Ldeltas_v = C_L_deltas(M);
C_mdeltas_v = C_m_deltas(M);

% parameters
AoA = 3;
delta_s = -1.5;
CL(i) = (7400*9.81)/(0.5*((M*c_sound)^2)*0.38*(18.22));
C_m_00 = 0.0506;

% system of equations
A = [C_Ldeltas_v C_Lalpha_v; 
     C_mdeltas_v C_malpha_v];
B = [CL(i); 
    -C_m_00];

x = A\B;

% storage
X(i,:) = x';

% keep count
i=i+1;
end

% plot
figure
plot(0.7:0.01:2 , rad2deg(X(:,2)),'g','LineWidth',2)
hold on
plot(0.7:0.01:2 , rad2deg(X(:,1)),'r','LineWidth',2)
ylabel('[Angle in degrees]')
xlabel('Mach')
legend('α','δs')
grid on

%% 2 M, -700 kg

alpha = -10:0.01:15;
X5 = zeros(length(alpha),2);
M = 2;
mass = 7400 - 700;
C_L_target = mass*g/(0.5*density*((M*c_sound)^2)*S);
C_m_00 = 0.0506;


delta_s = zeros(10,1);
AoA_current = zeros(10,1);
error = 1e-5;
for j = 2:10

% for L = W, find current α for defined δς
AoA_current(j) = (C_L_target - deg2rad(delta_s(j-1))*C_Ldeltas_v)/(C_Lalpha_v);

% calculate C_m
Cm_current = C_m_00 +  AoA_current(j)*C_malpha_v + C_mdeltas_v*(delta_s(j-1));

% calculate new δς to trim for current α
delta_s(j) = -(C_m_00 + AoA_current(j)*C_malpha_v)/(C_mdeltas_v);

% convergange check
Dif = (delta_s(j) - delta_s(j - 1));
if abs(Dif) <= error

    delta_s_optimal = delta_s(j);
break
end

end

% plotting
CL_plot = -0.8:0.01:0.3;
Cm_plot1 = C_m_00 + ((CL_plot - delta_s(1)*C_Ldeltas_v )/C_Lalpha_v)*C_malpha_v + delta_s(1)*C_mdeltas_v;
Cm_plot2 = C_m_00 + ((CL_plot - delta_s(2)*C_Ldeltas_v )/C_Lalpha_v)*C_malpha_v + delta_s(2)*C_mdeltas_v;


% Create a tiled layout with 1 row and 2 columns
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% First plot

nexttile;
Aaaa = -2:0.01:2;
plot(Aaaa ,C_Ldeltas_v*delta_s(2) + C_Lalpha_v*deg2rad(Aaaa),'r','LineWidth',2  )
hold on
plot(Aaaa ,C_Ldeltas_v*delta_s(1) + C_Lalpha_v*deg2rad(Aaaa), 'g', 'LineWidth',2 )
xline(0, 'k', 'LineWidth', 0.5); % Add vertical line at x=0
yline(0, 'k', 'LineWidth', 0.5); % Add horizontal line at y=0
yline(C_L_target, '--k', 'LineWidth', 0.5)
legend('trimmed', 'δς=0')
xlabel('α')
ylabel('CL')
hold off
axis tight
set(gca, 'PlotBoxAspectRatio', [1 2 1]); % Taller-than-wider aspect ratio

% Second plot
nexttile;
plot(Cm_plot1 , CL_plot , 'g','LineWidth',2)
hold on
plot(Cm_plot2 , CL_plot , 'r','LineWidth',2)
xline(0, 'k', 'LineWidth', 0.5); % Add vertical line at x=0
yline(0, 'k', 'LineWidth', 0.5); % Add horizontal line at y=0
ylabel('CL')
xlabel('Cm')
yline(C_L_target, '--k', 'LineWidth', 0.5)
set(gca, 'XDir', 'reverse'); % Reverse the direction of the x-axis
legend('δς=0', 'trimmed')
hold off
axis tight
set(gca, 'PlotBoxAspectRatio', [1 2 1]); % Taller-than-wider aspect ratio







va = (4*K/3/Cd0)^(0.25) *(W./R/S).^0.5;

vv = (P - 0.5*R.*vc.^3*S*Cd0 - K*W^2./(0.5*R.*vc.*S))/W;

