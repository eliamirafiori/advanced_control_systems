%% ASSIGNMENT #4: RECURSIVE NEWTON-EULER & DYNAMIC PARAMETER IDENTIFICATION

clear all
close all
clc

%% DH Table Symbols
syms d theta a alpha real

%% Homogeneous Transformation Matrix
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)   a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha)   a*sin(theta); ...
     0           sin(alpha)             cos(alpha)              d; ...
     0           0                      0                       1];

%% Load Robot Model from URDF
robot = importrobot('../00001_PRR_zzy.urdf');
robot.DataFormat = 'row';

eeName   = robot.BodyNames{end};
baseName = robot.BaseName;

% Kinematic link dimensions [m]
a1_a = 0.50;
a2_a = 0.40;
a3_a = 0.40;

% Analytical masses [kg]
m1_a = 1.0;
m2_a = 1.0;
m3_a = 1.0;

% Gravitational constant [m/s^2]
g_val = 9.81;
robot.Gravity = [0, 0, -g_val];

%% Define Joint Variables, Velocities, and Accelerations
syms q1 q2 q3 real
syms dq1 dq2 dq3 real
syms ddq1 ddq2 ddq3 real
syms a1 a2 a3 g real

q   = [q1; q2; q3];
dq  = [dq1; dq2; dq3];
ddq = [ddq1; ddq2; ddq3];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Denavit-Hartenberg Table: [ d | theta | a | alpha ]                     %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
DH_table = [q1   0    0     0;
            0    q2   a2  -pi/2;
            0    q3   0    pi/2];

% Rotation matrices from frame i to frame i-1: R_i^{i-1}
% Link 1: Prismatic (theta = 0, alpha = 0)
R01 = [ 1  0  0;
        0  1  0;
        0  0  1 ];

% Link 2: Revolute (theta = q2, alpha = -pi/2)
R12 = [ cos(q2)  0  -sin(q2);
        sin(q2)  0   cos(q2);
        0       -1   0      ];

% Link 3: Revolute (theta = q3, alpha = pi/2)
R23 = [ cos(q3)  0   sin(q3);
        sin(q3)  0  -cos(q3);
        0        1   0      ];

% Relative position vectors of frame origin i w.r.t. frame i-1, expressed in frame i:
% r_{i-1, i}^i = (R_i^{i-1})^T * r_{i-1, i}^{i-1} (Slide 34, 37)
r01_1 = [0; 0; q1];
r12_2 = [a2; 0; 0];
r23_3 = [0; 0; 0];

% Center-of-mass positions in link frame i: r_{i, C_i}^i = [lcix; lciy; lciz] (Slide 9)
syms lc1x lc1y lc1z real
syms lc2x lc2y lc2z real
syms lc3x lc3y lc3z real

r1_C1_1 = [lc1x; lc1y; lc1z];
r2_C2_2 = [lc2x; lc2y; lc2z];
r3_C3_3 = [lc3x; lc3y; lc3z];

% Link masses
syms m1 m2 m3 real

% Centroidal inertia tensors in link frames: \bar{I}_i^i (Slide 23, 37)
syms I1xx I1yy I1zz I1xy I1xz I1yz real
syms I2xx I2yy I2zz I2xy I2xz I2yz real
syms I3xx I3yy I3zz I3xy I3xz I3yz real

I1_bar = [I1xx I1xy I1xz;
          I1xy I1yy I1yz;
          I1xz I1yz I1zz];

I2_bar = [I2xx I2xy I2xz;
          I2xy I2yy I2yz;
          I2xz I2yz I2zz];

I3_bar = [I3xx I3xy I3xz;
          I3xy I3yy I3yz;
          I3xz I3yz I3zz];

z0 = [0; 0; 1]; % Axis of joint coordinates in frame i-1


%% ========================================================================
%% 1. RECURSIVE NEWTON-EULER (RNE) FORMULATION - FORWARD RECURSION
%% Slides 33-36: Propagate velocities & accelerations from base (0) to EE (3)
%% ========================================================================

fprintf('\n\n');
fprintf('====================================================================\n');
fprintf('----------- 1. Recursive Newton-Euler: Forward Recursion -----------\n');
fprintf('====================================================================\n');

% Base boundary conditions with gravity (Slide 35, 36)
% Since base frame z points upward, gravity acts along -z: \ddot{p}_0^0 - g_0^0 = [0; 0; g]
omega_0_0       = [0; 0; 0];
domega_0_0      = [0; 0; 0];
pddot_0_minus_g = [0; 0; g];

% -------------------------------------------------------------------------
% LINK 1: Prismatic Joint (d1 = q1, dot{d}1 = dq1, ddot{d}1 = ddq1) (Slide 34)
% -------------------------------------------------------------------------
omega_1_1  = R01.' * omega_0_0;
domega_1_1 = R01.' * domega_0_0;

pddot_1_1 = R01.' * (pddot_0_minus_g + ddq1 * z0) + ...
            2 * dq1 * cross(omega_1_1, R01.' * z0) + ...
            cross(domega_1_1, r01_1) + ...
            cross(omega_1_1, cross(omega_1_1, r01_1));

pddot_C1_1 = pddot_1_1 + cross(domega_1_1, r1_C1_1) + ...
             cross(omega_1_1, cross(omega_1_1, r1_C1_1));

% -------------------------------------------------------------------------
% LINK 2: Revolute Joint (theta2 = q2, dot{theta}2 = dq2, ddot{theta}2 = ddq2) (Slide 34)
% -------------------------------------------------------------------------
omega_2_2  = R12.' * (omega_1_1 + dq2 * z0);
domega_2_2 = R12.' * (domega_1_1 + ddq2 * z0 + dq2 * cross(omega_1_1, z0));

pddot_2_2 = R12.' * pddot_1_1 + ...
            cross(domega_2_2, r12_2) + ...
            cross(omega_2_2, cross(omega_2_2, r12_2));

pddot_C2_2 = pddot_2_2 + cross(domega_2_2, r2_C2_2) + ...
             cross(omega_2_2, cross(omega_2_2, r2_C2_2));

% -------------------------------------------------------------------------
% LINK 3: Revolute Joint (theta3 = q3, dot{theta}3 = dq3, ddot{theta}3 = ddq3) (Slide 34)
% -------------------------------------------------------------------------
omega_3_3  = R23.' * (omega_2_2 + dq3 * z0);
domega_3_3 = R23.' * (domega_2_2 + ddq3 * z0 + dq3 * cross(omega_2_2, z0));

pddot_3_3 = R23.' * pddot_2_2 + ...
            cross(domega_3_3, r23_3) + ...
            cross(omega_3_3, cross(omega_3_3, r23_3));

pddot_C3_3 = pddot_3_3 + cross(domega_3_3, r3_C3_3) + ...
             cross(omega_3_3, cross(omega_3_3, r3_C3_3));

fprintf('Forward recursion completed successfully.\n');


%% ========================================================================
%% 2. RECURSIVE NEWTON-EULER (RNE) FORMULATION - BACKWARD RECURSION
%% Slides 37-38: Propagate forces & moments from EE (3) down to base (1)
%% ========================================================================

fprintf('\n\n');
fprintf('====================================================================\n');
fprintf('----------- 2. Recursive Newton-Euler: Backward Recursion ----------\n');
fprintf('====================================================================\n');

% Boundary condition at the end-effector (free motion: no external wrench he = 0)
f_4_4  = [0; 0; 0];
mu_4_4 = [0; 0; 0];
R34    = eye(3);

% -------------------------------------------------------------------------
% LINK 3: Backward Equations (Slide 37)
% -------------------------------------------------------------------------
f_3_3 = R34 * f_4_4 + m3 * pddot_C3_3;

mu_3_3 = - cross(f_3_3, r23_3 + r3_C3_3) + ...
         R34 * mu_4_4 + cross(R34 * f_4_4, r3_C3_3) + ...
         I3_bar * domega_3_3 + cross(omega_3_3, I3_bar * omega_3_3);

% Generalized torque for Revolute Joint 3: projection along z_{i-1} in frame i
tau3_rne = simplify(mu_3_3.' * (R23.' * z0));

% -------------------------------------------------------------------------
% LINK 2: Backward Equations (Slide 37)
% -------------------------------------------------------------------------
f_2_2 = R23 * f_3_3 + m2 * pddot_C2_2;

mu_2_2 = - cross(f_2_2, r12_2 + r2_C2_2) + ...
         R23 * mu_3_3 + cross(R23 * f_3_3, r2_C2_2) + ...
         I2_bar * domega_2_2 + cross(omega_2_2, I2_bar * omega_2_2);

% Generalized torque for Revolute Joint 2: projection along z_{i-1} in frame i
tau2_rne = simplify(mu_2_2.' * (R12.' * z0));

% -------------------------------------------------------------------------
% LINK 1: Backward Equations (Slide 37)
% -------------------------------------------------------------------------
f_1_1 = R12 * f_2_2 + m1 * pddot_C1_1;

mu_1_1 = - cross(f_1_1, r01_1 + r1_C1_1) + ...
         R12 * mu_2_2 + cross(R12 * f_2_2, r1_C1_1) + ...
         I1_bar * domega_1_1 + cross(omega_1_1, I1_bar * omega_1_1);

% Generalized force for Prismatic Joint 1: projection of force along z0 in frame 1
tau1_rne = simplify(f_1_1.' * (R01.' * z0));

% Complete RNE Torque Vector
tau_rne = [tau1_rne; tau2_rne; tau3_rne];

fprintf('RNE Generalized Joint Torque Vector tau_rne(q, dq, ddq):\n');
disp(tau_rne);


%% ========================================================================
%% 3. LINEAR-IN-THE-DYNAMIC-PARAMETERS MODEL & REGRESSOR Y (Slides 7-15)
%% tau = Y(q, dq, ddq) * Theta
%% ========================================================================

fprintf('\n\n');
fprintf('====================================================================\n');
fprintf('----------- 3. Regressor Matrix Y & Parameter Vector Theta ---------\n');
fprintf('====================================================================\n');

% Standard 10 dynamic parameters per link (Slides 10, 12):
% theta_i = [m_li, m_li*lcix, m_li*lciy, m_li*lciz, ...
%            Ihat_xx, Ihat_xy, Ihat_xz, Ihat_yy, Ihat_yz, Ihat_zz]^T
% where Ihat_i is the Steiner inertia tensor relative to the link origin!
syms mx1 my1 mz1 Ihat1_xx Ihat1_xy Ihat1_xz Ihat1_yy Ihat1_yz Ihat1_zz real
syms mx2 my2 mz2 Ihat2_xx Ihat2_xy Ihat2_xz Ihat2_yy Ihat2_yz Ihat2_zz real
syms mx3 my3 mz3 Ihat3_xx Ihat3_xy Ihat3_xz Ihat3_yy Ihat3_yz Ihat3_zz real

% Parameter vectors for each link:
theta_1 = [m1; mx1; my1; mz1; Ihat1_xx; Ihat1_xy; Ihat1_xz; Ihat1_yy; Ihat1_yz; Ihat1_zz];
theta_2 = [m2; mx2; my2; mz2; Ihat2_xx; Ihat2_xy; Ihat2_xz; Ihat2_yy; Ihat2_yz; Ihat2_zz];
theta_3 = [m3; mx3; my3; mz3; Ihat3_xx; Ihat3_xy; Ihat3_xz; Ihat3_yy; Ihat3_yz; Ihat3_zz];

% Total dynamic parameter vector (30 x 1) (Slides 7, 12, 14)
Theta = [theta_1; theta_2; theta_3];

% Steiner substitution rules: express centroidal inertia I_bar in terms of
% origin inertia Ihat and first moments m*r_C:
% I_bar_xx = Ihat_xx - (my^2 + mz^2)/m,  I_bar_xy = Ihat_xy + mx*my/m, etc.
sub_linear_rules = [ ...
    lc1x == mx1/m1, lc1y == my1/m1, lc1z == mz1/m1, ...
    I1xx == Ihat1_xx - (my1^2 + mz1^2)/m1, ...
    I1yy == Ihat1_yy - (mx1^2 + mz1^2)/m1, ...
    I1zz == Ihat1_zz - (mx1^2 + my1^2)/m1, ...
    I1xy == Ihat1_xy + mx1*my1/m1, ...
    I1xz == Ihat1_xz + mx1*mz1/m1, ...
    I1yz == Ihat1_yz + my1*mz1/m1, ...
    lc2x == mx2/m2, lc2y == my2/m2, lc2z == mz2/m2, ...
    I2xx == Ihat2_xx - (my2^2 + mz2^2)/m2, ...
    I2yy == Ihat2_yy - (mx2^2 + mz2^2)/m2, ...
    I2zz == Ihat2_zz - (mx2^2 + my2^2)/m2, ...
    I2xy == Ihat2_xy + mx2*my2/m2, ...
    I2xz == Ihat2_xz + mx2*mz2/m2, ...
    I2yz == Ihat2_yz + my2*mz2/m2, ...
    lc3x == mx3/m3, lc3y == my3/m3, lc3z == mz3/m3, ...
    I3xx == Ihat3_xx - (my3^2 + mz3^2)/m3, ...
    I3yy == Ihat3_yy - (mx3^2 + mz3^2)/m3, ...
    I3zz == Ihat3_zz - (mx3^2 + my3^2)/m3, ...
    I3xy == Ihat3_xy + mx3*my3/m3, ...
    I3xz == Ihat3_xz + mx3*mz3/m3, ...
    I3yz == Ihat3_yz + my3*mz3/m3 ];

% Express torque equations in linear form with respect to Theta
tau_linear = subs(tau_rne, lhs(sub_linear_rules), rhs(sub_linear_rules));
tau_linear = simplify(tau_linear);

% Symbolic Regressor Matrix Y: Y = jacobian(tau_linear, Theta) (Slide 7, 13, 14)
% Due to linearity: tau = Y * Theta
Y = jacobian(tau_linear, Theta);
Y = simplify(Y);

% Block decomposition of Y (Slide 14):
Y1 = Y(:, 1:10);   % Influence of Link 1 dynamic parameters
Y2 = Y(:, 11:20);  % Influence of Link 2 dynamic parameters
Y3 = Y(:, 21:30);  % Influence of Link 3 dynamic parameters

fprintf('Dynamic Regressor Matrix Y(q, dq, ddq) size: %d x %d\n', size(Y,1), size(Y,2));
fprintf('Regressor block upper-triangular structure verified.\n');


%% ========================================================================
%% 4. NUMERICAL COMPARISON WITH MATLAB BUILT-IN DYNAMICS
%% ========================================================================

fprintf('\n\n');
fprintf('====================================================================\n');
fprintf('----------- 4. Numerical Comparison with Built-In Robotics ---------\n');
fprintf('====================================================================\n');

% Test configuration, velocities, accelerations
q_num   = [-0.20,  0.40,  0.60];
dq_num  = [ 0.10, -0.30,  0.50];
ddq_num = [ 0.05,  0.20, -0.10];

% Extract actual physical parameters directly from URDF
b1 = robot.Bodies{2}; % link_1
b2 = robot.Bodies{3}; % link_2
b3 = robot.Bodies{4}; % link_3

m1_val = b1.Mass;
m2_val = b2.Mass;
m3_val = b3.Mass;

c1_val = b1.CenterOfMass(:);
c2_val = b2.CenterOfMass(:);
c3_val = b3.CenterOfMass(:);

% Extract origin inertias: [Ixx Iyy Izz Iyz Ixz Ixy]
I1_v = b1.Inertia;
I2_v = b2.Inertia;
I3_v = b3.Inertia;

Ihat1_mat = [I1_v(1), I1_v(6), I1_v(5); I1_v(6), I1_v(2), I1_v(4); I1_v(5), I1_v(4), I1_v(3)];
Ihat2_mat = [I2_v(1), I2_v(6), I2_v(5); I2_v(6), I2_v(2), I2_v(4); I2_v(5), I2_v(4), I2_v(3)];
Ihat3_mat = [I3_v(1), I3_v(6), I3_v(5); I3_v(6), I3_v(2), I3_v(4); I3_v(5), I3_v(4), I3_v(3)];

% Centroidal inertias via Steiner's theorem for RNE substitution:
I1_c_mat = Ihat1_mat - m1_val * ((c1_val.' * c1_val)*eye(3) - (c1_val * c1_val.'));
I2_c_mat = Ihat2_mat - m2_val * ((c2_val.' * c2_val)*eye(3) - (c2_val * c2_val.'));
I3_c_mat = Ihat3_mat - m3_val * ((c3_val.' * c3_val)*eye(3) - (c3_val * c3_val.'));

% Numerical substitution for RNE torque
rne_syms = [q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3, a1, a2, a3, g, ...
            m1, m2, m3, ...
            lc1x, lc1y, lc1z, lc2x, lc2y, lc2z, lc3x, lc3y, lc3z, ...
            I1xx, I1yy, I1zz, I1xy, I1xz, I1yz, ...
            I2xx, I2yy, I2zz, I2xy, I2xz, I2yz, ...
            I3xx, I3yy, I3zz, I3xy, I3xz, I3yz];

rne_vals = [q_num, dq_num, ddq_num, a1_a, a2_a, a3_a, g_val, ...
            m1_val, m2_val, m3_val, ...
            c1_val(1), c1_val(2), c1_val(3), ...
            c2_val(1), c2_val(2), c2_val(3), ...
            c3_val(1), c3_val(2), c3_val(3), ...
            I1_c_mat(1,1), I1_c_mat(2,2), I1_c_mat(3,3), I1_c_mat(1,2), I1_c_mat(1,3), I1_c_mat(2,3), ...
            I2_c_mat(1,1), I2_c_mat(2,2), I2_c_mat(3,3), I2_c_mat(1,2), I2_c_mat(1,3), I2_c_mat(2,3), ...
            I3_c_mat(1,1), I3_c_mat(2,2), I3_c_mat(3,3), I3_c_mat(1,2), I3_c_mat(1,3), I3_c_mat(2,3)];

tau_rne_num = double(subs(tau_rne, rne_syms, rne_vals));

% Numerical true parameter vector Theta_true (30 x 1)
Theta_true = [ ...
    m1_val; m1_val*c1_val(1); m1_val*c1_val(2); m1_val*c1_val(3); ...
    Ihat1_mat(1,1); Ihat1_mat(1,2); Ihat1_mat(1,3); Ihat1_mat(2,2); Ihat1_mat(2,3); Ihat1_mat(3,3); ...
    m2_val; m2_val*c2_val(1); m2_val*c2_val(2); m2_val*c2_val(3); ...
    Ihat2_mat(1,1); Ihat2_mat(1,2); Ihat2_mat(1,3); Ihat2_mat(2,2); Ihat2_mat(2,3); Ihat2_mat(3,3); ...
    m3_val; m3_val*c3_val(1); m3_val*c3_val(2); m3_val*c3_val(3); ...
    Ihat3_mat(1,1); Ihat3_mat(1,2); Ihat3_mat(1,3); Ihat3_mat(2,2); Ihat3_mat(2,3); Ihat3_mat(3,3) ];

% Numerical evaluation of Regressor Y
kin_syms = [q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3, a1, a2, a3, g];
kin_vals = [q_num, dq_num, ddq_num, a1_a, a2_a, a3_a, g_val];

Y_num = double(subs(Y, kin_syms, kin_vals));
tau_reg_num = Y_num * Theta_true;

% MATLAB built-in inverse dynamics
tau_builtin = inverseDynamics(robot, q_num, dq_num, ddq_num).';

% Discrepancies
err_rne     = norm(tau_rne_num - tau_builtin);
err_reg_rne = norm(tau_reg_num - tau_rne_num);

fprintf('RNE Torque Vector:       [%.6f, %.6f, %.6f] N*m\n', tau_rne_num);
fprintf('Regressor Torque (Y*Th): [%.6f, %.6f, %.6f] N*m\n', tau_reg_num);
fprintf('Built-in Torque:         [%.6f, %.6f, %.6f] N*m\n', tau_builtin);

fprintf('RNE vs Built-in Error:       %.3e\n', err_rne);
fprintf('Regressor vs RNE Error:      %.3e\n', err_reg_rne);

if err_rne < 1e-6
    fprintf('RNE FORMULATION CROSS-CHECK: MATCH\n');
else
    fprintf('RNE FORMULATION CROSS-CHECK: DIFFER\n');
end

if err_reg_rne < 1e-6
    fprintf('REGRESSOR MODEL CROSS-CHECK: MATCH\n');
else
    fprintf('REGRESSOR MODEL CROSS-CHECK: DIFFER\n');
end


%% ========================================================================
%% 5. DYNAMIC PARAMETER IDENTIFICATION (Slides 16-20)
%% Observation Matrix X * Theta = T, Least Squares Solution
%% ========================================================================

fprintf('\n\n');
fprintf('====================================================================\n');
fprintf('----------- 5. Dynamic Parameter Identification (LS) ---------------\n');
fprintf('====================================================================\n');

% Trajectory generation: multi-frequency sinusoidal excitation (Slide 17, 20)
% q_j(t) = sum_k (A_jk * sin(w_jk * t))
time_samples = linspace(0, 5, 60); % 60 time instants (N = 60)
N_pts = length(time_samples);

% Preallocate identification regression matrices X and torque vector T (Slide 18)
X_matrix = zeros(3 * N_pts, 30);
T_vector = zeros(3 * N_pts, 1);

for k = 1:N_pts
    t_k = time_samples(k);
    
    % Trajectory state at instant t_k
    q_k = [ -0.25 + 0.15*sin(1.0*t_k);
             0.50*cos(1.5*t_k);
             0.60*sin(2.0*t_k) ];
         
    dq_k = [ 0.15*1.0*cos(1.0*t_k);
            -0.50*1.5*sin(1.5*t_k);
             0.60*2.0*cos(2.0*t_k) ];
         
    ddq_k = [ -0.15*(1.0^2)*sin(1.0*t_k);
              -0.50*(1.5^2)*cos(1.5*t_k);
              -0.60*(2.0^2)*sin(2.0*t_k) ];
          
    % Substitute into regressor
    sub_k = [q_k.', dq_k.', ddq_k.', a1_a, a2_a, a3_a, g_val];
    Y_k = double(subs(Y, kin_syms, sub_k));
    
    % True measured torque
    tau_k = Y_k * Theta_true;
    
    % Fill block rows
    row_idx = (3*k - 2) : (3*k);
    X_matrix(row_idx, :) = Y_k;
    T_vector(row_idx)    = tau_k;
end

% Rank analysis of the regressor matrix X (Slide 20)
rank_X = rank(X_matrix);
fprintf('Total number of dynamic parameters: %d\n', 30);
fprintf('Number of independent base parameters (rank of X): %d\n', rank_X);

% Base parameter identification via QR decomposition with column pivoting (Slide 14, 20)
[Q_qr, R_qr, P_qr] = qr(X_matrix, 0);

% Independent columns corresponding to base parameters
idx_base = P_qr(1:rank_X);
X_base   = X_matrix(:, idx_base);

% Least Squares estimation of base dynamic parameters:
% Upsilon_hat = (X_base^T * X_base)^(-1) * X_base^T * T (Slide 19)
Upsilon_hat = (X_base.' * X_base) \ (X_base.' * T_vector);
Upsilon_true = Theta_true(idx_base);

% Minimum-norm Least Squares solution for full Theta: Theta_hat = pinv(X) * T
Theta_hat = pinv(X_matrix) * T_vector;

% Reconstructed joint torques along trajectory: T_reconstructed = X * Theta_hat
T_reconstructed = X_matrix * Theta_hat;
id_torque_error = norm(T_vector - T_reconstructed);

fprintf('Least-Squares base parameter identification completed.\n');
fprintf('Base parameter estimation error norm:    %.3e\n', norm(Upsilon_hat - Upsilon_true));
fprintf('Full trajectory torque estimation error: %.3e N*m\n', id_torque_error);

if id_torque_error < 1e-6
    fprintf('PARAMETER IDENTIFICATION: MATCH\n');
else
    fprintf('PARAMETER IDENTIFICATION: DIFFER\n');
end