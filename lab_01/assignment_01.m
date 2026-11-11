clear all
close all
clc

%% DH Parameter Symbols (Standard Definition)
syms d theta a alpha real

%% General Homogeneous Transformation Matrix (Denavit-Hartenberg)
% Standard DH transformation from frame i-1 to frame i
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)   a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha)   a*sin(theta); ...
     0           sin(alpha)             cos(alpha)              d; ...
     0           0                      0                       1];

%% Load Robot Model from URDF
robot = importrobot('../00001_PRR_zzy.urdf');
robot.DataFormat = 'row';

showdetails(robot);

eeName   = robot.BodyNames{end}; % End-effector frame name ('ee')
baseName = robot.BaseName;        % Base frame name ('base_link')

% Fixed link structural dimensions [m] from URDF geometry
a1 = 0.50; % Length of link 1 (prismatic travel offset)
a2 = 0.40; % Offset of link 2 along x-axis
a3 = 0.40; % Length of link 3 along z-axis to the end-effector

%% Define Joint Variables and Kinematic Parameters
% Joint 1: Prismatic along z0 (variable d1 = q1)
% Joint 2: Revolute about z1  (variable theta2 = q2)
% Joint 3: Revolute about z2  (variable theta3 = q3)
syms q1 q2 q3 real

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
% Denavit-Hartenberg Table: [ d | theta | a | alpha ]                     %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
DH_table = [q1   0    0     0;
            0    q2   a2  -pi/2;
            0    q3   0    pi/2];

%% ASSIGNMENT #1: KINEMATICS AND JACOBIANS (MANUAL & BUILT-IN COMPARISON)

% 1. Direct Kinematics - Manual Symbolic Derivation

fprintf('\n########### Direct Kinematics - Manual ###########\n');

% Intermediate transformation matrices between consecutive DH frames
A01 = subs(A, [d theta a alpha], DH_table(1,:));
A12 = subs(A, [d theta a alpha], DH_table(2,:));
A23 = subs(A, [d theta a alpha], DH_table(3,:));

% Base frame offset transformation (URDF base_link to DH frame 0):
% 0.15 m pedestal + 0.50 m link_1 housing offset along z
Ab0 = [ 1 0 0  0;
        0 1 0  0;
        0 0 1  0.15 + a1;
        0 0 0  1];

% End-effector tool frame transformation (DH frame 3 to 'ee' body):
% 0.40 m along z3 to tool tip
A3e = [ 1 0 0  0;
        0 1 0  0;
        0 0 1  a3;
        0 0 0  1];

% Complete transformation from base_link to end-effector
T03 = A01 * A12 * A23;
Tbe = simplify(Ab0 * T03 * A3e);

pe  = Tbe(1:3, 4);   % End-effector position vector [px; py; pz]
Rbe = Tbe(1:3, 1:3); % End-effector rotation matrix

fprintf('End-effector Homogeneous Transformation Matrix Tbe(q):\n');
disp(Tbe);

fprintf('End-effector Position Vector pe(q):\n');
disp(pe);

% 2. Direct Kinematics - Comparison with Built-In Function

fprintf('\n########### Direct Kinematics - Comparison ###########\n');
fprintf('Base frame: %s\n', baseName);
fprintf('End-effector body: %s\n', eeName);

% Test joint configuration
q_check = [-0.20, 0.40, 0.60];

% Manual symbolic FK evaluated at q_check
T_manual_sym = subs(Tbe, [q1, q2, q3], q_check);
T_manual     = double(T_manual_sym);

% Built-in FK from MATLAB Robotics System Toolbox
T_builtin = getTransform(robot, q_check, eeName, baseName);

% Matrix and pose discrepancies
T_error        = T_manual - T_builtin;
max_error      = max(abs(T_error(:)));
position_error = norm(T_manual(1:3,4) - T_builtin(1:3,4));
rotation_error = norm(T_manual(1:3,1:3) - T_builtin(1:3,1:3), 'fro');

fprintf('\nManual FK Transformation Matrix:\n');
disp(T_manual);

fprintf('Built-in FK Transformation Matrix:\n');
disp(T_builtin);

fprintf('Maximum matrix error: %.3e\n', max_error);
fprintf('Position error:       %.3e m\n', position_error);
fprintf('Rotation error:       %.3e\n', rotation_error);

tol = 1e-9;
if max_error < tol
    fprintf('FK CROSS-CHECK: MATCH\n');
else
    fprintf('FK CROSS-CHECK: DIFFER\n');
end


% 3. Inverse Kinematics - Manual Analytical Derivation

fprintf('\n########### Inverse Kinematics - Manual ###########\n');

% Desired position variables in the base frame
syms px_d py_d pz_d real

% Radial projection on the xy plane
r = sqrt(px_d^2 + py_d^2);

% Joint 2: Rotation in the horizontal xy plane
q2Sol = simplify(atan2(py_d, px_d));

% Joint 3: Elevation angle from radial arm offset (two trigonometric branches)
sin_q3   = (r - a2) / a3;
q3Sol_1 = simplify(asin(sin_q3));        % Branch 1 (Elbow Up)
q3Sol_2 = simplify(pi - asin(sin_q3));   % Branch 2 (Elbow Down)

% Joint 1: Prismatic vertical translation compensating link 3 tilt and base offset
q1Sol_1 = simplify(pz_d - a3*cos(q3Sol_1) - 0.65);
q1Sol_2 = simplify(pz_d - a3*cos(q3Sol_2) - 0.65);

fprintf('Analytical Inverse Kinematics Solution - Branch 1:\n');
fprintf('q1 = '); disp(q1Sol_1);
fprintf('q2 = '); disp(q2Sol);
fprintf('q3 = '); disp(q3Sol_1);

fprintf('Analytical Inverse Kinematics Solution - Branch 2:\n');
fprintf('q1 = '); disp(q1Sol_2);
fprintf('q2 = '); disp(q2Sol);
fprintf('q3 = '); disp(q3Sol_2);


% 4. Inverse Kinematics - Comparison with Built-In IK Solver

fprintf('\n########### Inverse Kinematics - Comparison ###########\n');

% Reachable target end-effector position
px_t = 0.45;
py_t = 0.20;
pz_t = 0.80;

ik_vars = [px_d, py_d, pz_d];
ik_vals = [px_t, py_t, pz_t];

% Evaluate Branch 1 numerically
q1_val_1 = double(subs(q1Sol_1, ik_vars, ik_vals));
q2_val_1 = double(subs(q2Sol,   ik_vars, ik_vals));
q3_val_1 = double(subs(q3Sol_1, ik_vars, ik_vals));
q_analytical_1 = [q1_val_1, q2_val_1, q3_val_1];

% Evaluate Branch 2 numerically
q1_val_2 = double(subs(q1Sol_2, ik_vars, ik_vals));
q2_val_2 = double(subs(q2Sol,   ik_vars, ik_vals));
q3_val_2 = double(subs(q3Sol_2, ik_vars, ik_vals));
q_analytical_2 = [q1_val_2, q2_val_2, q3_val_2];

fprintf('Analytical IK Solution (Branch 1): [%.6f, %.6f, %.6f]\n', q_analytical_1);
fprintf('Analytical IK Solution (Branch 2): [%.6f, %.6f, %.6f]\n', q_analytical_2);

% Target pose for built-in inverse kinematics
T_target = eye(4);
T_target(1:3, 4) = [px_t; py_t; pz_t];

% Instantiate MATLAB inverse kinematics solver
ikSolver = inverseKinematics('RigidBodyTree', robot);
weights  = [0 0 0 1 1 1]; % Position tracking prioritized, orientation free

% --- Branch 1 Evaluation ---
[q_sol_1, solInfo_1] = ikSolver(eeName, T_target, weights, q_analytical_1);
T_analytical_1 = double(subs(Tbe, [q1, q2, q3], q_analytical_1));
T_solution_1   = getTransform(robot, q_sol_1, eeName, baseName);
err_analytical_1 = norm(T_analytical_1(1:3,4) - [px_t; py_t; pz_t]);
err_builtin_1    = norm(T_solution_1(1:3,4)   - [px_t; py_t; pz_t]);

fprintf('\n--- Branch 1 Validation ---\n');
fprintf('Analytical IK q: [%.6f, %.6f, %.6f]\n', q_analytical_1);
fprintf('Built-in   IK q: [%.6f, %.6f, %.6f]\n', q_sol_1);
fprintf('Analytical FK Position Error: %.3e m\n', err_analytical_1);
fprintf('Built-in   FK Position Error: %.3e m\n', err_builtin_1);
if err_analytical_1 < 1e-6
    fprintf('ANALYTICAL IK (Branch 1): MATCH\n');
else
    fprintf('ANALYTICAL IK (Branch 1): DIFFER\n');
end
if err_builtin_1 < 1e-6
    fprintf('BUILT-IN IK   (Branch 1): MATCH\n');
else
    fprintf('BUILT-IN IK   (Branch 1): DIFFER\n');
end

% --- Branch 2 Evaluation ---
[q_sol_2, solInfo_2] = ikSolver(eeName, T_target, weights, q_analytical_2);
T_analytical_2 = double(subs(Tbe, [q1, q2, q3], q_analytical_2));
T_solution_2   = getTransform(robot, q_sol_2, eeName, baseName);
err_analytical_2 = norm(T_analytical_2(1:3,4) - [px_t; py_t; pz_t]);
err_builtin_2    = norm(T_solution_2(1:3,4)   - [px_t; py_t; pz_t]);

fprintf('\n--- Branch 2 Validation ---\n');
fprintf('Analytical IK q: [%.6f, %.6f, %.6f]\n', q_analytical_2);
fprintf('Built-in   IK q: [%.6f, %.6f, %.6f]\n', q_sol_2);
fprintf('Analytical FK Position Error: %.3e m\n', err_analytical_2);
fprintf('Built-in   FK Position Error: %.3e m\n', err_builtin_2);
if err_analytical_2 < 1e-6
    fprintf('ANALYTICAL IK (Branch 2): MATCH\n');
else
    fprintf('ANALYTICAL IK (Branch 2): DIFFER\n');
end
if err_builtin_2 < 1e-6
    fprintf('BUILT-IN IK   (Branch 2): MATCH\n');
else
    fprintf('BUILT-IN IK   (Branch 2): DIFFER\n');
end


% 5. Visualization of the Two IK Branches

figure('Name', 'Assignment 1 - Inverse Kinematics Branches', 'NumberTitle', 'off');

% Plot Branch 1
ax1 = subplot(1, 2, 1);
show(robot, q_analytical_1, 'Frames', 'on', 'Parent', ax1);
hold(ax1, 'on');
scatter3(ax1, px_t, py_t, pz_t, 100, 'ro', 'filled', 'MarkerEdgeColor', 'k');
grid(ax1, 'on');
axis(ax1, 'equal');
xlabel(ax1, 'X [m]'); ylabel(ax1, 'Y [m]'); zlabel(ax1, 'Z [m]');
title(ax1, 'Analytical IK - Branch 1 (Elbow Up)');
view(ax1, 3);

% Plot Branch 2
ax2 = subplot(1, 2, 2);
show(robot, q_analytical_2, 'Frames', 'on', 'Parent', ax2);
hold(ax2, 'on');
scatter3(ax2, px_t, py_t, pz_t, 100, 'ro', 'filled', 'MarkerEdgeColor', 'k');
grid(ax2, 'on');
axis(ax2, 'equal');
xlabel(ax2, 'X [m]'); ylabel(ax2, 'Y [m]'); zlabel(ax2, 'Z [m]');
title(ax2, 'Analytical IK - Branch 2 (Elbow Down)');
view(ax2, 3);

drawnow;

% 6. Geometric Jacobian [P R R] - Manual Symbolic Derivation

fprintf('\n########### Geometric Jacobian [P R R] - Manual ###########\n');

% Joint 1: Prismatic joint along z0
z0 = Ab0(1:3, 3);
p0 = Ab0(1:3, 4);

% Joint 2: Revolute joint about z1
Ab1 = Ab0 * A01;
z1  = Ab1(1:3, 3);
p1  = Ab1(1:3, 4);

% Joint 3: Revolute joint about z2
Ab2 = Ab0 * A01 * A12;
z2  = Ab2(1:3, 3);
p2  = Ab2(1:3, 4);

% Columns of the linear velocity Jacobian J_P:
% Joint 1 (Prismatic): jp1 = z0
% Joint 2 (Revolute):  jp2 = z1 x (pe - p1)
% Joint 3 (Revolute):  jp3 = z2 x (pe - p2)
jp1 = z0;
jp2 = cross(z1, pe - p1);
jp3 = cross(z2, pe - p2);
JGv = simplify([jp1, jp2, jp3]);

% Columns of the angular velocity Jacobian J_O:
% Joint 1 (Prismatic): jo1 = [0; 0; 0]
% Joint 2 (Revolute):  jo2 = z1
% Joint 3 (Revolute):  jo3 = z2
jo1 = zeros(3, 1);
jo2 = z1;
jo3 = z2;
JGw = simplify([jo1, jo2, jo3]);

% Full Geometric Jacobian JG = [JGv; JGw]
JG = [JGv; JGw];

fprintf('Geometric Linear Jacobian JGv (3x3):\n');
disp(JGv);

fprintf('Geometric Angular Jacobian JGw (3x3):\n');
disp(JGw);

fprintf('Complete Geometric Jacobian JG (6x3):\n');
disp(JG);

% 7. Geometric Jacobian [P R R] - Comparison with Built-In Function

fprintf('\n########### Geometric Jacobian - Comparison ###########\n');

% Joint configuration chosen for evaluation
qJacobian = q_analytical_1;

% Numerical evaluation of manual geometric Jacobian
JG_manual_sym = subs(JG, [q1, q2, q3], qJacobian);
JG_manual     = double(JG_manual_sym);

% MATLAB built-in geometric Jacobian returns: [Jw; Jv]
J_builtin = geometricJacobian(robot, qJacobian, eeName);

% Reorder manual Jacobian to match MATLAB ordering: [JGw; JGv]
JG_manual_matlab_order = [JG_manual(4:6, :); JG_manual(1:3, :)];

% Evaluate numerical differences
JG_error      = JG_manual_matlab_order - J_builtin;
max_JG_error  = max(abs(JG_error(:)));
linear_error  = norm(JG_manual_matlab_order(4:6,:) - J_builtin(4:6,:), 'fro');
angular_error = norm(JG_manual_matlab_order(1:3,:) - J_builtin(1:3,:), 'fro');

fprintf('Manual Geometric Jacobian [Jw; Jv]:\n');
disp(JG_manual_matlab_order);

fprintf('Built-in Geometric Jacobian [Jw; Jv]:\n');
disp(J_builtin);

fprintf('Maximum absolute error: %.3e\n', max_JG_error);
fprintf('Linear Jacobian error:  %.3e\n', linear_error);
fprintf('Angular Jacobian error: %.3e\n', angular_error);

tol_JG = 1e-9;
if max_JG_error < tol_JG
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: MATCH\n');
else
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: DIFFER\n');
end


%% ========================================================================
%% 8. Analytical Jacobian - Manual Symbolic Derivation (ZYZ Euler Angles)
%% ========================================================================

fprintf('\n########### Analytical Jacobian - Manual (ZYZ) ###########\n');

% Extract ZYZ Euler angles analytically from end-effector rotation matrix Rbe
% R = Rz(phi) * Ry(theta) * Rz(psi)
phi_euler   = atan2(Rbe(2,3), Rbe(1,3));
theta_euler = atan2(sqrt(Rbe(1,3)^2 + Rbe(2,3)^2), Rbe(3,3));
psi_euler   = atan2(Rbe(3,2), -Rbe(3,1));

euler_angles = [phi_euler; theta_euler; psi_euler];

% Translational Analytical Jacobian: JAv = d(pe)/dq
JAv = jacobian(pe, [q1, q2, q3]);

% Euler-angle-rate Analytical Jacobian: JAw = d(euler_angles)/dq
JAw = jacobian(euler_angles, [q1, q2, q3]);

% Complete Analytical Jacobian: JA = [JAv; JAw]
JA = simplify([JAv; JAw]);

fprintf('ZYZ Euler Angles [phi; theta; psi]:\n');
disp(simplify(euler_angles));

fprintf('Translational Analytical Jacobian JAv (3x3):\n');
disp(simplify(JAv));

fprintf('Euler-Rate Analytical Jacobian JAw (3x3):\n');
disp(simplify(JAw));


%% ========================================================================
%% 9. Analytical Jacobian - Comparison with Built-In Function
%% ========================================================================

fprintf('\n########### Analytical Jacobian - Comparison (ZYZ) ###########\n');

% Substitute configuration into manual analytical Jacobian
JA_manual_sym = subs(JA, [q1, q2, q3], qJacobian);
JA_manual     = double(JA_manual_sym);

% Built-in geometric angular and linear velocities
Jw_builtin = J_builtin(1:3, :);
Jv_builtin = J_builtin(4:6, :);

% Evaluate numerical rotation matrix at current configuration
R_num = double(subs(Rbe, [q1, q2, q3], qJacobian));

% Extract numerical Euler angles for conversion matrix
phi_num   = atan2(R_num(2,3), R_num(1,3));
theta_num = atan2(sqrt(R_num(1,3)^2 + R_num(2,3)^2), R_num(3,3));

% Transformation matrix relating Euler angle rates to angular velocity:
% omega = E * [phi_dot; theta_dot; psi_dot]
E = [0, -sin(phi_num), cos(phi_num)*sin(theta_num);
     0,  cos(phi_num), sin(phi_num)*sin(theta_num);
     1,  0,            cos(theta_num)];

% Invert E to obtain Euler angle rate Jacobian: J_euler = E \ Jw
J_euler_builtin = E \ Jw_builtin;

% Full converted built-in analytical Jacobian: [Jv; J_euler]
JA_builtin = [Jv_builtin; J_euler_builtin];

% Evaluate discrepancies
JA_error        = JA_manual - JA_builtin;
max_JA_error    = max(abs(JA_error(:)));
linear_JA_err   = norm(JA_manual(1:3,:) - JA_builtin(1:3,:), 'fro');
angular_JA_err  = norm(JA_manual(4:6,:) - JA_builtin(4:6,:), 'fro');

fprintf('Manual Analytical Jacobian:\n');
disp(JA_manual);

fprintf('Built-in Jacobian converted to ZYZ Euler rates:\n');
disp(JA_builtin);

fprintf('Maximum absolute error:    %.3e\n', max_JA_error);
fprintf('Linear Jacobian error:     %.3e\n', linear_JA_err);
fprintf('Euler-rate Jacobian error: %.3e\n', angular_JA_err);

tol_JA = 1e-9;
if max_JA_error < tol_JA
    fprintf('ANALYTICAL JACOBIAN CROSS-CHECK: MATCH\n');
else
    fprintf('ANALYTICAL JACOBIAN CROSS-CHECK: DIFFER\n');
end
