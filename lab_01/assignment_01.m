clear all
close all
clc

%% DH table syms (DON'T EDIT)
syms d theta a alpha real

%% Homogeneous transformation matrix (DON'T EDIT)
A = [cos(theta) -sin(theta)*cos(alpha)  sin(theta)*sin(alpha)   a*cos(theta); ...
     sin(theta)  cos(theta)*cos(alpha) -cos(theta)*sin(alpha)   a*sin(theta); ...
     0           sin(alpha)             cos(alpha)              d; ...
     0           0                      0                       1];

%% Load robot
robot = importrobot('../00001_PRR_zzy.urdf');
robot.DataFormat = 'row';

showdetails(robot);

eeName = robot.BodyNames{end};
baseName = robot.BaseName;

% Fixed link dimensions
a1 = 0.5;
a2 = 0.4;
a3 = 0.4;

%% Define the joint variables and link variables
syms q1 q2 q3 real
 
% % % % % % % % % % % % % % % % % % % % % 
%                                       %
% DH_table = | d | theta | a | alpha |  %
%                                       %
% % % % % % % % % % % % % % % % % % % % %

DH_table = [q1      0   0     0;
            0       q2  a2  -pi/2;
            0       q3  0     pi/2];
         
%% Direct Kinematics - Manual

fprintf('\n########### Direct Kinematics - Manual ###########\n');

A01 = subs(A,[d theta a alpha],DH_table(1,:));
A12 = subs(A,[d theta a alpha],DH_table(2,:));
A23 = subs(A,[d theta a alpha],DH_table(3,:));
    
% Transform the T matrix according to base frame and EE frame
Ab0 = [ 1 0 0  0;
        0 1 0  0;
        0 0 1 0.65; % 0.15 from base + 0.50 from link 1 length
        0 0 0  1]; % Trasformazione frame base

A3e = [ 1   0  0  0;
        0   1  0  0;
        0   0  1  0.40; % 0.40 from link 3 length
        0   0  0  1]; % Trasformazione frame end-effector

T03 = A01 * A12 * A23;
Tbe = Ab0 * T03 * A3e;

fprintf("Tbe matrix:\n");
disp(simplify(Tbe));

%% Direct Kinematics - Comparison

fprintf('\n########### Direct Kinematics - Built-In ###########\n');
fprintf('Base frame: %s\n', baseName);
fprintf('End-effector body: %s\n', eeName);

% Test joint configuration
q_check = [-0.2, 0.4, 0.6];

% Manual symbolic FK
T_manual_sym = subs(Tbe, ...
    [q1, q2, q3], ...
    [q_check]);

T_manual = double(T_manual_sym);

% MATLAB FK
T_builtin = getTransform(robot, q_check, eeName, baseName);

% Compare transformations
T_error = T_manual - T_builtin;
max_error = max(abs(T_error(:)));
position_error = norm(T_manual(1:3,4) - T_builtin(1:3,4));
rotation_error = norm(T_manual(1:3,1:3) - T_builtin(1:3,1:3), 'fro');

fprintf('\nManual FK:\n');
disp(T_manual);

fprintf('MATLAB FK:\n');
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

%% Inverse Kinematics - Manual

fprintf('\n########### Inverse Kinematics - Manual ###########\n');

syms px_d py_d pz_d real

r = sqrt(px_d^2 + py_d^2);

% Joint 2
q2Sol = simplify(atan2(py_d, px_d));

% Two branches for joint 3
sin_q3 = (r - a2) / 0.4;

q3Sol_1 = simplify(asin(sin_q3));
q3Sol_2 = simplify(pi - asin(sin_q3));

% Corresponding solutions for joint 1
q1Sol_1 = simplify(pz_d - 0.4*cos(q3Sol_1) - 0.65);
q1Sol_2 = simplify(pz_d - 0.4*cos(q3Sol_2) - 0.65);

fprintf('Inverse Kinematics solutions:\n');

fprintf('\nBranch 1:\n');
disp(q1Sol_1);
disp(q2Sol);
disp(q3Sol_1);

fprintf('\nBranch 2:\n');
disp(q1Sol_2);
disp(q2Sol);
disp(q3Sol_2);

%% Inverse Kinematics - Comparison

fprintf('\n########### Inverse Kinematics - Comparison ###########\n');

% Desired end-effector position
% Choose a reachable target in the same frame as Tbe
px_t = 0.45;
py_t = 0.20;
pz_t = 0.80;

% Substitute the desired position and link parameters
ik_subs_vars = [px_d, py_d, pz_d];
ik_subs_vals = [px_t, py_t, pz_t];

% Evaluate the analytical IK solutions already defined above
q2v = double(subs(q2Sol, ik_subs_vars, ik_subs_vals));

q3v = [ ...
    double(subs(q3Sol_1, ik_subs_vars, ik_subs_vals)), ...
    double(subs(q3Sol_2, ik_subs_vars, ik_subs_vals))];

q1v = [ ...
    double(subs(q1Sol_1, ik_subs_vars, ik_subs_vals)), ...
    double(subs(q1Sol_2, ik_subs_vars, ik_subs_vals))];

% Store the two analytical solutions
q_analytical = [q1v(:), repmat(q2v, 2, 1), q3v(:)];

fprintf('\nAnalytical IK solutions:\n');
disp(q_analytical);

% Target pose for MATLAB IK.
% Position is constrained; orientation is left unconstrained.
T_target = eye(4);
T_target(1:3,4) = [px_t; py_t; pz_t];

% Create MATLAB inverse kinematics solver
ik = inverseKinematics('RigidBodyTree', robot);

% Weights: [orientation XYZ, position XYZ]
weights = [0 0 0 1 1 1];

% Store built-in IK results
q_builtin_ik = zeros(2,3);

for k = 1:2

    % Use the corresponding analytical solution as initial guess
    q0 = q_analytical(k,:);

    [q_sol, solInfo] = ik(eeName, T_target, weights, q0);

    q_builtin_ik(k,:) = q_sol;

    % Evaluate the analytical solution with manual FK
    T_analytical = double(subs(Tbe, ...
        [q1, q2, q3], ...
        [q_analytical(k,:)]));

    % Evaluate the built-in solution with MATLAB FK
    T_solution = getTransform(robot, q_sol, eeName, baseName);

    % Position errors relative to the target
    err_analytical = norm(T_analytical(1:3,4) ...
        - [px_t; py_t; pz_t]);

    err_builtin = norm(T_solution(1:3,4) ...
        - [px_t; py_t; pz_t]);

    fprintf('\nBranch %d:\n', k);

    fprintf('Analytical solution:\n');
    fprintf('q1 = %.6f, q2 = %.6f, q3 = %.6f\n', ...
        q_analytical(k,1), ...
        q_analytical(k,2), ...
        q_analytical(k,3));

    fprintf('Built-in solution:\n');
    fprintf('q1 = %.6f, q2 = %.6f, q3 = %.6f\n', ...
        q_sol(1), q_sol(2), q_sol(3));

    fprintf('MATLAB IK status: %s\n', solInfo.Status);
    fprintf('Analytical FK position error: %.3e\n', err_analytical);
    fprintf('Built-in IK position error:  %.3e\n', err_builtin);

    if err_analytical < 1e-6
        fprintf('ANALYTICAL IK: MATCH\n');
    else
        fprintf('ANALYTICAL IK: DIFFER\n');
    end

    if err_builtin < 1e-6
        fprintf('BUILT-IN IK: MATCH\n');
    else
        fprintf('BUILT-IN IK: DIFFER\n');
    end
end

%% Visualization of the two IK branches

figure;

for k = 1:2
    % Joint configuration for the current analytical IK branch
    q1Target = q1v(k);
    q2Target = q2v;
    q3Target = q3v(k);

    % DataFormat = 'row': configuration is a numeric row vector
    config = [q1Target, q2Target, q3Target];

    % Plot robot configuration
    ax = subplot(1, 2, k);

    show(robot, config, 'Frames', 'on', 'Parent', ax);
    hold(ax, 'on');

    % Plot desired end-effector position
    scatter3(ax, px_t, py_t, pz_t, ...
        100, 'ro', 'filled', 'MarkerEdgeColor', 'k');

    grid(ax, 'on');
    axis(ax, 'equal');
    xlabel(ax, 'X [m]');
    ylabel(ax, 'Y [m]');
    zlabel(ax, 'Z [m]');
    view(ax, 3);

    title(ax, sprintf('Analytical IK Branch %d', k));
end

drawnow;

%% Geometric Jacobian [P R R] - Manual

fprintf('\n########### Geometric Jacobian [P R R] - Manual ###########\n');

s = @(w) [ 0    -w(3)  w(2);
           w(3)  0    -w(1);
          -w(2)  w(1)  0];

pe = Tbe(1:3, 4);

% Joint 1 - Prismatic
z0 = Ab0(1:3, 3);
p0 = Ab0(1:3, 4);

% Joint 2 - Revolut
Ab1 = Ab0*A01;
z1 = Ab1(1:3, 3);
p1 = Ab1(1:3, 4);

% Joint 3 - Revolut
Ab2 = Ab0*A01*A12;
z2 = Ab2(1:3, 3);
p2 = Ab2(1:3, 4);

% Assemble the geometric Jacobian columns
JGv = [z0, cross(z1, pe - p1), cross(z2, pe - p2)];
JGv_a = [z0, s(z1)*(pe - p1), s(z2)*(pe - p2)]; % Alternative with Skew matrix
JGw = [zeros(3,1), z1, z2];
JG = [JGv; JGw];

disp('Geometric Jacobian:');
disp(simplify(JG));

%% Geometric Jacobian [P R R] - Comparison

fprintf('\n########### Geometric Jacobian Comparison ###########\n');

% Evaluate manual Jacobian at the first analytical IK solution
qJacobian = [q1v(1), q2v, q3v(1)];

% Substitute joint values into the symbolic manual Jacobian
JG_manual_sym = subs(JG, ...
    [q1, q2, q3], qJacobian);

JG_manual = double(JG_manual_sym);

% MATLAB built-in geometric Jacobian
J_builtin = geometricJacobian(robot, qJacobian, eeName);

% MATLAB ordering: [Jw; Jv]
% Manual ordering: [Jv; Jw]
JG_manual_matlab_order = [JG_manual(4:6,:); JG_manual(1:3,:)];

% Compare Jacobians
J_error = JG_manual_matlab_order - J_builtin;

max_J_error = max(abs(J_error(:)));
linear_error = norm(JG_manual_matlab_order(4:6,:) ...
    - J_builtin(4:6,:), 'fro');
angular_error = norm(JG_manual_matlab_order(1:3,:) ...
    - J_builtin(1:3,:), 'fro');

fprintf('\nManual Jacobian:\n');
disp(JG_manual_matlab_order);

fprintf('Built-in Jacobian:\n');
disp(J_builtin);

fprintf('Maximum absolute error: %.3e\n', max_J_error);
fprintf('Linear Jacobian error:  %.3e\n', linear_error);
fprintf('Angular Jacobian error: %.3e\n', angular_error);

tol_J = 1e-9;

if max_J_error < tol_J
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: MATCH\n');
else
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: DIFFER\n');
end
