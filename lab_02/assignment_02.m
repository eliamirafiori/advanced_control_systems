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

%Change gravity

% Fixed link analytical dimensions
a1_a = 0.5;
a2_a = 0.4;
a3_a = 0.4;

% Fixed link analytical weights
m1_a = 1;
m2_a = 1;
m3_a = 1;

%% Define the joint variables and link variables
syms q1 q2 q3 a1 a2 a3 real
 
% % % % % % % % % % % % % % % % % % % % % 
%                                       %
% DH_table = | d | theta | a | alpha |  %
%                                       %
% % % % % % % % % % % % % % % % % % % % %

DH_table = [q1  0   0    0;
            0   q2  a2  -pi/2;
            0   q3  0    pi/2];
         
%% Direct Kinematics - Manual

fprintf('\n\n');
fprintf('--------------------------------------------------\n');
fprintf('----------- Direct Kinematics - Manual -----------\n');
fprintf('--------------------------------------------------\n');

A01 = subs(A,[d theta a alpha],DH_table(1,:));
A12 = subs(A,[d theta a alpha],DH_table(2,:));
A23 = subs(A,[d theta a alpha],DH_table(3,:));
    
% Transform the T matrix according to base frame and EE frame
Ab0 = [ 1 0 0  0;
        0 1 0  0;
        0 0 1 0.15 + a1_a; % 0.15 from base + 0.50 from link 1 length
        0 0 0  1]; % Base frame transformation

A3e = [ 1   0  0  0;
        0   1  0  0;
        0   0  1  a3_a; % 0.40 from link 3 length
        0   0  0  1]; % EE frame transformation

T03 = A01 * A12 * A23;
Tbe = Ab0 * T03 * A3e;

fprintf("Tbe matrix:\n");
disp(simplify(Tbe));

%% Direct Kinematics - Comparison

fprintf('\n\n');
fprintf('------------------------------------------------------\n');
fprintf('----------- Direct Kinematics - Comparison -----------\n');
fprintf('------------------------------------------------------\n');
fprintf('Base frame: %s\n', baseName);
fprintf('End-effector body: %s\n', eeName);

% Test joint configuration
q_check = [-0.2, 0.4, 0.6];

% Manual symbolic FK
T_manual_sym = subs(Tbe, ...
    [q1, q2, q3, a1, a2, a3], ...
    [q_check, a1_a, a2_a, a3_a]);

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

fprintf('\n\n');
fprintf('---------------------------------------------------\n');
fprintf('----------- Inverse Kinematics - Manual -----------\n');
fprintf('---------------------------------------------------\n');

syms px py pz real

r = sqrt(px^2 + py^2);

% Joint 2
q2Sol = simplify(atan2(py, px));

% Two branches for joint 3
sin_q3 = (r - a2) / 0.4;

q3Sol_1 = simplify(asin(sin_q3));
q3Sol_2 = simplify(pi - asin(sin_q3));

% Corresponding solutions for joint 1
q1Sol_1 = simplify(pz - 0.4*cos(q3Sol_1) - 0.65);
q1Sol_2 = simplify(pz - 0.4*cos(q3Sol_2) - 0.65);

fprintf('\nBranch 1:\n');
disp(q1Sol_1);
disp(q2Sol);
disp(q3Sol_1);

fprintf('\nBranch 2:\n');
disp(q1Sol_2);
disp(q2Sol);
disp(q3Sol_2);

%% Inverse Kinematics - Comparison

fprintf('\n\n');
fprintf('-------------------------------------------------------\n');
fprintf('----------- Inverse Kinematics - Comparison -----------\n');
fprintf('-------------------------------------------------------\n');

% Target end-effector position
% Choose a reachable target in the same frame as Tbe
px_t = 0.45;
py_t = 0.20;
pz_t = 0.80;

% Substitute the desired position and link parameters
ik_syms = [px, py, pz, a1, a2, a3];
ik_vals = [px_t, py_t, pz_t, a1_a, a2_a, a3_a];

% Evaluate the analytical IK solutions already defined above
q2v = double(subs(q2Sol, ik_syms, ik_vals));

q3v = [ ...
    double(subs(q3Sol_1, ik_syms, ik_vals)), ...
    double(subs(q3Sol_2, ik_syms, ik_vals))];

q1v = [ ...
    double(subs(q1Sol_1, ik_syms, ik_vals)), ...
    double(subs(q1Sol_2, ik_syms, ik_vals))];

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
weights = [0 0 0 1 1 1]; % ignore orientation error entirely, care only about position error

% Store built-in IK results
q_builtin_ik = zeros(2,3);

% Get joint limits once before the loop
jointLimits = zeros(3, 2);
j = 0;
for i = 1:numel(robot.Bodies)
    joint = robot.Bodies{i}.Joint;
    if ~strcmp(joint.Type, "fixed")
        j = j + 1;
        jointLimits(j,:) = joint.PositionLimits;
    end
end

for k = 1:2

    % Use the corresponding analytical solution as initial guess
    q0 = q_analytical(k,:);

    withinLimits = all(q0 >= jointLimits(:,1).' & q0 <= jointLimits(:,2).');
    if ~withinLimits
        fprintf('Branch %d skipped: initial guess violates joint limits.\n', k);
        continue
    end

    [q_sol, solInfo] = ik(eeName, T_target, weights, q0);

    q_builtin_ik(k,:) = q_sol;

    % Evaluate the analytical solution with manual FK
    T_analytical = double(subs(Tbe, ...
        [q1, q2, q3, a1, a2, a3], ...
        [q_analytical(k,:), a1_a, a2_a, a3_a]));

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

    % fprintf('MATLAB IK status: %s\n', solInfo.Status);
    % fprintf('Analytical FK position error: %.3e\n', err_analytical);
    % fprintf('Built-in IK position error:  %.3e\n', err_builtin);

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

fprintf('\n\n');
fprintf('-----------------------------------------------------------\n');
fprintf('----------- Geometric Jacobian [P R R] - Manual -----------\n');
fprintf('-----------------------------------------------------------\n');

s = @(w) [ 0    -w(3)  w(2);
           w(3)  0    -w(1);
          -w(2)  w(1)  0];

pe = Tbe(1:3, 4); % End-Effector position

% Joint 1 - Prismatic
z0 = Ab0(1:3, 3);
p0 = Ab0(1:3, 4);

% Joint 2 - Revolute
Ab1 = Ab0*A01;
z1 = Ab1(1:3, 3);
p1 = Ab1(1:3, 4);

% Joint 3 - Revolute
Ab2 = Ab0*A01*A12;
z2 = Ab2(1:3, 3);
p2 = Ab2(1:3, 4);

% Assemble the geometric Jacobian columns
JGv = [z0, cross(z1, pe - p1), cross(z2, pe - p2)]; % Linear velocities
% JGv_a = [z0, s(z1)*(pe - p1), s(z2)*(pe - p2)]; % Alternative with Skew matrix
JGw = [zeros(3,1), z1, z2]; % Angular velocities
JG = [JGv; JGw]; % Complete Geometric Jacobian

disp('Geometric Jacobian:\n');
disp(simplify(JG));

%% Geometric Jacobian [P R R] - Comparison

fprintf('\n\n');
fprintf('-------------------------------------------------------------\n');
fprintf('----------- Geometric Jacobian [P R R] - Built-In -----------\n');
fprintf('-------------------------------------------------------------\n');

% Evaluate manual Jacobian at the first analytical IK solution
qJacobian = [q1v(1), q2v, q3v(1)];

% Substitute joint values into the symbolic manual Jacobian
JG_manual_sym = subs(JG, ...
    [q1, q2, q3, a1, a2, a3], ...
    [qJacobian, a1_a, a2_a, a3_a]);

JG_manual = double(JG_manual_sym);

% MATLAB built-in geometric Jacobian
JG_builtin = geometricJacobian(robot, qJacobian, eeName);
JG_builtin = round(JG_builtin, 4);

% MATLAB ordering: [Jw; Jv]
% Manual ordering: [Jv; Jw]
JG_manual_matlab_order = [JG_manual(4:6,:); JG_manual(1:3,:)];

% Compare Jacobians
JG_error = JG_manual_matlab_order - JG_builtin;

max_JG_error = max(abs(JG_error(:)));
linear_error = norm(JG_manual_matlab_order(4:6,:) ...
    - JG_builtin(4:6,:), 'fro');
angular_error = norm(JG_manual_matlab_order(1:3,:) ...
    - JG_builtin(1:3,:), 'fro');

fprintf('\nManual Jacobian:\n');
disp(JG_manual_matlab_order);

fprintf('Built-in Jacobian:\n');
disp(JG_builtin);

% fprintf('Maximum absolute error: %.3e\n', max_JG_error);
% fprintf('Linear Jacobian error:  %.3e\n', linear_error);
% fprintf('Angular Jacobian error: %.3e\n', angular_error);

tol_JG = 1e-9;

if max_JG_error < tol_JG
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: MATCH\n');
else
    fprintf('GEOMETRIC JACOBIAN CROSS-CHECK: DIFFER\n');
end

%Check why differ


%% Analytical Jacobian - Manual (ZYZ Euler angles)

fprintf('\n\n');
fprintf('----------------------------------------------------------\n');
fprintf('----------- Analytical Jacobian - Manual (ZYZ) -----------\n');
fprintf('----------------------------------------------------------\n');

pe = Tbe(1:3, 4); % End-Effector position
Rbe = Tbe(1:3, 1:3); % End-Effector orientation

% ZYZ Euler angles
phi_euler = atan2(Rbe(2,3), Rbe(1,3));

theta_euler = atan2( ...
    sqrt(Rbe(1,3)^2 + Rbe(2,3)^2), ...
    Rbe(3,3));

psi_euler = atan2(Rbe(3,2), -Rbe(3,1));

euler_angles = [phi_euler; theta_euler; psi_euler];

% Translational analytical Jacobian
JAv = [diff(pe, q1), diff(pe, q2), diff(pe, q3)];
% JAv_a = jacobian(pe, [q1, q2, q3]); % Alternative way to compute derivatives

% Euler-angle-rate Jacobian
JAw = [diff(euler_angles, q1), diff(euler_angles, q2), diff(euler_angles, q3)];
% JAw_a = jacobian(euler_angles, [q1, q2, q3]); % Alternative way to compute derivatives

% Complete analytical Jacobian
JA = simplify([JAv; JAw]); % Complete Analytical Jacobian

disp('ZYZ Euler angles:');
disp(simplify(euler_angles));

disp('Analytical Jacobian (ZYZ):');
disp(JA);

%% Analytical Jacobian [P R R] - Comparison (ZYZ Euler Angles)

fprintf('\n\n');
fprintf('------------------------------------------------------------\n');
fprintf('----------- Analytical Jacobian - Built-In (ZYZ) -----------\n');
fprintf('------------------------------------------------------------\n');

% Evaluate manual Jacobian at the first analytical IK solution
qJacobian = [q1v(1), q2v, q3v(1)];

% Substitute joint values into the symbolic analytical Jacobian
JA_manual_sym = subs(JA, ...
    [q1, q2, q3, a1, a2, a3], ...
    [qJacobian, a1_a, a2_a, a3_a]);
JA_manual = double(JA_manual_sym);

% MATLAB built-in geometric Jacobian
J_builtin = geometricJacobian(robot, qJacobian, eeName);

% MATLAB ordering: [Jw; Jv]
Jw_builtin = J_builtin(1:3,:);
Jv_builtin = J_builtin(4:6,:);

% Evaluate ZYZ Euler angles at this configuration
R_test = double(subs(Rbe, ...
    [q1, q2, q3, a1, a2, a3], ...
    [qJacobian, a1_a, a2_a, a3_a]));

phi_test = atan2(R_test(2,3), R_test(1,3));
theta_test = atan2(sqrt(R_test(1,3)^2 + R_test(2,3)^2), ...
    R_test(3,3));

% Convert spatial angular velocity to ZYZ Euler angle rates
% omega = E * [phi_dot; theta_dot; psi_dot]

E = [0, -sin(phi_test), cos(phi_test)*sin(theta_test);
    0,  cos(phi_test), sin(phi_test)*sin(theta_test);
    1,  0,             cos(theta_test)];

J_euler = E \ JGw;

% Reorder built-in Jacobian to [Jv; Euler angle rates]
JA_builtin = [JGv; J_euler];
JA_builtin = subs(JA_builtin, ...
    [q1, q2, q3, a1, a2, a3], ...
    [qJacobian, a1_a, a2_a, a3_a]);
JA_builtin = double(round(JA_builtin, 4));

% Compare Jacobians
J_error = JA_manual - JA_builtin;

max_J_error = max(abs(J_error(:)));
linear_error = norm(JA_manual(1:3,:) - JA_builtin(1:3,:), 'fro');
angular_error = norm(JA_manual(4:6,:) - JA_builtin(4:6,:), 'fro');

fprintf('\nManual Analytical Jacobian (ZYZ):\n');
disp(JA_manual);

fprintf('Built-in Jacobian converted to ZYZ Euler rates:\n');
disp(JA_builtin);

fprintf('Maximum absolute error: %.3e\n', max_J_error);
fprintf('Linear Jacobian error:  %.3e\n', linear_error);
fprintf('Euler-rate Jacobian error: %.3e\n', angular_error);

tol_J = 1e-9;

if max_J_error < tol_J
    fprintf('ANALYTICAL JACOBIAN CROSS-CHECK: MATCH\n');
else
    fprintf('ANALYTICAL JACOBIAN CROSS-CHECK: DIFFER\n');
end

%% ASSIGNMENT #2: KINETIC AND POTENTIAL ENERGY OF A RIGID BODY

fprintf('\n\n');
fprintf('-----------------------------------------------\n');
fprintf('----------- Kinetic and Potential Energy ------\n');
fprintf('-----------------------------------------------\n');

%% 1. Definition of Dynamic Variables and Gravity Vector
% Generalized coordinates and velocities
syms dq1 dq2 dq3 real
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% Analytical link masses (Slides 4, 20)
m_l1 = m1_a;
m_l2 = m2_a;
m_l3 = m3_a;

% Acceleration of gravity vector in base frame Sigma_0 (Slides 20, 25)
% Base z-axis points upward, gravity acts downwards along -z0
syms g real
g0 = [0; 0; -g];

% Local Center-of-Mass coordinates expressed in link frame Sigma_i (Slide 21)
syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
p_l1_1 = [xc1; yc1; zc1]; % CoM of link 1 in Frame 1
p_l2_2 = [xc2; yc2; zc2]; % CoM of link 2 in Frame 2
p_l3_3 = [xc3; yc3; zc3]; % CoM of link 3 in Frame 3

% Link inertia tensors relative to CoM expressed in link frame Sigma_i: I_li^i
% Symmetric matrices relative to central axes (Slides 10-15)
syms I1xx I1yy I1zz I1xy I1xz I1yz real
syms I2xx I2yy I2zz I2xy I2xz I2yz real
syms I3xx I3yy I3zz I3xy I3xz I3yz real

I_l1_1 = [I1xx I1xy I1xz;
          I1xy I1yy I1yz;
          I1xz I1yz I1zz];

I_l2_2 = [I2xx I2xy I2xz;
          I2xy I2yy I2yz;
          I2xz I2yz I2zz];

I_l3_3 = [I3xx I3xy I3xz;
          I3xy I3yy I3yz;
          I3xz I3yz I3zz];


%% 2. Steiner's Theorem (Parallel Axis Theorem) Formulation (Slides 14-15)
% The inertia tensor I with respect to a frame obtained by translating the CoM 
% frame Sigma_C by vector r is:
% I = I_C + m * S^T(r)*S(r) = I_C + m * (r^T*r*eye(3) - r*r^T)

% Skew-symmetric operator function S(r)
s = @(w) [  0    -w(3)  w(2);
           w(3)   0    -w(1);
          -w(2)  w(1)   0  ];

% Link 1: Shift from CoM to Frame 1 origin by vector r1 = p_l1_1
r1 = p_l1_1;
I_O1_1 = simplify(I_l1_1 + m_l1 * ((r1.' * r1)*eye(3) - (r1 * r1.')));

% Link 2: Shift from CoM to Frame 2 origin by vector r2 = p_l2_2
r2 = p_l2_2;
I_O2_2 = simplify(I_l2_2 + m_l2 * ((r2.' * r2)*eye(3) - (r2 * r2.')));

% Link 3: Shift from CoM to Frame 3 origin by vector r3 = p_l3_3
r3 = p_l3_3;
I_O3_3 = simplify(I_l3_3 + m_l3 * ((r3.' * r3)*eye(3) - (r3 * r3.')));

fprintf('\n--- Steiner Inertia Tensors about Link Frame Origins (Slide 14) ---\n');
fprintf('Link 1 Inertia about Frame 1 Origin:\n'); disp(I_O1_1);
fprintf('Link 2 Inertia about Frame 2 Origin:\n'); disp(I_O2_2);
fprintf('Link 3 Inertia about Frame 3 Origin:\n'); disp(I_O3_3);


%% 3. Forward Kinematics for Link Frames and CoM Positions (Slide 21)
% Base-to-link homogeneous transformations: T_0^i = Ab0 * A01 * ... * A(i-1)i
T_b1 = Ab0 * A01;
T_b2 = Ab0 * A01 * A12;
T_b3 = Ab0 * A01 * A12 * A23;

% Rotation matrices from base frame Sigma_0 to link frames Sigma_i
R1 = T_b1(1:3, 1:3);
R2 = T_b2(1:3, 1:3);
R3 = T_b3(1:3, 1:3);

% Origins of link frames in base frame
p_O1 = T_b1(1:3, 4);
p_O2 = T_b2(1:3, 4);
p_O3 = T_b3(1:3, 4);

% Joint axes and origins in base frame for Jacobian assembly (Slides 16-17)
z0 = Ab0(1:3, 3);        p0 = Ab0(1:3, 4);        % Joint 1 (Prismatic)
z1 = T_b1(1:3, 3);       p1 = T_b1(1:3, 4);       % Joint 2 (Revolute)
z2 = T_b2(1:3, 3);       p2 = T_b2(1:3, 4);       % Joint 3 (Revolute)

% Positions of Link Centers of Mass in base frame Sigma_0 (Slide 21):
% p_li = R_i * p_li^i + p_Oi
p_l1 = simplify(R1 * p_l1_1 + p_O1);
p_l2 = simplify(R2 * p_l2_2 + p_O2);
p_l3 = simplify(R3 * p_l3_3 + p_O3);

fprintf('\n--- Center of Mass Positions in Base Frame Sigma_0 (Slide 21) ---\n');
fprintf('p_l1 =\n'); disp(p_l1);
fprintf('p_l2 =\n'); disp(p_l2);
fprintf('p_l3 =\n'); disp(p_l3);


%% 4. Partial Jacobians for Link Centers of Mass (Slides 16-17)
% For link i, the partial Jacobians relate joint velocities to Cartesian velocities:
% dp_li = J_P_li * dq
% omega_i = J_O_li * dq

% Link 1 (Affected only by Joint 1 - Prismatic):
% Column 1: j_P1 = z0, j_O1 = 0
% Columns 2, 3: zero
J_P_l1 = [z0, zeros(3,1), zeros(3,1)];
J_O_l1 = [zeros(3,1), zeros(3,1), zeros(3,1)];

% Link 2 (Affected by Joint 1 - Prismatic, and Joint 2 - Revolute):
% Column 1: j_P1 = z0,                      j_O1 = 0
% Column 2: j_P2 = z1 x (p_l2 - p1),        j_O2 = z1
% Column 3: zero
J_P_l2 = [z0, cross(z1, p_l2 - p1), zeros(3,1)];
J_O_l2 = [zeros(3,1), z1, zeros(3,1)];

% Link 3 (Affected by Joints 1, 2, and 3):
% Column 1: j_P1 = z0,                      j_O1 = 0
% Column 2: j_P2 = z1 x (p_l3 - p1),        j_O2 = z1
% Column 3: j_P3 = z2 x (p_l3 - p2),        j_O3 = z2
J_P_l3 = [z0, cross(z1, p_l3 - p1), cross(z2, p_l3 - p2)];
J_O_l3 = [zeros(3,1), z1, z2];


%% 5. Kinetic Energy - Manual (Slides 6-18)
% König's Theorem (Slides 12, 16):
% T_li = 1/2 * m_li * dp_li^T * dp_li + 1/2 * omega_i^T * I_li * omega_i
% where I_li = R_i * I_li^i * R_i^T is the CoM inertia tensor in base frame.

% CoM Inertia tensors rotated into base frame Sigma_0 (Slides 10, 15):
I_l1_0 = simplify(R1 * I_l1_1 * R1.');
I_l2_0 = simplify(R2 * I_l2_2 * R2.');
I_l3_0 = simplify(R3 * I_l3_3 * R3.');

% Linear velocities of the CoMs
dp_l1 = J_P_l1 * dq;
dp_l2 = J_P_l2 * dq;
dp_l3 = J_P_l3 * dq;

% Angular velocities of the links
omega_1 = J_O_l1 * dq;
omega_2 = J_O_l2 * dq;
omega_3 = J_O_l3 * dq;

% Translational kinetic energy of each link (Slide 7)
T_l1_trans = simplify(1/2 * m_l1 * (dp_l1.' * dp_l1));
T_l2_trans = simplify(1/2 * m_l2 * (dp_l2.' * dp_l2));
T_l3_trans = simplify(1/2 * m_l3 * (dp_l3.' * dp_l3));

% Rotational kinetic energy of each link (Slide 10)
T_l1_rot = simplify(1/2 * (omega_1.' * I_l1_0 * omega_1));
T_l2_rot = simplify(1/2 * (omega_2.' * I_l2_0 * omega_2));
T_l3_rot = simplify(1/2 * (omega_3.' * I_l3_0 * omega_3));

% Total kinetic energy of each link (Slide 16)
T_l1 = simplify(T_l1_trans + T_l1_rot);
T_l2 = simplify(T_l2_trans + T_l2_rot);
T_l3 = simplify(T_l3_trans + T_l3_rot);

% Total kinetic energy of the manipulator (Slide 18): T = sum_i T_li
T_total = simplify(T_l1 + T_l2 + T_l3);

% Manipulator Inertia Matrix B(q) (Slides 18-19, 22-23):
% T = 1/2 * dq^T * B(q) * dq
% B(q) = sum_i [ m_li * (J_P^li)^T * J_P^li + (J_O^li)^T * R_i * I_li^i * R_i^T * J_O^li ]
B1 = simplify(m_l1 * (J_P_l1.' * J_P_l1) + J_O_l1.' * I_l1_0 * J_O_l1);
B2 = simplify(m_l2 * (J_P_l2.' * J_P_l2) + J_O_l2.' * I_l2_0 * J_O_l2);
B3 = simplify(m_l3 * (J_P_l3.' * J_P_l3) + J_O_l3.' * I_l3_0 * J_O_l3);

B = simplify(B1 + B2 + B3);

% Equivalent calculation from Hessian of kinetic energy:
B_hessian = simplify(jacobian(jacobian(T_total, dq).', dq));

fprintf('\n--- Kinetic Energy & Inertia Matrix B(q) (Slide 18) ---\n');
fprintf('Total Kinetic Energy T(q, dq) [J]:\n'); disp(T_total);
fprintf('Manipulator Inertia Matrix B(q):\n'); disp(B);

% Checks on Kinetic Energy and Inertia Matrix:
% Check 1: Symmetry B(q) = B^T(q) (Slide 19)
sym_check = simplify(B - B.');
fprintf('Inertia matrix symmetry check B - B^T (should be zero):\n');
disp(sym_check);

% Check 2: Quadratic form equivalence T = 1/2*dq'*B*dq (Slide 19)
quad_check = simplify(T_total - 1/2*dq.'*B*dq);
fprintf('Kinetic energy quadratic form check (should be zero):\n');
disp(quad_check);

% Check 3: Jacobian summation vs Hessian derivation
hess_check = simplify(B - B_hessian);
fprintf('Hessian vs Jacobian summation check (should be zero):\n');
disp(hess_check);


%% 6. Potential Energy - Manual (Slides 20-21)
% Gravitational potential energy of link i (Slide 20):
% U_li = - m_li * g0^T * p_li
U_l1 = simplify(- m_l1 * (g0.' * p_l1));
U_l2 = simplify(- m_l2 * (g0.' * p_l2));
U_l3 = simplify(- m_l3 * (g0.' * p_l3));

% Total gravitational potential energy (Slide 20): U = sum_i U_li
U_total = simplify(U_l1 + U_l2 + U_l3);

% Gravity generalized force vector g(q) = (dU/dq)^T (Slide 25)
g1 = simplify(diff(U_total, q1));
g2 = simplify(diff(U_total, q2));
g3 = simplify(diff(U_total, q3));
g_vec = [g1; g2; g3];

% Alternative calculation from partial linear Jacobians (Slide 26):
% g_i(q) = - sum_j m_lj * g0^T * j_Pi^lj(q)
g1_J = simplify(- (m_l1*g0.'*J_P_l1(:,1) + m_l2*g0.'*J_P_l2(:,1) + m_l3*g0.'*J_P_l3(:,1)));
g2_J = simplify(- (m_l1*g0.'*J_P_l1(:,2) + m_l2*g0.'*J_P_l2(:,2) + m_l3*g0.'*J_P_l3(:,2)));
g3_J = simplify(- (m_l1*g0.'*J_P_l1(:,3) + m_l2*g0.'*J_P_l2(:,3) + m_l3*g0.'*J_P_l3(:,3)));
g_vec_J = [g1_J; g2_J; g3_J];

fprintf('\n--- Potential Energy & Gravity Vector g(q) (Slides 20, 25) ---\n');
fprintf('Potential Energy of Link 1 U_l1 [J]:\n'); disp(U_l1);
fprintf('Potential Energy of Link 2 U_l2 [J]:\n'); disp(U_l2);
fprintf('Potential Energy of Link 3 U_l3 [J]:\n'); disp(U_l3);
fprintf('Total Potential Energy U(q) [J]:\n');     disp(U_total);
fprintf('Gravity Generalized Force Vector g(q):\n'); disp(g_vec);

% Check on potential energy derivations:
g_equiv_check = simplify(g_vec - g_vec_J);
fprintf('Gravity vector equivalence check (g_from_U - g_from_J == 0):\n');
disp(g_equiv_check);
