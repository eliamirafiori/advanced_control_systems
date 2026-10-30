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

%% Kinetic Energy - Manual

fprintf('\n\n');
fprintf('-----------------------------------------------\n');
fprintf('----------- Kinetic Energy - Manual -----------\n');
fprintf('-----------------------------------------------\n');

% Generalized coordinates and joint velocities.
% q1 is prismatic; q2 and q3 are revolute.
q = [q1; q2; q3];
syms dq1 dq2 dq3 real
dq = [dq1; dq2; dq3];

% Use the analytical link masses already defined above (in kg).
linkMass = [m1_a; m2_a; m3_a];

% Center-of-mass coordinates expressed in each link's own frame.
syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
comLocal = { [xc1; yc1; zc1], ...
             [xc2; yc2; zc2], ...
             [xc3; yc3; zc3] };

% Inertia tensors about each link's center of mass, expressed in that
% link's own frame. The tensors are symmetric but not assumed diagonal.
syms I1xx I1yy I1zz I1xy I1xz I1yz real
syms I2xx I2yy I2zz I2xy I2xz I2yz real
syms I3xx I3yy I3zz I3xy I3xz I3yz real

inertiaLocal = {
    [I1xx I1xy I1xz; I1xy I1yy I1yz; I1xz I1yz I1zz], ...
    [I2xx I2xy I2xz; I2xy I2yy I2yz; I2xz I2yz I2zz], ...
    [I3xx I3xy I3xz; I3xy I3yy I3yz; I3xz I3yz I3zz] };

% Transform each link frame into the base frame.
T_base_link = {Ab0*A01, Ab0*A01*A12, Ab0*A01*A12*A23};

% Joint axes and joint origins in the base frame, consistent with the
% geometric Jacobian derived in the previous block.
T_base_joint1 = Ab0;
T_base_joint2 = Ab0*A01;
T_base_joint3 = Ab0*A01*A12;

z0 = T_base_joint1(1:3,3);   % Prismatic joint axis

z1 = T_base_joint2(1:3,3);   % Revolute joint 2 axis
p1 = T_base_joint2(1:3,4);   % Origin of revolute joint 2

z2 = T_base_joint3(1:3,3);   % Revolute joint 3 axis
p2 = T_base_joint3(1:3,4);   % Origin of revolute joint 3

% Preallocate symbolic outputs for each link.
pCOM_base = cell(1,3);       % Center-of-mass position in base frame
Jv_COM = cell(1,3);          % Linear COM Jacobian
Jw_link = cell(1,3);         % Angular velocity Jacobian
inertiaCOM_base = cell(1,3); % COM inertia tensor rotated into base frame
linkKineticEnergy = sym(zeros(3,1));

for i = 1:3
    % Extract rotation and translation from base-to-link transform.
    R_base_link = T_base_link{i}(1:3,1:3);
    p_base_link = T_base_link{i}(1:3,4);

    % Express the center of mass in the base frame:
    % p_COM = R_base_link * p_COM_local + p_base_link.
    pCOM_base{i} = simplify(R_base_link*comLocal{i} + p_base_link);

    % Linear Jacobian, assembled using the same geometric-Jacobian rules
    % as the end-effector Jacobian above:
    %   prismatic joint:       Jv(:,1) = z0
    %   revolute joint j:      Jv(:,j) = z_j x (p_COM - p_j)
    % A link is affected only by joints preceding or attached to it.
    if i == 1
        Jv_COM{i} = [z0, zeros(3,1), zeros(3,1)];
        Jw_link{i} = [zeros(3,1), zeros(3,1), zeros(3,1)];
    elseif i == 2
        Jv_COM{i} = [z0, cross(z1, pCOM_base{i} - p1), zeros(3,1)];
        Jw_link{i} = [zeros(3,1), z1, zeros(3,1)];
    else
        Jv_COM{i} = [z0, ...
                     cross(z1, pCOM_base{i} - p1), ...
                     cross(z2, pCOM_base{i} - p2)];
        Jw_link{i} = [zeros(3,1), z1, z2];
    end

    % Rotate the COM inertia tensor from the local link frame into
    % the base frame: I_base = R * I_local * R^T.
    inertiaCOM_base{i} = simplify( ...
        R_base_link*inertiaLocal{i}*R_base_link.');

    % Obtain COM linear and link angular velocities from the Jacobians.
    vCOM = Jv_COM{i}*dq;
    omegaLink = Jw_link{i}*dq;

    % Rigid-body kinetic energy of this link:
    % T_i = 1/2*m_i*v_COM^T*v_COM + 1/2*omega^T*I_COM*omega.
    linkKineticEnergy(i) = simplify( ...
        1/2*linkMass(i)*(vCOM.'*vCOM) + ...
        1/2*(omegaLink.'*inertiaCOM_base{i}*omegaLink));
end

% Sum the kinetic energies of all three links.
totalKineticEnergy = simplify(sum(linkKineticEnergy));

% The joint-space inertia matrix B(q) satisfies:
% T = 1/2 * dq^T * B(q) * dq.
% Since T is quadratic in dq, its Hessian with respect to dq gives B.
manipulatorInertiaMatrix = simplify( ...
    jacobian(jacobian(totalKineticEnergy, dq).', dq));

% fprintf('\nCenter-of-mass positions in the base frame:\n');
% for i = 1:3
%     fprintf('pCOM_base{%d} =\n', i);
%     disp(pCOM_base{i});
% end
% 
% fprintf('\nLinear Jacobians of the link centers of mass:\n');
% for i = 1:3
%     fprintf('Jv_COM{%d} =\n', i);
%     disp(Jv_COM{i});
% end
% 
% fprintf('\nAngular Jacobians of the links:\n');
% for i = 1:3
%     fprintf('Jw_link{%d} =\n', i);
%     disp(Jw_link{i});
% end

% fprintf('\nKinetic energy of each link [J]:\n');
% disp(linkKineticEnergy);

fprintf('Total kinetic energy T(q,dq) [J]:\n');
disp(totalKineticEnergy);

fprintf('Joint-space inertia matrix B(q):\n');
disp(manipulatorInertiaMatrix);

% Checks: both expressions should simplify to zero.
fprintf('\nInertia matrix symmetry check (should be zero):\n');
disp(simplify(manipulatorInertiaMatrix - manipulatorInertiaMatrix.'));

fprintf('Kinetic-energy quadratic-form check (should be zero):\n');
disp(simplify(totalKineticEnergy - ...
    1/2*dq.'*manipulatorInertiaMatrix*dq));


%% Potential Energy - Manual

fprintf('\n\n');
fprintf('-------------------------------------------------\n');
fprintf('----------- Potential Energy - Manual -----------\n');
fprintf('-------------------------------------------------\n');

% Gravity expressed in the base frame.
% The base z-axis points upward, so gravity acts in the negative z direction.
syms g real
gravityBase = [0; 0; -g];

% Reuse the analytical link masses from the beginning of this script.
linkMass = [m1_a; m2_a; m3_a];

% Define the local center-of-mass coordinates again so this section can
% be understood and run independently of the kinetic-energy calculations.
syms xc1 yc1 zc1 xc2 yc2 zc2 xc3 yc3 zc3 real
comLocal = { [xc1; yc1; zc1], ...
             [xc2; yc2; zc2], ...
             [xc3; yc3; zc3] };

% Base-to-link transformations (repeated intentionally for clarity).
T_base_link = {Ab0*A01, Ab0*A01*A12, Ab0*A01*A12*A23};

pCOM_base_potential = cell(1,3);
linkPotentialEnergy = sym(zeros(3,1));

for i = 1:3
    % Convert each link's local COM position to the base frame.
    R_base_link = T_base_link{i}(1:3,1:3);
    p_base_link = T_base_link{i}(1:3,4);
    pCOM_base_potential{i} = simplify( ...
        R_base_link*comLocal{i} + p_base_link);

    % Gravitational potential energy relative to the base-frame origin:
    % U_i = -m_i * gravityBase^T * p_COM.
    linkPotentialEnergy(i) = simplify( ...
        -linkMass(i)*gravityBase.'*pCOM_base_potential{i});
end

% Total gravitational potential energy of the manipulator.
totalPotentialEnergy = simplify(sum(linkPotentialEnergy));

% fprintf('\nCenter-of-mass positions in the base frame:\n');
% for i = 1:3
%     fprintf('pCOM_base_potential{%d} =\n', i);
%     disp(pCOM_base_potential{i});
% end
% 
% fprintf('\nPotential energy of each link [J]:\n');
% disp(linkPotentialEnergy);

fprintf('Total potential energy U(q) [J]:\n');
disp(totalPotentialEnergy);

% Check: Potential-energy gradient should be zero

gravityGradient = simplify(jacobian(totalPotentialEnergy, q).');

% fprintf('\nGravity generalized-force vector G(q):\n');
% disp(gravityGradient);

fprintf('Potential-energy gradient check (should be zero):\n');
disp(simplify(jacobian(totalPotentialEnergy, q).' - gravityGradient));

%% Equations of Motion - Manual

fprintf('\n\n');
fprintf('---------------------------------------------------\n');
fprintf('----------- Equations of Motion - Manual ----------\n');
fprintf('---------------------------------------------------\n');

% Joint accelerations (Slide 25)
syms ddq1 ddq2 ddq3 real
ddq = [ddq1; ddq2; ddq3];

% Manipulator inertia matrix B(q) from kinetic energy (Slides 18, 25)
B = manipulatorInertiaMatrix;
n = length(q);

% -------------------------------------------------------------------------
% 1. Christoffel symbols of the first type and Matrix C(q, dq) (Slide 30)
%    c_ijk = 1/2 * ( d(b_ij)/d(q_k) + d(b_ik)/d(q_j) - d(b_jk)/d(q_i) )
%    c_ij  = sum_{k=1}^n c_ijk * dq_k
% -------------------------------------------------------------------------
c_ijk = sym(zeros(n, n, n));
for i = 1:n
    for j = 1:n
        for k = 1:n
            c_ijk(i,j,k) = simplify(1/2 * (diff(B(i,j), q(k)) + ...
                                           diff(B(i,k), q(j)) - ...
                                           diff(B(j,k), q(i))));
        end
    end
end

C = sym(zeros(n, n));
for i = 1:n
    for j = 1:n
        for k = 1:n
            C(i,j) = C(i,j) + c_ijk(i,j,k) * dq(k);
        end
    end
end
C = simplify(C);

% -------------------------------------------------------------------------
% 2. Centrifugal and Coriolis coefficients h_ijk (Slides 26-27)
%    h_ijk = d(b_ij)/d(q_k) - 1/2 * d(b_jk)/d(q_i)
%    h_ijj * dq_j^2 : centrifugal effect induced on Joint i by Joint j
%    h_ijk * dq_j * dq_k : Coriolis effect induced on Joint i by Joints j, k
% -------------------------------------------------------------------------
h_ijk = sym(zeros(n, n, n));
for i = 1:n
    for j = 1:n
        for k = 1:n
            h_ijk(i,j,k) = simplify(diff(B(i,j), q(k)) - 1/2 * diff(B(j,k), q(i)));
        end
    end
end

h_vec = sym(zeros(n, 1));
for i = 1:n
    for j = 1:n
        for k = 1:n
            h_vec(i) = h_vec(i) + h_ijk(i,j,k) * dq(k) * dq(j);
        end
    end
end
h_vec = simplify(h_vec);

% -------------------------------------------------------------------------
% 3. Gravity generalized force vector g(q) (Slides 20, 25-27)
%    Method A: Potential energy gradient: g(q) = (dU/dq)^T
%    Method B: Partial linear Jacobians: g_i(q) = - sum_j m_lj * g0^T * J_Pi^lj
% -------------------------------------------------------------------------
g_vec = gravityGradient; % Evaluated in Potential Energy section: jacobian(totalPotentialEnergy, q).'

% Gravity acceleration vector in base frame Sigma_0 (Slide 20)
g0 = gravityBase; 

g_from_J = sym(zeros(n, 1));
for i = 1:n
    for j = 1:n
        % Jv_COM{j}(:, i) is the i-th column of the j-th link COM Jacobian
        g_from_J(i) = g_from_J(i) - linkMass(j) * g0.' * Jv_COM{j}(:, i);
    end
end
g_from_J = simplify(g_from_J);

% -------------------------------------------------------------------------
% 4. Time derivative of inertia matrix B_dot and skew-symmetry check
%    N(q, dq) = B_dot - 2*C is skew-symmetric: N + N.' == 0
% -------------------------------------------------------------------------
B_dot = sym(zeros(n, n));
for i = 1:n
    for j = 1:n
        for k = 1:n
            B_dot(i,j) = B_dot(i,j) + diff(B(i,j), q(k)) * dq(k);
        end
    end
end
B_dot = simplify(B_dot);

N_skew = simplify(B_dot - 2*C);

% -------------------------------------------------------------------------
% 5. Standard Equations of Motion (Slide 31)
%    B(q)*ddq + C(q, dq)*dq + g(q) = tau
% -------------------------------------------------------------------------
tau = simplify(B*ddq + C*dq + g_vec);

% Optional: Extended model with friction and external wrench (Slide 31)
% B(q)*ddq + C(q, dq)*dq + Fv*dq + Fs*sign(dq) + g(q) = tau - J^T(q)*he
syms Fv1 Fv2 Fv3 real                 % Viscous friction coefficients
syms Fs1 Fs2 Fs3 real                 % Coulomb friction coefficients
syms fex fey fez mux muy muz real     % External wrench he = [fe; mue]

Fv = diag([Fv1, Fv2, Fv3]);
Fs = diag([Fs1, Fs2, Fs3]);
he = [fex; fey; fez; mux; muy; muz];

tau_friction = Fv*dq + Fs*sign(dq);
tau_external = JG.' * he;

tau_extended = simplify(B*ddq + C*dq + tau_friction + g_vec + tau_external);

% Display results
fprintf('Manipulator Inertia Matrix B(q):\n');
disp(B);

fprintf('Coriolis & Centrifugal Matrix C(q, dq):\n');
disp(C);

fprintf('Gravity Vector g(q):\n');
disp(g_vec);

fprintf('Equations of Motion: Generalized Torque Vector tau(q, dq, ddq):\n');
disp(tau);

% Analytical property checks (all must simplify to 0)
fprintf('\n--- Analytical Consistency Checks ---\n');

% Christoffel symmetry check: c_ijk = c_ikj
sym_c_err = sym(0);
for i = 1:n
    for j = 1:n
        for k = 1:n
            sym_c_err = sym_c_err + simplify(abs(c_ijk(i,j,k) - c_ijk(i,k,j)));
        end
    end
end
fprintf('Christoffel symmetry check (c_ijk - c_ikj == 0): %s\n', char(sym_c_err == 0));

% Centrifugal diagonal nullity check: h_iii == 0 (Slide 27)
h_iii_err = sym(0);
for i = 1:n
    h_iii_err = h_iii_err + abs(h_ijk(i,i,i));
end
fprintf('Centrifugal diagonal nullity check (h_iii == 0): %s\n', char(h_iii_err == 0));

% Equivalence between C(q, dq)*dq and sum(sum(h_ijk*dq_k*dq_j))
fprintf('C(q,dq)*dq vs h_vec equivalence (should be zero):\n');
disp(simplify(C*dq - h_vec));

% Gravity derivation equivalence
fprintf('Gravity vector equivalence check (g_from_U - g_from_J == 0):\n');
disp(simplify(g_vec - g_from_J));

% Skew-symmetry of B_dot - 2C (should be zero matrix)
fprintf('Skew-symmetry check (B_dot - 2*C) + (B_dot - 2*C).'' (should be zero):\n');
disp(simplify(N_skew + N_skew.'));


%% Equations of Motion - Comparison

fprintf('\n\n');
fprintf('-------------------------------------------------------\n');
fprintf('----------- Equations of Motion - Comparison ----------\n');
fprintf('-------------------------------------------------------\n');

% Ensure consistent gravity specification in the MATLAB rigid body tree
g_val = 9.81;
robot.Gravity = [0, 0, -g_val];

% Test configuration, velocities, and accelerations
q_check   = [-0.2,  0.4,  0.6];
dq_check  = [ 0.1, -0.3,  0.5];
ddq_check = [ 0.05, 0.2, -0.1];

% Extract non-fixed moving bodies corresponding to joints 1, 2, and 3
movingBodies = {};
for i = 1:numel(robot.Bodies)
    if ~strcmp(robot.Bodies{i}.Joint.Type, 'fixed')
        movingBodies{end+1} = robot.Bodies{i};
    end
end

% Extract actual physical parameters directly from the URDF model
m_val = zeros(1, 3);
com_val = cell(1, 3);
I_C_val = cell(1, 3);

for i = 1:3
    b = movingBodies{i};
    m_val(i) = b.Mass;
    c = b.CenterOfMass(:); % 3x1 vector [xc; yc; zc]
    com_val{i} = c;
    
    % MATLAB b.Inertia is [Ixx Iyy Izz Iyz Ixz Ixy] relative to body origin
    I_vec = b.Inertia;
    I_origin = [I_vec(1), I_vec(6), I_vec(5);
                I_vec(6), I_vec(2), I_vec(4);
                I_vec(5), I_vec(4), I_vec(3)];
            
    % Convert to inertia tensor about the Center of Mass via Steiner's theorem:
    % I_C = I_origin - m * ( (c'*c)*eye(3) - c*c' ) (Slides 14-15)
    I_C = I_origin - m_val(i) * ( (c'*c)*eye(3) - (c * c') );
    I_C_val{i} = I_C;
end

% Numerical substitution list for symbolic expressions
dyn_syms = [q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3, ...
            a1, a2, a3, g, ...
            xc1, yc1, zc1, xc2, yc2, zc2, xc3, yc3, zc3, ...
            I1xx, I1yy, I1zz, I1xy, I1xz, I1yz, ...
            I2xx, I2yy, I2zz, I2xy, I2xz, I2yz, ...
            I3xx, I3yy, I3zz, I3xy, I3xz, I3yz];

dyn_vals = [q_check, dq_check, ddq_check, ...
            a1_a, a2_a, a3_a, g_val, ...
            com_val{1}(1), com_val{1}(2), com_val{1}(3), ...
            com_val{2}(1), com_val{2}(2), com_val{2}(3), ...
            com_val{3}(1), com_val{3}(2), com_val{3}(3), ...
            I_C_val{1}(1,1), I_C_val{1}(2,2), I_C_val{1}(3,3), ...
            I_C_val{1}(1,2), I_C_val{1}(1,3), I_C_val{1}(2,3), ...
            I_C_val{2}(1,1), I_C_val{2}(2,2), I_C_val{2}(3,3), ...
            I_C_val{2}(1,2), I_C_val{2}(1,3), I_C_val{2}(2,3), ...
            I_C_val{3}(1,1), I_C_val{3}(2,2), I_C_val{3}(3,3), ...
            I_C_val{3}(1,2), I_C_val{3}(1,3), I_C_val{3}(2,3)];

% Evaluate manual analytical dynamic expressions
B_manual_num   = double(subs(B, dyn_syms, dyn_vals));
C_manual_num   = double(subs(C, dyn_syms, dyn_vals));
g_manual_num   = double(subs(g_vec, dyn_syms, dyn_vals));
tau_manual_num = double(subs(tau, dyn_syms, dyn_vals));

% Built-in MATLAB Robotics System Toolbox dynamic evaluations
B_builtin    = massMatrix(robot, q_check);
C_dq_builtin = velocityProduct(robot, q_check, dq_check).';
g_builtin    = gravityTorque(robot, q_check).';
tau_builtin  = inverseDynamics(robot, q_check, dq_check, ddq_check).';

% Cross-Check 1: Inertia Matrix B(q)
err_B = norm(B_manual_num - B_builtin, 'fro');
fprintf('Inertia Matrix B(q) Frobenius Error: %.3e\n', err_B);
if err_B < 1e-6
    fprintf('INERTIA MATRIX B(q) CROSS-CHECK: MATCH\n');
else
    fprintf('INERTIA MATRIX B(q) CROSS-CHECK: DIFFER\n');
end

% Cross-Check 2: Velocity Product Vector C(q, dq)*dq
err_C = norm(C_manual_num * dq_check.' - C_dq_builtin);
fprintf('Velocity Product C*dq Error:         %.3e\n', err_C);
if err_C < 1e-6
    fprintf('VELOCITY PRODUCT C*dq CROSS-CHECK: MATCH\n');
else
    fprintf('VELOCITY PRODUCT C*dq CROSS-CHECK: DIFFER\n');
end

% Cross-Check 3: Gravity Torque Vector g(q)
err_g = norm(g_manual_num - g_builtin);
fprintf('Gravity Vector g(q) Error:          %.3e\n', err_g);
if err_g < 1e-6
    fprintf('GRAVITY TORQUE g(q) CROSS-CHECK:   MATCH\n');
else
    fprintf('GRAVITY TORQUE g(q) CROSS-CHECK:   DIFFER\n');
end

% Cross-Check 4: Total Inverse Dynamics Torque Vector tau
err_tau = norm(tau_manual_num - tau_builtin);
fprintf('Inverse Dynamics Torque tau Error:   %.3e\n', err_tau);
if err_tau < 1e-6
    fprintf('INVERSE DYNAMICS tau CROSS-CHECK:  MATCH\n');
else
    fprintf('INVERSE DYNAMICS tau CROSS-CHECK:  DIFFER\n');
end
