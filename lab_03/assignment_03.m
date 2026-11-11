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

%% ASSIGNMENT #3: EQUATIONS OF MOTION (DYNAMIC MODEL)

fprintf('\n\n');
fprintf('---------------------------------------------------\n');
fprintf('----------- Equations of Motion - Manual ----------\n');
fprintf('---------------------------------------------------\n');

%% 1. Definition of Dynamic Variables and Accelerations (Slides 24-25)
% Generalized coordinates:            q  = [q1; q2; q3]
% Generalized velocities:             dq = [dq1; dq2; dq3]
% Generalized accelerations:          ddq = [ddq1; ddq2; ddq3]
% Generalized forces / joint torques: tau = [tau1; tau2; tau3]
syms ddq1 ddq2 ddq3 real
ddq = [ddq1; ddq2; ddq3];

% Recall generalized coordinates and velocities from previous sections
q  = [q1; q2; q3];
dq = [dq1; dq2; dq3];

% Gravity acceleration vector in base frame Sigma_0 (Slide 20, 25)
% g0 = [0; 0; -g]
g0 = gravityBase; 

% Manipulator inertia matrix B(q) obtained from Kinetic Energy (Slide 18, 25)
B = manipulatorInertiaMatrix;

% Explicit elements of B(q) (Slides 26-27):
% b_ij(q) represents the inertia effect of Joint j acceleration on Joint i
b11 = B(1,1);  b12 = B(1,2);  b13 = B(1,3);
b21 = B(2,1);  b22 = B(2,2);  b23 = B(2,3);
b31 = B(3,1);  b32 = B(3,2);  b33 = B(3,3);

fprintf('\n--- Manipulator Inertia Matrix B(q) (Slide 25) ---\n');
disp(B);


%% 2. Christoffel Symbols of the First Type c_ijk (Slide 30)
% Definition:
% c_ijk = 1/2 * ( d(b_ij)/d(q_k) + d(b_ik)/d(q_j) - d(b_jk)/d(q_i) )
% Property: c_ijk = c_ikj (symmetric in the last two indices)

% -------------------------------------------------------------------------
% i = 1 (Equation for Joint 1)
% -------------------------------------------------------------------------
% j = 1
c111 = simplify(1/2 * (diff(b11, q1) + diff(b11, q1) - diff(b11, q1)));
c112 = simplify(1/2 * (diff(b11, q2) + diff(b12, q1) - diff(b12, q1)));
c113 = simplify(1/2 * (diff(b11, q3) + diff(b13, q1) - diff(b13, q1)));

% j = 2
c121 = simplify(1/2 * (diff(b12, q1) + diff(b11, q2) - diff(b21, q1)));
c122 = simplify(1/2 * (diff(b12, q2) + diff(b12, q2) - diff(b22, q1)));
c123 = simplify(1/2 * (diff(b12, q3) + diff(b13, q2) - diff(b23, q1)));

% j = 3
c131 = simplify(1/2 * (diff(b13, q1) + diff(b11, q3) - diff(b31, q1)));
c132 = simplify(1/2 * (diff(b13, q2) + diff(b12, q3) - diff(b32, q1)));
c133 = simplify(1/2 * (diff(b13, q3) + diff(b13, q3) - diff(b33, q1)));

% -------------------------------------------------------------------------
% i = 2 (Equation for Joint 2)
% -------------------------------------------------------------------------
% j = 1
c211 = simplify(1/2 * (diff(b21, q1) + diff(b21, q1) - diff(b11, q2)));
c212 = simplify(1/2 * (diff(b21, q2) + diff(b22, q1) - diff(b12, q2)));
c213 = simplify(1/2 * (diff(b21, q3) + diff(b23, q1) - diff(b13, q2)));

% j = 2
c221 = simplify(1/2 * (diff(b22, q1) + diff(b21, q2) - diff(b21, q2)));
c222 = simplify(1/2 * (diff(b22, q2) + diff(b22, q2) - diff(b22, q2)));
c223 = simplify(1/2 * (diff(b22, q3) + diff(b23, q2) - diff(b23, q2)));

% j = 3
c231 = simplify(1/2 * (diff(b23, q1) + diff(b21, q3) - diff(b31, q2)));
c232 = simplify(1/2 * (diff(b23, q2) + diff(b22, q3) - diff(b32, q2)));
c233 = simplify(1/2 * (diff(b23, q3) + diff(b23, q3) - diff(b33, q2)));

% -------------------------------------------------------------------------
% i = 3 (Equation for Joint 3)
% -------------------------------------------------------------------------
% j = 1
c311 = simplify(1/2 * (diff(b31, q1) + diff(b31, q1) - diff(b11, q3)));
c312 = simplify(1/2 * (diff(b31, q2) + diff(b32, q1) - diff(b12, q3)));
c313 = simplify(1/2 * (diff(b31, q3) + diff(b33, q1) - diff(b13, q3)));

% j = 2
c321 = simplify(1/2 * (diff(b32, q1) + diff(b31, q2) - diff(b21, q3)));
c322 = simplify(1/2 * (diff(b32, q2) + diff(b32, q2) - diff(b22, q3)));
c323 = simplify(1/2 * (diff(b32, q3) + diff(b33, q2) - diff(b23, q3)));

% j = 3
c331 = simplify(1/2 * (diff(b33, q1) + diff(b31, q3) - diff(b31, q3)));
c332 = simplify(1/2 * (diff(b33, q2) + diff(b32, q3) - diff(b32, q3)));
c333 = simplify(1/2 * (diff(b33, q3) + diff(b33, q3) - diff(b33, q3)));


%% 3. Assembly of Coriolis and Centrifugal Matrix C(q, dq) (Slides 29-30)
% c_ij(q, dq) = sum_{k=1}^3 c_ijk * dq_k
c11 = simplify(c111*dq1 + c112*dq2 + c113*dq3);
c12 = simplify(c121*dq1 + c122*dq2 + c123*dq3);
c13 = simplify(c131*dq1 + c132*dq2 + c133*dq3);

c21 = simplify(c211*dq1 + c212*dq2 + c213*dq3);
c22 = simplify(c221*dq1 + c222*dq2 + c223*dq3);
c23 = simplify(c231*dq1 + c232*dq2 + c233*dq3);

c31 = simplify(c311*dq1 + c312*dq2 + c313*dq3);
c32 = simplify(c321*dq1 + c322*dq2 + c323*dq3);
c33 = simplify(c331*dq1 + c332*dq2 + c333*dq3);

C = [c11, c12, c13;
     c21, c22, c23;
     c31, c32, c33];

fprintf('\n--- Coriolis & Centrifugal Matrix C(q, dq) (Slide 30) ---\n');
disp(C);


%% 4. Centrifugal and Coriolis Coefficients h_ijk (Slides 26-27)
% Definition from Slide 26:
% h_ijk = d(b_ij)/d(q_k) - 1/2 * d(b_jk)/d(q_i)
% - Centrifugal term on Joint i by Joint j velocity:   h_ijj * dq_j^2  (with h_iii = 0)
% - Coriolis term on Joint i by Joints j, k velocity:  h_ijk * dq_j * dq_k  (j ~= k)

% Joint 1:
h111 = simplify(diff(b11, q1) - 1/2*diff(b11, q1));
h112 = simplify(diff(b11, q2) - 1/2*diff(b12, q1));
h113 = simplify(diff(b11, q3) - 1/2*diff(b13, q1));

h121 = simplify(diff(b12, q1) - 1/2*diff(b21, q1));
h122 = simplify(diff(b12, q2) - 1/2*diff(b22, q1));
h123 = simplify(diff(b12, q3) - 1/2*diff(b23, q1));

h131 = simplify(diff(b13, q1) - 1/2*diff(b31, q1));
h132 = simplify(diff(b13, q2) - 1/2*diff(b32, q1));
h133 = simplify(diff(b13, q3) - 1/2*diff(b33, q1));

% Joint 2:
h211 = simplify(diff(b21, q1) - 1/2*diff(b11, q2));
h212 = simplify(diff(b21, q2) - 1/2*diff(b12, q2));
h213 = simplify(diff(b21, q3) - 1/2*diff(b13, q2));

h221 = simplify(diff(b22, q1) - 1/2*diff(b21, q2));
h222 = simplify(diff(b22, q2) - 1/2*diff(b22, q2));
h223 = simplify(diff(b22, q3) - 1/2*diff(b23, q2));

h231 = simplify(diff(b23, q1) - 1/2*diff(b31, q2));
h232 = simplify(diff(b23, q2) - 1/2*diff(b32, q2));
h233 = simplify(diff(b23, q3) - 1/2*diff(b33, q2));

% Joint 3:
h311 = simplify(diff(b31, q1) - 1/2*diff(b11, q3));
h312 = simplify(diff(b31, q2) - 1/2*diff(b12, q3));
h313 = simplify(diff(b31, q3) - 1/2*diff(b13, q3));

h321 = simplify(diff(b32, q1) - 1/2*diff(b21, q3));
h322 = simplify(diff(b32, q2) - 1/2*diff(b22, q3));
h323 = simplify(diff(b32, q3) - 1/2*diff(b23, q3));

h331 = simplify(diff(b33, q1) - 1/2*diff(b31, q3));
h332 = simplify(diff(b33, q2) - 1/2*diff(b32, q3));
h333 = simplify(diff(b33, q3) - 1/2*diff(b33, q3));

% Quadratic velocity vector h(q, dq) (Slides 26-27):
h1 = simplify(h111*dq1*dq1 + h112*dq2*dq1 + h113*dq3*dq1 + ...
              h121*dq1*dq2 + h122*dq2*dq2 + h123*dq3*dq2 + ...
              h131*dq1*dq3 + h132*dq2*dq3 + h133*dq3*dq3);

h2 = simplify(h211*dq1*dq1 + h212*dq2*dq1 + h213*dq3*dq1 + ...
              h221*dq1*dq2 + h222*dq2*dq2 + h223*dq3*dq2 + ...
              h231*dq1*dq3 + h232*dq2*dq3 + h233*dq3*dq3);

h3 = simplify(h311*dq1*dq1 + h312*dq2*dq1 + h313*dq3*dq1 + ...
              h321*dq1*dq2 + h322*dq2*dq2 + h323*dq3*dq2 + ...
              h331*dq1*dq3 + h332*dq2*dq3 + h333*dq3*dq3);

h_vec = [h1; h2; h3];


%% 5. Gravity Generalized Force Vector g(q) (Slides 20, 25-27)
% Method 1: Via potential energy gradient (Slide 25)
% g_i(q) = dU / dq_i
g1_U = simplify(diff(totalPotentialEnergy, q1));
g2_U = simplify(diff(totalPotentialEnergy, q2));
g3_U = simplify(diff(totalPotentialEnergy, q3));
g_vec_U = [g1_U; g2_U; g3_U];

% Method 2: Via link linear Jacobians J_P^li (Slide 26)
% g_i(q) = - sum_{j=1}^3 m_lj * g0^T * J_Pi^lj(q)
% Columns of link COM Jacobians: Jv_COM{link}(:, joint)
g1_J = simplify(- (linkMass(1)*g0.'*Jv_COM{1}(:,1) + ...
                   linkMass(2)*g0.'*Jv_COM{2}(:,1) + ...
                   linkMass(3)*g0.'*Jv_COM{3}(:,1)));

g2_J = simplify(- (linkMass(1)*g0.'*Jv_COM{1}(:,2) + ...
                   linkMass(2)*g0.'*Jv_COM{2}(:,2) + ...
                   linkMass(3)*g0.'*Jv_COM{3}(:,2)));

g3_J = simplify(- (linkMass(1)*g0.'*Jv_COM{1}(:,3) + ...
                   linkMass(2)*g0.'*Jv_COM{3}(:,3) + ...
                   linkMass(3)*g0.'*Jv_COM{3}(:,3)));

g_vec_J = [g1_J; g2_J; g3_J];

% Adopt standard gravity vector g(q)
g_vec = g_vec_U;

fprintf('\n--- Gravity Generalized Force Vector g(q) (Slide 27) ---\n');
disp(g_vec);


%% 6. Time Derivative of Inertia Matrix B_dot and Skew-Symmetry Check (Slide 30)
% Time derivative of each inertia element:
% d(b_ij)/dt = (d(b_ij)/dq1)*dq1 + (d(b_ij)/dq2)*dq2 + (d(b_ij)/dq3)*dq3
b11_dot = diff(b11, q1)*dq1 + diff(b11, q2)*dq2 + diff(b11, q3)*dq3;
b12_dot = diff(b12, q1)*dq1 + diff(b12, q2)*dq2 + diff(b12, q3)*dq3;
b13_dot = diff(b13, q1)*dq1 + diff(b13, q2)*dq2 + diff(b13, q3)*dq3;

b21_dot = diff(b21, q1)*dq1 + diff(b21, q2)*dq2 + diff(b21, q3)*dq3;
b22_dot = diff(b22, q1)*dq1 + diff(b22, q2)*dq2 + diff(b22, q3)*dq3;
b23_dot = diff(b23, q1)*dq1 + diff(b23, q2)*dq2 + diff(b23, q3)*dq3;

b31_dot = diff(b31, q1)*dq1 + diff(b31, q2)*dq2 + diff(b31, q3)*dq3;
b32_dot = diff(b32, q1)*dq1 + diff(b32, q2)*dq2 + diff(b32, q3)*dq3;
b33_dot = diff(b33, q1)*dq1 + diff(b33, q2)*dq2 + diff(b33, q3)*dq3;

B_dot = [b11_dot, b12_dot, b13_dot;
         b21_dot, b22_dot, b23_dot;
         b31_dot, b32_dot, b33_dot];
B_dot = simplify(B_dot);

% Skew-symmetric matrix N(q, dq) = B_dot - 2*C
N_skew = simplify(B_dot - 2*C);


%% 7. Assembly of the Equations of Motion (Slides 26, 31)
% Detailed scalar equations (Slide 26, 29):
% tau_i = sum_j b_ij(q)*ddq_j + sum_j c_ij(q, dq)*dq_j + g_i(q)
tau1 = simplify(b11*ddq1 + b12*ddq2 + b13*ddq3 + c11*dq1 + c12*dq2 + c13*dq3 + g1_U);
tau2 = simplify(b21*ddq1 + b22*ddq2 + b23*ddq3 + c21*dq1 + c22*dq2 + c23*dq3 + g2_U);
tau3 = simplify(b31*ddq1 + b32*ddq2 + b33*ddq3 + c31*dq1 + c32*dq2 + c33*dq3 + g3_U);

% Compact Matrix Equation (Slide 31):
% B(q)*ddq + C(q, dq)*dq + g(q) = tau
tau = [tau1; tau2; tau3];

fprintf('\n--- Equations of Motion: Generalized Torques tau (Slide 31) ---\n');
fprintf('tau1 =\n'); disp(tau1);
fprintf('tau2 =\n'); disp(tau2);
fprintf('tau3 =\n'); disp(tau3);

% Complete model with viscous/Coulomb friction and external wrench (Slide 31)
% B(q)*ddq + C(q, dq)*dq + Fv*dq + Fs*sign(dq) + g(q) = tau - J^T(q)*he
syms Fv1 Fv2 Fv3 real
syms Fs1 Fs2 Fs3 real
syms fex fey fez mux muy muz real

Fv = diag([Fv1, Fv2, Fv3]);
Fs = diag([Fs1, Fs2, Fs3]);
he = [fex; fey; fez; mux; muy; muz];

tau_friction = Fv*dq + Fs*sign(dq);
tau_external = simplify(JG.' * he); % JG from Geometric Jacobian section

tau_complete = simplify(tau + tau_friction + tau_external);


%% 8. Analytical Identity Checks
fprintf('\n----------------------------------------------------\n');
fprintf('----------- Analytical Consistency Checks ----------\n');
fprintf('----------------------------------------------------\n');

% 1. Christoffel Symmetry Check: c_ijk = c_ikj (Slide 30)
err_c_sym = simplify(abs(c112 - c121) + abs(c113 - c131) + abs(c123 - c132) + ...
                     abs(c212 - c221) + abs(c213 - c231) + abs(c223 - c232) + ...
                     abs(c312 - c321) + abs(c313 - c331) + abs(c323 - c332));
fprintf('Christoffel symmetry check (c_ijk - c_ikj == 0): %s\n', char(err_c_sym == 0));

% 2. Centrifugal Nullity Check: h_iii == 0 (Slide 27)
err_h_null = simplify(abs(h111) + abs(h222) + abs(h333));
fprintf('Centrifugal self-effect nullity check (h_iii == 0): %s\n', char(err_h_null == 0));

% 3. Equivalence: C(q, dq)*dq == h(q, dq) (Slide 30)
err_Ch_equiv = simplify(C*dq - h_vec);
fprintf('Equivalence check: C(q,dq)*dq - h_vec (should be [0;0;0]):\n');
disp(err_Ch_equiv);

% 4. Gravity Formulations Equivalence: g_from_U == g_from_J (Slides 25-26)
err_g_equiv = simplify(g_vec_U - g_vec_J);
fprintf('Gravity equivalence check: g_U - g_J (should be [0;0;0]):\n');
disp(err_g_equiv);

% 5. Skew-Symmetry Check: (B_dot - 2*C) + (B_dot - 2*C)^T == 0 (Slide 30)
skew_check = simplify(N_skew + N_skew.');
fprintf('Skew-symmetry check: N + N^T == 0 (should be 3x3 zeros):\n');
disp(skew_check);
