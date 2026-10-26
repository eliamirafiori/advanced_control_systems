clear all
close 
clc

%% DH table syms (DON'T EDIT)

syms d theta a alpha real

%% Homogeneous transformation matrix (DON'T EDIT)

A = [cos(theta) -sin(theta)*cos(alpha) sin(theta)*sin(alpha) a*cos(theta); ...
    sin(theta) cos(theta)*cos(alpha) -cos(theta)*sin(alpha) a*sin(theta); ...
    0 sin(alpha) cos(alpha) d; ...
    0 0 0 1];

%% Define the joint variables and link variables
syms q1 q2 q3 q4 q5 a1 a2 d3
 
% % % % % % % % % % % % % % % % % % % % % 
%                                       %
% DH_table = | d | theta | a | alpha |  %
%                                       %
% % % % % % % % % % % % % % % % % % % % %  
DH_table = [0 0 0 0;
    0 0 0 0;
    0 0 0 0;
            0 0 0 0];
         

%% Compute the T matrix 

A01 = subs(A,[d theta a alpha],DH_table(1,:));
A12 = subs(A,[d theta a alpha],DH_table(2,:));
A23 = subs(A,[d theta a alpha],DH_table(3,:));
    
%% Eventually transform the T matrix according to a given base frame and EE frame

Ab0 = [0 0 0 0;
    0 0 0 0;
    0 0 0 0;
       0 0 0 0]; % base frame transform

A3e = [0 0 0 0;
    0 0 0 0;
    0 0 0 0;
       0 0 0 0]; % ee frame transform

T03 = A01 * A12 * A23;

Tbe = Ab0*T03*A3e;

fprintf("Tbe matrix:\n");
disp(simplify(Tbe));

%% Inverse Kinematics

%% Visualization

robot = importrobot('00001_PRR_zzy.urdf');
showdetails(robot)

figure;
config = homeConfiguration(robot);
show(robot,config);

% config(1).JointPosition = INSERT_VALUE; 
% config(2).JointPosition = INSERT_VALUE;
% config(3).JointPosition = INSERT_VALUE;

show(robot,config);