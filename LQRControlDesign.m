% Author: Antonio Sanchez
% Date: 27-07-2017
%
% Description:
% This script designs a LQR controller for an underactuated ship.
%
% References:
% Fossen, T. (2011). Handbook of marine craft hydrodynamics and motion control.
% Sanchez, A. (2017). Advanced Control of an Unmanned Surface Vehicle.
close all; clear variables; clc;

load('Halcyon.mat');    % Load first set of parameters.
load('HalcyonABC.mat'); % Load second set of paramenters.

%% 1. Determine state-space model matrices (6DOF).
% p175 (Fossen, 2011).
% where x = [eta_p, nu]. eta_p is vessel parallel coordinate system.
Minv = vesselABC.Minv;
G    = vesselABC.G;
D    = vessel.Bv;

A = [zeros(6,6) eye(6,6); -Minv*G -Minv*D];
B = [zeros(6,6); Minv];
E = [zeros(6,6); Minv];
C = eye(12,12);

%% 2. Account for vessel underactuation by including a matrix Kact.
Pw = vessel.main.rho;               % Water density [kg/m^3]
Ar = vessel.rudder.area;            % Area [m^2]
u0 = 5;                             % Surge velocity (u) in steady-state [m/s]
dCl = vessel.rudder.dCL;            % Partial derivative of non-dimensional lift [1/rad]
LCG = vessel.rudder.CG(1);          % Longitudinal distance between rudders and CoG [m]
VCG = vessel.rudder.CG(3);          % Vertical distance between rudders and CoG [m]

% Relation between control actuation forces and commanded thrust and rudder inputs.
Cy = 0.5*Pw*Ar*(u0^2)*dCl;          % Constant relating thrust and sway moment
Kprop = [2 0 0 0 0 0]';             % Matrix of propeller constants
Krudd = [0 Cy 0 -Cy*VCG 0 Cy*LCG]'; % Matrix of rudder constants
Kact = [Kprop Krudd];
sys = ss(A,B*Kact,C,0);             % Build state-space model

%% 3. Decomposition in Roll and Sway-Yaw Subsystems + Surge System.
% p160 (Fossen, 2011)
k = [8 12 6 10 4 7];    % This variable represents the order of states [v, r, psi, p, phi, u]'
An = sys.A(k,k);        % Reorder previous state space model to reveal roll-sway-yaw coupled dynamics + surge
Bn = sys.B(k,:);        % Corresponding input matrix with thrust and rudder
Cn = eye(6,6);
sysn = ss(An,Bn,Cn,0);  % Build state-space model

%% 4. Reference tracking for the Roll and Sway-Yaw subsystem + Surge System.
Cr = [0 0 0 0 0 1; 0 0 1 0 0 0]; % Tracking outputs.
% Construction of "zero" and "eye" matrices in kg.
% Row1: zeros(# of states, # of tracking outputs).
% Row2: eye(# of tracking outputs, # of tracking outputs).
kxku = [An Bn; Cr zeros(2,2)]\[zeros(6,2); eye(2,2)];
kx = kxku(1:6, :);
ku = kxku(7:end, :);

%% 5. LQR feedback controller for the Roll and Sway-Yaw subsystems.
Ag = [An zeros(6,2); -Cr zeros(2,2)];       % Augmented state matrix for offset-free tracking control with integral action.
Bg = [Bn; zeros(2,2)];                      % Augmented input matrix for offset-free tracking control with integral action.
Q = diag([1 500 500 1 1 1e4 15e4 1e-05]);   % Gains for velocity tracking.
R = eye(2,2);

[Kki,~,~] = lqr(Ag,Bg,Q,R);
K = Kki(:,1:6);                             % Optimal controller.
ki = -Kki(:,7:end);                         % Integral gains.
kr = ku + K*kx;                             % Feedforward controller.