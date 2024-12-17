clear; close all; clc

L = 1.2; % [m]
E = 68e9; % [Pa]
b = 40e-3; % [m]
h = 8e-3; % [m]
r = 2700; % [kg/m^3]
m = r*b*h;      % [kg/m]
J = 1/12*b*h^3; % [m^4]
A = b*h; % [m^2]
% EA = E*A; % [N]
EJ = E*J; % [Nm^2]

Omax = 100*2*pi;
a = 1.5;
Lmax = sqrt(pi^2/a/Omax * sqrt(EJ/m));

% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('ADMS_FEM_01');

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);