clear; close all; clc

m = 9.75;       % [kg/m]
EJ = 1.34e4;    % [Nm^2]
EA = 2.57e7;    % [N]
ml = 10;        % [kg] (lumped mass)
Jl = 1;         % [kg*m^2] (lumped moment of inertia)
kx = 2e6;       % [N/m] (Spring x)
ky = 3e6;       % [N/m] (Spring y)
k = 4e6;        % [N/m] (Spring)
h1 = 0.01;      % [-] (Damping ratio 1st mode)
h2 = 0.015;     % [-] (Damping ratio 2nd mode)
coef = 1.5;     % [-] (Safety factor)

Omax = 100*2*pi;
a = 1.5;
Lmax = sqrt(pi^2/a/Omax * sqrt(EJ/m));
disp(['Maximum element legnth recommended: ', num2str(Lmax)]);

% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('ADMS_FEM_03');

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

% adding lumped mass
Ml_local = [ml  0   0
            0   ml  0
            0   0   Jl];

Eml = zeros(3,size(M,1));
indices = idb(5,:);
Eml(:,indices) = eye(3);

Ml = Eml'*Ml_local*Eml;
M_tot = M + Ml;