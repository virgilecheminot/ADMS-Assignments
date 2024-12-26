clear; close all; clc

m = 9.75;       % [kg/m]
EJ = 1.34e4;    % [Nm^2]
ml = 10;        % [kg] (lumped mass)
Jl = 1;         % [kg*m^2] (lumped moment of inertia)
kx = 2e6;       % [N/m] (Spring x)
ky = 3e6;       % [N/m] (Spring y)
k = 4e6;        % [N/m] (Spring)
h1 = 0.01;      % [-] (Damping ratio 1st mode)
h2 = 0.015;     % [-] (Damping ratio 2nd mode)

Omax = 100*2*pi;
a = 1.5;
Lmax = sqrt(pi^2/a/Omax * sqrt(EJ/m));
disp(['Maximum element legnth recommended: ', num2str(Lmax)]);

% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('ADMS_FEM_03_mass');

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Concentrated elements contributions
% lumped mass
Ml_local = [ml  0   0
            0   ml  0
            0   0  Jl];

Eml = zeros(3,size(M,1)); % expansion matrix
indices = idb(5,:);
Eml(:,indices) = eye(3);
Eml_8 = zeros(3,size(M,1)); % expansion matrix for the 8th node
Eml_8(:,idb(8,:)) = eye(3);

Ml = Eml'*Ml_local*Eml + Eml_8'*Ml_local*Eml_8;
M_tot = M + Ml;

% concentrated springs (only at the 5th node)
Kl_local = [kx 0
            0 ky];

Ek = zeros(2,size(M,1));
Ek(:,idb(5,1:2)) = eye(2);
Kl = Ek'*Kl_local*Ek;

K_tot = K + Kl;

% concentrated springs between nodes 2 and 8
Kk_local = [1 0 0 -1 0 0]' * k * [1 0 0 -1 0 0];
g = 3*pi/2;
lambda_k = [cos(g)  sin(g) 0
            -sin(g) cos(g) 0
            0       0      1];
Lambda_k = [lambda_k zeros(3)
            zeros(3) lambda_k];
Kk_G = Lambda_k' * Kk_local * Lambda_k;
idofn2 = idb(2,:);
idofn8 = idb(8,:);
idofk = [idofn2 idofn8];

K_tot(idofk,idofk) = K_tot(idofk,idofk) + Kk_G;

%% Extract the free-free partition
[MFF, MFC, MCC] = freefree(M_tot,ndof);
[KFF, KFC, KCC] = freefree(K_tot,ndof);

%% Compute the modes
[x0, omega_squared] = eig(MFF\KFF);
omega = diag(sqrt(omega_squared));
[omega, ind] = sort(omega);
modes = x0(:,ind);

% Plot the first 3 modes
scale_factor = 1.5;
figure
for i = 1:3
    mode = modes(:,i);
    diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy);
end

% plot the mode frequencies
figure
stem(omega/(2*pi), 'filled', 'LineWidth', 1.5);
yscale('log');
xlabel('Mode number');
ylabel('Frequency [Hz]');
title('Mode frequencies');
grid on;

%% Damping matrix
B = [h1 h2]';
A = zeros(2,2);
for ii=1:2
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
ab = A\B;
C = ab(1)*M_tot + ab(2)*K_tot;
[CFF, CFC, CCC] = freefree(C,ndof);

%% Frequency response function
F0 = zeros(ndof,1);
index = idb(4,2);
F0(index) = 1;

om = (0:0.2:100)*2*pi;

X = zeros(ndof, length(om));
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

% Plot the FRF
figure
subplot(2,1,1);
semilogy(om/(2*pi), abs(X(idb(5,2),:)), 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('Amplitude');
title('Frequency response function');
grid on;
subplot(2,1,2);
plot(om/(2*pi), angle(X(idb(5,2),:)), 'LineWidth', 1.5);
xlabel('Frequency [Hz]');
ylabel('Phase [rad]');
grid on;

