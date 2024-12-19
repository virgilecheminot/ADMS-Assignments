clear; close all; clc

L = 1.2; % [m]
E = 68e9; % [Pa]
b = 40e-3; % [m]
h = 8e-3; % [m]
r = 2700; % [kg/m^3]
m = r*b*h;      % [kg/m]
J = 1/12*b*h^3; % [m^4]
% A = b*h; % [m^2]
% EA = E*A; % [N]
EJ = E*J; % [Nm^2]

Omax = 500*2*pi;
a = 1.5;
Lmax = sqrt(pi^2/a/Omax * sqrt(EJ/m));
disp(['Maximum element legnth recommended: ', num2str(Lmax)]);

% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('ADMS_FEM_01');

% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);

% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

% extract the free-free partition
[MFF, MFC, MCC] = freefree(M,ndof);
[KFF, KFC, KCC] = freefree(K,ndof);

%% Compute the modes
[modes, omega2] = eig(MFF\KFF);
omega = diag(sqrt(omega2));

% Sort frequencies in ascending order
[omega_sorted, omega_sorted_indices] = sort(omega);

% Sort mode shapes in ascending order
modes_sorted = modes(:,omega_sorted_indices);

% Plot the first 3 modes
scale_factor = 1.5;
figure
for i = 1:3
    mode = modes_sorted(:,i);
    diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy);
end

%% Damping matrix

B = [0.01 0.015 0.0098 0.012]';
A = zeros(4,2);
for ii=1:4
    A(ii,:) = [1/(2*omega_sorted(ii)) omega_sorted(ii)/2];
end
% ab = (A'*A)\A'*B;
ab = A\B;
CFF = ab(1)*MFF + ab(2)*KFF;

%% Frequency response function

F0 = zeros(ndof,1);
index = idb(7,2); % force applied at node 7 in the y direction
F0(index) = 1;

om = (0:1:500)*2*pi;

% Preallocate X
X = zeros(ndof, length(om));

for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

% Plot the FRF

index_interest = [idb(4,2), idb(4,3), idb(7,2)];

figure
for ii=1:length(index_interest)
    subplot(2,1,1)  
    semilogy(om/(2*pi), abs(X(index_interest(ii),:)), 'LineWidth', 1.5)
    xlabel('Frequency [Hz]')
    ylabel('Amplitude')
    title('Frequency Response Function')
    grid on
    hold on
    subplot(2,1,2)
    plot(om/(2*pi), angle(X(index_interest(ii),:)), 'LineWidth', 1.5)  
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    grid on
    hold on
end