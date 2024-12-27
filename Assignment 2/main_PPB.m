clear; close all; clc

L = 1.2;        % [m]
E = 68e9;       % [Pa]
b = 40e-3;      % [m]
h = 8e-3;       % [m]
r = 2700;       % [kg/m^3]
m = r*b*h;      % [kg/m]
J = 1/12*b*h^3; % [m^4]
A = b*h;        % [m^2]
EA = E*A;       % [N]
EJ = E*J;       % [Nm^2]

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
[modes, omega_squared] = eig(MFF\KFF);
omega = diag(sqrt(omega_squared));

% Sort frequencies in ascending order
[omega, ind] = sort(omega);

% Sort mode shapes in ascending order
modes = modes(:,ind);

% Plot the first 3 modes
scale_factor = 1.5;
figure
for i = 1:3
    mode = modes(:,i);
    subplot(3,1,i)
    diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy);
    title(['Mode ', num2str(i)])
end

%% Damping matrix

B = [0.01 0.015 0.0098 0.012]';
A = zeros(4,2);
for ii=1:4
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
% ab = (A'*A)\A'*B;
ab = A\B;
CFF = ab(1)*MFF + ab(2)*KFF;

%% Frequency response function

F0 = zeros(ndof,1);
index = idb(7,2); % force applied at node 7 in the y direction
F0(index) = 1;
index_2 = idb(12,2);
F0_2 = zeros(ndof,1);
F0_2([index, index_2]) = [1, 1];

om = (0:1:500)*2*pi; % radiants per second

% Preallocate X
X = zeros(ndof, length(om));

for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

X_2 = zeros(ndof, length(om));
for ii=1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X_2(:,ii) = A\F0_2;
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

%% Modal superposition approach
% Modal matrices
ii = 1:3; % consider the first 3 modes
Phi = modes(:,ii);
Mmod = Phi'*MFF*Phi;
Kmod = Phi'*KFF*Phi;
Cmod = Phi'*CFF*Phi;
Fmod = Phi'*F0;

% FRF in modal superposition approach
Xmod = zeros(length(ii), length(om));
for i=1:length(om)
    Xmod(:,i) = (Kmod - om(i)^2*Mmod + 1i*om(i)*Cmod)\Fmod;
end
X_2m = Phi*Xmod;

idx4y = idb(4,2);

figure
subplot(2,1,1)
semilogy(om/(2*pi), abs(X(idx4y,:)), 'LineWidth', 1.5)
hold on
semilogy(om/(2*pi), abs(X_2m(idx4y,:)), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('Frequency Response Function')
grid on
legend('Direct approach', 'Modal superposition approach')
subplot(2,1,2)
plot(om/(2*pi), angle(X(idx4y,:)), 'LineWidth', 1.5)
hold on
plot(om/(2*pi), angle(X_2m(idx4y,:)), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on

%% Modal superposition considering multiple forces
Fmod_2 = Phi'*F0_2;
Xmod_2 = zeros(length(ii), length(om));
for i = 1:length(om)
    Xmod_2(:,i) = (Kmod - om(i)^2*Mmod + 1i*om(i)*Cmod)\Fmod_2;
end

X_2m_2 = Phi*Xmod_2;

figure
subplot(2,1,1)
semilogy(om/(2*pi), abs(X_2(idx4y,:)), 'LineWidth', 1.5)
hold on
semilogy(om/(2*pi), abs(X_2m_2(idx4y,:)), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Amplitude')
title('Frequency Response Function')
grid on
legend('Direct approach', 'Modal superposition approach')
subplot(2,1,2)
plot(om/(2*pi), angle(X_2(idx4y,:)), 'LineWidth', 1.5)
hold on
plot(om/(2*pi), angle(X_2m_2(idx4y,:)), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on