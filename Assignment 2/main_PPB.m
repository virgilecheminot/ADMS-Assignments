clear; close all; clc

L = 1.2;        % [m]
E = 68e9;       % [Pa]
b = 40e-3;      % [m]
h = 8e-3;       % [m]
r = 2700;       % [kg/m^3]
m = r*b*h;      % [kg/m]
J = 1/12*b*h^3; % [m^4]
% A = b*h;        % [m^2]
% EA = E*A;       % [N]
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
diseg2(modes(:,1:3), scale_factor, incid, l, gamma, posit, idb, xy);

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

F0_7 = zeros(ndof,1);
F0_7(idb(7,2)) = 1; % force applied at node 7 in the y direction

om = (0:1:500)*2*pi; % radiants per second

function X = FRF(MFF, CFF, KFF, F0, om)
    ndof = size(MFF,1);
    X = zeros(ndof, length(om));
    for ii=1:length(om)
        A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
        X(:,ii) = A\F0;
    end
end

% Compute the FRF using the function
X = FRF(MFF, CFF, KFF, F0_7, om);

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
Fmod = Phi'*F0_7;

% FRF in modal superposition approach
function X = FRFmod(Kmod, Mmod, Cmod, Fmod, Phi, om)
    Xmod = zeros(length(Fmod), length(om));
    for i=1:length(om)
        Xmod(:,i) = (Kmod - om(i)^2*Mmod + 1i*om(i)*Cmod)\Fmod;
    end
    X = Phi*Xmod;
end

Xmod = FRFmod(Kmod, Mmod, Cmod, Fmod, Phi, om);

idx4y = idb(4,2);

function FRF_compare(om, X, Xmod, idx, title_str)
    figure
    subplot(2,1,1)
    semilogy(om/(2*pi), abs(X(idx,:)), 'LineWidth', 1.5)
    hold on
    semilogy(om/(2*pi), abs(Xmod(idx,:)), 'LineWidth', 1.5)
    xlabel('Frequency [Hz]')
    ylabel('Amplitude')
    title('Frequency Response Function')
    grid on
    title(sprintf(title_str))
    legend('Direct approach', 'Modal superposition approach')
    subplot(2,1,2)
    plot(om/(2*pi), angle(X(idx,:)), 'LineWidth', 1.5)
    hold on
    plot(om/(2*pi), angle(Xmod(idx,:)), 'LineWidth', 1.5)
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    grid on
end

FRF_compare(om, X, Xmod, idx4y, 'FRF at node 4 in the y direction considering force at node 7')

%% Modal superposition considering forces at nodes 7 and 12
F0_7_12 = zeros(ndof,1);
F0_7_12([idb(7,2) idb(12,2)]) = [1, 1]; % forces applied at nodes 7 and 12 in the y direction
Fmod_2 = Phi'*F0_7_12;

% FRF in direct approach
X_2 = FRF(MFF, CFF, KFF, F0_7_12, om);

% FRF in modal superposition approach
X_2mod = FRFmod(Kmod, Mmod, Cmod, Fmod_2, Phi, om);

FRF_compare(om, X_2, X_2mod, idx4y, 'FRF at node 4 in the y direction considering forces at nodes 7 and 12')

%% Modal superposition considering forces at nodes 2 and 12
F0_2_12 = zeros(ndof,1);
F0_2_12([idb(2,2) idb(12,2)]) = [1, 1];
Fmod_3 = Phi'*F0_2_12;

% FRF in direct approach
X_3 = FRF(MFF, CFF, KFF, F0_2_12, om);

% FRF in modal superposition approach
X_3mod = FRFmod(Kmod, Mmod, Cmod, Fmod_3, Phi, om);

FRF_compare(om, X_3, X_3mod, idx4y, 'FRF at node 4 in the y direction considering forces at nodes 2 and 12')

%% Modal superposition considering forces at nodes 2 and 12 in opposite directions
F0_2_12_opposite = zeros(ndof,1);
F0_2_12_opposite([idb(2,2) idb(12,2)]) = [-1, 1];
Fmod_4 = Phi'*F0_2_12_opposite;

% FRF in direct approach
X_4 = FRF(MFF, CFF, KFF, F0_2_12_opposite, om);

% FRF in modal superposition approach
X_4mod = FRFmod(Kmod, Mmod, Cmod, Fmod_4, Phi, om);

FRF_compare(om, X_4, X_4mod, idx4y, 'FRF at node 4 in the y direction considering forces at \nnodes 2 and 12 in opposite directions')

%% Distributed load on elements
p0 = 100; % [N/m]
p0G = [0 p0]';
FG = zeros(3*nnod,1);
elements = 2;

for ii = elements
    Lk = l(ii);
    g = gamma(ii);
    p0L = [cos(g) sin(g);
          -sin(g) cos(g)] * p0G;
    FkL = [Lk/2 ;  0   ; 0 ; Lk/2 ;  0   ; 0] * p0L(1) + ...
          [ 0   ; Lk/2 ; 0 ;  0   ; Lk/2 ; 0] * p0L(2);
    FkG = [cos(g) -sin(g) 0 0       0      0;
           sin(g)  cos(g) 0 0       0      0;
           0       0      1 0       0      0;
           0       0      0 cos(g) -sin(g) 0;
           0       0      0 sin(g)  cos(g) 0;
           0       0      0 0       0      1] * FkL;
    Ek = zeros(6,3*nnod);
    for jj = 1:6
        Ek(jj,incid(ii,jj)) = 1;
    end
end
FG = FG + Ek' * FkG;

%% Static deflection
xF = KFF\FG(1:ndof);

% Plot the static deflection
figure
diseg2(xF, 100, incid, l, gamma, posit, idb, xy); % scale factor 100
title('Static deflection')