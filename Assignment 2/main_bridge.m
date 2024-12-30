clear; close all; clc
set(groot, 'defaultFigureRenderer', 'painters');

%% Initialization
% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('TRUSS_BRIDGE');
% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);
% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Task 1
% Define maximum angular frequency
Omax = 7*2*pi;

% Compute ratio for each element
ratios = zeros(nbeam, 1);
for ii = 1:nbeam
    omega_1 = (pi / l(ii))^2 * sqrt(EJ(ii) / m(ii));
    ratios(ii) = omega_1 / Omax;
end
if all(ratios >= 1.5)
    disp('All elements satisfy the condition omega_1,k / Omega_max >= 1.5');
else
    failed_ind = find(ratios < 1.5);
    error('Some elements fail the condition omega_1,k / Omega_max >= 1.5.\nFailed element(s): %s\nElement ratios: %s', num2str(failed_ind'), num2str(ratios(failed_ind)'));
end
% adding a node in the middle of every element solved the problem

%% Task 2
% Extract the free-free partition
[MFF, MFC, MCC] = freefree(M, ndof);
[KFF, KFC, KCC] = freefree(K, ndof);

% Compute the modes
[x0, omega_squared] = eig(MFF\KFF);
omega = diag(sqrt(omega_squared));
[omega, ind] = sort(omega);
modes = x0(:,ind);

% Plot the first 6 modes
% scale_factor = 10;
% figure('Position', [300 578 1311 480]);
% t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% for i = 1:6
%     mode = modes(:,i);
%     nexttile;
%     diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy);
%     legend('off');
%     title(sprintf('Mode %d - f = %.2f Hz', i, omega(i)/(2*pi)));
% end

%% Task 3
B = [0.01 0.0075]';
A = zeros(2,2);
for ii=1:2
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
ab = A\B;
C = ab(1)*M + ab(2)*K;
[CFF, CFC, CCC] = freefree(C, ndof);

%% Task 4
F0 = zeros(ndof, 1);
F0(idb(6,2)) = 1; % force applied at node 6 in the y direction
om = (0:0.01:7)*2*pi; 

% FRF of vertical displacement and acceleration at nodes 6, 7 and 31
X = zeros(ndof, length(om));
for ii = 1:length(om)
    A = -om(ii)^2*MFF + 1i*om(ii)*CFF + KFF;
    X(:,ii) = A\F0;
end

function plotFRF(om, X, idx, idb, plot_title)
    figure;
    legend_list = cell(length(idx), 1);
    for i = 1:length(idx)
        [r, ~] = find(idb==idx(i), 1, "first");
        legend_list{i} = sprintf('Node %d', r);
    end
    subplot(2,1,1)
    for ii = idx
        semilogy(om/(2*pi), abs(X(ii,:)), 'LineWidth', 1.5)
        hold on
    end
    xlabel('Frequency [Hz]')
    ylabel('Displacement [m]')
    title(plot_title)
    legend(legend_list, 'Location', 'best')
    grid on
    subplot(2,1,2)
    for ii = idx
        plot(om/(2*pi), angle(X(ii,:)), 'LineWidth', 1.5)
        hold on
    end
    xlabel('Frequency [Hz]')
    ylabel('Phase [rad]')
    title('Phase')
    legend(legend_list, 'Location', 'southwest')
    grid on
end

% displacement FRF
plotFRF(om, X, [idb(6,2), idb(7,2), idb(31,2)], idb, 'Vertical displacement FRF');

% acceleration FRF
plotFRF(om, -om.^2 .* X, [idb(6,2), idb(7,2), idb(31,2)], idb, 'Vertical acceleration FRF');

% FRF of shear force, bending moment and axial force at nodes 31 and 19
% Shear force
function coeffs = internal_coeffs(n_el, nod_i, nod_j, idb, l, gamma, X)
    L_el = l(n_el);
    idof_i = idb(nod_i,:);
    idof_j = idb(nod_j,:);
    lambda = [cos(gamma(n_el)) sin(gamma(n_el)) 0
             -sin(gamma(n_el)) cos(gamma(n_el)) 0
              0                0                1];
    Xi = lambda * X(idof_i,:);
    Xj = lambda * X(idof_j,:);

    a = Xi(2,:);
    b = Xj(3,:);
    c = -3/L_el^2*Xi(2,:) -3/L_el^2*Xj(2,:) -2/L_el*Xi(3,:) -1/L_el*Xj(3,:);
    d = 2/L_el^3*Xi(2,:) -2/L_el^3*Xj(2,:) +1/L_el^2*Xi(3,:) +1/L_el^2*Xj(3,:);
    coeffs = [a; b; c; d];
end
coeffs_41 = internal_coeffs(41, 30, 31, idb, l, gamma, X);
coeffs_23 = internal_coeffs(23, 30, 19, idb, l, gamma, X);
T1 = EJ(41) * 6 * coeffs_41(4,:);
T2 = EJ(23) * 6 * coeffs_23(4,:);
figure
subplot(2,1,1)
semilogy(om/(2*pi), abs([T1; T2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Shear force [N]')
title('Shear force FRF')
legend('Node 31', 'Node 19', 'Location', 'best')
grid on
subplot(2,1,2)
plot(om/(2*pi), angle([T1; T2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on

% Bending moment
M1 = EJ(41) * (2*coeffs_41(3,:) + 6*coeffs_41(4,:)*l(41));
M2 = EJ(23) * (2*coeffs_23(3,:) + 6*coeffs_23(4,:)*l(23));
figure
subplot(2,1,1)
semilogy(om/(2*pi), abs([M1; M2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Bending moment [Nm]')
title('Bending moment FRF')
legend('Node 31', 'Node 19', 'Location', 'best')
grid on
subplot(2,1,2)
plot(om/(2*pi), angle([M1; M2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on

% Axial force
N1 = EA(41) * coeffs_41(2,:);
N2 = EA(23) * coeffs_23(2,:);
figure
subplot(2,1,1)
semilogy(om/(2*pi), abs([N1; N2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Axial force [N]')
title('Axial force FRF')
legend('Node 31', 'Node 19', 'Location', 'best')
grid on
subplot(2,1,2)
plot(om/(2*pi), angle([N1; N2]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on

% Missing part c of task 4

%% Task 5
% Modal matrices
ii = 1:3; % consider the first 3 modes
Phi = modes(:,ii);
Mmod = Phi'*MFF*Phi;
Kmod = Phi'*KFF*Phi;
Cmod = Phi'*CFF*Phi;
Fmod = Phi'*F0;

% FRF in modal superposition approach
function X = FRFmod(Kmod, Mmod, Cmod, Fmod, Phi, om)
    Xmod = zeros(length(Fmod), length(om));
    for i=1:length(om)
        Xmod(:,i) = (Kmod - om(i)^2*Mmod + 1i*om(i)*Cmod)\Fmod;
    end
    X = Phi*Xmod;
end

Xmod = FRFmod(Kmod, Mmod, Cmod, Fmod, Phi, om);

nodes_no = [6, 7, 31];
function no = index_convert(i,N,M)
    [r,c] = ind2sub([N,M], i);
    no = sub2ind([M,N], c, r);
end

for k = 0:1
    figure
    for ii=1:length(nodes_no)
        plot_no = index_convert(2*ii-1, 2, 3);
        subplot(2,3,plot_no)
        semilogy(om/(2*pi), abs((-om.^2).^k.* X(nodes_no(ii),:)), 'LineWidth', 1.5)
        hold on
        semilogy(om/(2*pi), abs((-om.^2).^k.* Xmod(nodes_no(ii),:)), 'LineWidth', 1.5)
        xlabel('Frequency [Hz]')
        ylabel('Amplitude')
        str = 'Displacement';
        if k == 1
            str = 'Acceleration';
        end
        title(sprintf('Node %d - %s', nodes_no(ii), str))
        grid on
        legend('Direct FRF', 'Modal FRF', 'Location', 'best')
        plot_no = index_convert(2*ii, 2, 3);
        subplot(2,3,plot_no)
        plot(om/(2*pi), angle((-om.^2).^k.* X(nodes_no(ii),:)), 'LineWidth', 1.5)
        hold on
        plot(om/(2*pi), angle((-om.^2).^k.* Xmod(nodes_no(ii),:)), 'LineWidth', 1.5)
        xlabel('Frequency [Hz]')
        ylabel('Phase [rad]')
        title(sprintf('Node %d - Phase', nodes_no(ii)))
        grid on
    end
end

%% Task 6
% load TMD input file and assemble the structure
[file_i,xy,nnod,sizew,idb_TMD,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('TRUSS_BRIDGE_TMD');
[M_TMD,K_TMD]=assem(incid,l,m,EA,EJ,gamma,idb_TMD);
dis_stru(posit,l,gamma,xy,pr,idb_TMD,ndof);

fprintf('Natural frequency of the first mode: %.3f Hz\n', omega(1)/(2*pi));
whole_mass = sum(m.*l);
fprintf('Maximum TDM mass: %.3f kg\n', 0.02*whole_mass);

ml = 500;
kl = omega(1)^2 * ml;
cl = 0.3;

m1 = Mmod(1,1);
k1 = Kmod(1,1);
c1 = Cmod(1,1);