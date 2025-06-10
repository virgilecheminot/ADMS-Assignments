clear; close all; clc
set(groot, 'defaultFigureRenderer', 'painters');

%% Initialization
% load the input file and assemble the structure
[file_i,xy,nnod,sizew,idb,ndof,incid,l,gamma,m,EA,EJ,posit,nbeam,pr]=loadstructure('TRUSS_BRIDGE');
% draw the structure
dis_stru(posit,l,gamma,xy,pr,idb,ndof);
% assemble mass and stiffness matrices
[M,K]=assem(incid,l,m,EA,EJ,gamma,idb);

%% Task 1 → OK
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

%% Task 2 → OK
% Extract the free-free partition
[MFF, MFC, MCC] = freefree(M, ndof);
[KFF, KFC, KCC] = freefree(K, ndof);

% Compute the modes
[x0, omega_squared] = eig(MFF\KFF);
omega = diag(sqrt(omega_squared));
[omega, ind] = sort(omega);
modes = x0(:,ind);

% Plot the first 6 modes
scale_factor = 10;
figure('Position', [300 578 1311 480]);
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:6
    mode = modes(:,i);
    nexttile;
    diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy, 0.05);
    legend('off');
    title(sprintf('Mode %d - f = %.2f Hz', i, omega(i)/(2*pi)));
end

%% Task 3 → OK
B = [0.01 0.0075]';
A = zeros(2,2);
for ii=1:2
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
ab = A\B;
C = ab(1)*M + ab(2)*K;
[CFF, CFC, CCC] = freefree(C, ndof);

%% Task 4 → OK
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

    a = Xi(1,:);
    b = (Xj(1,:)-Xi(1,:))/L_el;
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

% Reaction forces
R = zeros(3*nnod-ndof, length(om));
for ii = 1:length(om)
    R(:,ii) = (-MFC'*om(ii)^2 + 1i*CFC'*om(ii) + KFC')*X(:,ii);
end

R1x = R(idb(1,1)-ndof,:);
R1y = R(idb(1,2)-ndof,:);
R13y = R(idb(13,2)-ndof,:);

figure
subplot(2,1,1)
semilogy(om/(2*pi), abs([R1x; R1y; R13y]), 'LineWidth', 1.5)
ylim([1e-4, 1e2])
xlabel('Frequency [Hz]')
ylabel('Reaction force [N]')
title('Reaction force FRF')
legend('Node 1 - x', 'Node 1 - y', 'Node 13 - y', 'Location', 'best')
grid on
subplot(2,1,2)
plot(om/(2*pi), angle([R1x; R1y; R13y]), 'LineWidth', 1.5)
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
grid on

%% Task 5 → OK
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

%% Task 7 → OK
bridge_m = sum(m.*l);
fprintf('Total mass of the structure: %.0f kg\n', bridge_m);
deck_elements = 1:12;
deck_nodes = 1:13;
function FG = distributed_load(deck_elements, incid, l, gamma, m, nnod)
    FG = zeros(3*nnod, 1);
    for ii = deck_elements
        Lk = l(ii);
        g = gamma(ii);
        p0 = -m(ii)*9.81;
        p0G = [0 p0]';
        p0L = [cos(g) sin(g);
            -sin(g) cos(g)] * p0G;
        FkL = [Lk/2 ;  0   ;    0    ; Lk/2 ;  0   ;    0    ] * p0L(1) + ...
              [ 0   ; Lk/2 ; Lk^2/12 ;  0   ; Lk/2 ; -Lk^2/12] * p0L(2);
        FkG = [cos(g) -sin(g) 0 0       0      0;
            sin(g)  cos(g) 0 0       0      0;
            0       0      1 0       0      0;
            0       0      0 cos(g) -sin(g) 0;
            0       0      0 sin(g)  cos(g) 0;
            0       0      0 0       0      1] * FkL;
        Ek = zeros(6,3*nnod);
        for jj = 1:6
            Ek(jj, incid(ii,jj)) = 1;
        end
        FG = FG + Ek'*FkG;
    end
end

function FG = weight_load(nodes, M, nnod, idb)
    G = zeros(3*nnod, 1);
    for ii = idb(nodes,2)'
        G(ii) = -9.81;
    end
    FG = M*G;
end

FG_deck = distributed_load(deck_elements, incid, l, gamma, m, nnod);
xF_deck = KFF\FG_deck(1:ndof);

FG_deck2 = weight_load(deck_nodes, M, nnod, idb);
xF_deck2 = KFF\FG_deck2(1:ndof);

FG = distributed_load(1:nbeam, incid, l, gamma, m, nnod);
xF = KFF\FG(1:ndof);

FG2 = weight_load(1:nnod, M, nnod, idb);
xF2 = KFF\FG2(1:ndof);

ratio = bridge_m/sum(m(deck_elements).*l(deck_elements));

xF_comp = ratio * xF_deck;

figure
diseg2([xF_deck, xF], 75, incid, l, gamma, posit, idb, xy, 0.05);
title('Vertical displacement of the structure due to the distributed load');

%% Task 6
% Choosing the parameters for the TMD
Phi_TMD = modes(:,1);
Mmod_TMD = Phi_TMD'*MFF*Phi_TMD;
Kmod_TMD = Phi_TMD'*KFF*Phi_TMD;
Cmod_TMD = Phi_TMD'*CFF*Phi_TMD;
FRF_noTMD = 1./(-om.^2*Mmod_TMD + 1i*om*Cmod_TMD + Kmod_TMD);

ml = 333;
hl = 0.2;
reduction = zeros(length(ml), length(hl));

figure;
for i = 1:length(ml)
    kl = omega(1)^2 * ml(i);
    for j = 1:length(hl)
        cl = hl(j) * 2 * ml(i) * omega(1);
        MM = [Mmod_TMD 0; 0 ml(i)];
        CC = [Cmod_TMD+cl -cl; -cl cl];
        KK = [Kmod_TMD+kl -kl; -kl kl];
        FF = [1 0]';
        a = zeros(2, length(om));
        for ii = 1:length(om)
            a(:,ii) = (-om(ii)^2*MM + 1i*om(ii)*CC + KK) \ FF;
        end
        FRF_TMD = a(1,:);
        reduction(i,j) = (1 - max(abs(FRF_TMD))/max(abs(FRF_noTMD)))*100;
        subplot(length(ml), length(hl), (i-1)*length(hl)+j);
        semilogy(om/(2*pi), abs(FRF_noTMD), 'LineWidth', 1.1);
        hold on;
        semilogy(om/(2*pi), abs(FRF_TMD), 'LineWidth', 1.1);
        title_str = sprintf('m_2 = %d kg, h_2 = %.2f, reduction = %.2f', ml(i), hl(j), reduction(i,j));
        title(title_str);
    end
end

% figure;
% surf(hl, ml, reduction);
% hold on;
% surf(hl, ml, ones(size(reduction))*50, 'FaceAlpha', 0.5);
% xlabel('h_2 [m]');
% ylabel('m_2 [kg]');
% zlabel('Reduction [%]');
% title('Reduction in the maximum displacement due to the TMD');

% Choosing the optimal parameters
ml = [333, 999, 5982.2];
hl = [0.2, 0.14, 0.12036];
kl = ml * omega(1)^2;
cl = 2 * hl .* ml * omega(1);

% Applying the TMD
[file_i2,xy2,nnod2,sizew2,idb2,ndof2,incid2,l2,gamma2,m2,EA2,EJ2,posit2,nbeam2,pr2]=loadstructure('TRUSS_BRIDGE_TMD');
[M_TMD,K_TMD]=assem(incid2,l2,m2,EA2,EJ2,gamma2,idb2);
dis_stru(posit2,l2,gamma2,xy2,pr2,idb2,ndof2);

figure;
for i = 1:3
    E_ml = zeros(1,3*nnod2);
    E_ml(idb2(37,2)) = 1;
    M_TMD = M_TMD + E_ml'*ml(i)*E_ml;

    K_k_G = [0 -1 0 0 1 0]' * kl(i) * [0 -1 0 0 1 0];

    idofn37 = idb2(37,:);
    idofn7 = idb2(7,:);
    idofk = [idofn37 idofn7];

    K_TMD(idofk, idofk) = K_TMD(idofk, idofk) + K_k_G;

    C_TMD = ab(1)*M_TMD + ab(2)*K_TMD;
    C_c_G = [0 -1 0 0 1 0]' * cl(i) * [0 -1 0 0 1 0];
    C_TMD(idofk, idofk) = C_TMD(idofk, idofk) + C_c_G;

    MFF_TMD = freefree(M_TMD, ndof2);
    KFF_TMD = freefree(K_TMD, ndof2);
    CFF_TMD = freefree(C_TMD, ndof2);

    [x0_TMD, omega_squared_TMD] = eig(MFF_TMD\KFF_TMD);
    omega_TMD = diag(sqrt(omega_squared_TMD));
    [omega_TMD, ind_TMD] = sort(omega_TMD);
    modes_TMD = x0_TMD(:,ind_TMD);

    F0_TMD = zeros(ndof2,1);
    F0_TMD(idb2(6,2)) = 1;

    X_TMD = zeros(ndof2, length(om));
    for ii = 1:length(om)
        A = -om(ii)^2*MFF_TMD + 1i*om(ii)*CFF_TMD + KFF_TMD;
        X_TMD(:,ii) = A\F0_TMD;
    end

    reduction_TMD = (1 - max(abs(X_TMD(idb2(6,2),200:300)))/max(abs(X(idb2(6,2),200:300))))*100;
    semilogy(om/(2*pi), abs(X_TMD(idb2(6,2),:)), 'LineWidth', 1.5)
    hold on
end
semilogy(om/(2*pi), abs(X(idb2(6,2),:)), 'LineWidth', 1.5)
title_str = sprintf('Frequency Response Function');
title(title_str)
xlabel('Frequency [Hz]')
ylabel('Amplitude')
xlim([2 3])
legend('TMD 20%', 'TMD 50%', 'TMD 82.5%', 'No TMD', 'Location', 'best')

% Plot the first mode
figure
diseg2(modes_TMD(:,1), 10, incid2, l2, gamma2, posit2, idb2, xy2, 0.05);
title('First mode of the structure with TMD');

%% Task 8 → OK
d = 27.5; % m
for k = 1:3
    for i = 1:2
        Vik = omega(i)*d/(2*pi*k); % critical velocity
        fprintf('Critical velocity for mode %d and k = %d: %.1f km/h\n', i, k, Vik*3.6);
    end
end