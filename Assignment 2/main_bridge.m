clear; close all; clc

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
scale_factor = 10;
figure('Position', [300 578 1311 480], 'Renderer', 'painters'); %#ok<*FGREN>
t = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:6
    mode = modes(:,i);
    nexttile;
    diseg2(mode, scale_factor, incid, l, gamma, posit, idb, xy);
    legend('off');
    title(sprintf('Mode %d - f = %.2f Hz', i, omega(i)/(2*pi)));
end

%% Task 3
B = [0.01 0.0075]';
A = zeros(2,2);
for ii=1:2
    A(ii,:) = [1/(2*omega(ii)) omega(ii)/2];
end
ab = A\B;
C = ab(1)*M + ab(2)*K;
[CFF, CFC, CCC] = freefree(C, ndof);