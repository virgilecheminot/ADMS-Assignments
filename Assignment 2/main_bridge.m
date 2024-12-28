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