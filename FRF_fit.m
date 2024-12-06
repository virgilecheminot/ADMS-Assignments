close all; clc; clear;

% Load the experimental FRF data
load 'cantilever_FRF.mat'

f = freqs;
G = G_exp;
nFRF = size(G, 2);

% Find the peaks of the FRFs
[pks, locs] = findpeaks(abs(max(G_exp, [], 2)));

% Plot the experimental FRFs
figure;
subplot(2, 1, 1);
semilogy(f, abs(G));
hold on;
plot(f(locs), pks, 'ro', 'LineWidth', 1.5);
title(['Number of FRFs = ', num2str(size(G, 2))]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(2, 1, 2);
plot(f, angle(G));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;

%% Identification of the modal parameters

% Number of modes
omega_i = 2 * pi * f(locs);
mode_no = 1;

% Frequency range
freq_center = f(locs(mode_no));
fmin = 3;
fmax = 6;
f_range = f(f>=fmin & f<=fmax);
index_range = find(f>=fmin & f<=fmax);
disp(['Frequency range: ', num2str(fmin), ' Hz to ', num2str(fmax), ' Hz']);

%% First guess of the modal parameters
omega_i_guess = 2 * pi * freq_center;

[peak, r] = max(abs(G(locs(mode_no), :)));
f1_indices = find(f<freq_center & f>fmin);
G1_values = abs(G(f1_indices, r));
[~, idx1] = min(abs(G1_values - peak/sqrt(2)));
omega_1 = f(f1_indices(idx1)) * 2 * pi;
f2_indices = find(f>freq_center & f<fmax);
G2_values = abs(G(f2_indices, r));
[~, idx2] = min(abs(G2_values - peak/sqrt(2)));
omega_2 = f(f2_indices(idx2)) * 2 * pi;
damp_factor_guess = (omega_2^2 - omega_1^2) / (4 * omega_i_guess^2);

A_guess = zeros(nFRF, 1);
for i = 1:nFRF
    A_guess(i) = real(G(locs(mode_no), i) * (2j * damp_factor_guess * omega_i_guess^2));
end

RL_guess = zeros(nFRF, 1);
RH_guess = zeros(nFRF, 1);

function FRFs = numericalFRFs(freqs, nFRF, x)
    FRFs = zeros(length(freqs), nFRF);
    omega_vals = 2 * pi * freqs;
    omega_i = x(1);
    df_i = x(2);
    A_i = x(3:nFRF+2);
    RL_i = x(nFRF+3:2*nFRF+2);
    RH_i = x(2*nFRF+3:3*nFRF+2);
    % Main loop
    for i = 1:length(freqs)
        omega = omega_vals(i);
        for j = 1:nFRF
            den = -omega^2 + 2j * df_i * omega_i * omega + omega_i^2;
            FRFs(i, j) = FRFs(i, j) + A_i(j) / den + RL_i(j) / omega^2 + RH_i(j);
        end
    end
end

x0 = [omega_i_guess; damp_factor_guess; A_guess; RL_guess; RH_guess];

FRF_guess = numericalFRFs(f_range, nFRF, x0);

%% Optimization
% err = @(x) norm(G(index_range, :) - numericalFRFs(f_range, nFRF, x), 'fro')^2;
err = @(x) abs(reshape(G(index_range, :) - numericalFRFs(f_range, nFRF, x), [], 1));

options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
options.Display = "iter-detailed";
sol = lsqnonlin(@(x) err(x), x0, [], [], options);

FRF_opt = numericalFRFs(f_range, nFRF, sol);

figure;
subplot(2, 1, 1);
semilogy(f, abs(G), "LineWidth", 1.5);
hold on;
semilogy(f_range, abs(FRF_opt), 'ro');
title('Optimized FRFs');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([fmin-(fmax-fmin), fmax+(fmax-fmin)]);
subplot(2, 1, 2);
plot(f, angle(G), "LineWidth", 1.5);
hold on;
plot(f_range, angle(FRF_opt), 'ro');
xlabel('Frequency (Hz)');
ylabel('Phase');
xlim([fmin-(fmax-fmin), fmax+(fmax-fmin)]);

disp('Optimized parameters:');
disp(['omega_i = ', num2str(real(sol(1)) / (2 * pi)), ' Hz']);
disp(['damping factor = ', num2str(real(sol(2)) * 100), '%']);


%% Identified mode shapes

x_r = [0.2, 0.3, 0.45, 0.66, 0.9, 1.13];
A_r = sol(3:nFRF+2);

load 'cantilever_mode1.mat'

[~, idx_r] = arrayfun(@(x) min(abs(x_vals - x)), x_r);
err_shape = @(x) A_r * x - shape_1(idx_r).';
options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt', 'Display', 'none');
scale_factor = lsqnonlin(@(x) err_shape(x), 0, -100, 100, options);

figure;
plot(x_vals, shape_1, 'LineWidth', 1.5);
hold on;
for i = 1:length(x_r)
    plot(x_r(i), A_r(i)*scale_factor, 'ro', 'LineWidth', 1.5);
end
yline(0, '--', 'Color', 'k');
ylim([-max(abs(shape_1)), max(abs(shape_1))]);
title(['Mode shape n°1', newline, 'Scale factor = ', num2str(scale_factor)]);
xlabel('x (m)');
ylabel('Mode shape');
legend('True mode shape', 'Identified mode shape');
grid on;