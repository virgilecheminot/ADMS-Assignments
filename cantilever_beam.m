tic; % Start timing

close all; clc; clear;

% Data of the reference structure
L = 1.2; % m
h = 8e-3; % m
b = 40e-3; % m
rho = 2700; % kg/m^3
E = 68e9; % Pa
J = b * h^3 / 12; % m^4
m = rho * b * h; % kg/m

% Boundary conditions
H = @(gamma) [  1,  0,  1,   0;
                0,  1,  0,   1;
                -cos(gamma*L), -sin(gamma*L),  cosh(gamma*L),  sinh(gamma*L);
                sin(gamma*L), -cos(gamma*L),  sinh(gamma*L),  cosh(gamma*L)];

% Characteristic equation
chareq = @(gamma) det(H(gamma));

% Solve the characteristic equation for gamma
gamma_vals = linspace(0, 50, 10000);
omega_vals = gamma_vals.^2 * sqrt(E * J / (rho * b * h));
chareq_vals = arrayfun(chareq, gamma_vals);

% Numerically find the roots of the characteristic equation
gamma_roots = zeros(1, length(gamma_vals) - 1);
root_count = 0;
for i = 1:length(gamma_vals)-1
    if chareq_vals(i) * chareq_vals(i+1) < 0
        root_count = root_count + 1;
        gamma_roots(root_count) = fzero(chareq, [gamma_vals(i), gamma_vals(i+1)]);
    end
end
gamma_roots = gamma_roots(1:root_count);

% Convert gamma roots to frequencies
omega_i = gamma_roots.^2 * sqrt(E * J / (rho * b * h));
f_i = omega_i / (2 * pi);

disp(['First 4 natural frequencies: ', num2str(f_i(1:4)), ' Hz']);

%% Mode shapes

% Compute the mode shapes
x_vals = linspace(0, L, 2000);
mode_shapes = zeros(4, length(gamma_roots));

for i = 1:length(gamma_roots)
    g_i = gamma_roots(i);
    H_i = H(g_i);
    H_hat = H_i(2:4, 2:4);
    N = H_i(2:4, 1);
    H_hat_inv = pinv(H_hat);
    w_hat = -H_hat_inv * N;
    w_i = [1; w_hat];
    mode_shapes(:, i) = w_i; % Normalize the mode shape
end

function phi = phi_shape(gamma, x)
    phi = [cos(gamma * x); sin(gamma * x); cosh(gamma * x); sinh(gamma * x)];
end

% Plot the mode shapes
for i = 1:1
    figure;
    g_i = gamma_roots(i);
    shape_i = mode_shapes(:, i)' * phi_shape(g_i, x_vals);
    plot(x_vals, shape_i);
    hold on;
    line([0, L], [0, 0], 'LineStyle', '--', 'Color', 'k');
    title(['Mode Shape n°', num2str(i), ' - f = ', num2str(f_i(i)), ' Hz']);
    xlabel('x (m)');
    ylabel('Mode Shape');
    grid on;
    ylim([-max(abs(shape_i)), max(abs(shape_i))]);
end

% Save mode 1 shape
shape_1 = mode_shapes(:, 1)' * phi_shape(gamma_roots(1), x_vals);
save 'cantilever_mode1.mat' x_vals shape_1;

%% Frequency response function
x_j = 0.2; % m
x_k = 1.2; % m
damp_factor = 0.01;

freqs = linspace(0, 200, 10000);

function FRF = computeFRF(freqs, gamma_roots, mode_shapes, x_j, x_k, x_vals, damp_factor, E, J, rho, b, h, m)
    % Preallocation of FRF
    FRF = zeros(1, length(freqs));

    % Calculation of constants
    omega_vals = 2 * pi * freqs;
    omega_i_roots = gamma_roots.^2 * sqrt(E * J / (rho * b * h));

    % Pre-calculation of phi_shape for x_j, x_k and x_vals
    phi_j_vals = arrayfun(@(g_i) phi_shape(g_i, x_j), gamma_roots, 'UniformOutput', false);
    phi_k_vals = arrayfun(@(g_i) -phi_shape(g_i, x_k), gamma_roots, 'UniformOutput', false); % Negative direction for x_k because de force is applied downwards
    phi_i_vals = arrayfun(@(g_i) phi_shape(g_i, x_vals), gamma_roots, 'UniformOutput', false);

    % Calculation of modal masses
    m_i_vals = arrayfun(@(j) trapz(x_vals, m * (mode_shapes(:, j)' * phi_i_vals{j}).^2), 1:length(gamma_roots));

    % Main loop
    for i = 1:length(freqs)
        omega = omega_vals(i);
        for j = 1:length(gamma_roots)
            omega_i = omega_i_roots(j);
            w_i = mode_shapes(:, j);
            phi_j = w_i' * phi_j_vals{j};
            phi_k = w_i' * phi_k_vals{j};
            m_i = m_i_vals(j);
            num = (phi_j * phi_k / m_i);
            den = -omega^2 + 2j * damp_factor * omega_i * omega + omega_i^2;
            FRF(i) = FRF(i) + num / den;
        end
    end
end

FRF = computeFRF(freqs, gamma_roots, mode_shapes, x_j, x_k, x_vals, damp_factor, E, J, rho, b, h, m);

% Plot the FRF
figure;
subplot(2, 1, 1);
semilogy(freqs, abs(FRF));
title('FRF at x_j = 0.2 m and x_k = 1.2 m');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(2, 1, 2);
plot(freqs, angle(FRF));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;

%% Considering N FRFs

x_js = [0.2, 0.3, 0.45, 0.66, 0.9, 1.13];
x_ks = 1.2;

G_exp = zeros(length(freqs), length(x_js) * length(x_ks));

for i = 1:length(x_js)
    for j = 1:length(x_ks)
        x_j = x_js(i);
        x_k = x_ks(j);
        FRF = computeFRF(freqs, gamma_roots, mode_shapes, x_j, x_k, x_vals, damp_factor, E, J, rho, b, h, m);
        G_exp(:, (i-1)*length(x_ks) + j) = FRF.';
    end
end

% Plot the FRFs
figure;
subplot(2, 1, 1);
semilogy(freqs, abs(G_exp));
title(['Number of FRFs = ', num2str(size(G_exp, 2))]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(2, 1, 2);
plot(freqs, angle(G_exp));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;

% Save the "experimental" FRFs
save 'cantilever_FRF.mat' freqs G_exp;

toc; % End timing and display elapsed time