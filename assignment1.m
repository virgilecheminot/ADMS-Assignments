close all; clc;

% Data of the reference structure
L = 1.2; % m
h = 8e-3; % m
b = 40e-3; % m
rho = 2700; % kg/m^3
E = 68e9; % Pa
J = b * h^3 / 12; % m^4
m = rho * b * h * L; % kg

% Boundary conditions
syms gamma
H = [ 1 0 1 0
      0 1 0 1
      -cos(gamma*L) -sin(gamma*L) cosh(gamma*L) sinh(gamma*L)
      sin(gamma*L) -cos(gamma*L) sinh(gamma*L) cosh(gamma*L)];

% Characteristic equation
chareq = det(H);

% Solve the characteristic equation for gamma
gamma_vals = linspace(0, 50, 10000);
chareq_vals = double(subs(chareq, gamma, gamma_vals));

% Find the roots of the characteristic equation
gamma_roots = zeros(1, length(gamma_vals) - 1);
root_count = 0;
for i = 1:length(gamma_vals)-1
    if chareq_vals(i) * chareq_vals(i+1) < 0
        root_count = root_count + 1;
        gamma_roots(root_count) = fzero(@(gamma) double(subs(chareq, gamma)), [gamma_vals(i), gamma_vals(i+1)]);
    end
end
gamma_roots = gamma_roots(1:root_count);

% Convert gamma roots to frequencies
omega_i = gamma_roots.^2 * sqrt(E * J / (rho * b * h));
f_i = omega_i / (2 * pi);

%% Mode shapes

% Compute the mode shapes
x_vals = linspace(0, L, 1000);
mode_shapes = zeros(4, length(gamma_roots));

for i = 1:length(gamma_roots)
    g_i = gamma_roots(i);
    H_i = double(subs(H, gamma, g_i));
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
    line([0, L], [0, 0], 'LineStyle', '--', 'Color', 'k');
    plot(x_vals, shape_i);
    title(['Mode Shape for f = ', num2str(f_i(i)), ' Hz']);
    xlabel('x (m)');
    ylabel('Mode Shape');
end

%% Frequency response function
x_j = 0.2; % m
x_k = 1.2; % m
damp_factor = 0.01;

freqs = linspace(0, 200, 2000);

function FRF = computeFRF(freqs, gamma_roots, mode_shapes, x_j, x_k, x_vals, damp_factor, E, J, rho, b, h, m)
    % Preallocation of FRF
    FRF = zeros(1, length(freqs));
    
    % Calculation of constants
    omega_vals = 2 * pi * freqs;
    omega_i_vals = gamma_roots.^2 * sqrt(E * J / (rho * b * h));
    
    % Pre-calculation of phi_shape for x_j, x_k and x_vals
    phi_j_vals = arrayfun(@(g_i) phi_shape(g_i, x_j), gamma_roots, 'UniformOutput', false);
    phi_k_vals = arrayfun(@(g_i) phi_shape(g_i, x_k), gamma_roots, 'UniformOutput', false);
    phi_i_vals = arrayfun(@(g_i) phi_shape(g_i, x_vals), gamma_roots, 'UniformOutput', false);
    
    % Calculation of modal masses
    m_i_vals = arrayfun(@(j) trapz(x_vals, m * (mode_shapes(:, j)' * phi_i_vals{j}).^2), 1:length(gamma_roots));
    
    % Main loop
    for i = 1:length(freqs)
        omega = omega_vals(i);
        FRF(i) = 0;
        for j = 1:length(gamma_roots)
            omega_i = omega_i_vals(j);
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
title('Frequency Response Function');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(2, 1, 2);
plot(freqs, angle(FRF));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;

%% Considering N FRFs

x_js = [0.2, 0.4, 0.6, 0.8, 1.0];
x_ks = [1.2, 0.9, 0.6, 0.3];

G_exp = zeros(length(x_js) * length(x_ks), length(freqs));

for i = 1:length(x_js)
    for j = 1:length(x_ks)
        x_j = x_js(i);
        x_k = x_ks(j);
        FRF = computeFRF(freqs, gamma_roots, mode_shapes, x_j, x_k, x_vals, damp_factor, E, J, rho, b, h, m);
        G_exp((i-1)*length(x_ks) + j, :) = FRF;
    end
end

% Plot the FRFs
figure;
subplot(2, 1, 1);
semilogy(freqs, abs(G_exp));
title('Frequency Response Function');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(2, 1, 2);
plot(freqs, angle(G_exp));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;
