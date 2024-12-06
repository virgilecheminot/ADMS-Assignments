close all; clc; clear;

load 'LRV_wheel_FRF.mat'

f = freq;
FRF = frf;
N = size(FRF, 2);

[pks, locs] = findpeaks(abs(max(FRF, [], 2)), 'SortStr', 'descend', 'NPeaks', 4);
% sort peaks according from lower to higher frequency
[~, idx] = sort(f(locs));
locs = locs(idx);
pks = pks(idx);

figure;
subplot(3, 1, 1);
semilogy(f, abs(FRF));
hold on;
plot(f(locs), pks, 'ro', 'LineWidth', 1.5);
title(['Number of FRFs = ', num2str(N)]);
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;
subplot(3, 1, 2);
plot(f, angle(FRF));
xlabel('Frequency (Hz)');
ylabel('Phase');
grid on;
subplot(3, 1, 3);
plot(f, cohe);
xlabel('Frequency (Hz)');
ylabel('Coherence');
grid on;

omega_i = 2 * pi * f(locs);

modes = 1:2;

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

solutions = zeros(3*N+2, length(modes));

for mode_no = modes
    disp(['Mode number: ', num2str(mode_no)]);
    freq_center = f(locs(mode_no));
    disp(['Guess frequency: ', num2str(freq_center), ' Hz']);
    fmin = freq_center - 100;
    fmax = freq_center + 100;
    f_range = f(f>=fmin & f<=fmax);
    index_range = find(f>=fmin & f<=fmax);
    disp(['Frequency range: ', num2str(fmin), ' Hz to ', num2str(fmax), ' Hz']);
    
    omega_i_guess = 2 * pi * freq_center;
    
    [peak, r] = max(abs(FRF(locs(mode_no), :)));
    f1_indices = find(f<freq_center & f>fmin);
    G1_values = abs(FRF(f1_indices, r));
    [~, idx1] = min(abs(G1_values - peak/sqrt(2)));
    omega_1 = f(f1_indices(idx1)) * 2 * pi;
    f2_indices = find(f>freq_center & f<fmax);
    G2_values = abs(FRF(f2_indices, r));
    [~, idx2] = min(abs(G2_values - peak/sqrt(2)));
    omega_2 = f(f2_indices(idx2)) * 2 * pi;
    damp_factor_guess = (omega_2^2 - omega_1^2) / (4 * omega_i_guess^2);
    disp(['Guess damping factor: ', num2str(damp_factor_guess), '%']);
    
    A_guess = zeros(N, 1);
    for i = 1:N
        A_guess(i) = real(FRF(locs(mode_no), i) * (2j * damp_factor_guess * omega_i_guess^2));
    end

    RL_guess = zeros(N, 1);
    RH_guess = zeros(N, 1);

    x0 = [omega_i_guess; damp_factor_guess; A_guess; RL_guess; RH_guess];

    err = @(x) abs(reshape(FRF(index_range, :) - numericalFRFs(f_range, N, x), [], 1));
    % err = @(x) abs(reshape((FRF(index_range, :) - numericalFRFs(f_range, N, x)) .* conj(FRF(index_range, :) - numericalFRFs(f_range, N, x)), [], 1));
    % err = @(x) norm(FRF(index_range, :) - numericalFRFs(f_range, N, x), 'fro')^2;

    options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
    options.Display = 'final-detailed';
    solutions(:,mode_no) = lsqnonlin(@(x) err(x), x0, [], [], options);

    FRF_opt = numericalFRFs(f_range, N, solutions(:,mode_no));

    figure;
    subplot(2, 1, 1);
    semilogy(f_range, abs(FRF(index_range, :)));
    hold on;
    semilogy(f_range, abs(FRF_opt), 'r--');
    title(['Mode number: ', num2str(mode_no), newline, 'Frequency: ', num2str(solutions(1,mode_no)/(2*pi)), ' Hz, Damping factor: ', num2str(solutions(2,mode_no)*100), '%']);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude');
    grid on;
    subplot(2, 1, 2);
    plot(f_range, angle(FRF(index_range, :)));
    hold on;
    plot(f_range, angle(FRF_opt), 'r--');
    xlabel('Frequency (Hz)');
    ylabel('Phase');
    grid on;

    disp('------------------');
end

theta_r = -90:15:-90+15*(N-1);
figure;
for i = 1:length(modes)
    A_r = solutions(3:N+2, i)/max(abs(solutions(3:N+2, i)));
    polarplot(deg2rad(theta_r'), A_r/3+1, 'Marker', 'o', 'LineWidth', 1.5, 'DisplayName', ['Mode ', num2str(modes(i))]);
    hold on;
end
title('Mode shapes');
polarplot(deg2rad(0:1:360), ones(1, 361), 'k--', 'DisplayName', 'Original shape');
legend();
grid on;