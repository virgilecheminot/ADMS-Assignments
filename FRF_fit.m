%% Identification of the modal parameters

x0s = zeros(1, 5 * nO * R);

freq_range = [4, 5];
omega_range = 2 * pi * freq_range;


function x0 = init_x0(freqs, omega_vals, omega_i, damp_factor, G_exp, nO)
    x0 = zeros(1, 5 * nO);
    x0(1:nO) = omega_i;
    x0(nO+1:2*nO) = zeros(1, nO) + damp_factor;
    for i = 1:nO
        omega = omega_i(i);
        index = find(omega_vals > omega, 1);
        if ~isempty(index) && index <= length(freqs)
            x0(2*nO+i) = abs(G_exp(index));
        else
            x0(2*nO+i) = x0(2*nO+i-1);
        end
    end
end

% Linear array containing the initial guesses for the optimization
for r = 1:R
    x0s((r-1)*5*nO+1:r*5*nO) = init_x0(freqs, omega_vals, omega_i, damp_factor, G_exp(:, r), nO);
end

function FRF = numericalFRF(freqs, omega_i_roots, damp_factor, A, RL, RH)
    FRF = zeros(1, length(freqs));
    omega_vals = 2 * pi * freqs;
    % Main loop
    for i = 1:length(freqs)
        omega = omega_vals(i);
        for j = 1:length(omega_i_roots)
            omega_i = omega_i_roots(j);
            df = damp_factor(j);
            A_i = A(j);
            RL_i = RL(j);
            RH_i = RH(j);
            den = -omega^2 + 2j * df * omega_i * omega + omega_i^2;
            FRF(i) = FRF(i) + A_i / den + RL_i / omega^2 + RH_i;
        end
    end
end

function G_num = compute_G_num(freqs, nO, R, x)
    G_num = zeros(length(freqs), R);
    for r = 1:R
        idx = (r-1)*5*nO+1:r*5*nO;
        G_num(:, r) = numericalFRF(freqs, x(idx(1:nO)), x(idx(nO+1:2*nO)), x(idx(2*nO+1:3*nO)), x(idx(3*nO+1:4*nO)), x(idx(4*nO+1:5*nO)));
    end
end

G_num0 = compute_G_num(freqs, nO, R, x0s);

% Optimization
% err = @(x) (sum(sum((G_exp - compute_G_num(freqs, nO, R, x)) .* conj(G_exp - compute_G_num(freqs, nO, R, x)))))*1000;
% % err = @(x) norm(G_exp - compute_G_num(freqs, nO, R, x), 'fro');
% options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt');
% sol = lsqnonlin(@(x) err(x), x0s, [], [], options);

% G_sol = compute_G_num(freqs, nO, R, sol);